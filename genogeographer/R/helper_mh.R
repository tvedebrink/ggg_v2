## Helpers for the microhaplotype (MH) profile module.
##
## 'MH_profiles_probs.Rds' (shipped in inst/extdata) contains genotype probabilities per MH
## marker and meta-population, precomputed offline by the scripts in the 'microhaps' analysis
## project (see 'precompute_likelihood_deviance.R'). Fitting those models involves expensive
## graph-structure learning and simulation, so this module only ever performs fast table
## look-ups and products against the precomputed file - it never refits anything.
##
## 'MH_snp_map.Rds' maps each individual SNP locus (as genotyped, e.g. position '344229') to
## the MH marker it belongs to and its ordering within that marker's motif (built from the
## 'haps' sheet of the pedro reference workbook - see 'data.R'/'fit_structures.R').
##
## Accepted input formats (see inst/extdata/mh_sample_*.csv for examples):
##  - long, two columns for alleles:  sample_name;locus;A;B      (A/B = haplotype copy 1/2)
##  - long, one column for alleles:   sample_name;locus;allele   (two-character genotype, e.g. "TG": T = copy 1, G = copy 2)
##  - wide, two rows per sample:      sample_name;<locus1>;...   (one row per haplotype copy)
##  - wide, one row per sample:       sample_name;<locus1>;...   (two-character genotype per cell, first character = copy 1)
## All formats are phased: "TG" and "GT" are different genotypes (the first character is always the
## allele on the first haplotype copy), as the haplotype-based model requires known linkage phase.

mh_cache <- new.env(parent = emptyenv())
mh_cache$count_scale <- list()

#' Load the precomputed MH genotype probabilities (cached)
#' @noRd
mh_profiles_probs <- function(){
  if(is.null(mh_cache$probs)){
    path <- system.file("extdata", "MH_profiles_probs.Rds", package = ggg_package)
    if(!nzchar(path)) stop("Could not find 'MH_profiles_probs.Rds' in the package 'extdata' folder.")
    mh_cache$probs <- readRDS(path)
  }
  mh_cache$probs
}

#' Combined long table of MH genotype probabilities across all markers
#'
#' Columns: `MH`, `geno` (the two haplotype motifs, e.g. `"CTTGGC;CTTGGC"`), and one column
#' per meta-population with the corresponding genotype probability.
#' @param probs_list A list as returned by `mh_profiles_probs()` (or a user-supplied equivalent).
#' @noRd
mh_build_probs_table <- function(probs_list){
  probs_list %>% purrr::map_dfr("probs_df", .id = "MH")
}

#' Combined long table of MH genotype probabilities across all markers (cached, bundled data)
#' @noRd
mh_geno_probs_table <- function(){
  if(is.null(mh_cache$probs_table)){
    mh_cache$probs_table <- mh_build_probs_table(mh_profiles_probs())
  }
  mh_cache$probs_table
}

#' Combined long table of MH genotype deviances across all markers - same shape as
#' `mh_build_probs_table()`/`mh_geno_probs_table()`, but from each marker's `devi_df`.
#' @noRd
mh_build_devi_table <- function(probs_list){
  probs_list %>% purrr::map_dfr("devi_df", .id = "MH")
}

#' Combined long table of MH genotype deviances across all markers (cached, bundled data)
#' @noRd
mh_geno_devi_table <- function(){
  if(is.null(mh_cache$devi_table)){
    mh_cache$devi_table <- mh_build_devi_table(mh_profiles_probs())
  }
  mh_cache$devi_table
}

#' Ensure each per-meta outlier null model has a `cdf` (ecdf of `sims`), adding one if missing
#' @noRd
mh_prepare_out_metas <- function(x){
  lapply(x, function(m){ if(is.null(m$cdf)) m$cdf <- stats::ecdf(m$sims); m })
}

#' Load the precomputed per-meta-population outlier null models (cached)
#'
#' Each element is a lean version of a `molic` `outlier_model` (only the fields needed for a
#' p-value/critical-value/z-score: `sims`, `mu`, `sigma`, `cv`, `alpha`), fit offline on the
#' FULL joint SNP graphical model (see `outlier_detect.R`'s `out_metas`) - the summed per-MH
#' deviance (see `mh_profile_likelihood()`) is calibrated against these `sims`, exactly as
#' `sampling_profiles.R` does (`pval(x = out_metas[[meta]], dz = <summed devi>)`).
#' @noRd
mh_out_metas <- function(){
  if(is.null(mh_cache$out_metas)){
    path <- system.file("extdata", "MH_out_metas.Rds", package = ggg_package)
    if(!nzchar(path)) stop("Could not find 'MH_out_metas.Rds' in the package 'extdata' folder.")
    mh_cache$out_metas <- mh_prepare_out_metas(readRDS(path))
  }
  mh_cache$out_metas
}

#' Prepare the reference haplotype data used for the log-likelihood confidence interval
#'
#' `ci_data` is `MH_ci_data.Rds` (see `data-raw/MH_ci_data.R`): per meta-population and MH marker the
#' reference haplotype motifs, the motif of each reference haplotype (row), and the cliques (`sign` = 1)
#' and separators (`sign` = -1) of the fitted graph as positions within the motif. This derives, for every
#' margin, the cell of each reference haplotype and the cell counts.
#' @noRd
mh_ci_prepare <- function(ci_data){
  lapply(ci_data, lapply, function(e){
    margins <- Map(function(pos, sign){
      key <- vapply(strsplit(e$motifs, ""), function(ch) paste(ch[pos], collapse = ""), character(1))
      cells <- unique(key)
      rows <- match(key, cells)[e$row_code]
      list(pos = pos, sign = sign, cells = cells, rows = rows, n = tabulate(rows, length(cells)))
    }, e$pos, e$sign)
    list(N = length(e$row_code), margins = margins)
  })
}

#' Reference haplotype data for the log-likelihood confidence interval (bundled, cached)
#' @noRd
mh_ci_data <- function(){
  if(is.null(mh_cache$ci_data)){
    path <- system.file("extdata", "MH_ci_data.Rds", package = ggg_package)
    if(!nzchar(path)) stop("Could not find 'MH_ci_data.Rds' in the package 'extdata' folder.")
    mh_cache$ci_data <- mh_ci_prepare(readRDS(path))
  }
  mh_cache$ci_data
}

#' Names of the meta-populations available in the precomputed MH data
#' @noRd
mh_meta_populations <- function(){
  setdiff(names(mh_geno_probs_table()), c("MH", "geno"))
}

#' MH marker to SNP-locus map (cached): columns `MH`, `position` (order within the motif), `snp`
#' @noRd
mh_snp_map <- function(){
  if(is.null(mh_cache$snp_map)){
    path <- system.file("extdata", "MH_snp_map.Rds", package = ggg_package)
    if(!nzchar(path)) stop("Could not find 'MH_snp_map.Rds' in the package 'extdata' folder.")
    mh_cache$snp_map <- readRDS(path)
  }
  mh_cache$snp_map
}

#' Normalise a SNP locus identifier for robust matching (strips whitespace/decimal noise)
#' @noRd
normalize_locus_id <- function(x){
  x <- trimws(as.character(x))
  num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(num), x, format(num, scientific = FALSE, trim = TRUE))
}

#' Read a file's raw cell values with the header row kept as ordinary data (as text)
#'
#' `rio::import(header = FALSE)` only honours `header` for delimited text files - for Excel
#' files it always treats the first row as column names, mangling numeric/duplicate headers.
#' This dispatches to `readxl` directly (as text, so numbers/genotypes are never auto-converted)
#' for Excel files, and to `rio` for everything else.
#' @noRd
mh_read_raw <- function(path){
  ext <- tolower(rio::get_ext(path))
  if(ext %in% c("xlsx", "xls")) suppressMessages(readxl::read_excel(path, col_names = FALSE, col_types = "text"))
  else rio::import(path, header = FALSE)
}

#' Detect which of the four supported (phased) MH input formats a file uses
#'
#' In every format the alleles are phased: for a two-character genotype such as `"TG"` the first
#' character is the allele on haplotype copy 1 and the second the allele on copy 2.
#' @param path Path to the uploaded CSV/Excel file.
#' @return One of `"long_two_col"`, `"long_one_col"`, `"wide_two_row"`, `"wide_one_row"`.
#' @noRd
detect_mh_format <- function(path){
  raw <- mh_read_raw(path)
  if(nrow(raw) < 2) stop("The file must contain a header row and at least one data row.")
  header <- raw[1, ] %>% unlist(use.names = FALSE) %>% as.character() %>% trimws()
  header_lower <- tolower(header)
  if(length(header) >= 3 && header_lower[2] == "locus"){
    if(length(header) == 4 && all(header_lower[3:4] %in% c("a", "b"))) return("long_two_col")
    if(length(header) == 3 && header_lower[3] == "allele") return("long_one_col")
    stop("Unrecognised long-format header: expected 'sample_name;locus;A;B' or 'sample_name;locus;allele'.")
  }
  ## Wide format: distinguish by the width of a data cell (1 character = one haplotype copy per row;
  ## 2 characters = both haplotype copies in one cell)
  cell <- as.character(raw[2, 2])
  if(nchar(trimws(cell)) == 1) return("wide_two_row")
  if(nchar(trimws(cell)) == 2) return("wide_one_row")
  stop("Unrecognised wide-format data: expected single-character alleles (two rows per sample) or two-character genotypes (one row per sample).")
}

#' Read a phased MH profile file into canonical long format, regardless of the specific input format
#'
#' @param path Path to the uploaded CSV/Excel file.
#' @return A tibble with columns `sample_name`, `locus`, `allele1`, `allele2` (the alleles of
#'   haplotype copy 1 and 2).
#' @noRd
read_mh_profile_file <- function(path){
  format <- detect_mh_format(path)
  switch(format,
    long_two_col = parse_mh_long_two_col(path),
    long_one_col = parse_mh_long_one_col(path),
    wide_two_row = parse_mh_wide_two_row(path),
    wide_one_row = parse_mh_wide_one_row(path)
  )
}

#' @noRd
parse_mh_long_two_col <- function(path){
  raw <- rio::import(path, header = TRUE) %>% tibble::as_tibble()
  tibble::tibble(
    sample_name = as.character(raw[[1]]),
    locus = normalize_locus_id(raw[["locus"]]),
    allele1 = toupper(trimws(raw[["A"]])),
    allele2 = toupper(trimws(raw[["B"]]))
  )
}

#' @noRd
parse_mh_long_one_col <- function(path){
  raw <- rio::import(path, header = TRUE) %>% tibble::as_tibble()
  geno <- toupper(trimws(raw[["allele"]]))
  tibble::tibble(
    sample_name = as.character(raw[[1]]),
    locus = normalize_locus_id(raw[["locus"]]),
    allele1 = substr(geno, 1, 1),
    allele2 = substr(geno, 2, 2)
  )
}

#' @noRd
parse_mh_wide_two_row <- function(path){
  raw <- mh_read_raw(path)
  header <- raw[1, ] %>% unlist(use.names = FALSE) %>% as.character() %>% trimws()
  body <- raw[-1, , drop = FALSE]
  names(body) <- make.unique(header)
  loci <- normalize_locus_id(header[-1])
  sample_raw <- as.character(body[[1]])
  base_sample <- substr(sample_raw, 1, nchar(sample_raw) - 1L)
  copy <- ave(seq_along(base_sample), base_sample, FUN = seq_along)
  if(!all(table(base_sample) == 2)){
    stop("Wide two-row format expects exactly two rows (one per haplotype copy) per sample.")
  }
  purrr::map_dfr(seq_along(loci), function(j){
    col <- body[[j + 1L]]
    tibble::tibble(sample_name = base_sample, copy = copy, locus = loci[j], allele = toupper(trimws(col)))
  }) %>%
    tidyr::pivot_wider(names_from = copy, values_from = allele, names_prefix = "allele")
}

#' @noRd
parse_mh_wide_one_row <- function(path){
  raw <- mh_read_raw(path)
  header <- raw[1, ] %>% unlist(use.names = FALSE) %>% as.character() %>% trimws()
  body <- raw[-1, , drop = FALSE]
  names(body) <- make.unique(header)
  loci <- normalize_locus_id(header[-1])
  sample_name <- as.character(body[[1]])
  purrr::map_dfr(seq_along(loci), function(j){
    geno <- toupper(trimws(body[[j + 1L]]))
    tibble::tibble(
      sample_name = sample_name, locus = loci[j],
      allele1 = substr(geno, 1, 1), allele2 = substr(geno, 2, 2)
    )
  })
}

#' Haplotype-motif pair of a phased MH genotype
#' @param allele1,allele2 Character vectors of the alleles of haplotype copy 1 and 2 at each locus,
#'   ordered by `position`.
#' @return A two-element character vector `c(motif1, motif2)`.
#' @noRd
mh_motif_pair <- function(allele1, allele2){
  c(paste(allele1, collapse = ""), paste(allele2, collapse = ""))
}

#' Look up a motif pair in a combined `MH`/`geno` table, trying both allele orders
#' @noRd
mh_lookup_candidate <- function(pair, tbl, meta_cols){
  row <- tbl[tbl$geno == paste(pair[1], pair[2], sep = ";"), meta_cols, drop = FALSE]
  if(nrow(row) == 0) row <- tbl[tbl$geno == paste(pair[2], pair[1], sep = ";"), meta_cols, drop = FALSE]
  if(nrow(row) == 0) return(NULL)
  unlist(row[1, ])
}

#' Small-count correction factor of a cell's first-order standard error
#'
#' In the extended database a cell has count `n_ext = X + e`: `X` reference haplotypes (random) plus `e`
#' haplotypes of the profile (fixed). With `X ~ Binomial(N, n_ext / (N + 2))`, the first-order variance of
#' `log(n_ext)` is `N (1 - pi) / (n_ext (N + 2))`, which overestimates the exact variance of `log(X + e)` for
#' small counts (4x at `n_ext = 1`). Returns the square root of exact over first-order variance (1 for large counts).
#' @param n_ext,e Extended and profile counts of the cells (vectors of equal length).
#' @param N Number of reference haplotypes.
#' @noRd
mh_count_scale <- function(n_ext, e, N){
  vapply(seq_along(n_ext), function(i){
    if(n_ext[i] >= 100) return(1)
    key <- paste(n_ext[i], e[i], N)
    if(is.null(s <- mh_cache$count_scale[[key]])){
      pi_hat <- n_ext[i] / (N + 2); x <- 0:N
      w <- stats::dbinom(x, N, pi_hat); l <- log(x + e[i])
      s <- sqrt((sum(w * l^2) - sum(w * l)^2) / (N * (1 - pi_hat) / (n_ext[i] * (N + 2))))
      mh_cache$count_scale[[key]] <- s
    }
    s
  }, numeric(1))
}

#' First-order (delta method) variance of a profile's log10 likelihood in one meta-population
#'
#' The probability of an MH genotype is computed conditionally on the profile, i.e. with the profile's
#' two haplotypes added to the reference haplotypes (`x + x0`). The variance is derived in the same
#' framework: the extended database of `N + 2` haplotypes is the plug-in distribution, from which the
#' `N` reference haplotypes are drawn i.i.d. while the profile stays fixed. With
#' `g_i = sum_a c_a s_a 1[haplotype i of the extended database in cell a] / p_a` (`c_a` = +1/-1 times the
#' number of the profile's haplotypes in the clique/separator cell, `p_a` the cell proportion in the extended
#' database, `s_a` the small-count factor of `mh_count_scale()`), `Var(log P) = N / (N + 2)^2 * Var_i(g_i)`.
#' `g_i` is summed over all markers, which keeps the dependence between markers. Every cell of the profile
#' occurs in the extended database, so the variance is always defined.
#' @param prep The meta-population's slice of `mh_ci_prepare()`.
#' @param items One element per MH marker used: `list(MH, pair)` with the profile's motif pair.
#' @return The variance of `log10_likelihood` (`NA` if the marker data are not available).
#' @noRd
mh_log10_var <- function(prep, items){
  g <- NULL; N <- NULL
  for(it in items){
    e <- prep[[it$MH]]
    if(is.null(e)) return(NA_real_)
    if(is.null(g)){ N <- e$N; g <- numeric(N + 2) }
    haps <- strsplit(it$pair, "")
    for(M in e$margins){
      key <- vapply(haps, function(ch) paste(ch[M$pos], collapse = ""), character(1))
      cells <- c(M$cells, setdiff(key, M$cells))
      cell <- match(key, cells)
      rows <- c(M$rows, cell)
      n <- tabulate(rows, length(cells))
      coef <- numeric(length(cells))
      coef[cell[1]] <- coef[cell[1]] + M$sign
      coef[cell[2]] <- coef[cell[2]] + M$sign
      occ <- unique(cell)
      coef[occ] <- coef[occ] * mh_count_scale(n[occ], tabulate(cell, length(cells))[occ], N)
      g <- g + (coef / (n / (N + 2)))[rows]
    }
  }
  if(is.null(g)) return(NA_real_)
  N / (N + 2)^2 * mean((g - mean(g))^2) / log(10)^2
}

#' Standard error of the log10 likelihood for every sample and meta-population
#' @param per_mh,stats The per-(sample, MH) table and the matching list of per-marker results
#'   (`pair`) built in `mh_profile_likelihood()`.
#' @return A list with `se` (`sample_name`, `meta`, `se_log10`) and `unobserved`, an empty tibble (the
#'   variance is always defined, see `mh_log10_var()`; kept until the module is updated).
#' @noRd
mh_likelihood_ci <- function(per_mh, stats, metas, ci_data, snp_map){
  has_pair <- !vapply(stats, function(x) is.null(x$pair), logical(1))
  se <- list()
  for(s in unique(per_mh$sample_name)){
    idx <- which(per_mh$sample_name == s & has_pair)
    items <- lapply(idx, function(i) list(MH = per_mh$MH[i], pair = stats[[i]]$pair))
    for(meta in metas){
      se[[length(se) + 1]] <- data.frame(sample_name = s, meta = meta, se_log10 = sqrt(mh_log10_var(ci_data[[meta]], items)))
    }
  }
  list(se = tibble::as_tibble(do.call(rbind, se)),
       unobserved = tibble::tibble(sample_name = character(), meta = character(), MH = character(),
                                   haplotype = character(), snps = character()))
}

#' Compute per-sample, per-meta-population likelihoods and outlier test of a phased MH profile
#'
#' For each sample and MH marker with complete SNP coverage in the input, reconstructs the
#' haplotype motif pair and looks up its precomputed genotype probability/deviance. No
#' model fitting or simulation is performed - this is table look-ups and products only.
#' A microhaplotype marker is used only when ALL of its constituent SNPs are present for a
#' given sample - if any SNP is missing, the *entire* marker is dropped for that sample (a
#' missing allele only ever affects the single MH marker it belongs to).
#' @param canon_long Tibble as returned by `read_mh_profile_file()`.
#' @param snp_map The MH-to-SNP map, see `mh_snp_map()`.
#' @param probs_tbl The combined probability table, see `mh_geno_probs_table()`.
#' @param devi_tbl The combined deviance table, see `mh_geno_devi_table()`.
#' @param out_metas The per-meta outlier null models, see `mh_out_metas()`.
#' @param ci_data Prepared reference haplotype data for the log-likelihood confidence interval, see
#'   `mh_ci_data()`. If `NULL`, no standard errors are computed.
#' @param alpha Significance level for the outlier test (default `0.05`, mirrors the AIMs `CI` control).
#' @return `NULL` if nothing could be computed, otherwise a list with:
#'   `result`: tibble with `sample_name`, `meta`, `n_loci` (MH markers used), `likelihood`,
#'   `log10_likelihood`, `LR_best` (likelihood ratio vs. the most likely meta-population),
#'   `rel_prob` (likelihood normalised to sum to one across meta-populations, assuming a
#'   uniform prior - NOT a p-value), `devi` (summed deviance, calibrated against `out_metas`),
#'   `z_score` (`(devi - mu) / sqrt(2 * sigma)`), `p_value` (`1 - ecdf(sims)(devi)`, i.e. `pval()` from
#'   `sampling_profiles.R`, computed without depending on `molic` at query time), `critical_value`
#'   (deviance threshold at `alpha`) and `accept` (`p_value > alpha`).
#'   `se_log10` (first-order standard error of `log10_likelihood`, see `mh_log10_var()`; `NA` if no
#'   `ci_data` was supplied).
#'   `ci_unobserved`: always empty (the variance is defined for every profile), kept for the module.
#'   `unmatched_loci`: locus ids in the file not part of any known MH marker.
#'   `incomplete_mh`: tibble of (sample, MH) combinations dropped due to missing SNPs.
#'   `not_found`: tibble of (sample, MH) genotype combinations missing from the reference data.
#' @noRd
mh_profile_likelihood <- function(canon_long, snp_map = mh_snp_map(), probs_tbl = mh_geno_probs_table(),
                                   devi_tbl = mh_geno_devi_table(), out_metas = mh_out_metas(),
                                   ci_data = mh_ci_data(), alpha = 0.05){
  if(is.null(canon_long) || nrow(canon_long) == 0) return(NULL)
  meta_cols <- setdiff(names(probs_tbl), c("MH", "geno"))

  unmatched_loci <- canon_long %>% dplyr::anti_join(snp_map, by = c("locus" = "snp")) %>%
    dplyr::distinct(locus) %>% dplyr::pull(locus)

  joined <- canon_long %>% dplyr::inner_join(snp_map, by = c("locus" = "snp"))

  ## A microhaplotype marker is dropped ENTIRELY for a sample if any of its SNPs is missing -
  ## a missing allele only ever affects the single MH marker it belongs to.
  mh_snp_counts <- snp_map %>% dplyr::count(MH, name = "n_snps_expected")
  coverage <- joined %>% dplyr::distinct(sample_name, MH, locus) %>%
    dplyr::count(sample_name, MH, name = "n_snps_found") %>%
    dplyr::left_join(mh_snp_counts, by = "MH")
  incomplete_mh <- coverage %>% dplyr::filter(n_snps_found != n_snps_expected) %>%
    dplyr::select(sample_name, MH, n_snps_found, n_snps_expected)

  complete_keys <- coverage %>% dplyr::filter(n_snps_found == n_snps_expected) %>%
    dplyr::select(sample_name, MH)
  if(nrow(complete_keys) == 0){
    return(list(result = NULL, unmatched_loci = unmatched_loci, incomplete_mh = incomplete_mh, not_found = NULL))
  }

  per_mh <- joined %>%
    dplyr::semi_join(complete_keys, by = c("sample_name", "MH")) %>%
    dplyr::arrange(sample_name, MH, position) %>%
    dplyr::summarise(
      pair = list(mh_motif_pair(allele1, allele2)),
      .by = c(sample_name, MH)
    )

  probs_by_mh <- split(probs_tbl, probs_tbl$MH)
  devi_by_mh <- split(devi_tbl, devi_tbl$MH)
  mh_meta_stats <- purrr::pmap(per_mh, function(sample_name, MH, pair){
    na_out <- list(prob = stats::setNames(rep(NA_real_, length(meta_cols)), meta_cols),
                   devi = stats::setNames(rep(NA_real_, length(meta_cols)), meta_cols))
    ptbl <- probs_by_mh[[MH]]; dtbl <- devi_by_mh[[MH]]
    if(is.null(ptbl) || is.null(dtbl)) return(na_out)
    prob <- mh_lookup_candidate(pair, ptbl, meta_cols)
    devi <- mh_lookup_candidate(pair, dtbl, meta_cols)
    if(is.null(prob) || is.null(devi)) return(na_out)
    list(prob = prob, devi = devi, pair = pair)
  })
  not_found <- per_mh %>% dplyr::select(sample_name, MH) %>%
    dplyr::filter(vapply(mh_meta_stats, function(x) all(is.na(x$prob)), logical(1)))

  result_prob <- per_mh %>% dplyr::select(sample_name, MH) %>%
    dplyr::bind_cols(dplyr::bind_rows(purrr::map(mh_meta_stats, "prob")))
  result_devi <- per_mh %>% dplyr::select(sample_name, MH) %>%
    dplyr::bind_cols(dplyr::bind_rows(purrr::map(mh_meta_stats, "devi")))

  agg_prob <- result_prob %>%
    dplyr::filter(!is.na(.data[[meta_cols[1]]])) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(meta_cols), prod), n_loci = dplyr::n(), .by = sample_name)
  if(nrow(agg_prob) == 0){
    return(list(result = NULL, unmatched_loci = unmatched_loci, incomplete_mh = incomplete_mh, not_found = not_found))
  }
  agg_devi <- result_devi %>%
    dplyr::filter(!is.na(.data[[meta_cols[1]]])) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(meta_cols), sum), .by = sample_name)

  result <- agg_prob %>%
    tidyr::pivot_longer(cols = dplyr::all_of(meta_cols), names_to = "meta", values_to = "likelihood") %>%
    dplyr::left_join(
      agg_devi %>% tidyr::pivot_longer(cols = dplyr::all_of(meta_cols), names_to = "meta", values_to = "devi"),
      by = c("sample_name", "meta")
    ) %>%
    dplyr::mutate(
      log10_likelihood = log10(likelihood),
      LR_best = max(likelihood) / likelihood,
      rel_prob = likelihood / sum(likelihood),
      .by = sample_name
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      mu = out_metas[[meta]]$mu,
      sigma = out_metas[[meta]]$sigma,
      z_score = (devi - mu) / sqrt(2 * sigma),
      p_value = 1 - out_metas[[meta]]$cdf(devi),
      critical_value = stats::quantile(out_metas[[meta]]$sims, probs = 1 - alpha, names = FALSE),
      accept = p_value > alpha
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(sample_name, dplyr::desc(likelihood))
  ci <- if(is.null(ci_data)) NULL else mh_likelihood_ci(per_mh, mh_meta_stats, meta_cols, ci_data, snp_map)
  result <- if(is.null(ci)) dplyr::mutate(result, se_log10 = NA_real_) else dplyr::left_join(result, ci$se, by = c("sample_name", "meta"))
  list(result = result, unmatched_loci = unmatched_loci, incomplete_mh = incomplete_mh, not_found = not_found,
       ci_unobserved = ci$unobserved)
}

#' Pairwise likelihood ratios for an MH result table
#'
#' Mirrors `LR_table()` (see `api_LR_table.R`), but without the log-likelihood confidence
#' interval decomposition - MH likelihoods come from precomputed table look-ups, not a fitted
#' graphical model, so there is no `varlogP`-style variance available to derive a CI from.
#' @param result A single sample's slice of `mh_profile_likelihood()`'s `result` tibble.
#' @param who A vector of meta-population names to include. If `NULL`, all are used (ordered by likelihood).
#' @param only_accepted Restrict the ratios to pairs including at least one accepted population.
#' @return A tibble with `numerator`, `denominator`, `log10_LR`, `num_accept`, `den_accept`,
#'   `num_p_value`, `den_p_value` - or `NULL` if fewer than two meta-populations are available.
#' @noRd
mh_LR_table <- function(result, who = NULL, only_accepted = TRUE){
  if(is.null(result) || nrow(result) <= 1) return(NULL)
  if(is.null(who)) who <- result %>% dplyr::arrange(dplyr::desc(log10_likelihood)) %>% dplyr::pull(meta)
  lr_df <- result %>% dplyr::filter(meta %in% who) %>%
    dplyr::select(meta, log10_likelihood, p_value, accept)
  lr_list <- dplyr::cross_join(lr_df, lr_df, suffix = c("_num", "_den")) %>%
    dplyr::filter(meta_num != meta_den) %>%
    dplyr::mutate(log10_LR = log10_likelihood_num - log10_likelihood_den) %>%
    dplyr::rename(numerator = meta_num, denominator = meta_den, num_accept = accept_num, den_accept = accept_den,
                  num_p_value = p_value_num, den_p_value = p_value_den) %>%
    dplyr::select(numerator, denominator, log10_LR, num_accept, den_accept, num_p_value, den_p_value)
  rnk <- stats::setNames(seq_along(who), who)
  lr_list <- lr_list %>% dplyr::filter(rnk[numerator] < rnk[denominator])
  if(only_accepted) lr_list <- lr_list %>% dplyr::filter(num_accept | den_accept)
  lr_list %>% dplyr::arrange(numerator, denominator)
}

#' Plot of the log10 profile probability per meta-population with approximate confidence intervals
#'
#' Mirrors `error_bar_plot()`/`error_bar_plotly()` (see `helper_plots.R`) for the AIMs analysis: same
#' axis, ordering and accepted/rejected colour scale (`bar_colour()`). The intervals are
#' `log10_likelihood +/- z * se_log10` (first-order standard error, see `mh_log10_var()`); they are omitted
#' where `se_log10` is `NA`.
#' @param result_df A single sample's slice of `mh_profile_likelihood()`'s `result` tibble
#'   (must contain `meta`, `log10_likelihood`, `se_log10`, `n_loci`, `z_score`, `p_value`, `accept`).
#' @param ci Confidence level for the intervals (default `0.95`, mirrors the AIMs `CI` control).
#' @return A `ggplot` object.
#' @noRd
mh_error_bar_plot <- function(result_df, ci = 0.95){
  z <- stats::qnorm((1 - ci) / 2, lower.tail = FALSE)
  df <- result_df %>%
    dplyr::mutate(
      logP = log10_likelihood,
      logP_lwr = logP - z * se_log10,
      logP_upr = logP + z * se_log10,
      ci_text = ifelse(is.na(se_log10), "CI not available", sprintf("[%0.2f, %0.2f]", logP_lwr, logP_upr)),
      pop_label = sprintf(
        "<b>%s</b> (<i>%d MH markers</i>)<br><b><i>z</i>-score:</b> %0.2f (%s)<br><b>log<sub>10</sub> P</b>: %0.2f %s<br><b>p-value:</b> %0.4f",
        meta, n_loci, z_score, ifelse(accept, "accepted", "rejected"), logP, ci_text, p_value
      ),
      meta = forcats::fct_reorder(meta, logP)
    )
  df %>%
    ggplot(aes(y = meta, x = logP, xmin = logP_lwr, xmax = logP_upr, colour = meta, text = pop_label)) +
    labs(y = "", x = "log<sub>10</sub> P(Profile | Meta-population)") +
    guides(colour = "none") +
    coord_cartesian(xlim = rev(range(c(df$logP, df$logP_lwr, df$logP_upr), na.rm = TRUE))) +
    geom_point() + geom_errorbarh(na.rm = TRUE) +
    scale_colour_manual(values = bar_colour(df[, c("logP", "accept", "meta")])) +
    theme_bw()
}

#' Read a small preview (first rows/columns) of one of the bundled MH example files
#' @param filename File name within `inst/extdata` (e.g. `"mh_sample_long.csv"`).
#' @param n_row,n_col Maximum number of rows/columns to show.
#' @noRd
mh_example_preview <- function(filename, n_row = 6, n_col = 8){
  path <- system.file("extdata", filename, package = ggg_package)
  df <- rio::import(path, header = TRUE)
  df[seq_len(min(n_row, nrow(df))), seq_len(min(n_col, ncol(df))), drop = FALSE]
}
