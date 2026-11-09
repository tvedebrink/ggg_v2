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
##  - long, two columns for alleles:  sample_name;locus;A;B      (phased: A/B = haplotype copy 1/2)
##  - long, one column for alleles:   sample_name;locus;allele   (unphased genotype, e.g. "GT")
##  - wide, two rows per sample:      sample_name;<locus1>;...   (phased: one row per haplotype copy)
##  - wide, one row per sample:       sample_name;<locus1>;...   (unphased genotype per cell, e.g. "GT")

mh_cache <- new.env(parent = emptyenv())

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

#' Detect which of the four supported MH input formats a file uses
#'
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
  ## Wide format: distinguish by the width of a data cell (1 character = phased/one haplotype
  ## per row; 2 characters = unphased genotype per cell)
  cell <- as.character(raw[2, 2])
  if(nchar(trimws(cell)) == 1) return("wide_two_row")
  if(nchar(trimws(cell)) == 2) return("wide_one_row")
  stop("Unrecognised wide-format data: expected single-character alleles (two rows per sample) or two-character genotypes (one row per sample).")
}

#' Read an MH profile file into canonical long format, regardless of the specific input format
#'
#' @param path Path to the uploaded CSV/Excel file.
#' @return A tibble with columns `sample_name`, `locus`, `allele1`, `allele2`, `phased`
#'   (`TRUE` when `allele1`/`allele2` are known to come from a specific haplotype copy,
#'   `FALSE` when only the unordered pair of alleles at that locus is known).
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
    allele2 = toupper(trimws(raw[["B"]])),
    phased = TRUE
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
    allele2 = substr(geno, 2, 2),
    phased = FALSE
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
    tidyr::pivot_wider(names_from = copy, values_from = allele, names_prefix = "allele") %>%
    dplyr::mutate(phased = TRUE)
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
      allele1 = substr(geno, 1, 1), allele2 = substr(geno, 2, 2), phased = FALSE
    )
  })
}

#' Enumerate the distinct haplotype-motif pairs consistent with a (possibly unphased) MH genotype
#'
#' When `phased` is `TRUE`, or there is at most one heterozygous locus, the motif pair is
#' unambiguous. With two or more heterozygous loci, linkage phase is unknown, so every
#' distinct (unordered) phase resolution is returned - the caller sums their probabilities to
#' marginalise over the unknown phase.
#' @param allele1,allele2 Character vectors of alleles at each locus, ordered by `position`.
#' @param phased Logical; whether `allele1`/`allele2` are known to come from specific haplotype copies.
#' @return A list of two-element character vectors `c(motif1, motif2)`.
#' @noRd
mh_motif_candidates <- function(allele1, allele2, phased){
  het <- allele1 != allele2
  n_het <- sum(het)
  if(phased || n_het <= 1){
    return(list(c(paste(allele1, collapse = ""), paste(allele2, collapse = ""))))
  }
  het_idx <- which(het)
  n <- length(het_idx)
  swaps <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
  raw <- lapply(seq_len(nrow(swaps)), function(i){
    a1 <- allele1; a2 <- allele2
    flip <- het_idx[swaps[i, ]]
    if(length(flip) > 0){
      tmp <- a1[flip]; a1[flip] <- a2[flip]; a2[flip] <- tmp
    }
    c(paste(a1, collapse = ""), paste(a2, collapse = ""))
  })
  canon_keys <- vapply(raw, function(v) paste(sort(v), collapse = ";"), character(1))
  raw[!duplicated(canon_keys)]
}

#' Look up a candidate motif pair in a combined `MH`/`geno` table, trying both allele orders
#' @noRd
mh_lookup_candidate <- function(cand, tbl, meta_cols){
  row <- tbl[tbl$geno == paste(cand[1], cand[2], sep = ";"), meta_cols, drop = FALSE]
  if(nrow(row) == 0) row <- tbl[tbl$geno == paste(cand[2], cand[1], sep = ";"), meta_cols, drop = FALSE]
  if(nrow(row) == 0) return(NULL)
  unlist(row[1, ])
}

#' Compute per-sample, per-meta-population likelihoods and outlier test of an MH profile
#'
#' For each sample and MH marker with complete SNP coverage in the input, reconstructs the
#' haplotype motif pair (marginalising over unknown phase when needed, see
#' `mh_motif_candidates()`) and looks up its precomputed genotype probability/deviance. No
#' model fitting or simulation is performed - this is combinatorics and table look-ups only.
#' A microhaplotype marker is used only when ALL of its constituent SNPs are present for a
#' given sample - if any SNP is missing, the *entire* marker is dropped for that sample (a
#' missing allele only ever affects the single MH marker it belongs to).
#' When phase is unknown and a marker has 2+ heterozygous SNPs, the marker's deviance is the
#' probability-weighted average of its candidate phase resolutions' deviances (reducing to the
#' single value when there is no ambiguity).
#' @param canon_long Tibble as returned by `read_mh_profile_file()`.
#' @param snp_map The MH-to-SNP map, see `mh_snp_map()`.
#' @param probs_tbl The combined probability table, see `mh_geno_probs_table()`.
#' @param devi_tbl The combined deviance table, see `mh_geno_devi_table()`.
#' @param out_metas The per-meta outlier null models, see `mh_out_metas()`.
#' @param alpha Significance level for the outlier test (default `0.05`, mirrors the AIMs `CI` control).
#' @return `NULL` if nothing could be computed, otherwise a list with:
#'   `result`: tibble with `sample_name`, `meta`, `n_loci` (MH markers used), `likelihood`,
#'   `log10_likelihood`, `LR_best` (likelihood ratio vs. the most likely meta-population),
#'   `rel_prob` (likelihood normalised to sum to one across meta-populations, assuming a
#'   uniform prior - NOT a p-value), `devi` (summed deviance, calibrated against `out_metas`),
#'   `z_score` (`(devi - mu) / sigma`), `p_value` (`1 - ecdf(sims)(devi)`, i.e. `pval()` from
#'   `sampling_profiles.R`, computed without depending on `molic` at query time), `critical_value`
#'   (deviance threshold at `alpha`) and `accept` (`p_value > alpha`).
#'   `unmatched_loci`: locus ids in the file not part of any known MH marker.
#'   `incomplete_mh`: tibble of (sample, MH) combinations dropped due to missing SNPs.
#'   `not_found`: tibble of (sample, MH) genotype combinations missing from the reference data.
#' @noRd
mh_profile_likelihood <- function(canon_long, snp_map = mh_snp_map(), probs_tbl = mh_geno_probs_table(),
                                   devi_tbl = mh_geno_devi_table(), out_metas = mh_out_metas(), alpha = 0.05){
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
      candidates = list(mh_motif_candidates(allele1, allele2, phased[1])),
      .by = c(sample_name, MH)
    )

  probs_by_mh <- split(probs_tbl, probs_tbl$MH)
  devi_by_mh <- split(devi_tbl, devi_tbl$MH)
  mh_meta_stats <- purrr::pmap(per_mh, function(sample_name, MH, candidates){
    na_out <- list(prob = stats::setNames(rep(NA_real_, length(meta_cols)), meta_cols),
                   devi = stats::setNames(rep(NA_real_, length(meta_cols)), meta_cols))
    ptbl <- probs_by_mh[[MH]]; dtbl <- devi_by_mh[[MH]]
    if(is.null(ptbl) || is.null(dtbl)) return(na_out)
    hits_p <- purrr::compact(purrr::map(candidates, mh_lookup_candidate, tbl = ptbl, meta_cols = meta_cols))
    hits_d <- purrr::compact(purrr::map(candidates, mh_lookup_candidate, tbl = dtbl, meta_cols = meta_cols))
    if(length(hits_p) == 0 || length(hits_d) == 0) return(na_out)
    prob <- Reduce(`+`, hits_p)
    ## Probability-weighted average deviance across ambiguous phase candidates
    devi <- Reduce(`+`, purrr::map2(hits_p, hits_d, `*`)) / prob
    list(prob = prob, devi = devi)
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
      z_score = (devi - mu) / sigma,
      p_value = 1 - out_metas[[meta]]$cdf(devi),
      critical_value = stats::quantile(out_metas[[meta]]$sims, probs = 1 - alpha, names = FALSE),
      accept = p_value > alpha
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(sample_name, dplyr::desc(likelihood))
  list(result = result, unmatched_loci = unmatched_loci, incomplete_mh = incomplete_mh, not_found = not_found)
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

#' Error-bar plot of the outlier test statistic (deviance) per meta-population
#'
#' Mirrors `error_bar_plot()`/`error_bar_plotly()` (see `helper_plots.R`) for the AIMs analysis,
#' but plots the MH deviance statistic instead of `logP` - the error bars use `sigma`, the
#' standard error of the deviance estimated (offline) from the precomputed outlier null models
#' (see `mh_out_metas()`), exactly as `mu`/`sigma` are used for `z_score` in `mh_profile_likelihood()`.
#' @param result_df A single sample's slice of `mh_profile_likelihood()`'s `result` tibble
#'   (must contain `meta`, `devi`, `sigma`, `z_score`, `p_value`, `accept`).
#' @param ci Confidence level for the error bars (default `0.95`, mirrors the AIMs `CI` control).
#' @return A `ggplot` object.
#' @noRd
mh_error_bar_plot <- function(result_df, ci = 0.95){
  z <- stats::qnorm((1 - ci) / 2, lower.tail = FALSE)
  df <- result_df %>%
    dplyr::mutate(
      devi_lwr = devi - z * sigma,
      devi_upr = devi + z * sigma,
      pop_label = sprintf(
        "<b>%s</b><br><b><i>z</i>-score:</b> %0.2f (%s)<br><b>Deviance:</b> %0.1f [%0.1f, %0.1f]<br><b>p-value:</b> %0.4f",
        meta, z_score, ifelse(accept, "accepted", "rejected"), devi, devi_lwr, devi_upr, p_value
      ),
      meta = forcats::fct_reorder(meta, devi)
    )
  devi_range <- c(min(df$devi_lwr), max(df$devi_upr))
  df %>%
    ggplot(aes(y = meta, x = devi, xmin = devi_lwr, xmax = devi_upr, colour = accept, text = pop_label)) +
    labs(y = "", x = "Deviance") +
    guides(colour = "none") +
    coord_cartesian(xlim = devi_range) +
    geom_point() + geom_errorbarh() +
    scale_colour_manual(values = c(`TRUE` = "darkgreen", `FALSE` = "darkred")) +
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
