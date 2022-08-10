genotype_x <- function(profile, locus = "locus", genotype = "genotype", tol = 0.9){
  DNA_bases <- c("A", "C", "G", "T")
  profile <- na.omit(profile)
  if(is.null(genotype) || !(genotype[1] %in% names(profile))) genotype <- genotype_guess(profile, tol = tol)[1]
  if(!is.null(genotype)) genotype_ <- sym(genotype)
  else stop("Profile contains no genotype information")
  if(is.null(locus) || !(locus[1] %in% names(profile))) locus <- locus_guess(profile, tol = tol)[1]
  if(!is.null(locus)) locus_ <- sym(locus)
  else stop("Profile contains no locus information")
  profile %>% select(locus = !!locus_, genotype = !!genotype_) %>%
    tidyr::extract(genotype, into = c("g1", "g2"), regex = "(.{1})(.{1})") %>%
    pivot_longer(cols = c("g1", "g2"), values_to = "genotype", names_to = "_drop") %>%
    count(locus, genotype) %>%
    bind_rows(tibble(locus = "FIX", genotype = DNA_bases, n = 0L)) %>%
    mutate(genotype = factor(genotype, levels = DNA_bases)) %>%
    pivot_wider(names_from = "genotype", values_from = "n", values_fill = 0L) %>%
    filter(locus != "FIX") %>% select(locus, any_of(DNA_bases)) %>%
    filter(rowSums(select(., any_of(DNA_bases))) == 2L)
}

make_x <- function(df, groups, allele_list = ggg_allele_list, ...){
  DNA_bases <- c("A", "C", "G", "T")
  x0s <- replicate(length(DNA_bases), 0:2, simplify = FALSE) %>% set_names(tolower(DNA_bases)) %>%
    expand.grid() %>% filter(rowSums(.) == 2)
  dfc <- count_alleles(df, groups = groups, ...) ## returns: list(samples={}, alleles = {{groups}}, locus, allele, n)
  df_c <- dfc$alleles %>% bind_rows(tibble({{groups}} := "FIX", locus = "FIX", allele = DNA_bases, count = 0L)) %>%
    pivot_wider(names_from = "allele", values_from = "count", values_fill = 0L) %>%
    filter(locus != "FIX") %>%
    group_by(across(c({{groups}}, "locus"))) %>%
    mutate(n = sum(c_across(any_of(DNA_bases)))) %>%
    ungroup()
  df_c <- df_c %>% mutate(x0 = list(x0s)) %>%
    browser()
  select(-x2) %>%  ## ASSUMES ONLY BI-ALLELIC LOCI
    crossing(x0 = 0:2) %>%
    bind_cols(ss_moments(x0 = .$x0, x1 = .$x1, N = .$n)) %>%
    mutate(
      score = xlx(x1) + xlx(n-x1) - 2*log(2)*(x0==1), ## The unstandardised locus score
      freq = (x1 + x0)/(n + 2), ## freqs relevant for log P
      logP = log_P(x0 = x0, p = freq)/log(10), ## log_10 scale
      varlogP = varlog_P(x0 = x0, n = n, p = freq)/(log(10)^2), ## log_10 scale
      z_raw = score - z_exp
    ) %>% select(-score, -z_exp)
  ## Meta population
  df_c %>% nest(data = c(locus, x1, n, x0, z_raw, z_var, freq, logP, varlogP))
}

