score_add_df <- function(db){
  ## build fixes : start ##
  x0 <- NULL; meta <- NULL; x1 <- NULL
  lat <- NULL; lon <- NULL; n_ <- NULL
  locus <- NULL; freq <- NULL
  ## build fixes : end ##
  grouping <- match.arg(grouping, c("pop", "meta"))
  if(!(grouping %in% names(db))) return(NULL)
  if(grouping == "meta"){
    db <- db %>% filter(x0 == X0, x0==0) %>%
      select_(.dots = vars(meta:x1, grouping)) %>%
      mutate(n_ = ifelse(out_of_place, 0, n)) %>%
      group_by_(.dots = c(grouping, "locus", "main_allele", "other_allele")) %>%
      summarise(
        out_of_place = out_of_place[1],
        lat = weighted.mean(lat, w = n_, na.rm = TRUE),
        lon = weighted.mean(lon, w = n_, na.rm = TRUE),
        n = sum(n),
        x1 = sum(x1)
      ) %>% ungroup() %>%
      select_(.dots = vars(grouping, n, lat:lon, out_of_place, locus:x1))
  }
  ## X0 and X1 are relevant for simulations
  dd <- crossing(db, x0 = 0:2, X0 = 0:2) %>%
    mutate(X1 = x1 + x0 - X0)
  ## NB! if x0 == X0 we have that x1 == X1
  dd <- dd %>% bind_cols(ss_moments_(x0 = .$X0, x1 = .$X1, N = .$n))
  ## Fix correct simulation probabilities
  dd <- dd %>% mutate(
    score = xlx(X1) + xlx(n-X1) - 2*log(2)*(X0==1), ## The unstandardised locus score
    freq = (x1 + x0)/(n + 2), ## freqs relevant for log P
    logP = log_P_(x0 = x0, p = freq)/log(10), ## log_10 scale
    varlogP = varlog_P_(x0 = x0, n = n, p = freq)/(log(10)^2)) ## log_10 scale
  ## Meta population
  dd
}

#' Pre-compute the scores for a given reference database
#'
#' Convert the counts from each population over a range of AIMs SNPs q
#' to observed likelihood ratio test, its mean and variance.
#' Based on these pre-computed the evaluation of a specific profile is done
#' using \code{genogeo} with the resulting dataframe as \code{df}.
#' @param db A dataframe with columns similar to those of \code{simulate_pops()}.
#' If \code{db} contains information (recommended!) about "meta" (meta population)
#' and "lat"/"lon" (location) these are carried over into the calculations
#' @param ... Additional arguments passed to \code{score_add_df}
#' @return A tibble with population and locus specific score information
#' @export
#' @examples
#' df_ <- simulate_pops(pop_n = 4, aims_n = 50)
#' df_db <- pops_to_DB(df_)

pops_to_DB <- function(db, ...){
  ## build fixes : start ##
  meta <- NULL
  lon <- NULL
  lat <- NULL
  n_row <- NULL
  pop <- NULL
  population <- NULL
  metapopulation <- NULL
  cluster <- NULL
  ## build fixes : end ##
  dd_pop <- score_add_df(db, ...)
  if("meta" %in% names(dd_pop)){
    meta_hulls <- convex_hulls(dd_pop) %>%
      select(meta, lat, lon , n_row) %>% group_by(meta) %>% nest(.key = "hulls_meta")
    dd_meta <- score_add_df(db = dd_pop, grouping = "meta") %>%
      inner_join(meta_hulls, by = "meta") %>%
      group_by(meta) %>%
      nest(.key = "meta_data")
    dd <- inner_join(
      dd_pop %>% group_by(pop, population, meta, metapopulation) %>% nest(.key = "population_data"),
      dd_meta,
      by = "meta")
  }
  else dd <- dd_pop %>% group_by(pop, population) %>% nest(.key = "population_data")
  if("cluster" %in% names(db)){
    cluster_hulls <- convex_hulls(dd_pop, grouping = "cluster") %>%
      select(cluster, lat, lon , n_row) %>% group_by(cluster) %>% nest(.key = "hulls_cluster")
    dd_cluster <- score_add_df(db = dd_pop, grouping = "cluster") %>%
      inner_join(cluster_hulls, by = "cluster") %>%
      group_by(cluster) %>%
      nest(.key = "cluster_data")
    dd <- inner_join(
      dd_pop %>% group_by(pop, population, meta, metapopulation, cluster) %>% nest(.key = "population_data"),
      dd_meta, by = "meta") %>%
      inner_join(dd_cluster, by = "cluster")
  }
  dd
}

DB_locus_set <- function(db, locusset){
  db <- db %>%
    mutate(population_data = purrr::map(data, ~ semi_join(.x, locusset, by = "locus")),
           meta_data = purrr::map(data, ~ semi_join(.x, locusset, by = "locus")))
}

