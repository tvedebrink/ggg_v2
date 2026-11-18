
## Remember tri and tetra allelic loci! x1, x2, x3, x4 (at some point)

genotype_x0 <- function(profile, locus = "locus", genotype = "genotype", ggg = NULL){
  profile <- na.omit(profile)
  if(genotype %in% names(profile)) genotype_ <- sym(genotype)
  else genotype_ <- sym(genotype_guess(profile)[1])
  if(locus %in% names(profile)) locus_ <- sym(locus)
  else locus_ <- sym(locus_guess(profile)[1])
  ggg %>% inner_join(profile %>% select(locus = !!locus_, genotype = !!genotype_), by = "locus") %>%
    tidyr::extract(genotype, into = c("g1", "g2"), regex = "(.{1})(.{1})") %>%
    mutate(
      x0 = (g1 == x1) + (g2 == x1),
      n2_x0 = (g1 == x2) + (g2 == x2)
      ) %>%
    filter(x0 + n2_x0 == 2) %>%
    select(locus, x0)
}

count_alleles <- function(df, groups = NULL, snps = NULL, exclude = NULL, missing = c("NN", "99", "--", "-")){
  if(is.null(groups)) stop("Argument 'groups' must be specified")
  if(!is.null(exclude)) df <- df %>% select(-any_of(exclude))
  df_n <- names(df)
  if(is.null(snps)) snps <- setdiff(df_n, groups)
  else{
    exclude <- setdiff(df_n, c(groups, snps))
    df <- df %>% select(-any_of(exclude))
  }
  ##
  df <- df %>% group_by(across({{groups}}))
  group_n <- df %>% count() %>% ungroup()
  df <- df %>%
    pivot_longer(cols = {{snps}}, names_to = "locus", values_to = "genotype") %>%
    ungroup() %>%
    filter(!(genotype %in% missing)) %>% filter(!is.na(genotype)) %>%
    extract(col = genotype, into = c("g1", "g2"), regex = "(.{1})(.{1})") %>%
    pivot_longer(cols = c("g1", "g2"), names_to = "pair", values_to = "allele") %>%
    select(-pair)
  df_n <- names(df)
  df_c <- df %>% group_by(across({{df_n}})) %>% count(name = "count") %>%
    ungroup()
  list(samples = group_n, alleles = df_c)
  }

lat_lon <- function(df_latlon, df_n, groups, lat = "lat", lon = "lon", out_of_place = NULL){
  if(is.null(groups)) stop("Argument 'groups' must be specified")
  groups_ <- sym(groups)
  latlon_exists <- all(c(lat, lon) %in% names(df_latlon))
  if(!latlon_exists) stop(glue("Columns '{{lat}}' and/or '{{lon}}' not in df_latlon"))
  df_latlon <- df_latlon %>% rename(lat = !!sym(lat), lon = !!sym(lon))
  df <- df_latlon %>% inner_join(df_n, by = "pop") %>% mutate(n = paste(n))
  if(groups == "meta"){
    df <- df %>% filter(!(pop %in% out_of_place)) %>%
      select(-starts_with("pop")) %>% nest(latlon = -starts_with("meta")) %>%
      mutate(
        convex_hull = purrr::map(latlon, ~ .x[chull(x = .x[[lon]], y = .x[[lat]]),] %>% select(-n)),
        n = map_int(latlon, ~.x %>% pull(n) %>% as.integer() %>% sum()) %>% paste(),
        latlon = purrr::map(latlon, ~.x %>% summarise(across(-n, ~weighted.mean(.x, w = as.integer(n)))))
        ) %>%
      unnest(latlon)
  }
  else df <- df %>% select(-starts_with("meta"))
  df
}

incProgress_ <- function(amount = 0.1, message = NULL, detail = NULL, session = getDefaultReactiveDomain()){
  p <- session$progressStack$peek()
  p$inc(amount)
  invisible()
}

incProgress_trick <- function(trick = FALSE, message = NULL, detail = NULL, session = getDefaultReactiveDomain()){
  p <- session$progressStack$peek()
  if(trick) p$inc(1)
  invisible()
}


make_x1 <- function(df, latlon = NULL, groups, allele_list = NULL,
                    out_of_place = NULL, lat = "lat", lon = "lon", shiny = NULL, ...){
  dfc <- count_alleles(df, groups = groups, ...) ## returns: list(samples={}, alleles = {{groups}}, locus, allele, n)
  df_s <- dfc$samples
  if(groups == "meta") df_s <- count_alleles(df, groups = "pop", exclude = "meta")$samples
  if(!is.null(allele_list)) non_ggg <- dfc$alleles %>% anti_join(allele_list, by = c("locus"))
  else non_ggg <- dfc$alleles
  allele_list <- bind_rows(
    allele_list,
    non_ggg %>%
      group_by(locus, allele) %>% summarise(count = sum(count)) %>%
      slice_max(order_by = count, n = 2) %>% select(-count) %>%
      arrange(locus, allele) %>% mutate(ggg = paste0("x", row_number())) %>%
      spread(ggg, allele)
  ) %>% ungroup()
  df_c <- dfc$alleles %>% inner_join(allele_list, by = c("locus")) %>%
    mutate(across(starts_with("x"), ~count*(allele == .x))) %>%
    group_by(across(c({{groups}}, "locus"))) %>%
    summarise(
      across(starts_with("x"), ~sum(.x, na.rm = TRUE)),
      n = sum(c_across(starts_with("x"))),
    .groups = "drop")
  if(!is.null(shiny)){
    df_c <- df_c %>% split(.[[groups]])
    # print(length(df_c))
    df_c <- purrr::map(seq_along(df_c), ~ {
      incProgress_trick(trick = floor((shiny$start + .x)*100/shiny$n) > floor((shiny$start + .x - 1)*100/shiny$n), session = shiny$session)
      df_c[[.x]] %>% select(-x2) %>%  ## ASSUMES ONLY BI-ALLELIC LOCI
        crossing(x0 = 0:2) %>%
        bind_cols(ss_moments(x0 = .$x0, x1 = .$x1, N = .$n)) %>%
        mutate(
          score = xlx(x1) + xlx(n-x1) - 2*log(2)*(x0==1), ## The unstandardised locus score
          freq = (x1 + x0)/(n + 2), ## freqs relevant for log P
          logP = log_P(x0 = x0, p = freq)/log(10), ## log_10 scale
          varlogP = varlog_P(x0 = x0, n = n, p = freq)/(log(10)^2), ## log_10 scale
          z_raw = score - z_exp
        ) %>% select(-score, -z_exp)
    }) %>% bind_rows()
  }
  else df_c <- df_c %>%
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
  x1 <- df_c %>% nest(data = c(locus, x1, n, x0, z_raw, z_var, freq, logP, varlogP))
  if(!is.null(latlon)){
    latlon <- lat_lon(df_latlon = latlon, df_n = df_s, lat = lat, lon = lon, groups = groups, out_of_place = out_of_place)
    attr(x1, "info") <- latlon
  }
  attr(x1, "allele_list") <- allele_list
  x1
}

locus_set <- function(db, locusset){
  db <- db %>% mutate(data = purrr::map(data, ~ semi_join(.x, locusset, by = "locus")))
  attr(db, "allele_list") <- semi_join(attr(db, "allele_list"), locusset, by = "locus")
  db
}

db_set <- function(db, locusset){
  purrr::map(db, ~ .x %>% purrr::map(~ .x %>% purrr::map(~locus_set(.x, locusset = locusset))))[[1]]
}



admix_dbs <- function(dbs, tol = 1e-8, pops = NULL, no_cores = NULL, shiny = NULL){
  if("info" %in% names(attributes(dbs))) info <- attr(dbs, "info")
  if("allele_list" %in% names(attributes(dbs))) allele_list <- attr(dbs, "allele_list")
  ### FIX FIX for names to be used in tables and plots (set lat/lon to NA)
  if(is.null(no_cores)) no_cores <- ceiling(parallel::detectCores()/2L)
  else{
    if(no_cores > parallel::detectCores()) no_cores <- parallel::detectCores()
    if(no_cores < 1) no_cores <- 1
  }
  groups <- names(dbs)[1]
  groups_ <- sym(groups)
  groups1 <- paste0(groups, "1")
  groups1_ <- sym(groups1)
  groups2 <- paste0(groups, "2")
  groups2_ <- sym(groups2)
  #
  dbs <- dbs %>%
    mutate(data = purrr::map(data, ~.x %>% select(locus, x = x1, n, x0))) %>%
    unnest(cols = data)
  if(is.null(pops)){
    dbs1 <- dbs
    dbs2 <- dbs
  }
  else if(is.list(pops)){
    if(length(pops) == 1){
      dbs1 <- dbs %>% filter(!!groups_ %in% pops[[1]])
      dbs2 <- dbs
    }
    else {
      dbs1 <- dbs %>% filter(!!groups_ %in% pops[[1]])
      dbs2 <- dbs %>% filter(!!groups_ %in% pops[[2]])
    }
  }
  else{
    dbs1 <- dbs %>% filter(!!groups_ %in% pops)
    dbs2 <- dbs
  }
  dbs <- dbs1 %>% inner_join(dbs2, by = c("locus", "x0"), suffix = c("1", "2")) %>%
    filter(!!groups1_ != !!groups2_)
  dbs_pairs <- dbs %>% distinct(!!groups1_, !!groups2_) %>%
    mutate(
      P1 = ifelse(!!groups1_<!!groups2_, !!groups1_, !!groups2_),
      P2 = ifelse(!!groups1_<!!groups2_, !!groups2_, !!groups1_)
    ) %>% distinct(P1, P2, .keep_all = TRUE)
  grouping <- if(groups == "pop") "population" else "metapopulation"
  grouping_ <- sym(grouping)
  info_pairs <- dbs_pairs %>%
    inner_join(info %>% rename(!!groups1_ := groups_), by = groups1) %>%
    inner_join(info %>% rename(!!groups2_ := groups_), by = groups2, suffix = c("1", "2")) %>%
    unite(col = {{groups}}, !!groups1_, !!groups2_, sep = " & ") %>%
    unite(col = {{grouping}}, !!sym(paste0(grouping,"1")), !!sym(paste0(grouping,"2")), sep = " & ") %>%
    unite(col = n, n1, n2, sep = " & ") %>%
    select(!!groups_, !!grouping_, n)
  dbs <- dbs %>% semi_join(dbs_pairs, by = c(groups1, groups2)) %>%
    select(!!groups1_, !!groups2_, locus, x0, x1, n1, x2, n2)
  dbs <- dbs %>% nest(data = -c(!!groups1_, !!groups2_)) %>%
    unite(col = {{groups}}, !!groups1_, !!groups2_, sep = " & ")
  if(!is.null(shiny)){
    # print(length(dbs$data))
    dbs$data <- purrr::map(seq_along(dbs$data), ~{
      incProgress_trick(trick = floor((shiny$start + .x)*100/shiny$n) > floor((shiny$start + .x - 1)*100/shiny$n), session = shiny$session)
      x <- dbs$data[[.x]]
      bind_cols(x, with(x, markers_z(x0, x1, n1, x2, n2, tol = tol)) %>% bind_rows()) %>%
        # likelihood and variance
        mutate(
          f1 = (x1 + ifelse(x0==1, p1, x0/2))/(n1 + 1),
          f2 = (x2 + ifelse(x0==1, p2, x0/2))/(n2 + 1),
          logP = log_Padmix(x0, f1, f2)/log(10),
          varlogP = varlog_Padmix(x0, f1, f2, n1, n2)/log(10)^2
        )})
  }
  else dbs <- dbs %>%
    mutate(data = mclapply(data, function(x)
      bind_cols(x, with(x, markers_z(x0, x1, n1, x2, n2, tol = tol)) %>% bind_rows()) %>%
        # likelihood and variance
        mutate(
          f1 = (x1 + ifelse(x0==1, p1, x0/2))/(n1 + 1),
          f2 = (x2 + ifelse(x0==1, p2, x0/2))/(n2 + 1),
          logP = log_Padmix(x0, f1, f2)/log(10),
          varlogP = varlog_Padmix(x0, f1, f2, n1, n2)/log(10)^2
          ),
      mc.cores = no_cores)
      )
  #
  dbs <- dbs %>%
    mutate(
      data = purrr::map(data, ~.x %>%
                   unite(x1, x1, x2, sep = " & ") %>% unite(n, n1, n2, sep = " & ") %>%
                   mutate(across(c(f1,f2), ~ round(.x, 2))) %>% unite(freq, f1, f2, sep = " & ") %>%
                   select(-p1,-p2)
                 )
      )
  attr(dbs, "info") <- info_pairs
  attr(dbs, "allele_list") <- allele_list
  dbs
}

loo_update <- function(df_loo, x0){
  df_loo %>% unnest(data) %>% inner_join(x0, by = c("locus", "x0")) %>%
    mutate(x1 = x1 - x0, n = n-2L) %>%
    mutate(x1 = ifelse(x1 < 0, 0L, x1)) %>% ## Account for non observed alleles in ref
    select(-c(x0, z_raw, z_var, freq, logP, varlogP)) %>%
    crossing(x0 = 0:2) %>%
    bind_cols(ss_moments(x0 = .$x0, x1 = .$x1, N = .$n)) %>%
    mutate(
      score = xlx(x1) + xlx(n-x1) - 2*log(2)*(x0==1), ## The unstandardised locus score
      freq = (x1 + x0)/(n + 2), ## freqs relevant for log P
      logP = log_P(x0 = x0, p = freq)/log(10), ## log_10 scale
      varlogP = varlog_P(x0 = x0, n = n, p = freq)/(log(10)^2),  ## log_10 scale
      z_raw = score - z_exp
    ) %>% select(-score, -z_exp) %>%
    nest(data = c(locus, x1, n, x0, z_raw, z_var, freq, logP, varlogP))
}



