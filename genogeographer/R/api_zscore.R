non_list_cols <- function(df){
  names(df)[df %>% map_lgl(~ !any(class(.x) %in% "list"))]
}

#' @title Calculates z-scores based on genogeographer methodology
#' @param profile_x0 The AIMs profile to analyse
#' @param DB The AIMs database to use for the analysis
#' @param CI The level for used for confidence intervals
#' @param tilt Should exponential tilting be applied
#' @param LOO Leave one out analysis - provide the name of the (meta)population
#' @param min_n Minimum number of required samples
#' @param tilt_ctrl Controls for exponential tilting
#' @export
ggg_score <- function(profile_x0, DB, min_n = 75, CI = 0.95, tilt = FALSE, LOO = NULL, tilt_ctrl = list(n = 500, p_range = c(0.001, 0.1))){
  groups <- non_list_cols(DB)
  groups_ <- sym(groups)
  if("info" %in% names(attributes(DB))){
    DB_info <- attr(DB, "info") %>% mutate(n = map_dbl(n, ~eval(parse(text = (sub("\\&", "\\+", .x))))))
    DB_info <- DB_info %>% filter(n >= min_n)
    DB <- DB %>% semi_join(DB_info, by = groups)
  }
  if(nrow(DB) == 0) return(NULL)
  ## Fix for Leave-One-Out analysis!
  if(!is.null(LOO) && any(LOO == DB[[groups]])){
    DB_loo <- DB %>% filter(LOO == !!groups_) %>% loo_update(profile_x0)
    DB_ <- DB %>% filter(!(LOO == !!groups_))
    DB <- bind_rows(DB_loo, DB_) %>% arrange(!!groups_)
  }
  ##
  z_CI <- qnorm(1 - (1 - CI)/2)
  result <- DB %>%
    mutate(data = purrr::map(data, ~ .x %>% select(locus, x0, logP, varlogP, z_raw, z_var))) %>%
    unnest(cols = data) %>%
    inner_join(profile_x0, by = c("locus", "x0")) %>%
    group_by(across({{groups}})) %>%
    summarise(
      logP = sum(logP),
      varlogP = sum(varlogP),
      logP_upr = logP + z_CI*sqrt(varlogP),
      logP_lwr = logP - z_CI*sqrt(varlogP),
      z_score = sum(z_raw)/sqrt(sum(z_var)),
      p_value = pnorm(z_score, lower.tail = FALSE),
      .groups = "drop"
    ) %>%
    arrange(desc(logP)) %>%
    select({{groups}}, everything())
  if(tilt){
    result_tilt <- result %>% filter(between(p_value, tilt_ctrl$p_range[1], tilt_ctrl$p_range[2]))
    DB_tilt <- DB %>% semi_join(result_tilt, by = groups)
    ## Check for admixture - tilting not implemented
    admixed <- DB_tilt %>% mutate(admix = map_lgl(data, ~ class(.x$x1) == "character")) %>% filter(admix) %>% pull(!!groups_)
    if(length(admixed) > 0) warning("Exponential tilting is not implemented for admixed populations")
    result_tilt <- result_tilt %>% filter(!(!!groups_ %in% admixed))
    if(nrow(result_tilt) > 0){
      DB_tilt <- DB_tilt %>% filter(!(!!groups_ %in% admixed)) %>%
        mutate(data = purrr::map(data, ~.x %>% semi_join(profile_x0, by = c("locus", "x0")))) %>%
        mutate(p_value = map_dbl(data, ~ with(.x, exponent_tilt(x0 = x0, x1 = x1, n = n, B = tilt_ctrl$n, p_limit = tilt_ctrl$p_range[2]))))
      result_tilt <- result_tilt %>% select(-p_value) %>%
        left_join(DB_tilt %>% select(-data), by = groups)
      result <- result %>% anti_join(result_tilt, by = groups) %>%
        bind_rows(result_tilt) %>% arrange(desc(logP))
    }
  }
  result <- result %>% mutate(accept = p_value > 1-CI)
  if("info" %in% names(attributes(DB))) attr(result, "info") <- DB_info
  result
}

