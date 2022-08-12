#' Compute pairwise likelihood ratios
#'
#' For each pair of a specified vector of profiles the likelihood ratios are computed.
#' The list can include all populations in the data or only a subset.
#' We may for inference purposes restrict to ratios including at least one "accepted" population.

#' @param z_df The output from \code{genogeo}
#' @param who A vector of population names (index elements (e.g. `pop`) in \code{z_df}). If NULL all populations are used.
#' @param only_accepted Restrict the ratios to include minimum one accepted population.
#' @param CI The level of confidence interval to be computed
#' @author Torben Tvedebrink \email{tvede@math.aau.dk}
#' @return A tibble with numerator and denominator populations with their log10 LR and uncertainty.
#' @export

LR_table <- function(z_df, who = NULL, only_accepted = TRUE, CI = 0.95){
  z <- qnorm((1-CI)/2, lower.tail = FALSE)
  groups <- names(z_df)[1] ## First cols from z-result is grouping
  groups_ <- sym(groups)
  num_groups_ <- sym(paste0("num_", groups_))
  den_groups_ <- sym(paste0("den_", groups_))
  # If empty, make all LRs
  if(is.null(who)) who <- z_df %>% arrange(desc(logP)) %>% pull({{groups}})
  lr_df <- z_df %>%
    filter(!!groups_ %in% who) %>%
    select(!!groups_, logP, varlogP, z_score, p_value) %>%
    mutate(!!groups_ := factor(!!groups_, levels = who, ordered = TRUE))
  ##
  num <- names(lr_df)
  lr_list <- lr_df %>% rowwise() %>%
    mutate(den = list(.)) %>%
    ungroup() %>%
    set_names(c(paste0("num_",num), "den")) %>%
    unnest(cols = den, names_sep = "_") %>%
    filter(!!num_groups_ < !!den_groups_)
  lr_list <- lr_list %>%
    mutate(logLR = num_logP - den_logP,
           var_logLR = num_varlogP + den_varlogP,
           CI_lwr = logLR - z*sqrt(var_logLR),
           CI_upr = logLR + z*sqrt(var_logLR)
           )
  lr_list <- lr_list %>%
    mutate(num_accept = num_p_value > (1-CI),
           den_accept = den_p_value > (1-CI))
  if(only_accepted){
    lr_list <- lr_list %>% filter(num_accept | den_accept)
  }
  lr_list <- lr_list %>%
    arrange(!!num_groups_, !!den_groups_) %>%
    mutate(across(where(is.factor), paste)) %>%
    rename(numerator = !!num_groups_, denominator = !!den_groups_) %>%
    mutate(null_in_CI = sign(CI_lwr) != sign(CI_upr)) %>%
    select(numerator, denominator, logLR, var_logLR, CI_lwr, CI_upr, null_in_CI, num_accept, den_accept, everything())
  ## Status: One, Ambiguous, or rejected
  status <- lr_list %>%
    mutate(min_logP = ifelse(den_logP > num_logP, den_logP, num_logP)) %>%
    top_n(n = 1, wt = min_logP) %>%
    top_n(n = 1, wt = desc(logLR)) %>%
    mutate(
      status = case_when(
        xor(num_accept, den_accept) ~ "Accepted",
        !num_accept & !den_accept ~ "Rejected",
        num_accept & den_accept & !null_in_CI ~ "Accepted",
        null_in_CI ~ "Ambiguous",
        TRUE ~ "Check"
      ),
      status = case_when(
        status == "Accepted" ~ paste0("Accepted@", ifelse(num_accept, numerator, denominator)),
        status == "Ambiguous" ~ paste0("Ambiguous@", paste(numerator, denominator, sep = " = ")),
        TRUE ~ status
      )) %>% pull(status)
  #
  list(LR = lr_list, status = unlist(strsplit(status, "@")))
}


