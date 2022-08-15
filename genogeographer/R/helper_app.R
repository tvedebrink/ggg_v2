## app_helpers

## Formatting of result tables

result_table <- function(result, flat = FALSE, lr_listed = ""){
  groups <- names(result)[1]
  groups_ <- sym(groups)
  ## build fixes : start ##
  logP <- NULL; varlogP <- NULL; logP_lwr <- NULL; logP_upr <- NULL; z_score <- NULL
  p_value <- NULL; n <- NULL; meta <- NULL; pop <- NULL; . <- NULL; accept <- NULL
  lat <- NULL; lon <- NULL
  ## build fixes : end ##
  if (is.null(result)) return(NULL)
  row_colours_hex <- bar_colour(result[,c("logP","accept",groups), drop = FALSE], alpha = 0.1)
  row_colours <- rgba2rgb(row_colours_hex)
  ##
  result <- result %>% ## mutate(n = n/2) %>%
    mutate(
      across(where(is.numeric), ~round(.x,3)),
      logP_lwr = paste0("[", logP_lwr),
      logP_upr = paste0(logP_upr, "]")
      ) %>%
    unite(CI_logP, logP_lwr, logP_upr, sep = "; ") %>%
    select(!!groups_, logP, CI_logP, z_score, p_value) %>%
    rename(
      !!glue("CI[log10 P(G|{to_upper(groups_)})]") := CI_logP,
      !!glue("log10 P(G|{to_upper(groups_)})") := logP,
      !!glue("var[log10 P(G|{to_upper(groups_)})]") := varlogP,
      `z-score` = z_score,
      `p-value` = p_value
    )
  if(flat) return(kable(result))
  result %>% set_names(sub("10", "<sub>10</sub>", names(.))) %>%
    datatable(rownames=FALSE, filter = "bottom", escape = FALSE,
                  selection = list(mode = "multiple", target = "row"), #'none',
                  extensions = 'Buttons', options = list(
                    dom = 'Blfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                  )
    ) %>%
    formatStyle(columns = 1,
                target = "row",
                backgroundColor = styleEqual(result[[groups]], row_colours[result[[groups]]])) %>%
    formatStyle(columns = 1,
                target = "row",
                fontWeight = styleEqual(lr_listed, rep('bold', length(lr_listed))))
}

## Formatting of LR output

LR_list <- function(result = NULL, LR_tab = NULL, lr_pops = NULL, CI, accepted, flat = FALSE){
  groups <- names(result)[1]
  groups_ <- sym(groups)
  grouping_ <- sym(if(groups == "pop") "population" else "metapopulation")
  ## build fixes : start ##
  numerator <- NULL; Numerator <- NULL; denominator <- NULL; Denominator <- NULL; logLR <- NULL
  `log10 LR` <- NULL; var_logLR <- NULL; CI_lwr <- NULL; CI_upr <- NULL; `CI(log10 LR)` <- NULL
  null_in_CI <- NULL; `Null in CI` <- NULL;   z_score <- NULL;  accept <- NULL; . <- NULL
  ## build fixes : end ##
  if(is.null(LR_tab)){
    if(is.null(result)) return(NULL)
    if(is.null(lr_pops)) return(NULL)
    lr_list <- LR_table(z_df = result, who = lr_pops, CI = CI, only_accepted = !accepted)
    if(nrow(lr_list)==0) return(NULL)
  }
  else lr_list <- LR_tab
  lr_list <- lr_list %>% mutate(across(where(is.numeric), ~round(.x,3)))
  min_z_pop <- result %>% filter(accept) %>% slice_max(n = 1, order_by = desc(z_score)) %>% pull(!!groups_)
  if(length(min_z_pop)==0) min_z_pop <- ""
  db_info <- attr(result, "info")
  if(!is.null(db_info)){
    key_name <- db_info %>% select(!!groups_, !!grouping_) %>% deframe()
    if(min_z_pop != "") min_z_pop <- key_name[min_z_pop]
    lr_list <- lr_list %>%
      mutate(
        numerator = key_name[numerator],
        denominator = key_name[denominator]
        )
    result[[groups]] <- key_name[result[[groups]]]
  }
  ## lr_list <- format(lr_list, digits = 3, nsmall = 3)
  lr_list <- lr_list %>%
    mutate(`CI(log10 LR)` = paste0("[", CI_lwr, "; ", CI_upr,"]")) %>%
    rename(
      Numerator = numerator,
      Denominator = denominator,
      `log10 LR` = logLR,
      `var(log10 LR)` = var_logLR
    ) %>% mutate(`Significantly different` = ifelse(null_in_CI, "No", "Yes")) %>%
    select(Numerator, Denominator, `log10 LR`, `CI(log10 LR)`, `Significantly different`) %>%
    mutate(across(where(is.character), factor))
  ##
  if(flat) return(kable(lr_list))
  lr_list %>% set_names(sub("10", "<sub>10</sub>", names(.))) %>%
    datatable(rownames=FALSE, filter = 'bottom', extensions = 'Buttons', escape = FALSE,
                        selection = "none", # list(mode = "multiple", target = "row"), #'none',
                        options = list(
                          dom = 'Blfrtip',
                          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                          autoWidth = TRUE,
                          lengthMenu = list(c(10, 25, 50, 100, -1),
                                            c("10", "25", "50", "100", "All")),
                          columnDefs = list(list(width = '30%', targets = c(3)))
                          )) %>%
    formatStyle(
      columns = c('Numerator', 'Denominator'),
      fontWeight = styleEqual(min_z_pop, 'bold'),
      color = styleEqual(result[[groups]], bar_colour(result[,c("logP","accept",groups)]))
  )
}

## write variable names in italic font with , and 'and' separation
and_text <- function(x, anchor = "i", pre = "", post = ""){
  if(length(x) == 0) return("")
  if(length(x) == 1) return(paste0(pre,"<",anchor,">",x,"</",anchor,">",post))
  else paste0(pre,"<",anchor,">",paste0(x[-length(x)], collapse = paste0("</",anchor,">, <",anchor,">")), "</",anchor,"> and <",anchor,">",x[length(x)],"</",anchor,">",post)
}

## Same as subset but returns NULL in case of empty
nullset <- function(x, ...){
  x <- subset(x, ...)
  if(length(x) == 0) return(NULL)
  x
}
