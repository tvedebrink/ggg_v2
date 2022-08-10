if(FALSE){
  LR_list <- read_rds("LR_list.Rds")
LR_list_ <- LR_list %>% split(.$numerator) %>%
  imap_dfr(~.x %>%
             add_row(num_den = paste0("NA", .y), numerator = "", denominator = "",
                     logLR = -9999)) %>%
  slice(-nrow(.)) %>%
  mutate(num_den_ = as.numeric(fct_rev(fct_inorder(num_den))),
         tool_tip = paste0("<b>Numerator:</b> ", numerator, "<br>",
                  "<b>Denominator:</b> ", denominator, "<br>",
                  "<b>log<sub>10</sub> LR:</b> ", round(logLR, 2), "<br>",
                  "<b>CI(log<sub>10</sub> LR):</b> [", round(CI_lwr,2), "; ",round(CI_upr,2), "]<br>")
         )

tool_tip <-
"<b>Numerator:</b> %{numerator}
<b>Denominator:</b> %{denominator}
<b>log<sub>10</sub> LR:</b> {round(%{logLR}, 2)}
<b>CI(log<sub>10</sub> LR):</b> [{round(%{CI_lwr},2)};{round(%{CI_upr},2)}]"

LR_range <- LR_list_ %>% filter(logLR != -9999) %>%
  pull(logLR) %>% range() %>% {.*c(1/1.3, 1.05)}
if(LR_range[1] > 0) LR_range[1] <- 0
if(LR_range[2] < 0) LR_range[2] <- 0

LR_list_ %>% plot_ly(x = ~logLR, y = ~num_den_) %>%
  add_markers(x = ~logLR, yaxis = "y2", opacity = 0, hoverinfo = "none") %>%
  add_markers(error_x = ~list(array = sqrt(var_logLR)), text = ~tool_tip, hoverinfo = "text") %>%
  layout(
    xaxis = list(range = LR_range, title = "log<sub>10</sub> LR"),
    yaxis = list(
      side = "left",
      title = "Numerator",
      ticktext = as.list(LR_list_$numerator),
      tickvals = as.list(LR_list_$num_den_),
      tickmode = "array"),
    yaxis2 = list(
      side = "right",
      title = "Denominator",
      overlaying = "y",
      ticktext = as.list(LR_list_$denominator),
      tickvals = as.list(LR_list_$num_den_),
      tickmode = "array"),
    showlegend = FALSE,
    margin = list(l = 0, r = 150, b = 50, t = 10, pad = 0)
    ) %>%
  config(displayModeBar = FALSE)

}
