ggplotly_LR <- function(x, ...) {
  xh <- sum(x$h)
  x <- x$plot
  # browser()
  info <- shiny::getCurrentOutputInfo()
  height <- if (is.function(info$height)) info$height()
  width <- if (is.function(info$width)) info$width()
  plotly::ggplotly(x, width = width, height = width*0.8*xh, ...)
}


single_LRplot <- function(data, num, colours, x_lim, high = NULL){
  # data <- data %>% mutate(y = as.integer(denominator), denominator = paste(denominator))
  plot <- data %>%
    ggplot(aes(x = logLR, y = denominator, xmin = CI_lwr, xmax = CI_upr, color = denominator, text = tool_tip)) +
    labs(y = "Denominator", x = "log10 LR") + geom_vline(xintercept = 0, lty = 2, color = "#999999") +
    facet_wrap(~numerator, labeller = label_both) +
    scale_color_manual(values = colours$den) + guides(colour = "none") +
    coord_cartesian(xlim = x_lim) +
    theme_bw() + theme(strip.background = element_rect(fill = colours$num[num]))
  if(!is.null(high)){
    high <- high %>% filter(numerator == num)
    if(nrow(high) > 0){
      high <- data %>%
        select(denominator) %>%
        right_join(high, by = "denominator") %>%
        mutate(neg = -1000, pos = 1000,
               denominator = factor(denominator, levels = levels(data$denominator)))
      plot <- plot + geom_segment(inherit.aes = FALSE, data = high,
                                  mapping = aes(xend = neg, x = pos, yend = denominator, y = denominator),
                                  lwd = 5, color = "#999999", alpha = 0.3)
    }
  }
  plot + geom_errorbarh(height = 0.1) + geom_point()
}

#' LR plot
#'
#' Plots the LR of population likelihoods
#' @param result Result from \code{genogeo}. At least one of \code{result} and \code{LR_list} is needed.
#' @param LR_list Result from \code{LR_table}. At least one of \code{result} and \code{LR_list} is needed.
#' @param rows Which rows from LR list (or the computed) are used
#' @param theme_ The ggplot2 theme
#' @param ... Additional arguments passed to \code{LR_table}
#' @return A plot
#' @export
LR_plot <- function(result = NULL, LR_list = NULL, rows = NULL, plot_ly = TRUE, shiny = FALSE, which = NULL, ...){
  groups <- names(result)[1]
  groups_ <- sym(groups)
  grouping <- if(groups == "pop") "population" else "metapopulation"
  grouping_ <- sym(grouping)
  ##
  . <- NULL
  numerator <- NULL
  denominator <- NULL
  ##
  if(is.null(LR_list)) LR_list <- LR_table(z_df = result, ...)$LR
  if(is.null(result)) {
    colour_den <- LR_list %>% select(!!groups_) %>% mutate(colour_den = "#000000")
    colour_num <- LR_list %>% select(!!groups_) %>% mutate(colour_num = "#FFFFFF")
    }
  else{
    colour_den <- bar_colour(result[,c("logP", "accept", groups)])
    colour_num <- bar_colour(result[,c("logP", "accept", groups)], alpha = 0.5)
  }
  # if(!is.null(rows)) print(paste(rows, sep = ", "))
  ### LR_list defines order
  Num_order <- result %>% mutate(!!groups_ := fct_reorder(factor(!!groups_), desc(logP))) %>% pull(!!groups_) %>% levels()
  db_info <- attr(result, "info")
  if(!is.null(db_info)){
    Num_ordering <- db_info %>%
      mutate(
        !!groups_ := factor(!!groups_, levels = Num_order),
        !!grouping_ := fct_reorder(!!grouping_, !!groups_, .fun = as.numeric)
      ) %>% select(!!groups_, !!grouping_) %>% deframe()
    LR_list <- LR_list %>%
      mutate(
        Numerator = Num_ordering[numerator],
        numerator = paste(Numerator),
        denominator = paste(Num_ordering[denominator]),
      )
    names(colour_den) <- paste(Num_ordering[names(colour_den)])
    names(colour_num) <- paste(Num_ordering[names(colour_num)])
  }
  else LR_list <- LR_list %>% mutate(Numerator = factor(numerator, levels = Num_order))
  ##
  LR_list <- LR_list %>% # rename(Numerator = numerator) %>%
    mutate(
      # denominator_num = tidytext::reorder_within(denominator, desc(logLR), Numerator),
      tool_tip = glue("<b>Numerator:</b> {Numerator}<br>",
                      "<b>Denominator:</b> {denominator}<br>",
                      "<b>log<sub>10</sub> LR:</b> {round(logLR, 2)}<br>",
                      "<b>CI(log<sub>10</sub> LR):</b> [{round(CI_lwr,2)};{round(CI_upr,2)}]<br>"
                      )
      )
  ## Use DT table rows to select rows to plot?
  if(!is.null(which)) LR_highlight <- LR_list %>% select(numerator, denominator) %>% slice(which)
  else LR_highlight <- NULL
  LR_range <- c(min(LR_list$CI_lwr), max(LR_list$CI_upr))
  if(LR_range[1] > 0) LR_range[1] <- 0
  if(LR_range[2] < 0) LR_range[2] <- 0
  LR_plots <- LR_list %>%
    split(.$Numerator, drop = TRUE) %>%
    map(~.x %>% mutate(denominator = fct_reorder(denominator, logLR, .desc = TRUE))) %>%
    purrr::imap(~
      list(plot = single_LRplot(data = .x, num = .y, x_lim = LR_range,
                                colours = list(den = colour_den, num = colour_num), high = LR_highlight),
           height = .x %>% distinct(denominator) %>% nrow()
      )) # facet_col(~Numerator, scale = "free_y", space = "free", labeller = label_both)
  if(plot_ly){
    LR_plots_heights_ <- LR_plots %>% map_int(~.x$height)
    LR_plots_heights <- LR_plots_heights_ %>% {./sum(.)} %>% unname()
    sum_height <- sum(LR_plots_heights_)
    min_rel <- min(LR_plots_heights)
    LR_height <- max(500, sum_height*50)
    if(min_rel*LR_height < 100) LR_height <- 100/min_rel
    # if(shiny){
    #   info <- shiny::getCurrentOutputInfo()
    #   plot_height <- info$height()
    #   plot_width <- info$width()
    # }
    LR_plots <- LR_plots %>% purrr::map(~ ggplotly(.x$plot, tooltip = c("text"), height = LR_height))
    LR_plots <- subplot(LR_plots, nrows = length(LR_plots), shareX = TRUE, #shareY = TRUE,
                        heights = LR_plots_heights, margin = 0.025) %>%
      config(displayModeBar = FALSE) %>% layout(showlegend = FALSE)
  }
  else{
    LR_plots_heights <- LR_plots %>% map_int(~.x$height)
    LR_plots <- LR_plots %>% purrr::map(~ .x$plot) %>% wrap_plots(ncol = 1, heights = LR_plots_heights)
  }
  LR_plots
}


LR_plot_ly <- function(result = NULL, LR_list = NULL, ...){
  groups <- names(result)[1]
  groups_ <- sym(groups)
  grouping <- if(groups == "pop") "population" else "metapopulation"
  grouping_ <- sym(grouping)
  ##
  . <- NULL
  numerator <- NULL
  denominator <- NULL
  ##
  if(is.null(LR_list)){
    LR_list <- LR_table(z_df = result, ...)$LR
  }
  ##
  if(is.null(result)) {
    colour_den <- LR_list %>% select(!!groups_) %>% mutate(colour_den = "#000000")
    colour_num <- LR_list %>% select(!!groups_) %>% mutate(colour_num = "#FFFFFF")
  }
  else{
    Num_order <- result %>% mutate(!!groups_ := fct_reorder(factor(!!groups_), desc(logP))) %>% pull(!!groups_) %>% levels()
    colour_den <- bar_colour(result[,c("logP", "accept", groups)])
    colour_num <- bar_colour(result[,c("logP", "accept", groups)])
  }
  ### LR_list defines order
  db_info <- attr(result, "info")
  if(!is.null(db_info)){
    Num_ordering <- db_info %>%
      mutate(
        !!groups_ := factor(!!groups_, levels = Num_order),
        !!grouping_ := fct_reorder(!!grouping_, !!groups_, .fun = as.numeric)
      ) %>% select(!!groups_, !!grouping_) %>% deframe()
    LR_list <- LR_list %>%
      mutate(
        Numerator = Num_ordering[numerator],
        numerator = paste(Numerator),
        denominator = paste(Num_ordering[denominator]),
      )
    names(colour_den) <- paste(Num_ordering[names(colour_den)])
    names(colour_num) <- paste(Num_ordering[names(colour_num)])
  }
  else LR_list <- LR_list %>% mutate(numerator = factor(numerator, levels = Num_order))
  ##
  LR_list <- LR_list %>%
    mutate(
      tool_tip = glue("<b>Numerator:</b> {Numerator}<br>",
                      "<b>Denominator:</b> {denominator}<br>",
                      "<b>log<sub>10</sub> LR:</b> {round(logLR, 2)}<br>",
                      "<b>CI(log<sub>10</sub> LR):</b> [{round(CI_lwr,2)};{round(CI_upr,2)}]<br>"
      )
    )
  ## Use DT table rows to select rows to plot?
  LR_range <- c(min(LR_list$CI_lwr), max(LR_list$CI_upr))
  if(LR_range[1] > 0) LR_range[1] <- 0
  if(LR_range[2] < 0) LR_range[2] <- 0
  LR_range[1] <- if(LR_range[1] < 0) LR_range[1]*1.1 else LR_range[1]/1.1
  LR_range[2] <- if(LR_range[2] < 0) LR_range[2]/1.1 else LR_range[2]*1.1
  # browser()
  LR_list <- LR_list %>%
    unite(num_den, numerator, denominator, sep = " vs ", remove = FALSE) %>%
    mutate(num_den = fct_inorder(num_den),
           tool_tip = paste0("<b>Numerator:</b> ", numerator, "<br>",
                             "<b>Denominator:</b> ", denominator, "<br>",
                             "<b>log<sub>10</sub> LR:</b> ", round(logLR, 2), "<br>",
                             "<b>CI(log<sub>10</sub> LR):</b> [", round(CI_lwr,2), "; ",round(CI_upr,2), "]<br>")
    )
  LR_list %>% plot_ly(x = ~logLR, y = ~num_den, color = ~numerator, colors = colour_num) %>%
    add_markers(x = ~logLR, yaxis = "y2", opacity = 0, hoverinfo = "none") %>%
    add_markers(error_x = ~list(array = sqrt(var_logLR)), text = ~tool_tip, hoverinfo = "text") %>%
    layout(
      xaxis = list(
        range = LR_range,
        title = "log<sub>10</sub> LR",
        zerolinedash = "dashed"),
      yaxis = list(
        side = "left",
        title = "Numerator",
        autorange = "reversed",
        ticktext = as.list(LR_list$numerator),
        tickvals = as.list(LR_list$num_den),
        tickmode = "array"),
      yaxis2 = list(
        side = "right",
        title = "Denominator",
        overlaying = "y",
        ticktext = as.list(LR_list$denominator),
        tickvals = as.list(LR_list$num_den),
        tickmode = "array"),
    showlegend = FALSE,
    margin = list(l = 0, r = 150, b = 50, t = 10, pad = 0)
    ) %>%
    config(displayModeBar = FALSE)
}
