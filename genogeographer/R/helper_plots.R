to_upper <- function(x) paste0(toupper(substr(x,1,1)), tolower(substr(x, 2, nchar(x))))

ggplotly_ <- function(x, ...) {
  info <- shiny::getCurrentOutputInfo()
  height <- if (is.function(info$height)) info$height()
  width <- if (is.function(info$width)) info$width()
  plotly::ggplotly(x, width = width, height = width*0.8, ...)
}

#' bar_colour
#'
#' Creates the colour scale for the accepted and rejected populations based on z-score and the log likelihood (log P).
#' @param df A data.frame with at least three coloums. The first column is the logP, the second logical (z_score accept/reject), the third a unique naming column.
#' @param alpha Should the alpha opacity be applied? And what value, 1 = solid, 0 = transparent.
#' @export
bar_colour <- function(df, alpha = 1){
  if(alpha>1) alpha <- 1
  if(alpha<0) alpha <- 0
  stopifnot(ncol(df)==3)
  logP <- df$logP; accept <- df$accept; nm <- df[,-which(names(df) %in% c("logP", "accept")), drop = TRUE]
  if(any(accept)) col_accept <- leaflet::colorNumeric("Reds", c(logP[accept],min(logP[accept])*1.05)) ## Accepted
  if(any(!accept)) col_reject <- leaflet::colorNumeric("GnBu", c(logP[!accept], min(logP[!accept])*1.5)) ## Rejected
  bar_cols <- rep("#000000", length(accept))
  if(any(accept)){ bar_cols[accept] <- col_accept(logP[accept]) }
  if(any(!accept)){ bar_cols[!accept] <- col_reject(logP[!accept]) }
  bar_cols <- paste0(bar_cols, sprintf("%02X",as.hexmode(round(alpha*255))))
  setNames(bar_cols, paste(nm))
}

rgba2rgb <- function(hex_rgba){
  rgba <- col2rgb(hex_rgba, alpha = TRUE)
  rgb_ <- rgba[4,]*rgba[1:3,] + (255-rgba[4,])*col2rgb("#FFFFFF")[,rep(1,length(hex_rgba))]
  apply(rgb_, 2, function(RGB) rgb(RGB[1], RGB[2], RGB[3], maxColorValue = 255^2))
}

#' Plot log likelihoods of profiles with approximate confidence intervals
#'
#' Plots the estimated profile probabilities in each population.
#' The colour depends on the profiles likelihood and rejection/acceptance (blue/red) based on z-score
#'
#' @name error_bar_plot
#' @author Torben Tvedebrink, \email{tvede@@math.aau.dk}
#' @param result_df The output from the \code{genogeo} function
#' @param which The populations to highlight
#' @return A barplot of the log likelihoods for each population with confidence limits
#' @export


error_bar_plot <- function(result_df, which = NULL){
  groups <- names(result_df)[1]
  groups_ <- sym(groups)
  grouping <- if(groups == "pop") "population" else "metapopulation"
  grouping_ <- sym(grouping)
  ## build fixes : start ##
  logP <- NULL
  logP_lwr <- NULL
  logP_upr <- NULL
  ## build fixes : end ##
  db_info <- attr(result_df, "info")
  if(!is.null(db_info)){
    key_name <- db_info %>% select(!!groups_, !!grouping_) %>% deframe()
    result_df <- result_df %>% left_join(db_info %>% select(!!groups_, n), by = groups) %>%
      mutate(!!groups_ := key_name[!!groups_])
  }
  result_df <- result_df %>%
    mutate(
      pop_label = sprintf("<b>%s</b> (<i>%s samples</i>)
<b><i>z</i>-score:</b> %0.2f (%s)
<b>log<sub>10</sub> P</b>: %0.2f [%0.2f, %0.2f]",
                          !!groups_, n, z_score, ifelse(accept, "accepted", "rejected"), logP, logP_lwr, logP_upr),
      !!groups_ := fct_reorder(!!groups_, logP)
      )
  logP_range <- c(min(result_df$logP_lwr), max(result_df$logP_upr))
  p1 <- result_df %>%
    ggplot(aes(y=!!groups_,x=logP,xmin=logP_lwr,xmax=logP_upr, colour=!!groups_, text = pop_label)) +
    labs(y="",x=expression(log[10]~P(Genotype~"|"~Population))) +
    guides(colour="none") +
    coord_cartesian(xlim = rev(logP_range)) + geom_point() + geom_errorbarh() +
    scale_colour_manual(values = bar_colour(result_df[,c("logP","accept",groups)])) +
    # scale_x_reverse() +
    theme_bw()
  if(!is.null(which)){
    high <- result_df %>% slice(which) %>%
      select(!!groups_) %>%
      mutate(neg = -1000, pos = 1000, !!groups_ := factor(!!groups_, levels = levels(result_df[[groups]])))
    p1 <- p1 + geom_segment(inherit.aes = FALSE, data = high,
                            mapping = aes(xend = neg, x = pos, yend = !!groups_, y = !!groups_),
                            lwd = 5, color = "#999999", alpha = 0.3)
  }
  p1
}

  error_bar_plotly <- function(result_df, which = NULL){
  group <- ifelse(names(result_df)[1] == "pop", "Population", "Metapopulation")
  plot <- error_bar_plot(result_df = result_df, which = which)
  plot <- plot + labs(y="",x= paste0("log<sub>10</sub> P(Genotype | ", group, ")"))
  plot
  # plot %>% ggplotly_(tooltip = c("text")) %>%
  #   layout(showlegend = FALSE) %>%
  #   config(displayModeBar = FALSE)
}
