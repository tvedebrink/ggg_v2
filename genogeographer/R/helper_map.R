#' Plot LTR z-scores on map
#'
#' Plots the results from LRT on a map based on lat/lon info in the database.
#' If no location is found in the data (e.g. using \code{simulte_pops}) nothing is plotted.
#' @param z_df Result from `genogeo()`

leaflet_plot <- function(z_df = NULL){
  db_info <- attr(z_df, "info")
  if(!is.null(db_info)){
    group <- names(db_info)[1]
    group_ <- sym(group)
    groups <- if(group == "pop") "population" else "metapopulation"
    groups_ <- sym(groups)
    db_info <- db_info %>% rowwise() %>%
      mutate(lab = list(HTML(sprintf("<b>%s</b><br/>%s samples", !!groups_, n)))) %>%
      ungroup()
    if(!is.null(z_df)){
      db_info <- z_df %>% left_join(db_info, by = group)
      db_info$colour <- bar_colour(select(db_info, !!group_, logP, accept)) %>% unname()
      db_info <- db_info %>% mutate(log_P = logP*(-1)^(!accept))
      db_info <- db_info %>% rowwise() %>%
        mutate(lab = list(HTML(
          sprintf("
                  <b>%s</b> (<i>%s samples</i>)<br/>
                  <b><i>z</i>-score:</b> %0.2f (%s)<br/>
                  <b>log<sub>10</sub> P</b>: %0.2f [%0.2f, %0.2f]",
                  !!groups_, n, z_score, ifelse(accept, "accepted", "rejected"), logP, logP_lwr, logP_upr)))) %>%
        ungroup()
    }
    db_info <- db_info %>% filter(!(is.na(lat) | is.na(lon)))
  }
  ggg_map <- leaflet(db_info) %>%
    addTiles() %>%
    addProviderTiles(providers$CartoDB.Positron)
  if(is.null(db_info)) return(ggg_map)
  if("convex_hull" %in% names(db_info)){
    conv_df <- db_info %>%
      select(starts_with("meta"), convex_hull, lab, colour)
    for(m in conv_df$meta){
      conv_meta <- conv_df %>% filter(meta == m)
      conv_lab <- conv_meta %>% pull(lab)
      conv_col <- conv_meta %>% pull(colour)
      conv_meta <- conv_meta %>% pull(convex_hull)
      ggg_map <- ggg_map %>%
        addPolygons(data = conv_meta[[1]], lng = ~lon, lat = ~lat, opacity = 0.2, label = conv_lab,
                    fillColor = conv_col, stroke = NA)
    }
  }
  ggg_map <- ggg_map %>% addCircles(lng = ~lon, lat = ~lat, weight = ~sqrt(as.integer(n)), label = ~lab, color = ~colour)
  legend_cols <- db_info %>% select(logP, accept, colour) %>%
    slice_pretty(order_by = "logP", split_by = "accept") %>%
    mutate(labP = sprintf("%0.2f", logP))
  if(nrow(legend_cols) == 1) legend_cols <- bind_rows(legend_cols, legend_cols)
  ggg_map %>% addLegend(data = legend_cols, position = "bottomleft",
                        values = ~log_P, colors = ~colour, labels = ~labP, title = "<b>log<sub>10</sub> P</b>")
}

slice_pretty <- function(df, order_by, split_by = NULL){
  order_by_ <- sym(order_by)
  if(is.null(split_by)){
    df <- df %>% arrange(desc(!!order_by_))
    df_n <- nrow(df)
    if(df_n == 0) return(NULL)
    df_idx <- if(df_n <= 3) seq_len(df_n) else c(1, floor(df_n/2), df_n)
    df_ <- df %>% slice(df_idx)
    return(df_)
  }
  if(!is.null(split_by)){
    dfs <- split(df, df[[split_by]])
    dfs <- lapply(dfs, slice_pretty, order_by = order_by, split_by = NULL)
    dfs %>% bind_rows() %>% arrange(desc(!!order_by_))
  }
}

if(FALSE){

  leaflet_plot(z_df = dane_res_meta)

}

