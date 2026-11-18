#' Base leaflet map with the CARTO tile layer shared by the AIMs and MH maps
#' @noRd
ggg_base_map <- function(data = NULL){
  leaflet(data) %>%
    addTiles(
      url = "https://basemaps.cartocdn.com/rastertiles/light_all/{z}/{x}/{y}.png?key=cb1_3xnu_1_b9d34ccb5567a6125344645d",
      attribution = paste(
        '&copy; <a href="https://carto.com/attributions">CARTO</a>,',
        '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
      )
    )
}

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
  ggg_map <- ggg_base_map(db_info)
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

#' Map centres of the MH meta-populations, reusing the centres of the AIMs reference databases
#'
#' An MH meta-population is matched on the AIMs `meta` code first (so a database that uses the MH codes needs no
#' mapping) and then on the name of the AIMs `metapopulation`. Codes without a match get no centre.
#' @param db The named list of AIMs reference databases (`db_list`).
#' @param metas MH meta-population codes to look up.
#' @return Tibble with `meta` (MH code), `lat`, `lon` and `source` (AIMs meta-population and database).
#' @noRd
mh_meta_centres <- function(db, metas = mh_meta_populations()){
  name_pattern <- c(EUR = "^europe", AFR = "sub.?saharan|^africa", EAS = "^east.?asia",
                    SAS = "south.{0,4}central.{0,4}asia|south.?asia", ME = "middle.?east",
                    NAM = "native.?americ", AAM = "african.?americ")
  info <- purrr::imap(db, function(d, nm){
    i <- attr(d$meta$db, "info")
    if(is.null(i) || !all(c("meta", "metapopulation", "lat", "lon") %in% names(i))) return(NULL)
    i <- i[!is.na(i$lat) & !is.na(i$lon), c("meta", "metapopulation", "lat", "lon")]
    i$db <- nm
    i
  }) %>% purrr::compact() %>% bind_rows()
  out <- if(nrow(info) == 0) NULL else bind_rows(lapply(metas, function(m){
    hit <- info[tolower(info$meta) == tolower(m), ]
    if(nrow(hit) == 0 && m %in% names(name_pattern)) hit <- info[grepl(name_pattern[[m]], info$metapopulation, ignore.case = TRUE), ]
    if(nrow(hit) == 0) return(NULL)
    tibble(meta = m, lat = hit$lat[1], lon = hit$lon[1], source = paste0(hit$metapopulation[1], " (", hit$db[1], ")"))
  }))
  if(is.null(out) || nrow(out) == 0) tibble(meta = character(), lat = numeric(), lon = numeric(), source = character()) else out
}

#' Plot the MH profile probabilities on a map
#'
#' Like `leaflet_plot()` for the AIMs analysis, but with one circle at the centre of each meta-population
#' (no convex hulls). Meta-populations without a centre in `centres` are not drawn.
#' @param result_df A single sample's slice of `mh_profile_likelihood()`'s `result` tibble.
#' @param centres Centres of the meta-populations, see `mh_meta_centres()`.
#' @param ci Confidence level of the intervals shown in the labels.
#' @param n_ind Optional named vector with the number of reference individuals per meta-population (circle size).
#' @noRd
leaflet_mh_plot <- function(result_df, centres, ci = 0.95, n_ind = NULL){
  z <- stats::qnorm((1 - ci) / 2, lower.tail = FALSE)
  df <- result_df %>% mutate(logP = log10_likelihood)
  df$colour <- unname(bar_colour(df[, c("logP", "accept", "meta")]))
  df$log_P <- df$logP * (-1)^(!df$accept)
  df$individuals <- if(is.null(n_ind)) NA_real_ else unname(n_ind[as.character(df$meta)])
  df$lab <- lapply(seq_len(nrow(df)), function(i){
    HTML(sprintf("<b>%s</b>%s<br/><b><i>z</i>-score:</b> %0.2f (%s)<br/><b>log<sub>10</sub> P</b>: %0.2f %s",
                 df$meta[i], if(is.na(df$individuals[i])) "" else sprintf(" (<i>%d individuals</i>)", round(df$individuals[i])),
                 df$z_score[i], ifelse(df$accept[i], "accepted", "rejected"), df$logP[i],
                 if(is.na(df$se_log10[i])) "" else sprintf("[%0.2f, %0.2f]", df$logP[i] - z * df$se_log10[i], df$logP[i] + z * df$se_log10[i])))
  })
  mapped <- inner_join(df, centres[, c("meta", "lat", "lon")], by = "meta")
  ggg_map <- ggg_base_map(mapped)
  if(nrow(mapped) == 0) return(ggg_map)
  weight <- ifelse(is.na(mapped$individuals), 10, sqrt(mapped$individuals))
  ggg_map <- ggg_map %>% addCircles(lng = ~lon, lat = ~lat, weight = weight, label = ~lab, color = ~colour)
  legend_cols <- mapped %>% select(logP, accept, colour, log_P) %>%
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

