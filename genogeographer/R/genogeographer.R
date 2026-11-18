#'
#' genogeographer.devel: Methods for analysing forensic Ancestry Informative Markers
#'
#' The genogeographer.devel package provides:
#' genogeo() app_genegeo()
#'
#' @section genogeo functions:
#' See ?genogeo
#'
#' @docType package
#' @name genogeographer.devel
#' @importFrom dplyr starts_with distinct filter select mutate rename row_number
#' @importFrom dplyr top_n desc pull funs rowwise ungroup arrange case_when bind_rows slice_sample bind_cols
#' @importFrom dplyr vars group_by summarise n everything between group_vars slice count
#' @importFrom dplyr full_join inner_join right_join anti_join semi_join left_join across c_across slice_max
#' @importFrom forcats fct_reorder fct_inorder fct_rev
#' @importFrom purrr map_lgl map_int map2 set_names map_dbl map_chr
#' @importFrom tidyr unnest nest crossing extract unite pivot_longer spread
#' @importFrom tidyselect any_of all_of
#' @importFrom magrittr "%>%"
#' @importFrom readr write_csv read_csv
#' @importFrom leaflet colorNumeric leaflet addTiles addProviderTiles addPolygons addCircles providers addLegend leafletOutput renderLeaflet providerTileOptions
#' @importFrom rlang sym quo ":="
#' @importFrom knitr kable
#' @importFrom glue glue
#' @importFrom htmltools h5 includeHTML p img
#' @importFrom waiter withProgressWaitress useWaitress
#' @importFrom rio import get_ext
#' @importFrom rmarkdown render pdf_document html_document word_document
#' @importFrom stats optimize pnorm qnorm rbeta rbinom rpois sd setNames weighted.mean median na.omit
#' @importFrom utils data packageDescription packageVersion read.csv head globalVariables
#' @importFrom DT formatStyle styleEqual datatable DTOutput renderDT
#' @importFrom tibble tibble as_tibble enframe deframe
#' @importFrom shiny tags HTML observeEvent renderUI h3 span updateSelectInput
#' @importFrom shiny uiOutput div downloadHandler selectInput reactive outputOptions modalButton
#' @importFrom shiny downloadLink reactiveValues eventReactive checkboxGroupInput selectizeInput modalDialog
#' @importFrom shiny updateSelectizeInput textInput showModal
#' @importFrom shiny radioButtons downloadButton hr withProgress h4 fileInput sliderInput actionButton icon helpText actionLink
#' @importFrom shiny getCurrentOutputInfo getDefaultReactiveDomain req isolate tagList shinyApp
#' @importFrom shiny renderText textOutput moduleServer NS showNotification addResourcePath
#' @importFrom shinyjs runjs useShinyjs hidden disable enable show hide delay html
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyWidgets sliderTextInput
#' @importFrom ggplot2 ggplot aes labs geom_point guides geom_polygon scale_colour_manual scale_shape_manual scale_fill_manual label_both
#' @importFrom ggplot2 geom_errorbarh scale_x_reverse map_data coord_cartesian theme_bw fortify theme element_blank geom_vline facet_grid
#' @importFrom ggplot2 xlim geom_segment element_rect scale_color_manual facet_wrap
#' @importFrom maps map
#' @importFrom plotly plot_ly subplot config layout add_markers plotlyOutput renderPlotly ggplotly
#' @importFrom patchwork wrap_plots
#' @importFrom parallel mclapply
#' @importFrom grDevices col2rgb rgb
#' @importFrom bslib page_navbar bs_theme bs_theme_update bs_add_rules nav_panel nav_menu sidebar layout_sidebar card card_header layout_columns
utils::globalVariables("where")
utils::globalVariables("ggg_package")
utils::globalVariables(c(".", "locus", "x0", "genotype", "pop_label", "neg", "pos", "g1", "g2", "n2_x0",
                         "lat", "lon", "convex_hull", "freq", "score", "z_exp", "ggg", "ggg_allele_list", "CI_logP",
                         "logP", "varlogP", "z_raw", "z_var", "z_score", "p_value", "logP_lwr", "logP_upr", "min_logP",
                         "logLR", "var_logLR", "null_in_CI", "num_den", "CI_lwr", "CI_upr", "p_text",
                         "num", "Numerator", "numerator", "num_logP", "num_varlogP", "num_p_value", "num_accept",
                         "den", "denominator", "den_p_value", "den_logP", "den_varlogP", "den_accept", "tool_tip",
                         "meta", "accept", "p1", "p2", "pop", "latlon", "colour", "allele",
                         "P1", "P2", "n1", "n2", "x1", "x2", "f1", "f2", "pair", "admix", "Kidd", "Seldin",
                         "Significantly different", "lab", "population", "metapopulation", "db_x1",
                         "MH", "MH_raw", "sample_name", "geno", "likelihood",
                         "phased", "allele1", "allele2", "position", "snp", "n_snps_found", "n_snps_expected", "allele",
                         "devi", "critical_value", "log10_likelihood", "log10_LR",
                         "meta_num", "meta_den", "accept_num", "accept_den", "p_value_num", "p_value_den",
                         "neg_log10_p", "accept_chr", "sigma", "mu", "devi_lwr", "devi_upr"))
NULL
