#'
#' genogeographer: Methods for analysing forensic Ancestry Informative Markers
#'
#' The genogeographer package provides:
#' genogeo() app_genegeo()
#'
#' @section genogeo functions:
#' See ?genogeo
#'
#' @docType package
#' @name genogeographer
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
#' @importFrom leaflet colorNumeric leaflet addTiles addProviderTiles addPolygons addCircles providers addLegend leafletOutput renderLeaflet
#' @importFrom rlang sym quo ":="
#' @importFrom knitr kable
#' @importFrom glue glue
#' @importFrom htmltools h5 includeHTML p img
#' @importFrom waiter withProgressWaitress useWaitress
#' @importFrom shinythemes shinytheme
#' @importFrom rio import get_ext
#' @importFrom rmarkdown render pdf_document html_document word_document
#' @importFrom stats optimize pnorm qnorm rbeta rbinom rpois sd setNames weighted.mean median na.omit
#' @importFrom utils data packageDescription packageVersion read.csv head globalVariables
#' @importFrom DT formatStyle styleEqual datatable DTOutput renderDT
#' @importFrom tibble tibble as_tibble enframe deframe
#' @importFrom shiny shinyApp tags wellPanel HTML sidebarLayout sidebarPanel fluidPage observeEvent renderUI h3 mainPanel span updateSelectInput
#' @importFrom shiny uiOutput column div downloadHandler renderTable observe selectInput reactive outputOptions verticalLayout modalButton
#' @importFrom shiny downloadLink reactiveValues eventReactive checkboxGroupInput conditionalPanel selectizeInput renderPlot modalDialog
#' @importFrom shiny fluidRow plotOutput brushOpts clickOpts hoverOpts nearPoints updateSelectizeInput brushedPoints textInput showModal
#' @importFrom shiny radioButtons downloadButton hr withProgress h4 titlePanel fileInput sliderInput actionButton icon helpText actionLink
#' @importFrom shiny getCurrentOutputInfo getDefaultReactiveDomain req bootstrapPage navbarPage tabPanel isolate tagList navbarMenu
#' @importFrom shiny renderText textOutput
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
                         "Significantly different", "lab", "population", "metapopulation", "db_x1"))
NULL
