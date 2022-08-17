# library(shiny)
# library(shinythemes)
# library(shinycssloaders)
# library(tidyverse)
# library(DT)
# library(shinyjs)
# library(glue)
# library(plotly)
# library(ggforce)
# library(patchwork)
# library(shinyWidgets)
# library(parallel)
#
# library(leaflet)
# library(leaflet.providers)
# library(waiter)

ggg_package <- "genogeographerDEVEL"

#' Shiny application for GenoGeoGrapher
#' @param db_list A named list of databases of reference populations.
#' Each component is expected to be returned from \code{make_x1}.
#' @param reporting_panel Logical. Should report generate and download be available after sample analysis.
#' @export
app_genogeo <- function(db_list = NULL, reporting_panel = TRUE){
  ui <- ui_api()
  server <- function(input, output, session) server_api(input = input, output = output, session = session)
  shinyApp(ui, server)
}

# source("R/app_server.R")
# source("R/app_ui.R")
# app_genogeo()
