library(shiny)
library(genogeographer.devel)

options(shiny.maxRequestSize = 500 * 1024^2) ## 500 Mb

#' The object `db_list` should be initialised prior to lauching
#' For the settings at Department of Forensic Medicine, Section of Forensic Genetics,
#' University of Copenhagen, Copenhagen, Denmark we have a file "db_list.Rds"

#' This loads three databases generated using `make_x1()`
#' on a spreadsheet data containing the necessary columns and information.

#' global.R must define `db_list` and `reporting_panel` (optionally `mh_data`).
#' www/ must contain meta_structure.png and meta_pca.png (served from the app folder).
source("global.R", encoding = "UTF-8")

shinyApp(
  ui = genogeographer.devel::ui_api(),
  server = function(input, output, session){
    genogeographer.devel::server_api(input = input, output = output, session = session)
    }
)

