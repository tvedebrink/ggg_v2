library(shiny)

#' The object `db_list` should be initialised prior to lauching
#' For the settings at Department of Forensic Medicine, Section of Forensic Genetics,
#' University of Copenhagen, Copenhagen, Denmark we have a file "db_list.Rds"

#' This loads three databases generated using `make_x1()`
#' on a spreadsheet data containing the necessary columns and information.

source("global.R", encoding = "UTF-8")

shinyApp(
  ui = genogeographer::ui_api(),
  server = function(input, output, session){
    genogeographer::server_api(input = input, output = output, session = session)
    }
)

