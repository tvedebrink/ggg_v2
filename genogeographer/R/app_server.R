#' @title Server part for Shiny app
#' @description server function to be used in the `app_genogeo` shiny app.
#' Wires up the shiny modules that make up the app: analysis, reference
#' database information, and adding new reference populations.
#' @param input The list of inputs (leave empty)
#' @param output The list of outputs (leave empty)
#' @param session The shiny session - to refer to session elements (leave empty)
#' @export
server_api <- function(input, output, session){
  if(!exists("db_list")) db_list <- get("db_list", envir = -2)
  if(!exists("reporting_panel")) reporting_panel <- get("reporting_panel", envir = -2)
  if(!exists("mh_data")) mh_data <- tryCatch(get("mh_data", envir = -2), error = function(e) NULL)

  reactive_db_list <- reactiveValues(db = db_list)

  mod_analyse_server("analyse", db_list = reactive_db_list, reporting_panel = reporting_panel)
  mod_mh_analyse_server("mh", mh_data = mh_data, db_list = reactive_db_list)
  mod_reference_databases_server("refdb", db_list = reactive(reactive_db_list$db))
  mod_add_reference_server("addref")
}
