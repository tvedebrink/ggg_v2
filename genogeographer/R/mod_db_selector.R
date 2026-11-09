#' Reference database selector module - UI
#'
#' Allows uploading a user-provided reference database (.rds) and selecting
#' which loaded database to use for the analysis.
#' @param id Module id
#' @noRd
mod_db_selector_ui <- function(id){
  ns <- NS(id)
  uiOutput(ns("dbs"))
}

#' Reference database selector module - server
#'
#' @param id Module id
#' @param db_list A `reactiveValues` object with a `db` element holding the named list of reference databases.
#' @return A list of reactives: `snp_set` (name of the selected database)
#' @noRd
mod_db_selector_server <- function(id, db_list){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observeEvent(input$db_add, {
      req(input$db_add)
      user_db <- readRDS(file = input$db_add$datapath)
      db_names <- user_db %>% purrr::map(names) %>% unlist() %>% unname() %>% unique()
      if(length(setdiff(db_names, c("pop", "meta"))) == 0){
        db_list$db <- c(db_list$db, user_db)
        updateSelectInput(session, inputId = "snp_set", selected = names(user_db))
      }
    })

    output$dbs <- renderUI({
      div(
        fileInput(inputId = ns("db_add"), "Upload own reference database:", accept = ".rds", width = "100%"),
        helpText("Ensure that the '.rds'-file is created using genogeographer, e.g. by following the steps under tab 'Add referece population'"),
        selectInput(inputId = ns("snp_set"), label = "Select reference database:", choices = names(db_list$db))
      )
    })

    list(
      snp_set = reactive(input$snp_set),
      has_db = reactive(!is.null(db_list$db) && length(db_list$db) > 0)
    )
  })
}
