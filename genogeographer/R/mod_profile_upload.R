#' Profile upload module - UI
#'
#' Handles uploading an AIMs profile file and selecting the locus/genotype columns.
#' @param id Module id
#' @noRd
mod_profile_upload_ui <- function(id){
  ns <- NS(id)
  tagList(
    fileInput(ns('profile_file'), 'Choose CSV File', width = "100%",
              accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx', '.xls')),
    uiOutput(ns("column_panel"))
  )
}

#' Profile upload module - server
#'
#' @param id Module id
#' @return A list of reactives: `file`, `data`, `locus_col`, `genotype_col`, `locus_guess`, `genotype_guess`
#' @noRd
mod_profile_upload_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    reactive_columns <- reactiveValues(columns = NULL, locus_guess = NULL, genotype_guess = NULL)

    reactive_read_profile <- reactive({
      if(is.null(input$profile_file)) return(NULL) # User has not uploaded a file yet
      ext <- rio::get_ext(input$profile_file$datapath)
      if(grepl("xls", ext)) data <- rio::import(input$profile_file$datapath, header = TRUE) %>% as_tibble()
      else data <- rio::import(input$profile_file$datapath, header = TRUE) %>% as_tibble()
      ## Remove NA columns
      data <- data %>% select(which(map_lgl(.x = ., .f = ~ !all(is.na(.x)))))
      ndata <- names(data)
      ## Update columns:
      reactive_columns$columns <- ndata
      reactive_columns$locus_guess <- locus_guess(data)
      reactive_columns$genotype_guess <- genotype_guess(data)
      ##
      updateSelectInput(session, inputId = "col_locus", selected = reactive_columns$locus_guess[1])
      updateSelectInput(session, inputId = "col_genotype", selected = reactive_columns$genotype_guess[1])
      data
    })

    output$column_panel <- renderUI({
      sel_locus <- input$col_locus
      sel_genotype <- input$col_genotype
      opt_columns <- reactive_columns$columns

      div(
        selectInput(inputId = ns("col_locus"), label = "Locus column:", choices = opt_columns, multiple = FALSE, selected = sel_locus),
        selectInput(inputId = ns("col_genotype"), label = "Genotype column:", choices = opt_columns, multiple = FALSE, selected = sel_genotype)
      )
    })

    list(
      file = reactive(input$profile_file),
      data = reactive_read_profile,
      locus_col = reactive(input$col_locus),
      genotype_col = reactive(input$col_genotype),
      locus_guess = reactive(reactive_columns$locus_guess),
      genotype_guess = reactive(reactive_columns$genotype_guess)
    )
  })
}
