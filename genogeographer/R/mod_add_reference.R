#' Add reference populations module - UI
#'
#' Upload a dataset of genotypes plus an accompanying information file to build
#' new reference database(s) that can be used for analysis.
#' @param id Module id
#' @noRd
mod_add_reference_ui <- function(id){
  ns <- NS(id)
  nav_panel(
    title = "Add reference populations",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 350,
        h3("Dataset file"),
        fileInput(ns('dataset_file'), 'Choose CSV or Excel File', width = "100%",
                  accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx', '.xls')),
        uiOutput(ns("dataset_panel")),
        h3("Information file"),
        fileInput(ns('info_file'), 'Choose CSV or Excel File', width = "100%",
                  accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx', '.xls')),
        uiOutput(ns("info_panel"))
      ),
      card(
        fill = FALSE,
        card_header("Instructions"),
        card_body(
          fill = FALSE,
          HTML("Download the database and information files to see the required columns and data structure here:"),
          downloadLink(ns("download_reference_db"), "Reference database example"),
          HTML(" and "),
          downloadLink(ns("download_reference_info"), "Information example")
        )
      ),
      uiOutput(ns("upload_data"))
    )
  )
}

#' Add reference populations module - server
#' @param id Module id
#' @noRd
mod_add_reference_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    ## build fixes : start ##
    meta <- NULL; pop <- NULL; . <- NULL
    lat <- NULL; lon <- NULL
    ## build fixes : end ##

    output$download_reference_db <- downloadHandler(
      filename <- function(){
        paste("EUROFORGEN_global_db", "csv", sep = ".")
      },
      contentType = "text/csv",
      content = function(file) {
        db_example <- read_csv(system.file("deployable_app", "EUROFORGEN_global_db.csv", package = ggg_package), show_col_types = FALSE)
        write_csv(db_example, file = file)
      }
    )

    output$download_reference_info <- downloadHandler(
      filename <- function(){
        paste("EUROFORGEN_global_info", "csv", sep = ".")
      },
      contentType = "text/csv",
      content = function(file) {
        info_example <- read_csv(system.file("deployable_app", "EUROFORGEN_global_info.csv", package = ggg_package), show_col_types = FALSE)
        write_csv(info_example, file = file)
      }
    )

    ## UPLOAD DATA FOR REFERENCE x1

    output$dataUploaded <- reactive({
      return(nrow(reactive_read_data()) > 0)
    })
    outputOptions(output, 'dataUploaded', suspendWhenHidden = FALSE)

    output$infoUploaded <- reactive({
      return(nrow(reactive_read_info()) > 0)
    })
    outputOptions(output, 'infoUploaded', suspendWhenHidden = FALSE)

    reactive_read_data <- reactive({
      if(is.null(input$dataset_file)) { # User has not uploaded a file yet
        reactive_dataset$is_excel <- NULL
        return(tibble())
      }
      ext <- rio::get_ext(input$dataset_file$datapath)
      if(grepl("xls", ext)){
        data <- rio::import(input$dataset_file$datapath) %>% as_tibble()
        reactive_dataset$is_excel <- TRUE
      }
      else{
        data <- rio::import(input$dataset_file$datapath, header = TRUE) %>% as_tibble()
        reactive_dataset$is_excel <- FALSE
      }
      ## Remove NA columns
      data <- data %>% select(which(map_lgl(.x = ., .f = ~ !all(is.na(.x)))))
      ndata <- names(data)
      ## Update columns:
      reactive_dataset$columns <- ndata
      rs_count <- sum(grepl(input$col_rs, ndata))
      reactive_dataset$rs <- paste0("Matches ", rs_count, " column", ifelse(rs_count == 1, "", "s"), " in the dataset.")
      ## return read data
      data
    })

    reactive_selected_data <- eventReactive(list(input$dataset_file, input$process_dataset), {
      dat <- reactive_read_data()
      ncol_prior <- ncol(dat)
      dataset_cols <- c("col_sample", "col_pop", "col_meta", "col_rs")
      dataset_cols %>% map_chr(~ isolate(input[[.x]]))
      dataset_inputs <- dataset_cols %>% map_chr(~ input[[.x]])
      if(length(dataset_inputs) == length(unique(dataset_inputs))){
        dat <- dat %>% select(all_of(dataset_inputs[-4]), starts_with(dataset_inputs[4]))
        ncol_post <- ncol(dat)
        reactive_dataset$success <- paste("The table above shows the", ncol_post, "out of",
                                          ncol_prior, "columns that are used in the process.")
      }
      dat
    }, ignoreInit = TRUE)

    reactive_read_info <- reactive({
      if(is.null(input$info_file)) { # User has not uploaded a file yet
        reactive_info$is_excel <- NULL
        return(tibble())
      }
      ext <- rio::get_ext(input$info_file$datapath)
      if(grepl("xls", ext)){
        data <- rio::import(input$info_file$datapath) %>% as_tibble()
        reactive_info$is_excel <- TRUE
      }
      else{
        data <- rio::import(input$info_file$datapath, header = TRUE) %>% as_tibble()
        reactive_info$is_excel <- FALSE
      }
      ## Remove NA columns
      data <- data %>% select(which(map_lgl(.x = ., .f = ~ !all(is.na(.x)))))
      ndata <- names(data)
      ## Update columns:
      reactive_info$columns <- ndata
      ## return read data
      data
    })

    reactive_selected_info <- eventReactive(list(input$info_file, input$process_info), {
      dat <- reactive_read_info()
      ncol_prior <- ncol(dat)
      info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon")
      info_cols %>% map_chr(~ isolate(input[[.x]]))
      info_inputs <- info_cols %>% map_chr(~ input[[.x]])
      if(length(info_inputs) == length(unique(info_inputs))){
        dat <- dat %>% select(all_of(info_inputs))
        ncol_post <- ncol(dat)
        reactive_info$success <- paste("The table above shows the", ncol_post, "out of",
                                        ncol_prior, "columns that are used in the process.")
      }
      dat
    })

    reactive_dataset <- reactiveValues(columns = NULL, rs = NULL, success = NULL, is_excel = NULL)
    reactive_info <- reactiveValues(columns = NULL, outofplace = NULL, success = NULL, is_excel = NULL)
    reactive_x1 <- reactiveValues(db = NULL)

    observeEvent(input$col_POP, {
      reactive_info$outofplace <- unique(reactive_selected_info()[[input$col_POP]])
    })

    observeEvent(input$process_dataset, {
      dataset_cols <- c("col_sample", "col_pop", "col_meta", "col_rs")
      dataset_inputs <- dataset_cols %>% map_chr(~ input[[.x]])
      if(length(unique(dataset_inputs)) == length(dataset_inputs)) dataset_cols %>% purrr::map(~ disable(.x))
    })

    output$dataset_check <- renderUI({
      dataset_cols <- c("col_sample", "col_pop", "col_meta", "col_rs")
      dataset_cols_ <- dataset_cols %>% purrr::map(~sym(.x))
      if(all(map_lgl(dataset_cols_, ~!is.null(.x)))){
        dataset_inputs <- dataset_cols %>% map_chr(~ input[[.x]])
        if(input$process_dataset && (length(unique(dataset_inputs)) != length(dataset_inputs)))
          tags$p("Error: The selected columns must be different", style = "color: red;")
      }
    })

    observeEvent(input$unlock_dataset, {
      dataset_cols <- c("col_sample", "col_pop", "col_meta", "col_rs")
      reactive_dataset$success <- NULL
      dataset_cols %>% purrr::map(~ enable(.x))
    })

    observeEvent(input$process_info, {
      info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon")
      info_inputs <- info_cols %>% map_chr(~ input[[.x]])
      if(length(unique(info_inputs)) == length(info_cols)) c(info_cols, "col_outofplace") %>% purrr::map(~ disable(.x))
    })

    output$info_check <- renderUI({
      info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon")
      info_inputs <- info_cols %>% map_chr(~ input[[.x]])
      if(input$process_info && (length(unique(info_inputs)) != length(info_inputs)))
        tags$p("Error: The selected columns must be different", style = "color: red;")
    })

    observeEvent(input$unlock_info, {
      info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon", "col_outofplace")
      info_cols %>% purrr::map(~ enable(.x))
    })

    output$dataset_panel <- renderUI({
      sel_sample <- input$col_sample
      sel_pop <- input$col_pop
      sel_meta <- input$col_meta
      sel_rs <- if(is.null(input$col_rs)) "rs" else input$col_rs
      data_columns <- reactive_dataset$columns
      data_rs <- reactive_dataset$rs

      div(
        uiOutput(ns("data_excel")),
        helpText("Select columns containing the specified information below:"),
        selectInput(inputId = ns("col_sample"), label = "Sample", choices = data_columns, multiple = FALSE, selected = sel_sample),
        selectInput(inputId = ns("col_pop"), label = "Population identifier", choices = data_columns, multiple = FALSE, selected = sel_pop),
        selectInput(inputId = ns("col_meta"), label = "Metapopulation identifier", choices = data_columns, multiple = FALSE, selected = sel_meta),
        textInput(inputId = ns("col_rs"), label = 'SNP columns indicator (typically "rs")', value = sel_rs, placeholder = "rs"),
        helpText(data_rs),
        uiOutput(ns("dataset_check")),
        div(
          actionButton(inputId = ns("process_dataset"), label = "Process dataset", icon = icon("forward-step", verify_fa = FALSE)),
          actionButton(inputId = ns("unlock_dataset"), label = "Update selections", icon = icon("backward-step", verify_fa = FALSE))
        )
      )
    })

    output$info_panel <- renderUI({
      info_columns <- reactive_info$columns
      info_outofplace <- reactive_info$outofplace

      sel_POP <- input$col_POP
      sel_population <- input$col_population
      sel_META <- input$col_META
      sel_metapopulation <- input$col_metapopulation
      sel_lat <- input$col_lat
      sel_lon <- input$col_lon
      sel_outofplace <- input$col_outofplace

      div(
        helpText("Select columns containing the specified information below:"),
        h5("Population"),
        selectInput(inputId = ns("col_POP"), label = "Identifier (typically shorter names/abbreviations)", choices = info_columns, multiple = FALSE, selected = sel_POP),
        selectInput(inputId = ns("col_population"), label = "Name/description", choices = info_columns, multiple = FALSE, selected = sel_population),
        h5("Metapopulation"),
        selectInput(inputId = ns("col_META"), label = "Identifier (typically shorter names/abbreviations)", choices = info_columns, multiple = FALSE, selected = sel_META),
        selectInput(inputId = ns("col_metapopulation"), label = "Name/description", choices = info_columns, multiple = FALSE, selected = sel_metapopulation),
        h5("Location"),
        selectInput(inputId = ns("col_lat"), label = "Latitude", choices = info_columns, multiple = FALSE, selected = sel_lat),
        selectInput(inputId = ns("col_lon"), label = "Longitude", choices = info_columns, multiple = FALSE, selected = sel_lon),
        helpText(paste("If some of the populations are located 'out of region', e.g. Central Europeans samples in the US,
                       specify them below for prettier map plots. The options below are based on the population identifier set above.",
                       ifelse(is.null(sel_POP) || sel_POP == "", "", paste0('Currently given by column "', sel_POP, '"')))),
        selectInput(inputId = ns("col_outofplace"), label = "Out of region population identifiers", choices = info_outofplace, multiple = TRUE, selected = sel_outofplace),
        uiOutput(ns("info_check")),
        div(
          actionButton(inputId = ns("process_info"), label = "Process information", icon = icon("step-forward", verify_fa = FALSE)),
          actionButton(inputId = ns("unlock_info"), label = "Update selections", icon = icon("step-backward", verify_fa = FALSE))
        )
      )
    })

    reactive_info_rows <- reactive({
      data_ <- reactive_selected_data() %>%
        rename(sample = input$col_sample, pop = input$col_pop, meta = input$col_meta)
      info_ <- reactive_selected_info() %>%
        rename(pop = input$col_POP, population = input$col_population,
               meta = input$col_META, metapopulation = input$col_metapopulation,
               lat = input$col_lat, lon = input$col_lon)
      info_ %>% select(pop, population, meta, metapopulation) %>%
        semi_join(data_, by = c("pop", "meta"))
    })

    output$selected_populations <- renderDT({
      dt_table(reactive_info_rows())
    })

    output$info_data <- renderUI({
      if(!is.null(reactive_dataset$success) & !is.null(reactive_info$success)){
        tagList(
          card(
            fill = FALSE,
            card_header("Populations and metapopulations selected from the data"),
            card_body(fill = FALSE, DTOutput(ns("selected_populations"), fill = FALSE))
          ),
          card(
            fill = FALSE,
            card_header("Compute"),
            card_body(
              fill = FALSE,
              textInput(inputId = ns("db_name"), label = "Specify database name (will appear in app)", value = sub("\\..*$", "", input$dataset_file$name)),
              div(tags$b("Click on the button below to execute the calculations... "), tags$text("(May be slow)")),
              div(
                actionButton(inputId = ns("comp_x1"), label = "Compute databases", icon = icon("database", verify_fa = FALSE)),
                shinyjs::disabled(downloadButton(outputId = ns("download_x1"), label = "Download created databases", class = "btn-primary"))
              )
            )
          )
        )
      }
      else NULL
    })

    observeEvent(input$comp_x1, {
      data_ <- reactive_selected_data() %>%
        rename(sample = isolate(input$col_sample), pop = isolate(input$col_pop), meta = isolate(input$col_meta))
      n_pop <- length(unique(data_$pop))
      n_meta <- length(unique(data_$meta))
      info_ <- reactive_selected_info() %>%
        rename(pop = isolate(input$col_POP), population = isolate(input$col_population),
               meta = isolate(input$col_META), metapopulation = isolate(input$col_metapopulation),
               lat = isolate(input$col_lat), lon = isolate(input$col_lon))
      n_comp <- c(0, cumsum(c(n_pop, n_meta, choose(n_pop, 2), choose(n_meta, 2))))
      shiny_progress <- list(n = n_comp[length(n_comp)], session = session)
      withProgressWaitress({
        pop_DB <- make_x1(df = data_, latlon = info_, groups = "pop", exclude = c("sample", "meta"),
                          allele_list = ggg_allele_list, shiny = c(shiny_progress, list(start = n_comp[1])))
        meta_DB <- make_x1(df = data_, latlon = info_, groups = "meta", exclude = c("sample", "pop"),
                           allele_list = ggg_allele_list, out_of_place = isolate(input$col_outofplace), shiny = c(shiny_progress, list(start = n_comp[2])))
        pop_DB_1admix <- admix_dbs(pop_DB, shiny = c(shiny_progress, list(start = n_comp[3])))
        meta_DB_1admix <- admix_dbs(meta_DB, shiny = c(shiny_progress, list(start = n_comp[4])))
      }, selector = paste0("#", ns("download_x1")), max = 100, theme = "overlay-percent")
      db_list_new <- list(list(pop = list(db = pop_DB, admix = pop_DB_1admix),
                          meta = list(db = meta_DB, admix = meta_DB_1admix)))
      if(!is.null(input$db_name) || input$db_name != "") names(db_list_new) <- input$db_name
      else names(db_list_new) <- sub("\\..*$", "", isolate(input$dataset_file$name))
      reactive_x1$db <- db_list_new
      shinyjs::enable("download_x1")
    })

    output$download_x1 <- downloadHandler(
      filename = function() {
        fil <- if(!is.null(input$db_name) || input$db_name != "") input$db_name else sub("\\..*$", "", input$dataset_file$name)
        gsub("\\s", "_", paste0(fil, ".rds"))
      },
      content = function(file) {
        saveRDS(reactive_x1$db, file)
      }
    )

    output$data_preview <- renderDT({
      data_ <- reactive_selected_data()
      rs_cols <- if(is.null(input$col_rs)) FALSE else grepl(input$col_rs, names(data_))
      data_ %>% head(n = min(nrow(data_), 10)) %>% dt_simple() %>%
        formatStyle(columns = rs_cols, backgroundColor = "#ecf1f2")
    })

    output$info_preview <- renderDT({
      info_ <- reactive_selected_info()
      info_ %>% head(n = min(nrow(info_), 10)) %>% dt_simple()
    })

    output$upload_data <- renderUI({
      data_n <- nrow(reactive_selected_data())
      n_data <- if(data_n > 10) 10 else data_n
      info_n <- nrow(reactive_selected_info())
      n_info <- if(info_n > 10) 10 else info_n
      tagList(
        card(
          fill = FALSE,
          card_header("Uploaded dataset file"),
          card_body(
            fill = FALSE,
            tags$p(paste0('The first ', n_data, ' out of ', data_n, ' rows  the dataset file are shown below.
                          Columns matching SNP column indicator ("', input$col_rs, '") are highlighted.')),
            DTOutput(ns("data_preview"), fill = FALSE),
            helpText(reactive_dataset$success)
          )
        ),
        card(
          fill = FALSE,
          card_header("Uploaded information file"),
          card_body(
            fill = FALSE,
            tags$p(paste0('The first ', n_info, ' out of ', info_n, ' rows of the information file are shown below.')),
            DTOutput(ns("info_preview"), fill = FALSE),
            helpText(reactive_info$success)
          )
        ),
        uiOutput(ns("info_data"))
      )
    })
  })
}
