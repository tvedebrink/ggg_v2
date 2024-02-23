dt_simple <- function(tab, ...){
  tab %>%
    DT::datatable(...,
                  rownames = FALSE,
                  options = list(
                    scrollX = TRUE,
                    dom = 't'
                  ), class = 'white-space: nowrap')
}

dt_table <- function(x, ...){
  x %>%
    DT::datatable(data = ., ...,
                  rownames=FALSE, filter = "bottom", selection = 'none',
                  extensions = 'Buttons',
                  options = list(
                    dom = 'Blfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                    lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All"))
                  )
    )
}

#' @title Server part for Shiny app
#' @description server function to be used in the `app_genogeo` shiny app
#' @param input The list of inputs (leave empty)
#' @param output The list of outputs (leave empty)
#' @param session The shiny session - to refer to session elements (leave empty)
#' @export
server_api <- function(input, output, session){
  ## build fixes : start ##
  Target.ID <- NULL; Genotype <- NULL; locus <- NULL; genotype <- NULL
  accept <- NULL; p_value <- NULL; selected_ <- NULL
  meta <- NULL; pop <- NULL; . <- NULL
  lat <- NULL; lon <- NULL; aims_example <- NULL
  ## build fixes : end ##
  if(!exists("db_list")) db_list <- get("db_list", envir = -2)
  if(is.null(db_list)) disable("analyse")
  reactive_db_list <- reactiveValues(db = db_list)
  if(!exists("reporting_panel")) reporting_panel <- get("reporting_panel", envir = -2)

  observeEvent(input$analyse, {
    # When the button is clicked, wrap the code in a call to `withBusyIndicatorServer()`
    withBusyIndicatorServer("analyse", reactive_result())
  })

  output$pop_info <- renderUI({
    print(names(reactive_db_list$db))
    pop_info <- attr(reactive_db_list$db[[1]]$pop$db, "info")
    print(dim(pop_info))
    meta_info <- attr(reactive_db_list$db[[1]]$meta$db, "info")
    print(dim(meta_info))
    fluidPage(
    fluidRow(h3("Populations")),
    fluidRow(pop_info %>% select(population, n) %>% dt_table()),
    fluidRow(h3("Metapopulations")),
    fluidRow(meta_info %>% select(metapopulation, n) %>% dt_table())
    )
  })

  output$allele_info <- renderDT({
  ggg_allele_list %>%
    mutate(locus = paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/",locus,"' target='_blank'>", locus, "</a>")) %>%
    rename('AIMs SNP' = locus, 'Reference allele' = x1, 'Alternative allele' = x2) %>%
    mutate(across(Kidd:Seldin, ~ifelse(.x, "Yes", "No"))) %>%
    dt_table(escape = FALSE)
  })

  ## USER INTERFACE
  output$analysis <- renderUI({
    res <- reactive_result()
    dat <- reactive_read_profile()
    if(is.null(res)){
      if(is.null(dat) || nrow(dat) == 0){
        fluidPage(
          h3("Upload valid AIMs profile"),
          HTML("<p>Upload valid file: It must contain a column containing marker (locus) names, and a column with the genotypes as e.g. <tt>AA</tt> or <tt>AC</tt></p>"),
          HTML("You can download a sample file showing the necessary columns and data structure here:"),
          downloadLink("download_sample", "Download")
        )
      }
      else{
        ## Make column suggestions
        guess_locus_text <- and_text(pre = "<p>Suggestions for <b>locus column</b>: ", x = reactive_columns$locus_guess, post = "</p>")
        guess_genotype_text <- and_text(pre = "<p>Suggestions for <b>genotype column</b>: ", x = reactive_columns$genotype_guess, post = "</p>")
        #
        fluidPage(
          fluidRow(h3("Uploaded file")),
          fluidRow(
            HTML("<p>The first few rows of the uploaded file are shown below</p>"),
            HTML(guess_locus_text),
            HTML(guess_genotype_text)
          ),
          fluidRow(
            renderDT({
              datatable(dat, options(dom = "t"))
            })
          )
        )
      }
    }
    else{
      fluidPage(
        tags$head(tags$style(paste0(".modal-lg{ width: ", 2*session$clientData$output_barplot_width,"px}"))),
        fluidRow(HTML(paste0("<p><b>Analysis of file:</b> ",input$profile_file$name,":")),
                 actionLink("show_profile", label = "Show profile"),
                 icon("new-window", lib = "glyphicon", verify_fa = FALSE),
                 HTML("</p>")
                 ),
        fluidRow(tags$h2("Graphics")),
        fluidRow(
          column(width = 6, uiOutput("barplot_panel") ),
          column(width = 6, uiOutput("map_panel") )
        ),
        fluidRow(tags$h2("Tables")),
        fluidRow(DTOutput("result_table")),
        fluidRow(tags$h2("Likelihood ratios")),
        fluidRow(textOutput("LR_one_pop")),
        fluidRow(
          column(width = 6, DTOutput("lr_list")),
          column(width = 6, uiOutput("LRplot_panel") )
        )
      )
    }
  })

  ##

  output$download_sample <- downloadHandler(
    filename <- function(){
      paste("aims_example", "csv", sep = ".")
    },
    contentType = "text/csv",
    content = function(file) {
      aims_example <- read_csv(system.file("deployable_app", "aims_example.csv", package = ggg_package), col_types = "cc")
      write_csv(aims_example, file = file)
    }
  )

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

  ## REACTIVES

  observeEvent(input$show_profile, {
    A1 <- NULL
    A2 <- NULL
    raw_profile <- reactive_read_profile() %>% rename(locus = !!sym(input$col_locus))
    profile <- reactive_profile()
    profile_drop <- raw_profile %>% anti_join(profile, by = "locus") %>%
      rename(!!sym(input$col_locus) := locus)
    profile_x0 <- raw_profile %>% semi_join(profile, by = "locus")
    showModal(modalDialog(
      title = paste("Uploaded profile:",input$profile_file$name),
      size = "m",
      h3("Analysed loci"),
      helpText("The loci below has been included in the analysis."),
      renderDT(profile_x0 %>% dt_table),
      h3("Dropped or unused loci"),
      helpText(paste0("The loci below has been excluded from the analysis.\n
               Either because of state 'NN', locus not in '",input$snp_set,"'
                      or other typing error (e.g. different reference allele).")),
      renderDT(profile_drop %>% dt_table),
      footer = modalButton("Close"),
      easyClose = TRUE
      ))
    })

  ### ANALYSIS

  observeEvent(input$reset,{
    runjs("history.go(0)")
  })

  output$uploaded_profile <- renderTable({
    if (is.null(input$profile_file)) return(NULL)
    reactive_profile()
  })

  reactive_columns <- reactiveValues(columns = NULL, locus_guess = NULL, genotype_guess = NULL)

  output$column_panel <- renderUI({
    sel_locus <- input$col_locus
    sel_genotype <- input$col_genotype
    opt_columns <- reactive_columns$columns

    verticalLayout(
      selectInput(inputId = "col_locus", label = "Locus column:", choices = opt_columns, multiple = FALSE, selected = sel_locus),
      selectInput(inputId = "col_genotype", label = "Genotype column:", choices = opt_columns, multiple = FALSE, selected = sel_genotype)
    )
  })

  observeEvent(input$db_add, {
    req(input$db_add)
    user_db <- readRDS(file = input$db_add$datapath)
    db_names <- user_db %>% purrr::map(names) %>% unlist() %>% unname() %>% unique()
    if(length(setdiff(db_names, c("pop", "meta"))) == 0){
      reactive_db_list$db <- c(reactive_db_list$db, user_db)
      updateSelectInput(session, inputId = "snp_set", selected = names(user_db))
    }
    if(!is.null(reactive_db_list$db)) enable("analyse")
    })

  output$dbs <- renderUI({
    verticalLayout(
      fileInput(inputId = "db_add", "Upload own reference database:", accept = ".rds", width = "100%"),
      helpText("Ensure that the '.rds'-file is created using genogeographer, e.g. by following the steps under tab 'Add referece population'"),
      selectInput(inputId = "snp_set", label = "Select reference database:", choices = names(reactive_db_list$db))
    )
  })

  output$fileUploaded <- reactive({
    return(nrow(reactive_read_profile())>0)
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

  reactive_read_profile <- reactive({
    if(is.null(input$profile_file)) { # User has not uploaded a file yet
      return(NULL)
    }
    ## Import profile
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
    ## return read data
    # print("reactive_read_profile")
    # browser()
    data
  })

  reactive_profile <- reactive({
    profile <- reactive_read_profile()
    # req(profile)
    if(is.null(profile) || nrow(profile) == 0) return(NULL) ## No profile
    nprofile <- names(profile)
    col_return <- FALSE
    # req(input$col_locus)
    # req(input$col_genotype)
    if(!(input$col_locus %in% nprofile)){
      updateSelectInput(session, inputId = "col_locus", selected = NULL)
      col_return <- TRUE
      }
    if(!(input$col_genotype %in% nprofile)){
      updateSelectInput(session, inputId = "col_genotype", selected = NULL)
      col_return <- TRUE
    }
    if(col_return || (input$col_locus == input$col_genotype)) return(NULL)
    # browser()
    ##
    ggg_db <- attr(reactive_db_list$db[[input$snp_set]]$meta$db, "allele_list")
    ##
    genotype_x0(profile = profile, locus = input$col_locus, genotype = input$col_genotype, ggg = ggg_db)
  })

  reactive_result <- eventReactive(list(input$profile_file,input$analyse),{
    profile <- reactive_profile()
    if(is.null(profile)) return(NULL)
    ### TILTING
    tilt_control <- if(is.null(input$tilt)) FALSE else (input$tilt == "adjust")
    ## ADMIXTURE
    admix_control <- if(is.null(input$admix)) FALSE else (input$admix == "admix")
    ## COMPUTE
    res <- ggg_score(profile_x0 = profile, DB = reactive_db_list$db[[input$snp_set]][[input$meta]]$db,
                     min_n = input$min_n, CI = input$CI/100, tilt = tilt_control)
    if(admix_control){
      res_admix <- ggg_score(profile_x0 = profile, DB = reactive_db_list$db[[input$snp_set]][[input$meta]]$admix,
                             min_n = input$min_n, CI = input$CI/100, tilt = FALSE)
      info <- bind_rows(attr(res, "info"), attr(res_admix, "info"))
      res <- bind_rows(res, res_admix) %>% arrange(desc(logP))
      attr(res, "info") <- info
    }
    res
  })

  output$side_pvalue <- renderUI({
    res <- reactive_result()
    if(is.null(res)) return(NULL)
    input_meta <- isolate(input$meta)
    groups <- if(is.null(input_meta)) "meta" else input_meta
    groups_ <- sym(groups)
    grouping <- if(groups == "meta") "metapopulation" else "population"
    grouping_ <- sym(grouping)
    db_info <- attr(res, "info")
    if(!is.null(db_info)){
      key_name <- db_info %>% select(!!groups_, !!grouping_) %>% deframe()
      res <- res %>% mutate(!!groups_ := key_name[!!groups_])
    }
    largest_p_value <- res %>% filter(accept)
    if(nrow(largest_p_value)==0) largest_p_value_text <- paste0("All ", grouping, "s are rejected")
    else{
      largest_p_value <- largest_p_value %>% slice_max(n = 1, order_by = p_value) %>%
        select(p_value, !!groups_) %>%
        mutate(p_value = round(p_value, 3)) %>%
        unite(p_text, p_value, !!groups_, sep  = " (")
      largest_p_value_text <- paste0("<b>DB-score (largest p-value):</b><br/>", largest_p_value, ")")
    }
    verticalLayout(
      h4("Reporting"),
      HTML(paste0("<p>",largest_p_value_text,"</p>"))
    )
  })

  ### MAP

  plot_map <- reactive({
    res <- reactive_result()
    if(is.null(res)) return(NULL)
    leaflet_plot(z_df = res)
  })

  output$map <- renderLeaflet({ plot_map() })

  output$map_panel <- renderUI({
    # div(style = paste0("position:relative; ",
    #                    "height: ",0.8*session$clientData$output_barplot_width,"px;"),
        withSpinner(leafletOutput("map", width = "100%"), type = 4) #)
  })

  ### BARPLOT

  output$barplot_panel <- renderUI({
    withSpinner(plotlyOutput("barplot", width = "100%"), type = 4)
  })

  output$barplot <- renderPlotly({
    res <- reactive_result()
    if(is.null(res)) return(NULL)
    ebp <- error_bar_plotly(result_df = res, which = isolate(input$result_table_rows_selected))
    height <- session$clientData$output_barplot_height
    width <- session$clientData$output_barplot_width
    ebp %>% ggplotly(tooltip = c("text")) %>% # , height = width*0.8, width = width) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = FALSE)
    })

  ## TABLES
  # https://stackoverflow.com/questions/53768488/how-to-display-html-in-dt-header

  z_table <- reactive({
    res <- reactive_result()
    if(is.null(res)) return(NULL)
    db_info <- attr(res, "info")
    if(!is.null(db_info)){
      input_meta <- isolate(input$meta)
      groups <- if(is.null(input_meta)) "meta" else input_meta
      groups_ <- sym(groups)
      grouping <- if(groups == "meta") "metapopulation" else "population"
      grouping_ <- sym(grouping)
      key_name <- db_info %>% select(!!groups_, !!grouping_) %>% deframe()
      res <- res %>% mutate(!!groups_ := key_name[!!groups_]) %>%
        rename(!!grouping_ := !!groups_)
    }
    result_table(res)
  })

  output$result_table <- renderDT({
    z_table()
  })

  ## LR calculations and controls

  output$LR_select <- renderUI({
    res <- reactive_result()
    if(is.null(res)) return(NULL)
    input_meta <- isolate(input$meta)
    groups <- if(is.null(input_meta)) "meta" else input_meta
    groups_ <- sym(groups)
    grouping <- if(groups == "meta") "metapopulation" else "population"
    grouping_ <- sym(grouping)
    accepted <- res %>% filter(accept) %>% pull(!!groups_)
    LR_choices <- res %>% pull(!!groups_) %>% paste()
    db_info <- attr(res, "info")
    if(!is.null(db_info)){
      key_name <- db_info %>% select(!!groups_, !!grouping_) %>% deframe()
      LR_choices <- LR_choices %>% set_names(key_name[.])
    }
    verticalLayout(
      checkboxGroupInput(inputId = "LR_accept",
                         label = "LRs with two rejected populations:",
                         choiceNames = list("Allow (LRs may be misleading)"), choiceValues = list("allow")),
      selectizeInput(inputId = "LR_selected",
                     label =  "Include in LR calculations:",
                     choices = LR_choices,
                     selected = accepted,
                     multiple = TRUE,
                     options = list(plugins = list("remove_button", "drag_drop")),
                     width = "100%")
      )
  })

  output$LR_one_pop <- renderText({
    res <- reactive_result()
    input_meta <- isolate(input$meta)
    mm <- if(is.null(input_meta)) "meta" else input_meta
    MM <- if(mm == "meta") "metapopulation" else "population"
    if(nrow(res) <= 1) paste("Likelihood ratios are not computed, since there is just a single", MM, "in the results")
    else return(NULL)
  })

  reactive_LR_table <- reactive({
    accepted_control <- if(is.null(input$LR_accept)) FALSE else (input$LR_accept == "allow")
    res <- reactive_result()
    if(nrow(res) <= 1) return(NULL)
    list(
      res = res,
      LR_Tab = LR_table(z_df = res, who = input$LR_selected, CI = input$CI/100, only_accepted = !accepted_control)$LR
      )
  })

  output$lr_list <- renderDT({
    LR_Tab <- reactive_LR_table()
    req(LR_Tab)
    LR_list(result = LR_Tab$res, LR_tab = LR_Tab$LR_Tab)
  })

  observeEvent(input$LR_selected,{ ## Can't delete most probable (based on z_score)
    res <- reactive_result()
    if(is.null(res)) return(NULL)
    groups <- names(res)[1]
    groups_ <- sym(groups)
    max_score_pop <- res %>% filter(accept) %>% slice_max(n = 1, order_by = p_value) %>% pull(!!groups_) ## Highest p-value
    LR_selected <- unique(c(input$LR_selected, max_score_pop))
    updateSelectizeInput(session, "LR_selected", selected = LR_selected)
  })

  ## LR PLOTS

  output$LRplot <- renderPlotly({
    LR_Tab <- reactive_LR_table()
    req(LR_Tab)
    if(is.null(LR_Tab$LR_Tab) | nrow(LR_Tab$LR_Tab) == 0) return(NULL)
    LR_plot_ly(result = LR_Tab$res, LR_list = LR_Tab$LR_Tab)
    })

  output$LRplot_panel <- renderUI({
        withSpinner(plotlyOutput("LRplot", width = "100%"), type = 4)
  })

  ## UPLOAD DATA FOR REFERENCE x1

  output$dataUploaded <- reactive({
    return(nrow(reactive_read_data())>0)
  })
  outputOptions(output, 'dataUploaded', suspendWhenHidden=FALSE)

  output$infoUploaded <- reactive({
    return(nrow(reactive_read_info())>0)
  })
  outputOptions(output, 'infoUploaded', suspendWhenHidden=FALSE)

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

  reactive_selected_data <- eventReactive(list(input$dataset_file, input$process_dataset),{
    dat <- reactive_read_data()
    ncol_prior <- ncol(dat)
    dataset_cols <- c("col_sample", "col_pop", "col_meta", "col_rs")
    ## dataset_cols %>% map_chr(~ req(input[[.x]]))
    dataset_cols %>% map_chr(~ isolate(input[[.x]]))
    dataset_inputs <- dataset_cols %>% map_chr(~ input[[.x]])
    if(length(dataset_inputs) == length(unique(dataset_inputs))){
      dat <- dat %>% select(all_of(dataset_inputs[-4]), starts_with(dataset_inputs[4]))
      ncol_post <- ncol(dat)
      reactive_dataset$success <- paste("The table above shows the", ncol_post, "out of",
                                         ncol_prior, "columns that are used in the process.")
    }
    dat
  }, ignoreInit = TRUE) ##

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

  reactive_selected_info <- eventReactive(list(input$info_file, input$process_info),{
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
  }) #, ignoreInit = TRUE) ##

  reactive_dataset <- reactiveValues(columns = NULL, rs = NULL, success = NULL, is_excel = NULL)
  reactive_info <- reactiveValues(columns = NULL, outofplace = NULL, success = NULL, is_excel = NULL)
  reactive_x1 <- reactiveValues(db = NULL)

  observeEvent(input$col_POP,{
    reactive_info$outofplace <- unique(reactive_selected_info()[[input$col_POP]])
    })

  #

  observeEvent(input$process_dataset,{
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

  observeEvent(input$unlock_dataset,{
    dataset_cols <- c("col_sample", "col_pop", "col_meta", "col_rs")
    reactive_dataset$success <- NULL
    dataset_cols %>% purrr::map(~ enable(.x))
  })

  ##

  observeEvent(input$process_info,{
    info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon")
    info_inputs <- info_cols %>% map_chr(~ input[[.x]])
    if(length(unique(info_inputs)) == length(info_inputs)) info_cols %>% purrr::map(~ disable(.x))
  })

  output$info_check <- renderUI({
    info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon")
    info_cols_ <- info_cols %>% purrr::map(~sym(.x))
    if(all(map_lgl(info_cols_, ~!is.null(.x)))){
      info_inputs <- info_cols %>% map_chr(~ input[[.x]])
      if(input$process_info && (length(unique(info_inputs)) != length(info_inputs)))
        tags$p("Error: The selected columns must be different", style = "color: red;")
    }
  })

  observeEvent(input$unlock_info,{
    info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon")
    reactive_info$success <- NULL
    info_cols %>% purrr::map(~ enable(.x))
  })

  output$dataset_panel <- renderUI({
    ##
    sel_sample <- input$col_sample
    sel_pop <- input$col_pop
    sel_meta <- input$col_meta
    sel_rs <- if(is.null(input$col_rs)) "rs" else input$col_rs
    data_columns <- reactive_dataset$columns
    data_rs <- reactive_dataset$rs

    verticalLayout(
      uiOutput("data_excel"),
      helpText("Select columns containing the specified information below:"),
      selectInput(inputId = "col_sample", label = "Sample", choices = data_columns, multiple = FALSE, selected = sel_sample),
      selectInput(inputId = "col_pop", label = "Population identifier", choices = data_columns, multiple = FALSE, selected = sel_pop),
      selectInput(inputId = "col_meta", label = "Metapopulation identifier", choices = data_columns, multiple = FALSE, selected = sel_meta),
      textInput(inputId = "col_rs", label = 'SNP columns indicator (typically "rs")', value = sel_rs, placeholder = "rs"),
      helpText(data_rs),
      uiOutput("dataset_check"),
      div(
        (actionButton(inputId = "process_dataset", label = "Process dataset", icon = icon("forward-step", verify_fa = FALSE))),
        (actionButton(inputId = "unlock_dataset", label = "Update selections", icon = icon("backward-step", verify_fa = FALSE)))
      )
    )
  })

  observeEvent(input$process_info,{
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

  observeEvent(input$unlock_info,{
    info_cols <- c("col_POP", "col_population", "col_META", "col_metapopulation", "col_lat", "col_lon", "col_outofplace")
    info_cols %>% purrr::map(~ enable(.x))
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

    verticalLayout(
      helpText("Select columns containing the specified information below:"),
      wellPanel(style = "border-color:#c9d7e8;", h5("Population"),
      selectInput(inputId = "col_POP", label = "Identifier (typically shorter names/abbreviations)", choices = info_columns, multiple = FALSE, selected = sel_POP),
      selectInput(inputId = "col_population", label = "Name/description", choices = info_columns, multiple = FALSE, selected = sel_population),
      ),wellPanel(style = "border-color:#c9d7e8;", h5("Metapopulation"),
      selectInput(inputId = "col_META", label = "Identifier (typically shorter names/abbreviations)", choices = info_columns, multiple = FALSE, selected = sel_META),
      selectInput(inputId = "col_metapopulation", label = "Name/description", choices = info_columns, multiple = FALSE, selected = sel_metapopulation),
      ),wellPanel(style = "border-color:#c9d7e8;", h5("Location"),
      selectInput(inputId = "col_lat", label = "Latitude", choices = info_columns, multiple = FALSE, selected = sel_lat),
      selectInput(inputId = "col_lon", label = "Longitude", choices = info_columns, multiple = FALSE, selected = sel_lon),
      helpText(paste("If some of the populations are located 'out of region', e.g. Central Europeans samples in the US,
                     specify them below for prettier map plots. The options below are based on the population identifier set above.",
                     ifelse(is.null(sel_POP) || sel_POP == "", "", paste0('Currently given by column "', sel_POP, '"')))),
      selectInput(inputId = "col_outofplace", label = "Out of region population identifiers", choices = info_outofplace, multiple = TRUE, selected = sel_outofplace),
      ),
      uiOutput("info_check"),
      div(
        actionButton(inputId = "process_info", label = "Process information", icon = icon("step-forward", verify_fa = FALSE)),
        actionButton(inputId = "unlock_info", label = "Update selections", icon = icon("step-backward", verify_fa = FALSE))
      )    )
  })

  output$info_data <- renderUI({
    if(!is.null(reactive_dataset$success) & !is.null(reactive_info$success)){
      data_ <- reactive_selected_data() %>%
        rename(sample = input$col_sample, pop = input$col_pop, meta = input$col_meta)
      info_ <- reactive_selected_info() %>%
        rename(pop = input$col_POP, population = input$col_population,
               meta = input$col_META, metapopulation = input$col_metapopulation,
               lat = input$col_lat, lon = input$col_lon)
      info_rows <- info_ %>% select(pop, population, meta, metapopulation) %>%
        semi_join(data_, by = c("pop", "meta"))
      # n_pop <- length(unique(data_$pop))
      # n_meta <- length(unique(data_$meta))
      # n_comp <- n_pop + n_meta + choose(n_pop, 2) + choose(n_meta, 2)
      ## browser()
      verticalLayout(
        h3("Populations and metapopulations selected from the data"),
        dt_table(info_rows),
        h3("Compute"),
        textInput(inputId = "db_name", label = "Specify database name (will appear in app)", value = sub("\\..*$", "", input$dataset_file$name)),
        div(tags$b("Click on the button below to execute the calculations... "), tags$text("(May be slow)")),
        div(
          actionButton(inputId = "comp_x1", label = "Compute databases", icon = icon("database", verify_fa = FALSE)),
          shinyjs::disabled(downloadButton(outputId = "download_x1", label = "Download created databases", class = "btn-primary"))
        ),
        tags$p("")
      )
    }
    else tags$p("")
  })

  observeEvent(input$comp_x1, {
  ## reactive_db_x1 <- reactive({
    #req(input$comp_x1)
    data_ <- reactive_selected_data() %>%
      rename(sample = isolate(input$col_sample), pop = isolate(input$col_pop), meta = isolate(input$col_meta))
    n_pop <- length(unique(data_$pop))
    n_meta <- length(unique(data_$meta))
    info_ <- reactive_selected_info() %>%
      rename(pop = isolate(input$col_POP), population = isolate(input$col_population),
             meta = isolate(input$col_META), metapopulation = isolate(input$col_metapopulation),
             lat = isolate(input$col_lat), lon = isolate(input$col_lon))
    n_comp <- c(0,cumsum(c(n_pop, n_meta, choose(n_pop, 2), choose(n_meta, 2))))
    # print(n_comp)
    shiny_progress <- list(n = n_comp[length(n_comp)], session = session)
    # print(shiny_progress$step)
    ##
    ## browser()
    withProgressWaitress({
    pop_DB <- make_x1(df = data_, latlon = info_, groups = "pop", exclude = c("sample", "meta"),
                      allele_list = ggg_allele_list, shiny = c(shiny_progress, list(start = n_comp[1])))
    meta_DB <- make_x1(df = data_, latlon = info_, groups = "meta", exclude = c("sample", "pop"),
                       allele_list = ggg_allele_list, out_of_place = isolate(input$col_outofplace), shiny = c(shiny_progress, list(start = n_comp[2])))
    pop_DB_1admix <- admix_dbs(pop_DB, shiny = c(shiny_progress, list(start = n_comp[3])))
    meta_DB_1admix <- admix_dbs(meta_DB, shiny = c(shiny_progress, list(start = n_comp[4])))
    }, selector = "#download_x1", max = 100, theme = "overlay-percent")
    db_list <- list(list(pop = list(db = pop_DB, admix = pop_DB_1admix),
                    meta = list(db = meta_DB, admix = meta_DB_1admix)))
    if(!is.null(input$db_name) || input$db_name != "") names(db_list) <- input$db_name
    else names(db_list) <- sub("\\..*$", "", isolate(input$dataset_file$name))
    reactive_x1$db <- db_list
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

  output$upload_data <- renderUI({
    data_ <- reactive_selected_data() ## reactive_read_data()
    data_n <- nrow(data_)
    n_data <- if(data_n > 10) 10 else data_n
    info_ <- reactive_selected_info() ## reactive_read_info() ##
    info_n <- nrow(info_)
    n_info <- if(info_n > 10) 10 else info_n
    ##
    rs_cols <- if(is.null(input$col_rs)) FALSE else grepl(input$col_rs, names(data_))
    ## Make column suggestions
    fluidPage(
      h3("Instructions"),
      HTML("Download the database and information files to see the required columns and data structure here:"),
      downloadLink("download_reference_db", "Reference database example"),
      HTML(" and "),
      downloadLink("download_reference_info", "Information example"),
      ## Which columns are used for joins
      h3("Uploaded dataset file"),
      tags$p(paste0('The first ', n_data ,' out of ', data_n , ' rows  the dataset file are shown below.
                    Columns matching SNP column indicator ("', input$col_rs, '") are highlighted.')),
      renderDT({
        data_ %>% head(n = n_data) %>% dt_simple() %>%
          formatStyle(columns = rs_cols, backgroundColor = "#ecf1f2")
        }),
      helpText(reactive_dataset$success),
      h3("Uploaded information file"),
      tags$p(paste0('The first ', n_info ,' out of ', info_n , ' rows of the information file are shown below.')),
      renderDT({
        info_ %>% head(n = n_info) %>% dt_simple() }),
      helpText(reactive_info$success),
      uiOutput("info_data")
    )
  })

  ### RETURN pdf REPORT
  #
  output$report_panel <- renderUI({
    res <- reactive_result()
    if(!reporting_panel) return(verticalLayout())
    if(is.null(res)) return(verticalLayout())
    verticalLayout(
      h4("Report"),
      textInput(inputId = "name", label = "Name of analyst", width = "100%",
                placeholder = "Name as to appear in report", value = input$name),
      # radioButtons(inputId = 'format', label = 'Report format', choices = c('PDF', 'HTML', 'Word'), inline = TRUE, selected = input$format),
      # withBusyIndicatorUI(downloadButton(outputId = "report_download", label = paste0("Download (",input$format,")"), class = "btn-primary")),
      withBusyIndicatorUI(downloadButton(outputId = "report_download", label = "Download HTML report", class = "btn-primary")),
      hr()
    )
  })
  #

  output$report_download = downloadHandler(
    filename = function(){
      paste(sub("\\.[[:alnum:]]*.$","",input$profile_file$name),
            "html", #switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'), #
            sep ="."
            )
    },
    content = function(file) {
      src <- list(rmd_file = normalizePath(system.file("deployable_app", "aims_report.Rmd", package = ggg_package)))
      owd <- setwd(tempdir())
      cat(file = stderr(), owd, "\n")
      on.exit(setwd(owd))
      file.copy(src$rmd_file, 'aims_report.Rmd', overwrite = TRUE)
      out <- rmarkdown::render(input = 'aims_report.Rmd', clean = TRUE,
                               output_format = rmarkdown::html_document(), #switch(input$format,PDF = rmarkdown::pdf_document(), HTML = rmarkdown::html_document(), Word = rmarkdown::word_document()), #
        params = list(
          set_file = input$profile_file$name,
          set_author = input$name,
          set_output = input$format,
          set_locus = input$col_locus,
          set_genotype = input$col_genotype,
          set_db = input$snp_set
          )
        )
      file.rename(out, file)
    }
    )
  # pdf REPORT
}
