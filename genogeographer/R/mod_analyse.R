#' Analyse AIMs profile module - UI
#'
#' The main analysis tab: upload a profile, configure settings, and inspect results.
#' @param id Module id
#' @noRd
mod_analyse_ui <- function(id){
  ns <- NS(id)
  nav_panel(
    title = "Analyse AIMs profile",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 350, open = TRUE,
        h4("Input file"),
        mod_profile_upload_ui(ns("profile")),
        h4("Settings"),
        radioButtons(inputId = ns("meta"), label = "Grouping type:", choices = c("Meta-populations" = "meta", "Populations" = "pop")),
        sliderInput(inputId = ns("min_n"), label = "Minimum sample size:", min = 5, max = 200, step = 5, value = 75),
        checkboxGroupInput(inputId = ns("admix"), label = "Analyse 1st order admixture:",
                           choiceNames = list("Admixture"), choiceValues = list("admix")),
        shinyWidgets::sliderTextInput(inputId = ns("CI"), label = "Confidence level:", choices = c(95, 97.5, 99, 99.9, 99.99), selected = 95, post = "%"),
        mod_db_selector_ui(ns("db")),
        checkboxGroupInput(inputId = ns("tilt"), label = "Adjust p-values by exponential tilting:",
                           choiceNames = list("Adjust (may take some time)"), choiceValues = list("adjust")),
        uiOutput(ns("side_pvalue")),
        uiOutput(ns("LR_select")),
        tags$style(type = 'text/css', paste0(
          "#", ns("analyse"), ", #", ns("reset"), ", #", ns("report_download"), " { width:100%; margin-top: 25px;}"
        )),
        withBusyIndicatorUI(
          actionButton(inputId = ns("analyse"), label = "Analyse!", icon = icon("calculator", verify_fa = FALSE), class = "btn-primary")
        ),
        actionButton(inputId = ns("reset"), label = "Reset", icon = icon("trash", verify_fa = FALSE)),
        tags$hr(),
        uiOutput(ns("report_panel")),
        div(helpText(paste0("Version: ", ggg_package, " (", packageVersion(ggg_package), ")"))),
        div(helpText(paste0("Developer: ", packageDescription(ggg_package, fields = "Maintainer"))))
      ),
      uiOutput(ns("analysis"))
    )
  )
}

#' Analyse AIMs profile module - server
#'
#' @param id Module id
#' @param db_list A `reactiveValues` object with a `db` element holding the named list of reference databases.
#' @param reporting_panel Logical. Should the report generation/download be available.
#' @noRd
mod_analyse_server <- function(id, db_list, reporting_panel = TRUE){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    ## build fixes : start ##
    Target.ID <- NULL; Genotype <- NULL; locus <- NULL; genotype <- NULL
    accept <- NULL; p_value <- NULL; selected_ <- NULL
    meta <- NULL; pop <- NULL; . <- NULL
    lat <- NULL; lon <- NULL; aims_example <- NULL
    ## build fixes : end ##

    if(is.null(isolate(db_list$db))) shinyjs::disable("analyse")

    profile_upload <- mod_profile_upload_server("profile")
    db_selector <- mod_db_selector_server("db", db_list = db_list)

    observeEvent(input$analyse, {
      withBusyIndicatorServer(buttonId = "analyse", session = session, expr = reactive_result())
    })

    observeEvent(db_selector$has_db(), {
      if(db_selector$has_db()) shinyjs::enable("analyse")
    })

    ## USER INTERFACE
    output$analysis <- renderUI({
      res <- reactive_result()
      dat <- profile_upload$data()
      if(is.null(res)){
        if(is.null(dat) || nrow(dat) == 0){
          card(
            fill = FALSE,
            card_header("Upload valid AIMs profile"),
            card_body(
              fill = FALSE,
              HTML("<p>Upload valid file: It must contain a column containing marker (locus) names, and a column with the genotypes as e.g. <tt>AA</tt> or <tt>AC</tt></p>"),
              HTML("You can download a sample file showing the necessary columns and data structure here:"),
              downloadLink(ns("download_sample"), "Download")
            )
          )
        }
        else{
          ## Make column suggestions
          guess_locus_text <- and_text(pre = "<p>Suggestions for <b>locus column</b>: ", x = profile_upload$locus_guess(), post = "</p>")
          guess_genotype_text <- and_text(pre = "<p>Suggestions for <b>genotype column</b>: ", x = profile_upload$genotype_guess(), post = "</p>")
          ## ggg_score() returns NULL when no group reaches the minimum sample size
          min_n <- isolate(input$min_n)
          n_max <- NA_real_
          if(!is.null(reactive_profile()) && !is.null(min_n)){
            db_info <- attr(db_list$db[[db_selector$snp_set()]][[isolate(input$meta)]]$db, "info")
            if(!is.null(db_info)) n_max <- suppressWarnings(max(as.numeric(db_info$n), na.rm = TRUE))
          }
          no_group_msg <- if(is.finite(n_max) && n_max < min_n){
            div(class = "alert alert-warning",
                paste0("No group in the selected reference database has at least ", min_n,
                       " samples (the largest has ", n_max, "). Lower 'Minimum sample size' and click 'Analyse!'."))
          }
          tagList(
            no_group_msg,
            card(
              fill = FALSE,
              card_header("Uploaded file"),
              card_body(
                fill = FALSE,
                HTML("<p>The first few rows of the uploaded file are shown below</p>"),
                HTML(guess_locus_text),
                HTML(guess_genotype_text),
                renderDT({
                  datatable(dat, options(dom = "t"))
                })
              )
            )
          )
        }
      }
      else{
        tagList(
          tags$head(tags$style(paste0(".modal-lg{ width: ", 2*session$clientData$output_barplot_width, "px}"))),
          div(HTML(paste0("<p><b>Analysis of file:</b> ", profile_upload$file()$name, ":")),
              actionLink(ns("show_profile"), label = "Show profile"),
              icon("new-window", lib = "glyphicon", verify_fa = FALSE),
              HTML("</p>")
          ),
          card(
            fill = FALSE,
            card_header("Graphics"),
            card_body(
              fill = FALSE,
              layout_columns(
                fill = FALSE,
                col_widths = c(6, 6),
                uiOutput(ns("barplot_panel")),
                uiOutput(ns("map_panel"))
              )
            )
          ),
          card(
            fill = FALSE,
            card_header("Tables"),
            card_body(fill = FALSE, DTOutput(ns("result_table"), fill = FALSE))
          ),
          card(
            fill = FALSE,
            card_header("Likelihood ratios"),
            card_body(
              fill = FALSE,
              textOutput(ns("LR_one_pop")),
              layout_columns(
                fill = FALSE,
                col_widths = c(6, 6),
                DTOutput(ns("lr_list"), fill = FALSE),
                uiOutput(ns("LRplot_panel"))
              )
            )
          )
        )
      }
    })

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

    observeEvent(input$show_profile, {
      A1 <- NULL
      A2 <- NULL
      raw_profile <- profile_upload$data() %>% rename(locus = !!sym(profile_upload$locus_col()))
      profile <- reactive_profile()
      profile_drop <- raw_profile %>% anti_join(profile, by = "locus") %>%
        rename(!!sym(profile_upload$locus_col()) := locus)
      profile_x0 <- raw_profile %>% semi_join(profile, by = "locus")
      showModal(modalDialog(
        title = paste("Uploaded profile:", profile_upload$file()$name),
        size = "m",
        h3("Analysed loci"),
        helpText("The loci below has been included in the analysis."),
        renderDT(profile_x0 %>% dt_table),
        h3("Dropped or unused loci"),
        helpText(paste0("The loci below has been excluded from the analysis.\n
               Either because of state 'NN', locus not in '", db_selector$snp_set(), "'
                      or other typing error (e.g. different reference allele).")),
        renderDT(profile_drop %>% dt_table),
        footer = modalButton("Close"),
        easyClose = TRUE
      ))
    })

    ### ANALYSIS

    observeEvent(input$reset, {
      runjs("history.go(0)")
    })

    output$fileUploaded <- reactive({
      dat <- profile_upload$data()
      !is.null(dat) && nrow(dat) > 0
    })
    outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)

    reactive_profile <- reactive({
      profile <- profile_upload$data()
      if(is.null(profile) || nrow(profile) == 0) return(NULL) ## No profile
      nprofile <- names(profile)
      col_locus <- profile_upload$locus_col()
      col_genotype <- profile_upload$genotype_col()
      if(is.null(col_locus) || is.null(col_genotype)) return(NULL)
      if(!(col_locus %in% nprofile) || !(col_genotype %in% nprofile)) return(NULL)
      if(col_locus == col_genotype) return(NULL)
      ##
      req(db_selector$snp_set())
      ggg_db <- attr(db_list$db[[db_selector$snp_set()]]$meta$db, "allele_list")
      ##
      genotype_x0(profile = profile, locus = col_locus, genotype = col_genotype, ggg = ggg_db)
    })

    reactive_result <- eventReactive(list(profile_upload$file(), input$analyse), {
      profile <- reactive_profile()
      if(is.null(profile)) return(NULL)
      ### TILTING
      tilt_control <- if(is.null(input$tilt)) FALSE else (input$tilt == "adjust")
      ## ADMIXTURE
      admix_control <- if(is.null(input$admix)) FALSE else (input$admix == "admix")
      ## COMPUTE
      res <- ggg_score(profile_x0 = profile, DB = db_list$db[[db_selector$snp_set()]][[input$meta]]$db,
                       min_n = input$min_n, CI = input$CI/100, tilt = tilt_control)
      if(admix_control){
        res_admix <- ggg_score(profile_x0 = profile, DB = db_list$db[[db_selector$snp_set()]][[input$meta]]$admix,
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
      if(nrow(largest_p_value) == 0) largest_p_value_text <- paste0("All ", grouping, "s are rejected")
      else{
        largest_p_value <- largest_p_value %>% slice_max(n = 1, order_by = p_value) %>%
          select(p_value, !!groups_) %>%
          mutate(p_value = round(p_value, 3)) %>%
          unite(p_text, p_value, !!groups_, sep = " (")
        largest_p_value_text <- paste0("<b>DB-score (largest p-value):</b><br/>", largest_p_value, ")")
      }
      div(
        h4("Reporting"),
        HTML(paste0("<p>", largest_p_value_text, "</p>"))
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
      withSpinner(leafletOutput(ns("map"), width = "100%"), type = 4)
    })

    ### BARPLOT

    output$barplot_panel <- renderUI({
      withSpinner(plotlyOutput(ns("barplot"), width = "100%", fill = FALSE), type = 4)
    })

    output$barplot <- renderPlotly({
      res <- reactive_result()
      if(is.null(res)) return(NULL)
      ebp <- error_bar_plotly(result_df = res, which = isolate(input$result_table_rows_selected))
      ebp %>% ggplotly(tooltip = c("text")) %>%
        layout(showlegend = FALSE) %>%
        config(displayModeBar = FALSE)
    })

    ## TABLES

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
      div(
        checkboxGroupInput(inputId = ns("LR_accept"),
                           label = "LRs with two rejected populations:",
                           choiceNames = list("Allow (LRs may be misleading)"), choiceValues = list("allow")),
        selectizeInput(inputId = ns("LR_selected"),
                       label = "Include in LR calculations:",
                       choices = LR_choices,
                       selected = accepted,
                       multiple = TRUE,
                       options = list(plugins = list("remove_button", "drag_drop")),
                       width = "100%")
      )
    })

    output$LR_one_pop <- renderText({
      res <- reactive_result()
      if(is.null(res)) return(NULL)
      input_meta <- isolate(input$meta)
      mm <- if(is.null(input_meta)) "meta" else input_meta
      MM <- if(mm == "meta") "metapopulation" else "population"
      if(nrow(res) <= 1) paste("Likelihood ratios are not computed, since there is just a single", MM, "in the results")
      else return(NULL)
    })

    reactive_LR_table <- reactive({
      accepted_control <- if(is.null(input$LR_accept)) FALSE else (input$LR_accept == "allow")
      res <- reactive_result()
      if(is.null(res) || nrow(res) <= 1) return(NULL)
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

    observeEvent(input$LR_selected, { ## Can't delete most probable (based on z_score)
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
      withSpinner(plotlyOutput(ns("LRplot"), width = "100%", fill = FALSE), type = 4)
    })

    ### RETURN pdf REPORT

    output$report_panel <- renderUI({
      res <- reactive_result()
      if(!reporting_panel) return(div())
      if(is.null(res)) return(div())
      div(
        h4("Report"),
        textInput(inputId = ns("name"), label = "Name of analyst", width = "100%",
                 placeholder = "Name as to appear in report", value = input$name),
        withBusyIndicatorUI(downloadButton(outputId = ns("report_download"), label = "Download HTML report", class = "btn-primary")),
        hr()
      )
    })

    output$report_download <- downloadHandler(
      filename = function(){
        paste(sub("\\.[[:alnum:]]*.$", "", profile_upload$file()$name),
              "html",
              sep = "."
        )
      },
      content = function(file) {
        src <- list(rmd_file = normalizePath(system.file("deployable_app", "aims_report.Rmd", package = ggg_package)))
        owd <- setwd(tempdir())
        cat(file = stderr(), owd, "\n")
        on.exit(setwd(owd))
        file.copy(src$rmd_file, 'aims_report.Rmd', overwrite = TRUE)
        out <- rmarkdown::render(input = 'aims_report.Rmd', clean = TRUE,
                                 output_format = rmarkdown::html_document(),
          params = list(
            set_file = profile_upload$file()$name,
            set_author = input$name,
            set_output = "HTML",
            set_locus = profile_upload$locus_col(),
            set_genotype = profile_upload$genotype_col(),
            set_db = db_selector$snp_set()
          )
        )
        file.rename(out, file)
      }
    )
  })
}
