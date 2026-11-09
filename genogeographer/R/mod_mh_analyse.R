#' Microhaplotype (MH) profile module - UI
#'
#' Upload a microhaplotype (MH) profile (CSV/Excel) and compare it against the
#' precomputed MH genotype probabilities per meta-population. Unlike the AIMs
#' analysis, this uses a precomputed lookup table (`MH_profiles_probs.Rds`) rather
#' than fitting/simulating a model on the fly, since that fitting step is expensive.
#' @param id Module id
#' @noRd
mod_mh_analyse_ui <- function(id){
  ns <- NS(id)
  nav_panel(
    title = "Analyse MH profile",
    layout_sidebar(
      sidebar = sidebar(
        width = 350, open = TRUE,
        h4("Input file"),
        fileInput(ns("mh_file"), "Choose CSV/Excel file", width = "100%",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".xlsx", ".xls")),
        h4("Settings"),
        uiOutput(ns("sample_select")),
        shinyWidgets::sliderTextInput(inputId = ns("CI"), label = "Confidence level:",
                                      choices = c(95, 97.5, 99, 99.9, 99.99), selected = 95, post = "%"),
        uiOutput(ns("side_pvalue")),
        uiOutput(ns("LR_select")),
        tags$style(type = 'text/css', paste0(
          "#", ns("analyse"), ", #", ns("reset"), " { width:100%; margin-top: 25px;}"
        )),
        withBusyIndicatorUI(
          actionButton(inputId = ns("analyse"), label = "Analyse!", icon = icon("calculator", verify_fa = FALSE), class = "btn-primary")
        ),
        actionButton(inputId = ns("reset"), label = "Reset", icon = icon("trash", verify_fa = FALSE)),
        tags$hr(),
        div(helpText(paste0("Version: ", ggg_package, " (", packageVersion(ggg_package), ")"))),
        div(helpText(paste0("Developer: ", packageDescription(ggg_package, fields = "Maintainer"))))
      ),
      uiOutput(ns("analysis"))
    )
  )
}

#' Microhaplotype (MH) profile module - server
#' @param id Module id
#' @param mh_data A list with elements `probs`, `snp_map` and `out_metas` overriding the bundled
#' precomputed MH genotype probabilities / MH-to-SNP map / outlier null models (see
#' `mh_profiles_probs()`/`mh_snp_map()`/`mh_out_metas()`). If `NULL`, the data shipped with the
#' package is used.
#' @noRd
mod_mh_analyse_server <- function(id, mh_data = NULL){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    ## build fixes : start ##
    sample_name <- NULL; meta <- NULL
    likelihood <- NULL; log10_likelihood <- NULL; LR_best <- NULL; rel_prob <- NULL; n_loci <- NULL
    devi <- NULL; z_score <- NULL; p_value <- NULL; critical_value <- NULL; accept <- NULL
    ## build fixes : end ##

    snp_map_used <- if(!is.null(mh_data$snp_map)) mh_data$snp_map else mh_snp_map()
    probs_tbl_used <- if(!is.null(mh_data$probs)) mh_build_probs_table(mh_data$probs) else mh_geno_probs_table()
    devi_tbl_used <- if(!is.null(mh_data$probs)) mh_build_devi_table(mh_data$probs) else mh_geno_devi_table()
    out_metas_used <- if(!is.null(mh_data$out_metas)) mh_prepare_out_metas(mh_data$out_metas) else mh_out_metas()

    observeEvent(input$analyse, {
      withBusyIndicatorServer(buttonId = "analyse", session = session, expr = reactive_mh_likelihood())
    })

    observeEvent(input$reset, {
      runjs("history.go(0)")
    })

    reactive_mh_canon <- reactive({
      if(is.null(input$mh_file)) return(NULL)
      tryCatch(read_mh_profile_file(input$mh_file$datapath), error = function(e){
        shiny::showNotification(paste("Could not read file:", conditionMessage(e)), type = "error")
        NULL
      })
    })

    reactive_mh_likelihood <- eventReactive(list(input$mh_file, input$analyse), {
      canon <- reactive_mh_canon()
      if(is.null(canon)) return(NULL)
      alpha <- 1 - isolate(input$CI) / 100
      mh_profile_likelihood(canon, snp_map = snp_map_used, probs_tbl = probs_tbl_used,
                            devi_tbl = devi_tbl_used, out_metas = out_metas_used, alpha = alpha)
    })

    output$sample_select <- renderUI({
      res <- reactive_mh_likelihood()
      req(res$result)
      selectInput(ns("sample_id"), "Sample:",
                 choices = unique(res$result$sample_name))
    })

    sample_download <- function(filename){
      downloadHandler(
        filename = function() filename,
        contentType = "text/csv",
        content = function(file) file.copy(system.file("extdata", filename, package = ggg_package), file)
      )
    }
    output$dl_long_two <- sample_download("mh_sample_long.csv")
    output$dl_long_one <- sample_download("mh_sample_long2.csv")
    output$dl_wide_two <- sample_download("mh_sample_wide.csv")
    output$dl_wide_one <- sample_download("mh_sample_wide2.csv")

    example_preview_dt <- function(filename){
      renderDT({
        datatable(mh_example_preview(filename), rownames = FALSE,
                 options = list(dom = "t", scrollX = TRUE))
      })
    }
    output$preview_long_two <- example_preview_dt("mh_sample_long.csv")
    output$preview_long_one <- example_preview_dt("mh_sample_long2.csv")
    output$preview_wide_two <- example_preview_dt("mh_sample_wide.csv")
    output$preview_wide_one <- example_preview_dt("mh_sample_wide2.csv")

    ## USER INTERFACE
    output$analysis <- renderUI({
      lik <- reactive_mh_likelihood()
      dat <- reactive_mh_canon()
      if(is.null(lik) || is.null(lik$result)){
        if(is.null(dat) || nrow(dat) == 0){
          tagList(
            h3("Upload a valid MH profile"),
            HTML("<p>Upload a file with individual SNP genotypes for each sample, in one of the following formats
                 (detected automatically). Each is illustrated below with the first rows/columns of the
                 corresponding sample file, downloadable via the link underneath.</p>"),
            layout_columns(
              card(
                card_header("Long format"),
                h5("Two columns (phased)"), tags$code("sample_name;locus;A;B"),
                p("One row per sample/locus, alleles phased into haplotype copies A/B."),
                DTOutput(ns("preview_long_two")), downloadLink(ns("dl_long_two"), "Download"), tags$hr(),
                h5("One column (unphased)"), tags$code("sample_name;locus;allele"),
                p("One row per sample/locus, unphased genotype (e.g. \"GT\")."),
                DTOutput(ns("preview_long_one")), downloadLink(ns("dl_long_one"), "Download")
              ),
              card(
                card_header("Wide format"),
                h5("Two rows (phased)"),
                p("One column per locus, two rows per sample - one per haplotype copy."),
                DTOutput(ns("preview_wide_two")), downloadLink(ns("dl_wide_two"), "Download"), tags$hr(),
                h5("One row (unphased)"),
                p("One column per locus, one row per sample with two-character genotypes."),
                DTOutput(ns("preview_wide_one")), downloadLink(ns("dl_wide_one"), "Download")
              )
            )
          )
        }
        else{
          tagList(
            h3("Uploaded file"),
            p(paste0("The file was recognised as ", if(isTRUE(dat$phased[1])) "phased" else "unphased",
                     " input with ", dplyr::n_distinct(dat$sample_name), " sample(s) and ",
                     dplyr::n_distinct(dat$locus), " locus/loci.")),
            p(HTML("Click <b>Analyse!</b> in the sidebar to compute the meta-population likelihoods.")),
            card(
              renderDT({
                datatable(head(dat, 20), options = list(dom = "t"))
              })
            )
          )
        }
      }
      else{
        panels <- list()
        if(length(lik$unmatched_loci) > 0){
          panels <- c(panels, list(
            div(class = "alert alert-warning",
                HTML(paste0("<b>", length(lik$unmatched_loci), " locus/loci in the file were not recognised</b> ",
                            "and were excluded from the analysis: ", and_text(lik$unmatched_loci), "."))
            )
          ))
        }
        if(nrow(lik$incomplete_mh) > 0){
          panels <- c(panels, list(
            div(class = "alert alert-warning",
                HTML(paste0("<b>", nrow(lik$incomplete_mh), " sample/MH marker combination(s)</b> had missing SNPs ",
                            "and were excluded (a microhaplotype marker is only used when all its SNPs are present)."))
            )
          ))
        }
        if(!is.null(lik$not_found) && nrow(lik$not_found) > 0){
          panels <- c(panels, list(
            div(class = "alert alert-warning",
                HTML(paste0("<b>", nrow(lik$not_found), " sample/MH marker genotype combination(s)</b> were not ",
                            "found in the reference data and were excluded from the corresponding sample's likelihood."))
            )
          ))
        }
        panels <- c(panels, list(
          layout_columns(
            col_widths = 8,
            card(
              card_header("Outlier test"),
              uiOutput(ns("barplot_panel"))
            )
          ),
          card(
            card_header("Results"),
            HTML("<p>The table below shows, for the selected sample, the likelihood of observing its MH
                 profile in each meta-population (based on precomputed genotype probabilities), the
                 likelihood ratio versus the most likely meta-population (<code>LR_best</code>) and the
                 relative probability across meta-populations assuming a uniform prior
                 (<code>rel_prob</code>, NOT a p-value). The outlier test (<code>z_score</code>,
                 <code>p_value</code>, <code>accept</code>) is calibrated against precomputed null models
                 fit on the full microhaplotype panel. When linkage phase is unknown (unphased input with
                 2+ heterozygous SNPs in a marker), probabilities are summed and deviances
                 probability-weighted over all phase resolutions consistent with the observed genotypes.</p>"),
            DTOutput(ns("result_table"))
          ),
          card(
            card_header("Likelihood ratios"),
            textOutput(ns("LR_one_pop")),
            DTOutput(ns("lr_list"))
          )
        ))
        do.call(tagList, panels)
      }
    })

    reactive_sample_result <- reactive({
      res <- reactive_mh_likelihood()
      req(res$result, input$sample_id)
      res$result %>% dplyr::filter(sample_name == input$sample_id)
    })

    output$side_pvalue <- renderUI({
      df <- reactive_sample_result()
      if(is.null(df) || nrow(df) == 0) return(NULL)
      accepted <- df %>% dplyr::filter(accept)
      status_text <- if(nrow(accepted) == 0) "All meta-populations are rejected" else{
        best <- accepted %>% dplyr::slice_max(n = 1, order_by = p_value)
        paste0(best$meta, " (p = ", round(best$p_value, 3), ")")
      }
      div(
        h4("Reporting"),
        HTML(paste0("<p><b>Outlier test (largest p-value):</b><br/>", status_text, "</p>"))
      )
    })

    output$LR_select <- renderUI({
      df <- reactive_sample_result()
      if(is.null(df) || nrow(df) == 0) return(NULL)
      accepted <- df %>% dplyr::filter(accept) %>% dplyr::pull(meta)
      div(
        checkboxGroupInput(inputId = ns("LR_accept"),
                           label = "LRs with two rejected populations:",
                           choiceNames = list("Allow (LRs may be misleading)"), choiceValues = list("allow")),
        selectizeInput(inputId = ns("LR_selected"),
                       label = "Include in LR calculations:",
                       choices = df$meta,
                       selected = accepted,
                       multiple = TRUE,
                       options = list(plugins = list("remove_button", "drag_drop")),
                       width = "100%")
      )
    })

    output$LR_one_pop <- renderText({
      df <- reactive_sample_result()
      if(is.null(df)) return(NULL)
      if(nrow(df) <= 1) "Likelihood ratios are not computed, since there is just a single meta-population in the results"
      else return(NULL)
    })

    reactive_LR_table <- reactive({
      df <- reactive_sample_result()
      req(df)
      if(nrow(df) <= 1) return(NULL)
      accepted_control <- if(is.null(input$LR_accept)) FALSE else (input$LR_accept == "allow")
      mh_LR_table(df, who = input$LR_selected, only_accepted = !accepted_control)
    })

    output$lr_list <- renderDT({
      lr_tab <- reactive_LR_table()
      req(lr_tab)
      lr_tab %>%
        dplyr::mutate(log10_LR = signif(log10_LR, 4)) %>%
        dplyr::rename(Numerator = numerator, Denominator = denominator, `log10 LR` = log10_LR) %>%
        dplyr::select(Numerator, Denominator, `log10 LR`) %>%
        dt_table()
    })

    output$result_table <- renderDT({
      df <- reactive_sample_result()
      req(df)
      df %>%
        dplyr::mutate(dplyr::across(c(log10_likelihood, z_score, p_value, LR_best, rel_prob), signif, digits = 4)) %>%
        dplyr::select(meta, n_loci, log10_likelihood, z_score, p_value, accept, LR_best, rel_prob) %>%
        dplyr::arrange(dplyr::desc(p_value)) %>%
        dt_table()
    })

    output$barplot_panel <- renderUI({
      withSpinner(plotlyOutput(ns("barplot"), width = "100%"), type = 4)
    })

    output$barplot <- renderPlotly({
      df <- reactive_sample_result()
      req(df)
      ci <- input$CI / 100
      p <- mh_error_bar_plot(df, ci = ci)
      ggplotly(p, tooltip = "text") %>%
        layout(showlegend = FALSE) %>%
        config(displayModeBar = FALSE)
    })
  })
}

