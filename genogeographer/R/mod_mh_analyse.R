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
      fillable = FALSE,
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
#' @param mh_data A list with elements `probs`, `snp_map`, `out_metas` and `ci_data` overriding the bundled
#' precomputed MH genotype probabilities / MH-to-SNP map / outlier null models / reference haplotype data for
#' the confidence intervals (see `mh_profiles_probs()`/`mh_snp_map()`/`mh_out_metas()`/`mh_ci_data()`). If
#' `NULL`, the data shipped with the package is used. If `probs` is supplied without `ci_data`, no confidence
#' intervals are computed (the bundled reference haplotypes would not match).
#' @param db_list Optional `reactiveValues` with a `db` element holding the AIMs reference databases; the map
#' centres of the meta-populations are taken from these (see `mh_meta_centres()`).
#' @noRd
mod_mh_analyse_server <- function(id, mh_data = NULL, db_list = NULL){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    ## build fixes : start ##
    sample_name <- NULL; meta <- NULL
    likelihood <- NULL; log10_likelihood <- NULL; LR_best <- NULL; rel_prob <- NULL; n_loci <- NULL
    devi <- NULL; z_score <- NULL; p_value <- NULL; critical_value <- NULL; accept <- NULL; se_log10 <- NULL
    ## build fixes : end ##

    snp_map_used <- if(!is.null(mh_data$snp_map)) mh_data$snp_map else mh_snp_map()
    probs_tbl_used <- if(!is.null(mh_data$probs)) mh_build_probs_table(mh_data$probs) else mh_geno_probs_table()
    devi_tbl_used <- if(!is.null(mh_data$probs)) mh_build_devi_table(mh_data$probs) else mh_geno_devi_table()
    out_metas_used <- if(!is.null(mh_data$out_metas)) mh_prepare_out_metas(mh_data$out_metas) else mh_out_metas()
    ci_data_used <- if(!is.null(mh_data$ci_data)) mh_ci_prepare(mh_data$ci_data) else if(is.null(mh_data$probs)) mh_ci_data() else NULL

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
                            devi_tbl = devi_tbl_used, out_metas = out_metas_used, ci_data = ci_data_used, alpha = alpha)
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
            HTML("<p>Upload a file with phased SNP genotypes for each sample, in one of the following formats
                 (detected automatically). The alleles are always phased: in a two-character genotype such as
                 <code>TG</code> the first character is the allele on the first haplotype copy and the second
                 character the allele on the second, so <code>TG</code> and <code>GT</code> are different.
                 Each format is illustrated below with the first rows/columns of the corresponding sample file,
                 downloadable via the link underneath.</p>"),
            layout_columns(
              fill = FALSE,
              card(
                fill = FALSE,
                card_header("Long format"),
                card_body(
                  fill = FALSE,
                  h5("Two allele columns"), tags$code("sample_name;locus;A;B"),
                  p("One row per sample/locus, with the alleles of haplotype copies A and B."),
                  DTOutput(ns("preview_long_two"), fill = FALSE), downloadLink(ns("dl_long_two"), "Download"), tags$hr(),
                  h5("One genotype column"), tags$code("sample_name;locus;allele"),
                  p("One row per sample/locus, with the genotype of the two haplotype copies (e.g. \"TG\": T on copy 1, G on copy 2)."),
                  DTOutput(ns("preview_long_one"), fill = FALSE), downloadLink(ns("dl_long_one"), "Download")
                )
              ),
              card(
                fill = FALSE,
                card_header("Wide format"),
                card_body(
                  fill = FALSE,
                  h5("Two rows per sample"),
                  p("One column per locus, two rows per sample - one per haplotype copy."),
                  DTOutput(ns("preview_wide_two"), fill = FALSE), downloadLink(ns("dl_wide_two"), "Download"), tags$hr(),
                  h5("One row per sample"),
                  p("One column per locus, one row per sample with the two-character genotype of the haplotype copies (first character = copy 1)."),
                  DTOutput(ns("preview_wide_one"), fill = FALSE), downloadLink(ns("dl_wide_one"), "Download")
                )
              )
            )
          )
        }
        else{
          tagList(
            h3("Uploaded file"),
            p(paste0("The file contains ", dplyr::n_distinct(dat$sample_name), " sample(s) and ",
                     dplyr::n_distinct(dat$locus), " locus/loci.")),
            p(HTML("Click <b>Analyse!</b> in the sidebar to compute the meta-population likelihoods.")),
            card(
              fill = FALSE,
              card_body(
                fill = FALSE,
                renderDT({
                  datatable(head(dat, 20), options = list(dom = "t"))
                })
              )
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
        if(!is.null(lik$ci_unobserved) && nrow(lik$ci_unobserved) > 0){
          n_comb <- nrow(dplyr::distinct(lik$ci_unobserved, sample_name, meta))
          panels <- c(panels, list(
            div(class = "alert alert-warning",
                HTML(paste0("<b>", n_comb, " sample/meta-population combination(s) contain haplotypes that have not ",
                            "previously been observed</b> in the reference data of the meta-population (see the table below the plot). ",
                            "The likelihood of these rests on the profile's own haplotype alone, and the confidence interval is ",
                            "correspondingly uncertain."))
            )
          ))
        }
        if(!is.null(lik$ci_failed) && nrow(lik$ci_failed) > 0){
          panels <- c(panels, list(
            div(class = "alert alert-danger",
                HTML(paste0("<b>The confidence interval of the log<sub>10</sub> likelihood could not be calculated for ",
                            nrow(lik$ci_failed), " sample/meta-population combination(s)</b>: ", and_text(unique(lik$ci_failed$meta)),
                            ". The likelihood itself is not affected."))
            )
          ))
        }
        panels <- c(panels, list(
          layout_columns(
            fill = FALSE,
            col_widths = c(6, 6),
            card(
              fill = FALSE,
              card_header("Profile probability"),
              card_body(fill = FALSE, uiOutput(ns("barplot_panel")))
            ),
            card(
              fill = FALSE,
              card_header("Map"),
              card_body(fill = FALSE, uiOutput(ns("map_panel")), uiOutput(ns("map_note")))
            )
          ),
          uiOutput(ns("ci_unobserved_panel")),
          card(
            fill = FALSE,
            card_header("Results"),
            card_body(
              fill = FALSE,
              HTML("<p>The table below shows, for the selected sample, the likelihood of observing its MH
                 profile in each meta-population (based on precomputed genotype probabilities) and its
                 first-order standard error (<code>se_log10</code>, on the log<sub>10</sub> scale; missing if it could not be calculated), the likelihood ratio versus the most likely meta-population (<code>LR_best</code>) and the
                 relative probability across meta-populations assuming a uniform prior
                 (<code>rel_prob</code>, NOT a p-value). The outlier test (<code>z_score</code>,
                 <code>p_value</code>, <code>accept</code>) is calibrated against precomputed null models
                 fit on the full microhaplotype panel.</p>"),
              DTOutput(ns("result_table"), fill = FALSE)
            )
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
      ## Without standard errors (no reference haplotype data) only point estimates are available
      if(anyNA(df$se_log10)) return(list(res = NULL, LR_Tab = mh_LR_table(df, who = input$LR_selected, only_accepted = !accepted_control)))
      res <- mh_as_aims_result(df)
      list(res = res,
           LR_Tab = LR_table(z_df = res, who = input$LR_selected, CI = input$CI / 100, only_accepted = !accepted_control)$LR)
    })

    output$lr_list <- renderDT({
      LR_Tab <- reactive_LR_table()
      req(LR_Tab)
      if(is.null(LR_Tab$res)){
        return(LR_Tab$LR_Tab %>%
          dplyr::rename(Numerator = numerator, Denominator = denominator, `log10 LR` = log10_LR) %>%
          dplyr::select(Numerator, Denominator, `log10 LR`) %>%
          dt_table() %>%
          DT::formatSignif("log10 LR", digits = 4))
      }
      LR_list(result = LR_Tab$res, LR_tab = LR_Tab$LR_Tab)
    })

    observeEvent(input$LR_selected, { ## The most probable accepted meta-population cannot be removed
      df <- reactive_sample_result()
      if(is.null(df) || !any(df$accept)) return(NULL)
      max_p <- df$meta[df$accept][which.max(df$p_value[df$accept])]
      if(!max_p %in% input$LR_selected) updateSelectizeInput(session, "LR_selected", selected = unique(c(input$LR_selected, max_p)))
    })

    output$LRplot <- renderPlotly({
      LR_Tab <- reactive_LR_table()
      req(LR_Tab, LR_Tab$res)
      if(is.null(LR_Tab$LR_Tab) || nrow(LR_Tab$LR_Tab) == 0) return(NULL)
      LR_plot_ly(result = LR_Tab$res, LR_list = LR_Tab$LR_Tab)
    })

    output$LRplot_panel <- renderUI({
      withSpinner(plotlyOutput(ns("LRplot"), width = "100%", fill = FALSE), type = 4)
    })

    output$result_table <- renderDT({
      df <- reactive_sample_result()
      req(df)
      df %>%
        dplyr::select(meta, n_loci, log10_likelihood, se_log10, z_score, p_value, accept, LR_best, rel_prob) %>%
        dplyr::arrange(dplyr::desc(p_value)) %>%
        dt_table() %>%
        DT::formatSignif(c("log10_likelihood", "se_log10", "z_score", "LR_best", "rel_prob"), digits = 4) %>%
        DT::formatRound("p_value", digits = 4)
    })

    output$barplot_panel <- renderUI({
      withSpinner(plotlyOutput(ns("barplot"), width = "100%", fill = FALSE), type = 4)
    })

    output$barplot <- renderPlotly({
      df <- reactive_sample_result()
      req(df)
      p <- mh_error_bar_plot(df, ci = input$CI / 100)
      ggplotly(p, tooltip = "text") %>%
        layout(showlegend = FALSE) %>%
        config(displayModeBar = FALSE)
    })

    reactive_centres <- reactive({
      mh_meta_centres(if(is.null(db_list)) NULL else db_list$db)
    })

    ind_per_meta <- if(is.null(ci_data_used)) NULL else vapply(ci_data_used, function(x) x[[1]]$N / 2, numeric(1))

    output$map_panel <- renderUI({
      withSpinner(leafletOutput(ns("map"), width = "100%"), type = 4)
    })

    output$map <- renderLeaflet({
      df <- reactive_sample_result()
      req(df)
      leaflet_mh_plot(df, reactive_centres(), ci = input$CI / 100, n_ind = ind_per_meta)
    })

    output$map_note <- renderUI({
      df <- reactive_sample_result()
      req(df)
      no_centre <- setdiff(as.character(df$meta), reactive_centres()$meta)
      if(length(no_centre) == 0) return(NULL)
      p(class = "small text-muted", "Not shown on the map (no centre in the reference databases): ", paste(no_centre, collapse = ", "), ".")
    })

    reactive_ci_unobserved <- reactive({
      lik <- reactive_mh_likelihood()
      req(lik, input$sample_id)
      if(is.null(lik$ci_unobserved)) return(NULL)
      lik$ci_unobserved %>% dplyr::filter(sample_name == input$sample_id)
    })

    output$ci_unobserved_panel <- renderUI({
      unobs <- reactive_ci_unobserved()
      if(is.null(unobs) || nrow(unobs) == 0) return(NULL)
      card(
        fill = FALSE,
        card_header("Haplotypes not previously observed"),
        card_body(
          fill = FALSE,
          p("The profile contains a haplotype with a combination of SNP alleles that is not observed in the reference data of ",
            "these meta-populations. Its probability rests on the profile's own haplotype alone, so the confidence interval of ",
            "the log likelihood is correspondingly uncertain."),
          DTOutput(ns("ci_unobserved_table"), fill = FALSE)
        )
      )
    })

    output$ci_unobserved_table <- renderDT({
      unobs <- reactive_ci_unobserved()
      req(unobs)
      unobs %>%
        dplyr::select(`Meta-population` = meta, `MH marker` = MH, Haplotype = haplotype, `Unobserved SNP alleles` = snps) %>%
        dt_table()
    })
  })
}

