
#' @title UI part for Shiny app
#' @description The shiny front-end
#' @export
ui_api <- function(){
  ggg_references <- list(
    list(author = "T Tvedebrink, PS Eriksen, HS Mogensen, N Morling",
         title = "GenoGeographer - A tool for genogeographic inference.",
         journal = "Forensic Science International: Genetics Supplement Series 6, e463-e465.",
         year = 2017),
    list(author = "T Tvedebrink, PS Eriksen, HS Mogensen, N Morling",
         title = "Weight of the evidence of genetic investigations of ancestry informative markers.",
         journal = "Theoretical Population Biology 120, 1-10.",
         year = 2018),
    list(author = "T Tvedebrink, PS Eriksen",
         title = "Inference of admixed ancestry with Ancestry Informative Markers.",
         journal = "Forensic Science International: Genetics 42, 147-153.",
         year = 2019),
    list(author = "HS Mogensen, T Tvedebrink, C Borsting, V Pereira, N Morling",
         title = "Ancestry prediction efficiency of the software GenoGeographer using a z-score method and the ancestry informative markers in the Precision ID Ancestry Panel.",
         journal = "Forensic Science International: Genetics 44, 102154.",
         year = 2020)
  ) #  %>% bind_rows()

  ggg_ref <- function(x){
    tags$li(paste0(x$author, " (", x$year,")."), #tags$br(),
            tags$i(x$title), #tags$br(),
            x$journal)
  }

  bootstrapPage(
    tags$head(includeHTML("gtag.html")),
    useWaitress(),
    tags$style(appCSS),
    shinyjs::useShinyjs(),  # Set up shinyjs
    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Genogeographer</a>'), id="nav",
               windowTitle = "Genogeographer",

               ## ANALYSIS TAB
               tabPanel("Analyse AIMs profile",
                        sidebarLayout(
                          sidebarPanel =
                            sidebarPanel(width = 3,
                                         h4("Input file"),
                                         fileInput('profile_file', 'Choose CSV File', width = "100%",
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx', '.xls')),
                                         uiOutput("column_panel"),
                                         ##
                                         h4("Settings"),
                                         radioButtons(inputId = "meta", label = "Grouping type:", choices = c("Meta-populations" = "meta", "Populations" = "pop")),
                                         sliderInput(inputId = "min_n", label = "Minimum sample size:", min = 5, max = 200, step = 5, value = 75),
                                         checkboxGroupInput(inputId = "admix", label = "Analyse 1st order admixture:",
                                                            choiceNames = list("Admixture (may take some time)"), choiceValues = list("admix")),
                                         shinyWidgets::sliderTextInput(inputId = "CI", label = "Confidence level:", choices = c(95, 97.5, 99, 99.9, 99.99), selected = 95, post = "%"),
                                         uiOutput("dbs"),
                                         checkboxGroupInput(inputId = "tilt", label = "Adjust p-values by exponential tilting:",
                                                            choiceNames = list("Adjust (may take some time)"), choiceValues = list("adjust")),
                                         uiOutput("side_pvalue"),
                                         uiOutput("LR_select"),
                                         tags$style(type='text/css', "#analyse, #reset, #report_download { width:100%; margin-top: 25px;}"),
                                         withBusyIndicatorUI(
                                           actionButton(inputId = "analyse", label = "Analyse!", icon = icon("calculator"), class = "btn-primary")
                                         ),
                                         actionButton(inputId = "reset", label = "Reset", icon = icon("trash")),
                                         tags$hr(),
                                         uiOutput("report_panel"),
                                         div(helpText(paste0("Version: ", ggg_package, " (", packageVersion(ggg_package),")"))),
                                         div(helpText(paste0("Developer: ", packageDescription(ggg_package, fields = "Maintainer"))))
                            ),
                          mainPanel = mainPanel(width = 9, uiOutput("analysis"))
                        )
               ),

               ## DATA TAB
               tabPanel("Data",
                        img(src = "meta_structure.png"),
                        img(src = "meta_pca.png"),
               ),

               ## CONSTRUCT OWN DATA
               tabPanel("Add reference populations",
                        sidebarLayout(
                          sidebarPanel =
                            sidebarPanel(width = 3,
                                         h3("Dataset file"),
                                         fileInput('dataset_file', 'Choose CSV or Excel File', width = "100%",
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx', '.xls')),
                                         uiOutput("dataset_panel"),
                                         h3("Information file"),
                                         fileInput('info_file', 'Choose CSV or Excel File', width = "100%",
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx', '.xls')),
                                         uiOutput("info_panel"),
                            ),
                          mainPanel = mainPanel(width = 9, uiOutput("upload_data"))
                        )
               ),

               ## ABOUT TAB
               tabPanel("About",
                        fluidPage(
                          p("The Genogeograher is an implementation of the methodology described by Tvedebrink et al. (2017, 2018, 2019)
                          and benchmarked in Mogensen et al. (2020). The aim is to classify an individual into a list of reference
                          populations based on the individual's genotype. To this purpose ancestry informative markers (AIMs) are used,
                          which (typically) are biallelic SNPs with pronounced observed variation across geography, ethnicity or culture"),
                          p("The genogeographer methodology is similar to an outlier test, where a genotype is tested for being an outlier
                          in each of the references populations. Hence, an individual can be rejected in ",tags$i("all"), "populations or
                          accepted in one or more populations. The key point is that a genotyped profile can be rejected in all
                          populations if it is too different from the typed reference populations. In an ordinary classification approach,
                          the genotyped profile would be assigned to the least unlikely (or most probable population). That is, the
                          reference populations are exclusive, but not exhaustive, which may cause the conclusions to be wrong."),
                          h4("Version"),
                          p("This Genogeograher online app is implemented in R with a frontend in Shiny.
                            Care is taken in the implemenation, but the app comes with absolutely no warranty."),
                          p("The app is implemented by Torben Tvedebrink <",
                            tags$a("genogeographer@tvedebrink.dk", href = "mailto:genogeographer@tvedebrink.dk", .noWS = "outside"), ">."),
                          p("The app is based on ", tags$b(paste("genogeographer version", packageVersion(ggg_package))),
                            "and associated reference databases."),
                          h4("References"),
                          tags$ul(ggg_references %>% lapply(ggg_ref))
                        )
               )
    )
  )
}

