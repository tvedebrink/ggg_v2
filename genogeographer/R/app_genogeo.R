# library(shiny)
# library(shinythemes)
# library(shinycssloaders)
# library(tidyverse)
# library(DT)
# library(shinyjs)
# library(glue)
# library(plotly)
# library(ggforce)
# library(patchwork)
# library(shinyWidgets)
# library(parallel)
#
# library(leaflet)
# library(leaflet.providers)
# library(waiter)

ggg_package <- "genogeographer.devel"

#' Shiny application for GenoGeoGrapher
#' @param db_list A named list of databases of reference populations.
#' Each component is expected to be returned from \code{make_x1}.
#' @param reporting_panel Logical. Should report generate and download be available after sample analysis.
#' @param mh_data A list with elements `probs`, `snp_map`, `out_metas` and `ci_data` holding the precomputed
#' microhaplotype (MH) genotype probabilities/deviances, MH-to-SNP map, per-meta-population
#' outlier null models, and reference haplotype data for the log-likelihood confidence intervals, respectively
#' (see `precompute_likelihood_deviance.R`, `fit_structures.R` and `outlier_detect.R` in the 'microhaps'
#' analysis project, and `data-raw/MH_ci_data.R`). If `NULL` (the default), the data
#' bundled with the package (`inst/extdata/MH_profiles_probs.Rds`, `MH_snp_map.Rds`,
#' `MH_out_metas.Rds` and `MH_ci_data.Rds`) is used.
#' @export
app_genogeo <- function(db_list = NULL, reporting_panel = TRUE, mh_data = NULL){
  ui <- ui_api()
  server <- function(input, output, session) server_api(input = input, output = output, session = session)
  shinyApp(ui, server)
}

# source("R/app_server.R")
# source("R/app_ui.R")
# app_genogeo()
