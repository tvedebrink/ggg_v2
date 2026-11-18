#' Reference databases module - UI
#'
#' Static information about the reference databases: population samples,
#' AIMs information, and the STRUCTURE/PCA analyses used to derive metapopulations.
#' @param id Module id
#' @noRd
mod_reference_databases_ui <- function(id){
  ns <- NS(id)
  nav_menu(
    title = "Reference Databases",
    nav_panel(
      title = "Population samples",
      HTML(paste(paste0("<p><b>", ggg_package, "</b>"), "comes with three reference databases genotyped with
      Applied Biosystems Precision ID Ancestry Panel, which contains 165 Ancestry Informative Markers (AIMs)
      SNPs (AISNPs) selected by the Kidd and Seldin labs, respectively. One locus (rs10954737) was excluded
      the databases as it is underrepresented in several of the publicly available population samples used
      to form the reference databases.</p>")),
      p("Below are tables of the included populations and the induced metapopulations. The latter groups of
      populations identified using STRUCTURE (see STRUCTURE analysis) and other population genetic analyses
      (see Mogensen et. al, 2020, for details)."),
      p("All 164 AIMs can be used for analyses, but also the so-called 'Kidd' (55 markers) and
      'Seldin' (122 markers) subsets of markers can be used. For all three sets of markers analses can be
      performed using the individual populations or metapopulations, and in both cases either as non-admixed
      or admixed (first-order admixture of two distinct (meta-)populations)."),
      p("See AIMs information for specific details on the AISNPs and their allelic variants."),
      uiOutput(ns("pop_info"))
    ),
    nav_panel(
      title = "AIMs information",
      p("The table below contains an overview of the AIMs SNPs included in the reference databases.
        For more information on the individual SNP a link to the National Library of Medicine' dbSNP is provided."),
      DTOutput(ns("allele_info"))
    ),
    nav_panel(
      title = "STRUCTURE analysis",
      p("Based on a STRUCTURE analysis (and other population genetic analyses, see Mogensen et. al, 2020,
      for details) the individual populations were grouped into metapopulations. The metapopulations
      were identified by similar STRUCTURE cluster membership components and by no or a limited number
      of significant allele frequency differences among the included populations."),
      p("The STRUCTURE plot below originates from Mogensen et. al (2020)."),
      img(src = "meta_structure.png", class = "responsive")
    ),
    nav_panel(
      title = "PCA analysis",
      p("The plot below shows the outcome of a principal components analysis of the genotypes included in the
        reference databases. The colours are determined by the metapopulations. An analysis similar to
        EIGENSTRAT identified that the first five eigenvalues were significant in order to infer population
        structure. Thus the pairwise combinations of the first five PCs are plotted below to show the
        population structure in the reference databases."),
      img(src = "meta_pca.png", class = "responsive")
    )
  )
}

#' Reference databases module - server
#'
#' @param id Module id
#' @param db_list A reactive returning the named list of reference databases.
#' @noRd
mod_reference_databases_server <- function(id, db_list){
  moduleServer(id, function(input, output, session){
    output$pop_info <- renderUI({
      req(db_list())
      pop_info <- attr(db_list()[[1]]$pop$db, "info")
      meta_info <- attr(db_list()[[1]]$meta$db, "info")
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Populations"), pop_info %>% select(population, n) %>% dt_table()),
        card(card_header("Metapopulations"), meta_info %>% select(metapopulation, n) %>% dt_table())
      )
    })

    output$allele_info <- renderDT({
      ggg_allele_list %>%
        mutate(locus = paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/", locus, "' target='_blank'>", locus, "</a>")) %>%
        rename('AIMs SNP' = locus, 'Reference allele' = x1, 'Alternative allele' = x2) %>%
        mutate(across(Kidd:Seldin, ~ifelse(.x, "Yes", "No"))) %>%
        dt_table(escape = FALSE)
    })
  })
}
