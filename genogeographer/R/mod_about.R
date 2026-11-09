#' About tab - UI
#'
#' Static content describing the GenoGeographer methodology and references.
#' No server-side logic is required for this panel.
#' @noRd
mod_about_ui <- function(){
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
  )

  ggg_ref <- function(x){
    tags$li(paste0(x$author, " (", x$year, ")."),
            tags$i(x$title),
            x$journal)
  }

  nav_panel(
    title = "About",
    card(
      p("The Genogeograher is an implementation of the methodology described by Tvedebrink et al. (2017, 2018, 2019)
      and benchmarked in Mogensen et al. (2020). The aim is to classify an individual into a list of reference
      populations based on the individual's genotype. To this purpose ancestry informative markers (AIMs) are used,
      which (typically) are biallelic SNPs with pronounced observed variation across geography, ethnicity or culture"),
      p("The genogeographer methodology is similar to an outlier test, where a genotype is tested for being an outlier
      in each of the references populations. Hence, an individual can be rejected in ", tags$i("all"), "populations or
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
}
