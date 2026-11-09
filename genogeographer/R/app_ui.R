

#' @title UI part for Shiny app
#' @description The shiny front-end, composed of a bslib navbar and shiny modules.
#' @export
ui_api <- function(){

  ggg_theme <- bs_theme(font_scale = 0.85) #, spacer = "0.2rem", version = 5, bootswatch = "flatly") %>%
  #   bs_theme_update("navbar-padding-y" = "10px") %>%
  #   bs_add_rules(".wide_popover { max-width: none; }")

  page_navbar(
    title = "Genogeographer",
    id = "nav",
    window_title = "Genogeographer",
    theme = ggg_theme,
    header = tagList(
      useShinyjs(),
      useWaitress(),
      tags$style(appCSS)
    ),
    mod_analyse_ui("analyse"),
    mod_mh_analyse_ui("mh"),
    mod_reference_databases_ui("refdb"),
    mod_add_reference_ui("addref"),
    mod_about_ui()
  )
}

