#' User guide tab - UI
#'
#' Static tab linking to the PDF user guide shipped in `inst/www`.
#' The `ggg_docs` resource path is registered in `ui_api()`.
#' @noRd
mod_user_guide_ui <- function(){
  nav_panel(
    title = "User guide",
    card(
      card_header("User guide"),
      p("A user guide is available ",
        tags$a(href = "ggg_docs/GenoGeographer_UserGuide.pdf", target = "_blank", "here"),
        ".")
    )
  )
}
