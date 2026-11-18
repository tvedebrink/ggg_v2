#' User guide tab - UI
#'
#' The user guide rendered as HTML (`inst/www/userguide/user_guide.html`),
#' with links to the PDF version. The `ggg_docs` resource path is registered in `ui_api()`.
#' @noRd
mod_user_guide_ui <- function(){
  nav_panel(
    title = "User guide",
    div(
      class = "mx-auto", style = "max-width: 60rem;",
      card(
        fill = FALSE,
        card_header("User guide"),
        card_body(
          fill = FALSE,
          p("The user guide is also available as a PDF: ",
            tags$a(href = "ggg_docs/GenoGeographer_UserGuide.pdf", target = "_blank", rel = "noopener",
                   icon("file-pdf", verify_fa = FALSE), "open"),
            " or ",
            tags$a(href = "ggg_docs/GenoGeographer_UserGuide.pdf", download = "GenoGeographer_UserGuide.pdf",
                   icon("download", verify_fa = FALSE), "download as PDF"),
            ".")
        )
      ),
      card(
        fill = FALSE,
        card_body(
          fill = FALSE,
          includeHTML(system.file("www", "userguide", "user_guide.html", package = ggg_package))
        )
      )
    )
  )
}
