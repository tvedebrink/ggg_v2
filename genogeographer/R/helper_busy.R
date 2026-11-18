# All the code in this file needs to be copied to your Shiny app, and you need
# to call `withBusyIndicatorUI()` and `withBusyIndicatorServer()` in your app.
# You can also include the `appCSS` in your UI, as the example app shows.

# =============================================

# Set up a button to have an animated loading indicator and a checkmark
# for better user experience
# Need to use with the corresponding `withBusyIndicator` server function
withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      hidden(
        icon("rotate", class = "btn-loading-indicator", verify_fa = FALSE),
        icon("check", class = "btn-done-indicator", verify_fa = FALSE)
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("circle-exclamation", verify_fa = FALSE),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked. `session` is required so the
# CSS selectors below can be resolved to the fully namespaced id when used
# inside a shiny module (selector-based shinyjs calls are not auto-namespaced).
withBusyIndicatorServer <- function(buttonId, session, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  ns_id <- session$ns(buttonId)
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", ns_id)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", ns_id)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", ns_id)
  disable(id = ns_id, asis = TRUE)
  show(selector = loadingEl)
  hide(selector = doneEl)
  hide(selector = errEl)
  on.exit({
    enable(id = ns_id, asis = TRUE)
    hide(selector = loadingEl)
  })

  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    show(selector = doneEl)
    delay(2000, hide(selector = doneEl, anim = TRUE, animType = "fade", time = 0.5))
    value
  }, error = function(err) { errorFunc(err, ns_id) })
}

# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  html(html = errMessage, selector = errElMsg)
  show(selector = errEl, anim = TRUE, animType = "fade")
}

appCSS <- "
.btn-loading-container {
  margin-left: 10px;
  font-size: 1.2em;
}
.btn-done-indicator {
  color: green;
}
.btn-err {
  margin-top: 10px;
  color: red;
}
.responsive {
  width: 100%;
  max-width: 800px;
  height: auto;
}
"
