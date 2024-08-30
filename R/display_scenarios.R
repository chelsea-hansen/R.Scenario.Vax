#' Launch Shiny App
#'
#' This function launches the Shiny app included in the package.
#' @export
display_scenarios <- function() {
  app_dir <- system.file("shiny", package = "R.Scenario.Vax")
  if (app_dir == "") {
    stop("Could not find the Shiny app directory. Try re-installing `R.Scenario.Vax`.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
