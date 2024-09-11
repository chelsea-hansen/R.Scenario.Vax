#' Launch Shiny App to display scenario projections (function)
#'
#' This function launches the Shiny App included in the package.
#' The Shiny App allows users to upload the saved results from the package functions and display scenario projections.
#' In order for the Shiny App to work, make sure that one of the scenarios is called "Counterfactual" and that projection_intervals = TRUE
#' @export
#' @import stringr
#' @import ggplot2
display_scenarios <- function() {
  app_dir <- system.file("shiny", package = "R.Scenario.Vax")
  if (app_dir == "") {
    stop("Could not find the Shiny app directory. Try re-installing `R.Scenario.Vax`.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
