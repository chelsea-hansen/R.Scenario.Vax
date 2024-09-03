#' Fitted parameter values (list)
#'
#' Sample output from the fit_model() function
#'
#'
#'
#' @format ## `fitLL`
#' A list with 9 items:
#' \describe{
#'   \item{list item 1}{a numeric value representing the fitted baseline transmission rate.}
#'   \item{list item 2}{a numeric value representing the fitted amplitude of seasonal forcing.}
#'   \item{list item 3}{a numeric value representing the fitted phase of seasonal forcing.}
#'   \item{list item 4}{a numeric value representing the fitted value for the proportion of RSV infections leading to reported hospitalizations in infants <2 months.}
#'   \item{list item 5}{a numeric value representing the fitted value for the proportion of RSV infections leading to reported hospitalizations in children 1-4 years.}
#'   \item{list item 6}{a numeric value representing the fitted value for the proportion of RSV infections leading to reported hospitalizations in ages 5-64 years.}
#'   \item{list item 7}{a numeric value representing the fitted value for the proportion of RSV infections leading to reported hospitalizations in adults 65 - 74 years.}
#'   \item{list item 8}{a numeric value representing the fitted value for the proportion of RSV infections leading to reported hospitalizations in adults 75+ years.}
#'   \item{list item 9}{a matix with 3 columns and 100 rows with the parameters values (baseline transmission +/- 5%, amplitude +/- 5%, phase +/- 1.5%) resampled using latin hypercube sampling.}
#'   ...
#' }
#'
"fitLL"
