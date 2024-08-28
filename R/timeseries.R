#' Time series of RSV hospitalizations
#'
#' Weekly RSV hospitalizations for 12 states participating in RSVNet
#' Data from May 2020 through April 2021 have been removed.
#' Date for data prior to May 2020 have been increased by 1 year to create a continuous dataset
#' Ex: (May1,2020  = May1,2021)
#'
#'
#' @format ## `timeseries`
#' A data frame with 4,053 rows and 3 columns:
#' \describe{
#'   \item{date}{MMWR week}
#'   \item{state}{state name}
#'   \item{value}{estimated number of RSV hospitalizations}
#'   ...
#' }
#' @source <https://data.cdc.gov/Public-Health-Surveillance/Rates-of-Laboratory-Confirmed-RSV-COVID-19-and-Flu/kvib-3txy/about_data>
"timeseries"
