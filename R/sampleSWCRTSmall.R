#' @name sampleSWCRTSmall
#' @title simulated small SW-CRT data
#' @description
#'   Simulated cross-sectional individual-level SW-CRT data with 12 clusters and 4 periods. The cluster-period size is uniformly distributed between 5 and 10.
#'   The correlated binary and continous outcomes are used for analysis as examples.
#' @docType data
#' @format A data frame with 373 rows and 9 variables:
#' \describe{
#'   \item{period1}{indicator of being at period 1}
#'   \item{period2}{indicator of being at period 2}
#'   \item{period3}{indicator of being at period 3}
#'   \item{period4}{indicator of being at period 4}
#'   \item{treatment}{indicator of being treated}
#'   \item{id}{cluster identification number}
#'   \item{period}{period order number}
#'   \item{y_bin}{binary outcome variable}
#'   \item{y_con}{continous outcome variable}
#' }
NULL
