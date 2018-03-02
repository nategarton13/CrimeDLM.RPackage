#' Function to perform analysis of crime data using dynamic linear models
#'
#' Add more detailed information here,
#' see http://r-pkgs.had.co.nz/man.html for details on using roxygen2 for
#' documentation.
#'
#' @param data a \code{data.frame} containing aggregated crime data
#' @param prior a list with elements ... containing a prior specification,
#'   see details
#' @param ... passed to \code{sampling} in rstan
#' @return a \code{data.frame} containing posterior samples for model parameters
#' @details
#' Provide model and prior specification here.
#'
#'
#' @export
#' @examples
#' \dontrun{
#' dlm_analysis(chicago, list(a=1))
#' }
#'
dlm_analysis <- function(data, prior, ...) {
  data.frame(a=1)
}
