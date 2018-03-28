#' Function to compute the estimated smoothed time series means.
#' 
#' \code{get_smoothmeans} computes the time series means absent of seasonal effects.
#' 
#' @param states A three dimensional array of sampled states such as that given by \code{get_states}.
#' @param crime_types Character vector specifying the types of crimes to include in the analysis.
#' The function requires at least two crime types.
#' 
#' @return Returns a list with two objects: smoothmeans and crime_types.
#' smoothmeans is a three dimensional array where dimension one corresponds to the sample/iteration,
#' dimension two corresponds to the number of time points, and dimension three corresponds to the crime type.
#' crime_types returns the crime_types argument specified in the function. This serves as a reminder of
#' the ordering of crime types, which is important for interpreting the smoothmeans output 
#' as well as making sure that the order of crime types is consistent between other functions.
#' 
#' @export
#' @examples
#'  \dontrun{
#'  cov_samples <- run_mcmc(data = chicago, chains = 2, adapt_delta = 0.8)
#'  state_samples <- get_states(mcmc_samples = cov_samples$samples, data = chicago)
#'  smoothmeans <- get_smoothmeans(states = state_samples$state_samples)
#'  }


get_smoothmeans <- function(states, crime_types = c("burglary","robbery"))
{
  ## this is just pulling off the first p parameters
  p <- length(crime_types)
  return(list("smoothmeans" = states[,,1:p], "crime_types" = crime_types))
}
