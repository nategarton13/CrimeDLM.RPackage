#' Function to compute the estimated seasonal time series components.
#'
#' \code{get_seasonality} computes the seasonal component in the observation equation of the program narrative.
#'
#' @param states A three dimensional array of sampled states such as that given by \code{get_states}.
#' @param harmonics An integer in 1 to 6 specifiying the number of harmonics to model seasonality.
#'  Must be consistent with specification in \code{run_mcmc()}.
#' @param crime_types Character vector specifying the types of crimes to include in the analysis.
#' The function requires at least two crime types.
#'
#' @return Returns a list with two objects: seasonal_effects and crime_types.
#' seasonal_effects is a three dimensional array where dimension one corresponds to the sample/iteration,
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
#'  seasonality <- get_seasonality(states = state_samples$state_samples)
#'  }


get_seasonality <- function(states, harmonics = 4, crime_types = c("burglary","robbery"))
{
  ## Use dlm package and the appropriate portion of the FF matrix to get the sum of appropriate seasonal parameters
  p <- length(crime_types)
  q <- harmonics
  period <- 12
  mod.uni <- dlm::dlmModPoly(order = 2) +
    dlm::dlmModTrig(s = period, q = q)
  mod.multi <- mod.uni
  mod.multi$FF <- mod.uni$FF %x% diag(p)

  get_seas_fun <- function(FF, states, p)
  {
    FF[1:p, 1:(p*2)] = 0
    return(states %*% t(FF))
  }

  seasonal_effects <- array(dim = c(dim(states)[1:2], p))
  for(i in 1:dim(states)[1])
  {
    seasonal_effects[i,,] <- get_seas_fun(FF = mod.multi$FF, states = states[i,,], p = p)
  }

  return(list("seasonal_effects" = seasonal_effects, "crime_types" = crime_types))
}
