#' Function to compute the estimated time series means.
#' 
#' \code{get_means} computes the linear combination of parameters defined to be the data model mean for a given set of sampled states.
#' 
#' @param states A three dimensional array of sampled states such as that given by \code{get_states}.
#' @param harmonics An integer in 1 to 6 specifiying the number of harmonics to model seasonality. 
#'  Must be consistent with specification in \code{run_mcmc()}.
#' @param crime_types Character vector specifying the types of crimes to include in the analysis.
#' The function requires at least two crime types.
#' 
#' @return Returns a list with two objects: ts_means and crime_types. 
#' ts_means is a three dimensional array where dimension one corresponds to the sample/iteration,
#' dimension two corresponds to the number of time points, and dimension three corresponds to the crime type.
#' crime_types returns the crime_types argument specified in the function. This serves as a reminder of
#' the ordering of crime types, which is important for interpreting the ts_means output 
#' as well as making sure that the order of crime types is consistent between other functions.
#' 
#' @export
#' @examples
#'  \dontrun{
#'  cov_samples <- run_mcmc(data = chicago, chains = 2, adapt_delta = 0.8)
#'  state_samples <- get_states(mcmc_samples = cov_samples$samples, data = chicago)
#'  ts_means <- get_means(states = state_samples$state_samples)
#'  }

get_means <- function(states, harmonics = 4, crime_types = c("burglary","robbery"))
{
  p <- length(crime_types)
  period <- 12
  q <- harmonics
  mod.uni <- dlm::dlmModPoly(order = 2) +
    dlm::dlmModTrig(s = period, q = q)
  
  ## define multivariate model
  mod.multi <- mod.uni
  mod.multi$FF <- mod.uni$FF %x% diag(p)
  mod.multi$GG <- mod.uni$GG %x% diag(p)
  
  ## Create function to get time series means for a given set of sampled states
  ymean.fun <- function(FF,states)
  {
    return(states %*% t(FF))
  }
  
  ## get sampled time series means
  y.mean.samples.stan <- array(dim = c(dim(states)[c(1,2)],p))
  for(i in 1:dim(states)[1])
  {
    y.mean.samples.stan[i,,] <- ymean.fun(FF = mod.multi$FF, states = states[i,,])
  }
  
  dimnames(x = y.mean.samples.stan)[[3]] <- crime_types
  return(list("ts_means" = y.mean.samples.stan, "crime_types" = crime_types))
  
  ## returns array with 
    # dim 1 = mcmc samples
    # dim 2 = time (month)
    # dim 3 = crime type
  
}
