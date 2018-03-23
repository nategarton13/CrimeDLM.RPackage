#' Function to sample the latent states given sampled covariance matrices.
#' 
#' \code{get_states} uses the \code{dlm} package to sample the latent states via the
#' forward filtering-backward sampling algorithm as described in the research narrative.
#' The same arguments should be supplied for parameters shared by \code{get_states} and \code{run_mcmc}.
#' 
#' @param mcmc_samples A list of parameter samples with the same names and classes as those produced by \code{run_mcmc}.
#' @param data The full Chicago crime dataset as given by \code{chicago}.
#' @param harmonics An integer in 1 to 6 specifiying the number of harmonics to model seasonality.
#' @param crime_types Character vector specifying the types of crimes to include in the analysis.
#' The function requires at least two crime types.
#' @param initial_year Integer in 2007 to 2015 specifying the first year in the analysis.
#' @param final_year Integer in 2008 to 2016 specifying the last year in the analysis.
#' 
#' @return Returns a three dimensional array where dimension one corresponds to the sample/iteration,
#' dimension two corresponds to the number of time points, and dimension three corresponds to the parameter.
#' Parameters are in the following order: \eqn{\mu_1,\mu_2,...,\beta_1, \beta_2,...}. The ordering of a given set of parameters,
#' say \eqn{\mu}, is the ordering of the types of crimes specified by the \code{crime_types} argument. 
#' After \eqn{\beta} parameters, the remainder are needed for the construction of seasonal effects.
#' 
#' @export
#' @examples
#'  \dontrun{
#'  cov_samples <- run_mcmc(data = chicago, chains = 2, adapt_delta = 0.8)
#'  state_samples <- get_states(mcmc_samples = cov_samples$samples, data = chicago)
#'  }

get_states <- function(mcmc_samples, data, harmonics = 4, crime_types = c("burglary","robbery"), initial_year = 2012, final_year = 2016)
{
  ## function to take standard deviation diagonal matrix and correlation matrix and create covariance matrix
  LDL <- function(sdmat, cormat)
  {
    return(sdmat %*% cormat %*% sdmat)
  }
  
  ## get sampled covariance matrices for both the errors and evolutions
  sigma.beta.samples.stan <- array(dim = dim(mcmc_samples$omega_evo))
  sigma.samples.stan <- array(dim = dim(mcmc_samples$omega_error))
  for(i in 1:dim(sigma.beta.samples.stan)[1])
  {
    sigma.beta.samples.stan[i,,] <- LDL(sdmat = diag(mcmc_samples$sigma_evo[i,]), cormat = mcmc_samples$omega_evo[i,,])
    sigma.samples.stan[i,,] <- LDL(sdmat = diag(mcmc_samples$sigma_error[i,]), cormat = mcmc_samples$omega_error[i,,])
  }
  
  ## parameters for the internal get_states function
  q <- harmonics 
  ## modify the data by taking logs and turning it into a matrix
  y <- matrix(nrow = 12*(final_year - initial_year + 1), ncol = length(crime_types))
  for(i in 1:length(crime_types))
  {
    y[,i] <- log(data[data$type == crime_types[i] & data$year >= initial_year & data$year <= final_year,]$count)
  }
  p <- ncol(y)
  
  
  ## function to sample the states given both an evolution and an error covariance matrix
  sub_get_states <- function(sigma_error_samples, sigma_evo_samples, p = p, q = q, y = y)
  {
    states <- array(dim = c(dim(mcmc_samples[[1]])[1], nrow(y) + 1, ifelse(q < 6, p*(2 + 2*q), p*(2 + 2*q - 1))))
    period <- 12
    mod.uni <- dlm::dlmModPoly(order = 2) +
      dlm::dlmModTrig(s = period, q = q)
    
    ## define multivariate model
    mod.multi <- mod.uni
    mod.multi$FF <- mod.uni$FF %x% diag(p)
    mod.multi$GG <- mod.uni$GG %x% diag(p)
    m0 <- rep(0, times = ifelse(q < 6, p*(2 + 2*q), p*(2 + 2*q - 1)))
    C0 <- 1e7 * diag(ifelse(q < 6, p*(2 + 2*q), p*(2 + 2*q - 1)))
    mod.multi$m0 <- m0
    mod.multi$C0 <- C0
    
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    for(i in 1:dim(states)[1])
    {
      ## Parameter covariance matrix 
      W.mu <- matrix(0, nrow = p, ncol = p)
      W <- dlm::bdiag(W.mu, sigma.beta.samples.stan[i,,], diag(ifelse(q < 6,2*q, 2*q - 1)) %x% (0*diag(p)))
      mod.multi$W <- W
      
      ## Data covariance matrix
      mod.multi$V <- sigma.samples.stan[i,,]
      
      ## get filter distribution
      filter <- dlm::dlmFilter(y = y, mod = mod.multi)
      
      ## sample parameters
      states[i,,] <- dlm::dlmBSample(modFilt = filter)
      
      # update progress
      setTxtProgressBar(pb, 100*i/dim(states)[1])
    }
    return(states)
  }
  
  ## sample the states
  print("This may take a while... See progress bar.")
  state.samples <- sub_get_states(sigma_error_samples = sigma.samples.stan, sigma_evo_samples = sigma.beta.samples.stan, p = p,  q = q, y = y)
  return(state.samples)
  
  ## every 6 columns correspond to the same single parameter from each crime type
    ## first 6 columns are the smoothed means, mu
    ## second 6 columns are the linear trends, beta
    ## the remaining columns correspond to the harmonics and conjugate harmonics such that 
    ## the first six are the first six harmonics and the second six are the conjugate harmonics
    ## to the first six harmonics
    ## and so on...
}
