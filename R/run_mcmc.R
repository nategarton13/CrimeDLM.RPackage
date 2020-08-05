#' Function to sample observation and state covariance parameters.
#'
#' \code{run_mcmc} uses Stan to run a Hamiltonian Monte Carlo algorithm to sample the observation
#' and latent state covariance parameters in the dynamic linear model specificied in the research narrative.
#' The full chicago crime dataset must be passed to the \code{data} argument. Other parameter arguments dictate
#' the subset of the data on which the analysis is performed as well as details of the MCMC to be passed to \code{rstan::sampling}.
#'
#' @param data The full Chicago crime dataset as given by \code{chicago}.
#' @param harmonics An integer in 1 to 6 specifiying the number of harmonics to model seasonality.
#' @param chains An integer specifying the number of MCMC chains for Stan to run.
#' @param iter An integer specifying the total number of samples for each chain, including warmup.
#' @param warmup An integer specifying the number of burn-in samples in Stan.
#' @param adapt_delta See Stan language manual at \url{http://mc-stan.org/users/documentation/}.
#' @param crime_types Character vector specifying the types of crimes to include in the analysis.
#' The function requires at least two crime types.
#' @param initial_year Integer in 2007 to 2015 specifying the first year in the analysis.
#' @param final_year Integer in 2008 to 2016 specifying the last year in the analysis.
#'
#' @return Returns a list of three objects: samples, summary, and crime_types. summary gives a summary of the
#' samples. crime_types returns the crime_types argument specified in the function. This serves as a reminder of
#'  the ordering of crime types, which is important for interpreting covariance matrix output
#'  as well as making sure that the order of crime types is consistent between other functions.
#'  samples is itself a list where each element is either a two or three dimensional array
#'  containing all of the samples of a given parameter. The first dimension of the array always
#'  corresponds to the sample/iteration. Parameter names work as follows:
#'  \code{omega_evo} corresponds to evolution correlation matrices,
#'  \code{omega_error} corresponds to error correlation matrices,
#'  \code{sigma_evo} corresponds to a vector of evolution standard deviations,
#'  \code{sigma_error} corresponds to a vector of the error standard deviations,
#'  \code{sigma_evo_mat} corresponds to evolution covariance matrices,
#'  \code{sigma_error_mat} corresponds to error covariance matrices.
#'  See Stan documentation for \code{lp__}.
#'
#' @export
#' @examples
#'  \dontrun{
#'  cov_samples <- run_mcmc(data = chicago, chains = 2, adapt_delta = 0.8)
#'  }
#'

## Function that takes in raw data, number of harmonics, and stan options and runs the entire MCMC from stan
run_mcmc <- function(data, harmonics = 4, chains = 3, iter = 2000, warmup = 1000, adapt_delta = 0.9, crime_types = c("burglary","robbery"), initial_year = 2012, final_year = 2016)
{
  ## load necessary packages
  #require(zoo)
  # require(rstan)
  # require(dlm)

  ## modify the data by taking logs and turning it into a matrix
  y <- matrix(nrow = 12*(final_year - initial_year + 1), ncol = length(crime_types))
  for(i in 1:length(crime_types))
  {
    y[,i] <- log(data[data$type == crime_types[i] & data$year >= initial_year & data$year <= final_year,]$count)
  }

  ## arguments to be passed to the stan model
  p <- ncol(y)
  q <- harmonics
  period <- 12
  mod.uni <- dlm::dlmModPoly(order = 2) +
    dlm::dlmModTrig(s = period, q = q)
  mod.multi <- mod.uni
  mod.multi$FF <- mod.uni$FF %x% diag(p)
  mod.multi$GG <- mod.uni$GG %x% diag(p)
  m0 <- rep(0, times = ifelse(q < 6, p*(2 + 2*q), p*(2 + 2*q - 1)))
  C0 <- 1e7 * diag(ifelse(q < 6, p*(2 + 2*q), p*(2 + 2*q - 1)))
  mod.multi$m0 <- m0
  mod.multi$C0 <- C0

  ## Stan model
  # n: number of observations
  # q: number of harmonics
  # p: number of crime types
  # y: matrix of data of dimension p x n
  # m0: mean of initial state vector
  # C0: covariance matrix of initial state vector
  # FF:
  model <- "
  data {
  int n;
  int q;
  int p;
  matrix[p, n] y;
  vector[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1)] m0;
  matrix[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1), (q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1)] C0;
  matrix[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1),p] FF;
  matrix[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1), (q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1)] GG;
  }
  transformed data{
  matrix[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1),p] F;
  matrix[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1), (q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1)] G;
  F = FF;
  G = GG;
  }
  parameters {
  corr_matrix[p] omega_evo;
  vector<lower=0>[p] sigma_evo;
  corr_matrix[p] omega_error;
  vector<lower=0>[p] sigma_error;
  }
  transformed parameters {
  cov_matrix[p] sigma_evo_mat;
  cov_matrix[p] sigma_error_mat;
  sigma_evo_mat = quad_form_diag(omega_evo, sigma_evo);
  sigma_error_mat = quad_form_diag(omega_error, sigma_error);
  }
  model {
  int counter;
  matrix[(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1),(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1)] W;
  W = rep_matrix(0,(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1),(q < 6) ? p*(2*q + 2) : p*(2*q + 2 - 1));


  sigma_evo ~ cauchy(0,1);
  sigma_error ~ cauchy(0,1);

  omega_evo ~ lkj_corr(1);
  omega_error ~ lkj_corr(1);

  counter = 0;
  for(i in 1:p)
  {
  for(j in 1:p)
  {
  W[i + p, j + p] = sigma_evo_mat[i,j];
  }
  }

  y ~ gaussian_dlm_obs(F, G, sigma_error_mat, W, m0, C0);
  }
  "
  ## compile Stan model
  m <- rstan::stan_model(model_code = model)

  ## sample from Stan model
  ## Notes:
  # sampling takes approximately 16 hours to run for 3 chains, 5000 iterations each for all crime types and time periods
  stan_samples = rstan::sampling(object = m, data = list(p = p, n=nrow(y), y=t(y), q = q, m0 = m0, C0 = C0, FF = t(mod.multi$FF), GG = mod.multi$GG), iter = iter, warmup = warmup, chains = chains, control = list(adapt_delta = adapt_delta))

  ## return extracted stan samples
  return(list("samples" = rstan::extract(stan_samples), "summary" = summary(stan_samples), "crime_types" = crime_types))
}

