#' @title Compute PDF or CDF of delay
#' @description Compute PDF or CDF of lognormal delays 
#' 
#' @param lmu mean on the logscale
#' @param lsigma sd on the logscale
#' 
#' @return vector of either probability densities or probabilities
#' 
computeDelay <- function(lmu, lsigma, nmax, dt = 1, start = dt/2, pdf = T) {
  if (pdf) {
    res <- dlnorm(seq(start, nmax, by = dt), meanlog = lmu, sdlog = lsigma)
  } else {
    res <- plnorm(seq(start, nmax, by = dt), meanlog = lmu, sdlog = lsigma)
  }
  return(res)
}

#' @title Get time delay PDF or CDF
#' @description Get the PDF or CDF of the delay distributions assuming a log-normal
#' 
#' @param delay_params Parameters of the log-normal distribution
#' @param delay Name of the delay
#' @param nmax Max number of timepoints, here using timestep dt
#' @param dt timestep in days along which to compute the densities
#' @param start firts timepoint to use
#' @param pdf boolean, whether to return the PDF or the CDF
#' @param rnd Whether to draw distribution parameter values
#' 
#' @return vector of either probability or probability density values
#' 
getDelay <- function(delay_params = NULL, delay = NULL, nmax, dt = 1, start = dt/2, pdf = T, rnd = F) {
  
  if (!(delay %in% delay_params$delay))
    stop("Delay name not known")
  
  ind <- which(delay_params$delay == delay) 
  
  if (!rnd) {
    # Set lognorm parameters mu and sigma on the logscale
    lmu <- delay_params$logmu[ind]
    lsigma <- delay_params$logsigma[ind]
  } else {
    # Set lognorm parameters mu and sigma on the logscale
    lmu <- truncnorm::rtruncnorm(1, a = 0, mean = delay_params$logmu[ind], sd = delay_params$logmu.sd[ind])
    lsigma <- truncnorm::rtruncnorm(1, a = 0, mean = delay_params$logsigma[ind], sd = delay_params$logsigma.sd[ind])
  }
  
  res <- computeDelay(lmu, lsigma, nmax, dt = dt, start = start, pdf = pdf)
  return(res)
}

#' @title Set Delay PMF
#' @description Set the PMF of the delay distributions assuming a log-normal
#' 
#' @param delay_params Parameters of the log-normal distribution
#' @param delay Name of the delay
#' @param nmax Max number of timepoints, here using a 1 day dt
#' @param rnd Whether to draw distribution parameter values
#' 
#' @details The daily PMF values are computed by taking the difference between the
#' delay distribution's CDF at times t and t-1: PMF(t) = CDF(t+1) - CDF(t).
#' 
#' @return vector of probability values
#' 
setDelayPMF <- function(delay_params, delay, nmax, dt = 1, rnd = F) {
  
  if (!(delay %in% delay_params$delay))
    stop("Delay name not known")
  
  ntot <- nmax/dt + 1 # total number of timesteps
  
  cdf <- getDelay(delay_params, delay, nmax, start = 0, pdf = F, rnd = rnd)
  
  prob_values <- cdf[2:ntot] - cdf[1:(ntot-1)]
  return(prob_values)
}

#' @title Compute Istar
#' @description Function to compute the deconvolution of reported cases and the 
#' incubation period accounting for reporting delay. Returns an estimate of the
#' cumulative number of infections up to a proportionality constant (Istar).
#'
#' @param cases time series of daily cumulative reported cases
#' @param pdf_inc PMF of the incubation period
#' @param pdf_report PMF of the reporting delay
#' @param gamma the threshold value for bounding the inverse fourier transform values
#'
#' @return A vector of the values of Istar
#' 
computeIstar <- function(cases, pdf_inc, pdf_report,  gamma = .05) {
  
  nc <- length(cases)
  ni <- length(pdf_inc)
  nr <- length(pdf_report)
  
  # FFT convolution from `convolve` function in R
  # Pad pdfs
  ntot <- ni + nr - 1
  pdf_inc2 <- c(pdf_inc, rep.int(0, ntot - ni))
  pdf_report2 <- c(pdf_report, rep.int(0, ntot - nr))
  ntot <- length(pdf_report2)
  F_pdf_comb <- fft(pdf_inc2) * fft(pdf_report2)
  pdf_comb <- Re(fft(F_pdf_comb, inverse  = T))/ntot
  # Preserve sum(prob) = 1
  pdf_comb <- pdf_comb/sum(pdf_comb)
  
  # Pad cases
  ntot2 <- nc + ntot + - 1
  cases2 <- c(cases, rep.int(0, ntot2 - nc))
  eps <- 1e-10
  pdf_comb2 <- c(pdf_comb, rep(0, ntot2 - ntot))
  # fourier transform of convolution pdf
  F_cases <- fft(cases2)
  F_pdf_comb2 <- fft(pdf_comb2)
  
  # Water level regularization to prevent numerical instability
  # From https://cnx.org/resources/22c9f37591a06c51a0e632cc790ec83bcb853aa5/inverseFilter.m
  R <- F_pdf_comb2
  R1 <- F_pdf_comb2*(abs(F_pdf_comb2)>0)+1/gamma*(abs(F_pdf_comb2)==0)
  iR <- 1/F_pdf_comb2
  # invert the filter using threshold gamma
  G <- iR *((abs(R1)*gamma)>1)+gamma*abs(R1)*iR*((abs(R1)*gamma)<=1);
  Istar <- Re(fft(F_cases*G, inverse = T))[1:nc]
  # Sanity check for cumulative function
  Istar2 <- cummax(Istar)
  return(Istar2)
}

#' @title Compute Bayesian p-values
#' @description Compute the Bayesian p-value that two vectors of posterior draws 
#' have different means.
#' 
#' @param x vector of posterior draws of parameter to compare
#' @param y vector of posterior draws of reference parameter
#' @details the hypothesis that is tested is x = y using y - x = 0.
#' 
#' @return p-value
#' 
computePval <- function(x, y) {
  if (length(x) != length(y))
    stop("Vectors need to be the same length")
  
  return(min(sum((y-x) > 0), sum((y-x) < 0))*2/length(x))
}

#' @title Random initial values
#' @description Produces random initial values for Stan's HMC
#'
#' @return list with initial parameter values
#' 
rndInit <- function() {
  # Draw hyperprior parameters from priors 
  phi <- rbeta(1, 1, 6.5) # median of prior ~ 0.1
  lambda <- actuar::rpareto1(1, shape = 1.5, min = 0.1)
  # Draw IFR
  IFR <- runif(1, 1e-6, 1e-1)
  list(phi = phi, lambda = lambda, IFR = IFR)
}


#' @title Get SD
#' @description Computes the sd assuming a normal based on the 95% CI
#' 
#' @param mu mean
#' @param q025 0.025 quantile
#' @param q975 0.975 quantile
#' 
#' @return the standard deviation
#' 
getSD <- function(mu, q025, q975) {
  sd1 <- (mu - q025)/2
  sd2 <- (q975 - mu)/2
  return(mean(c(sd1, sd2)))
}

#' @title Log-norm mean
#' @description Computes the mean of a lognormal distribution
#' 
#' @param logmu mean on the logscale
#' @param logsigma sd on the logscale
#' 
#' @return the mean
#' 
lnormMean <- function(logmu, logsigma) {
  exp(logmu + .5 * logsigma^2)
}

#' @title Log-norm sd
#' @description Computes the sd of a lognormal distribution
#' @param logmu mean on the logscale
#' @param logsigma sd on the logscale
#' @return the sd
lnormSD <- function(logmu, logsigma) {
  sqrt((exp(logsigma^2) - 1) * exp(2 * logmu + logsigma^2))
}

#' @title Get quantiles
#' @description Computes the quantiles of a matrix of values
#' 
#' @param mat matrix over which to compute quantiles
#' 
#' @return dataframe of quantiles
#' 
getQuantiles <- function(mat) {
  mat %>% 
    as.data.frame() %>%
    mutate(sim = row_number()) %>%
    gather(var, value, -sim) %>%
    mutate(time = as.numeric(str_replace_all(var, "V", ""))) %>%
    group_by(sim) %>%
    arrange(time) %>%
    mutate(cumvalue = cumsum(value)) %>%
    group_by(time) %>%
    summarise(q025 = quantile(value, .025),
              q975 = quantile(value, .975),
              median = quantile(value, 0.5),
              mean = mean(value),
              cdf.q025 = quantile(cumvalue, .025),
              cdf.q975 = quantile(cumvalue, .975),
              cdf.median = quantile(cumvalue, 0.5),
              cdf.mean = mean(cumvalue))
}
