# Packages
library(spatstat)
library(fields)
library(FastGP)

# Simulate Covariate
win <- owin(xrange = c(0, 10), yrange = c(0, 10))
grid_res <- 30
x_seq <- seq(win$xrange[1], win$xrange[2], length.out = grid_res)
y_seq <- seq(win$yrange[1], win$yrange[2], length.out = grid_res)
coords <- as.matrix(expand.grid(x_seq, y_seq))
dists <- as.matrix(dist(coords))
n <- nrow(coords)

S <- 0.2 * exp(-dists/1.5)

mu <- rep(0, n)

covariate <- rcpp_rmvnorm(1,S,mu)

cov_field <- matrix(covariate, nrow = grid_res, ncol = grid_res)

cov_im <- im(mat = t(cov_field), xcol = x_seq, yrow = y_seq)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated Exponential Covariate")

# Simulate NHPP
win <- owin(x = c(0,10), y = c(0,10))

b_0 <- 2
b_1 <- 0.5

lambda <- exp(b_0 + b_1*(cov_im))

nhpp_sim <- rpoispp(lambda)

plot(nhpp_sim)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
  points(nhpp_sim)

# MCMC 
  
# Creating dataframe
fit_model <- function(points, cov_im){
  # Observed Points
  obs_point_df <- as.data.frame(points)
  obs_point_df$pt <- 1
  obs_point_df$wt <- 1e-6
  
  #Quadrature Points
  quad_points <- quadscheme(points, method = "grid")
  
  quad_df <- data.frame(
    x = quad_points$dummy$x, 
    y = quad_points$dummy$y,
    pt = 0,
    wt = area(win)/(quad_points$dummy$n)
  )
  
  # Merging and interpolating covariate
  data_merge <- rbind(obs_point_df, quad_df)
  
  data_merge$covariate <- interp.im(cov_im, data_merge$x, data_merge$y)
  
  return(data_merge)
}

update_betas <- function(data, params, priors, proposal_sd){
 
  beta_curr <- params$beta
  beta_prop <- rnorm(length(beta_curr), mean = beta_curr, sd = proposal_sd)

  #Log-Likelihood Function
  loglike <- function(beta) {
    log_lambda <- beta[1] + beta[2] * data$covariate
    lambda <- exp(log_lambda)
    
    likelihood <- sum(data$pt * log_lambda - lambda * data$wt)
    return(likelihood)
  }
  
  #Log-Likelihood for Prior
  logprior <- function(beta) {
    sum(dnorm(beta, mean = priors$mean, sd = priors$sd, log = TRUE))
  }
  
  # Posteriors
  post_curr <- loglike(beta_curr) + logprior(beta_curr)
  post_prop <- loglike(beta_prop) + logprior(beta_prop)
  
  # Metropolis Hastings Ratio
  log_acceptance_ratio <- post_prop - post_curr
  
  if(log(runif(1)) < log_acceptance_ratio) {
    params$beta <- beta_prop
  }
  
  return(params)
}

# Wrapper function
driver <- function(data, params, priors,  proposal_sd, iters){
  out=list()
  out$data=data
  out$params=params
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(params$beta),ncol = iters)
  
  for(k in 1:iters){
    params=update_betas(data,params, priors, proposal_sd)
    out$beta[,k]=params$beta 
  }
  return(out)
}

# Attempting

data <- fit_model(nhpp_sim, cov_im)
params <- list(beta = c(0,0))
proposal_sd <- c(0.05, 0.05)
priors <- list(mean = c(0,0), sd = c(10,10))
iters <- 5000

sim <- driver(data, params, priors, proposal_sd, iters)

apply(sim$beta, 1, mean)
apply(sim$beta, 1, sd)
