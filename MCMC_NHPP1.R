# Packages
library(spatstat)
library(fields)
library(FastGP)

# Simulate Covariate
win <- owin(xrange = c(0, 10), yrange = c(0, 10))

grid_res <- 30

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)

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

b_0 <- 1
b_1 <- 3

lambda <- exp(b_0 + b_1*(cov_im))

nhpp_sim <- rpoispp(lambda)

plot(nhpp_sim)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
  points(nhpp_sim)

## MCMC 
  
# Creating X and X_grid data frames
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")

X <- as.data.frame(nhpp_sim)

X_grid$covariate <- interp.im(cov_im, X_grid$x, X_grid$y)

X$covariate <- interp.im(cov_im, X$x, X$y)

#Log-Likelihood Function
loglike <- function(beta, X, X_grid, cell_area) {
  log_lambda_points <- beta[1] + beta[2] * X$covariate
  term1 <- sum(log_lambda_points)
  
  log_lambda_grid <- beta[1] + beta[2] * X_grid$covariate
  lambda_grid <- exp(log_lambda_grid)
  term2 <- sum(lambda_grid * cell_area)
  
  likelihood <- sum(term1 - term2)
  return(likelihood)
}

#Log-Likelihood for Prior
logprior <- function(beta) {
  sum(dnorm(beta, mean = priors$mean, sd = priors$sd, log = TRUE))
}

update_betas <- function(params, priors, proposal_sd, X, X_grid, cell_area){
 
  beta_curr <- params$beta
  beta_prop <- rnorm(length(beta_curr), mean = beta_curr, sd = proposal_sd)
  
  # Posteriors
  post_curr <- loglike(beta_curr, X, X_grid, cell_area) + logprior(beta_curr) #change to top
  post_prop <- loglike(beta_prop, X, X_grid, cell_area) + logprior(beta_prop) #change to bottom
  
  # Metropolis Hastings Ratio
  log_acceptance_ratio <- post_prop - post_curr
  
  if(log(runif(1)) < log_acceptance_ratio) {
    params$beta <- beta_prop
  }
  
  return(params)
}

# Wrapper function
driver <- function(X, X_grid, params, cell_area, priors,  proposal_sd, iters){
  out=list()
  out$params=params
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(params$beta),ncol = iters)
  
  for(k in 1:iters){
    params=update_betas(params, priors, proposal_sd, X, X_grid, cell_area)
    out$beta[,k]=params$beta 
  }
  return(out)
}

# Attempting

params <- list(beta = c(0,0))
cell_area <- (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res)
proposal_sd <- c(0.05, 0.05)
priors <- list(mean = c(0,0), sd = c(10,10))
iters <- 7000

sim <- driver(X,X_grid, params, cell_area,priors, proposal_sd, iters)

apply(sim$beta, 1, mean)
apply(sim$beta, 1, sd)


# Adding another covariate-------------------------------------------------------

# Simulate Covariate
win <- owin(xrange = c(0, 10), yrange = c(0, 10))

grid_res <- 30

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)

coords <- as.matrix(expand.grid(x_seq, y_seq))

dists <- as.matrix(dist(coords))
n <- nrow(coords)

S_1 <- 0.2 * exp(-dists/1.5)
mu_1 <- rep(0, n)
covariate_1 <- rcpp_rmvnorm(1,S_1,mu_1)

S_2 <- 1 * exp(-dists/3)
mu_2 <- rep(0, n)
covariate_2 <- rcpp_rmvnorm(1,S_2,mu_2)

cov_field_1 <- matrix(covariate_1, nrow = grid_res, ncol = grid_res)
cov_field_2 <- matrix(covariate_2, nrow = grid_res, ncol = grid_res)

cov_im_1 <- im(mat = t(cov_field_1), xcol = x_seq, yrow = y_seq)
cov_im_2 <- im(mat = t(cov_field_2), xcol = x_seq, yrow = y_seq)

image.plot(x_seq, y_seq, cov_field_1, col = terrain.colors(100), main = "Covariate 1")
image.plot(x_seq, y_seq, cov_field_2, col = terrain.colors(100), main = "Covariate 2")

# Simulate NHPP
win <- owin(x = c(0,10), y = c(0,10))

b_0 <- 1
b_1 <- 3
b_2 <- 2

lambda <- exp(b_0 + b_1*(cov_im_1) + b_2*(cov_im_2))

nhpp_sim <- rpoispp(lambda)

plot(nhpp_sim)

## MCMC 

# Creating X and X_grid data frames
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")

X <- as.data.frame(nhpp_sim)

X_grid$covariate_1 <- interp.im(cov_im_1, X_grid$x, X_grid$y)
X_grid$covariate_2 <- interp.im(cov_im_2, X_grid$x, X_grid$y)

X$covariate_1 <- interp.im(cov_im_1, X$x, X$y)
X$covariate_2 <- interp.im(cov_im_2, X$x, X$y)

#Log-Likelihood Function
loglike <- function(beta, X, X_grid, cell_area) {
  log_lambda_points <- beta[1] + beta[2] * X$covariate_1 + beta[3] * X$covariate_2
  term1 <- sum(log_lambda_points)
  
  log_lambda_grid <- beta[1] + beta[2] * X_grid$covariate_1 + beta[3] * X_grid$covariate_2
  lambda_grid <- exp(log_lambda_grid)
  term2 <- sum(lambda_grid * cell_area)
  
  likelihood <- sum(term1 - term2)
  return(likelihood)
}

#Log-Likelihood for Prior
logprior <- function(beta) {
  sum(dnorm(beta, mean = priors$mean, sd = priors$sd, log = TRUE))
}

update_betas <- function(params, priors, proposal_sd, X, X_grid, cell_area){
  
  beta_top <- params$beta
  beta_bottom <- rnorm(length(beta_top), mean = beta_top, sd = proposal_sd)
  
  # Posteriors
  post_top <- loglike(beta_top, X, X_grid, cell_area) + logprior(beta_top) 
  post_bottom <- loglike(beta_bottom, X, X_grid, cell_area) + logprior(beta_bottom) 
  
  # Metropolis Hastings Ratio
  log_acceptance_ratio <- post_bottom - post_top
  
  if(log(runif(1)) < log_acceptance_ratio) {
    params$beta <- beta_bottom
  }
  
  return(params)
}

# Wrapper function
driver <- function(X, X_grid, params, cell_area, priors,  proposal_sd, iters){
  out=list()
  out$params=params
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(params$beta),ncol = iters)
  
  for(k in 1:iters){
    params=update_betas(params, priors, proposal_sd, X, X_grid, cell_area)
    out$beta[,k]=params$beta 
  }
  return(out)
}


# Attempting

params <- list(beta = c(0,0,0))
cell_area <- (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res)
proposal_sd <- c(0.05, 0.05, 0.05)
priors <- list(mean = c(0,0,0), sd = c(10,10,10))
iters <- 7000

sim <- driver(X,X_grid, params, cell_area,priors, proposal_sd, iters)

apply(sim$beta, 1, mean)
apply(sim$beta, 1, sd)