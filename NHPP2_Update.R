# Packages
library(spatstat)
library(fields)
library(FastGP)
library(MASS)

# True NHPP --------------------------------------------------------------------
# Simulate Covariate
win <- owin(xrange = c(0, 10), yrange = c(0, 10))

grid_res <- 10

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)

coords <- as.matrix(expand.grid(x_seq, y_seq))
grid_coords <- expand.grid(x = x_seq, y = y_seq)

dists <- as.matrix(dist(coords))
n <- nrow(coords)

S <- 0.2 * exp(-dists/1.5)

mu <- rep(0, n)

covariate <- as.vector(rcpp_rmvnorm(1,S,mu))

cov_field <- matrix(covariate, 
                    nrow = grid_res, 
                    ncol = grid_res,
                    byrow = TRUE)

cov_field_ppp <- ppp(x = grid_coords$x, 
                     y = grid_coords$y,
                     window = win, 
                     marks = covariate)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated Exponential Covariate")


# Simulate NHPP
b_0 <- 1
b_1 <- 3

lambda <- exp(b_0 + b_1*(cov_field))
lambda_im <- im(t(lambda), xcol = x_seq, yrow = y_seq)

nhpp_sim <- rpoispp(lambda_im)

# Discretize using spatstat


par(mfrow = c(1,2))
image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
nhpp_discretize <- pixellate(nhpp_sim, eps = 1)
points(nhpp_sim)
plot(nhpp_discretize)
par(mfrow = c(1,1))


# Source 1  --------------------------------------------------------------------

set.seed(111)
sigma_2 <- 0.01
nrow <- length(nhpp_discretize$yrow)
ncol <- length(nhpp_discretize$xcol)

mu_1 <- -0.5 * sigma_2^2
f <- rnorm(nrow * ncol, mean = mu_1, sd = sigma_2)
exp_f <- exp(f)

lambda_1 <- nhpp_discretize * exp_f

nhpp_1 <- rpoispp(lambda_1)

par(mfrow=(c(1,2)))
plot(nhpp_discretize)
plot(lambda_1)
par(mfrow=(c(1,1)))

par(mfrow=(c(1,2)))
plot(nhpp_sim)
plot(nhpp_1)
par(mfrow=(c(1,1)))

# Source 2 ----------------------------------------------------------------------

tau_2 <- 1
alpha <- 0.75

g <- rnorm(nrow * ncol, mean = alpha, sd = tau_2)
exp_g <- exp(g)

lambda_2 <- nhpp_discretize * exp_g

nhpp_2 <- rpoispp(lambda_2)

par(mfrow=(c(1,2)))
plot(nhpp_discretize)
plot(lambda_2)
par(mfrow=(c(1,1)))

par(mfrow=(c(1,2)))
plot(nhpp_sim)
plot(nhpp_2)
par(mfrow=(c(1,1)))

par(mfrow = (c(1,3)))
plot(nhpp_sim)
plot(nhpp_1)
plot(nhpp_2)
par(mfrow=(c(1,1)))

# Data frame creation ----------------------------------------------------------

## X grid
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")
X_grid$covariate <- covariate

## Source 1 (small var)
X_1 <- as.data.frame(nhpp_1)
nn_X_1 <- nncross(nhpp_1, cov_field_ppp)
X_1$covariate <- cov_field_ppp$marks[nn_X_1$which]

## Source 2 (large var)
X_2 <- as.data.frame(nhpp_2)
nn_X_2 <- nncross(nhpp_2, cov_field_ppp)
X_2$covariate <- cov_field_ppp$marks[nn_X_2$which]


par(mfrow = (c(1,2)))
quilt.plot(X_1$x, X_1$y, X_1$covariate)
plot(nhpp_1)
par(mfrow = (c(1,1)))


par(mfrow = (c(1,2)))
quilt.plot(X_2$x, X_2$y, X_2$covariate)
plot(nhpp_2)
par(mfrow = (c(1,1)))

quilt.plot(X_grid$x, X_grid$y, X_grid$covariate)


# MCMC --------------------------------------------------------------------------

## log likelihood function
loglike <- function(parameters, data) {
  
  log_lambda_points1 <- parameters$beta[1] + parameters$beta[2] * data$X_1$covariate + parameters$f[data$nn_index_1]
  term1 <- sum(log_lambda_points1)
  
  log_lambda_grid1 <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$f
  lambda_grid1 <- exp(log_lambda_grid1)
  term2 <- sum(lambda_grid1 * data$cell_area)
  
  log_lambda_points2 <- parameters$beta[1] + parameters$beta[2] * data$X_2$covariate + parameters$g[data$nn_index_2]
  term3 <- sum(log_lambda_points2)
  
  log_lambda_grid2 <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$g
  lambda_grid2 <- exp(log_lambda_grid2)
  term4 <- sum(lambda_grid2 * data$cell_area)
  
  likelihood <- (term1 - term2) + (term3 - term4)
  return(likelihood)
}

## Updating slope estimates using Metropolis Hastings MCMC
update_betas<- function(parameters, priors, data){
  
  beta_cand <- rnorm(length(parameters$beta), mean = parameters$beta, sd = priors$beta_prop_sd) 
  
  params_top <- list(beta = beta_cand, f = parameters$f, g = parameters$g)
  params_bottom <- list(beta = parameters$beta, f = parameters$f, g = parameters$g)
  
  # Posteriors
  post_top <- loglike(params_top, data) + sum(dnorm(beta_cand, mean = priors$beta_mean, sd = priors$beta_sd, log = TRUE))
  post_bottom <- loglike(params_bottom, data) + sum(dnorm(parameters$beta, mean = priors$beta_mean, sd = priors$beta_sd, log = TRUE))
  
  # Metropolis Hastings Ratio
  log_acceptance_ratio <- post_top - post_bottom
  
  if(log(runif(1)) < log_acceptance_ratio) {
    parameters$beta <- beta_cand
  }
  
  return(parameters)
}

## Updating f - Source 1 measurement error
update_f <- function(parameters, priors, data){
  
  # Choosing ellipse (nu) from prior (f)
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$f)), Sigma = diag(parameters$sigma_2, length(parameters$f))))
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike(parameters, data) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
  repeat {
    # Calculate f'
    f_prime <- as.vector(parameters$f*cos(theta) + nu*sin(theta))
    
    params_prime <- parameters
    params_prime$f <- f_prime
    
    # Shrinking bracket
    if(loglike(params_prime, data) > log_y){
      parameters$f <- f_prime
      return(parameters) 
    } else {
      if(theta < 0) {
        theta_min <- theta
      } else {
        theta_max <- theta
      }
      theta <- runif(1, theta_min, theta_max)
    }
  }
}

## Updating g - Source 2 measurement error
update_g <- function(parameters, priors, data){
  
  # Choosing ellipse (nu) from prior (g)
## Is this correct mean and standard deviation?
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$g)), Sigma = diag(parameters$tau_2, length(parameters$g))))
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike(parameters, data) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
  repeat {
    # Calculate g'
    g_prime <- as.vector(parameters$g*cos(theta) + nu*sin(theta))
    
    params_prime <- parameters
    params_prime$g <- g_prime
    
    # Shrinking bracket
    if(loglike(params_prime, data) > log_y){
      parameters$g <- g_prime
      return(parameters) 
    } else {
      if(theta < 0) {
        theta_min <- theta
      } else {
        theta_max <- theta
      }
      theta <- runif(1, theta_min, theta_max)
    }
  }
}

## Updating source 1 variance using Gibbs Sampling- sigma_2

update_sigma_2 <- function(parameters, priors, data){
  n <- length(parameters$f)
  
  alpha_post <- priors$a_0_sigma + n/2
  beta_post  <- priors$b_0_sigma  + 0.5 * sum((parameters$f - priors$f_mean)^2)
  
  # Draw samples from Inverse-Gamma
  parameters$sigma_2 <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)
  
  return(parameters)
  
}

## Updating source 2 variance using Gibbs Sampling - tau_2

update_tau_2 <- function(parameters, priors, data){
  n <- length(parameters$g)
  
  alpha_post <- priors$a_0_tau + n/2
  beta_post  <- priors$b_0_tau  + 0.5 * sum((parameters$g - parameters$alpha)^2)
  
  # Draw samples from Inverse-Gamma
  parameters$tau_2 <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)
  
  return(parameters)
  
}

## Updating source 2 mean using Gibbs Sampling - alpha
  
update_alpha <- function(parameters, priors, data){
  n <- length(parameters$g)
  
  denom <- (n / parameters$tau_2) + (1 / priors$phi)
  mu_post <- ( sum(parameters$g) / parameters$tau_2 + priors$gamma / priors$phi ) / denom
  sigma_post <- sqrt(1 / denom)
  
  parameters$alpha <- rnorm(1, mean = mu_post, sd = sigma_post)
  
  return(parameters)
}

## Wrapper function

driver <- function(parameters, priors, data, iters){
  out=list()
  out$params=parameters
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(parameters$beta),ncol = iters)
  out$f=matrix(NA, nrow = length(parameters$f), ncol = iters)
  out$g=matrix(NA, nrow = length(parameters$g), ncol = iters) 
  out$sigma_2=matrix(NA, nrow = 1, ncol = iters)
  out$tau_2=matrix(NA, nrow = 1, ncol = iters)
  out$alpha=matrix(NA, nrow = 1, ncol = iters)
  
  
  for(k in 1:iters){
    parameters <- update_betas(parameters, priors, data)
    out$beta[,k]=parameters$beta 
    
    parameters <- update_f(parameters, priors, data)
    out$f[,k] <- parameters$f
    
    parameters <- update_g(parameters, priors, data)
    out$g[,k] <- parameters$g
    
    parameters <- update_sigma_2(parameters, priors, data)
    out$sigma_2[,k] <- parameters$sigma_2
    
    parameters <- update_alpha(parameters, priors, data)
    out$alpha[,k] <- parameters$alpha
    
    parameters <- update_tau_2(parameters, priors, data)
    out$tau_2[,k] <- parameters$tau_2
  }
  return(out)
}

# Attempting --------------------------------------------------------------------


data <- list(X_grid = X_grid, 
             X_1 = X_1, 
             X_2 = X_2,
             cell_area  = (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res),
             nn_index_1 = nn_X_1$which,
             nn_index_2 = nn_X_2$which)

parameters <- list(beta = c(0,0),
                   f = f,
                   g = g,
                   sigma_2 = 0.2,
                   alpha = 0.75,
                   tau_2 = 1) 

priors <- list(beta_mean = c(0,0), 
               beta_sd = c(10,10), 
               beta_prop_sd = c(0.05, 0.05),
               f_mean = 0,
               a_0_sigma = 2,
               b_0_sigma = 1,
               a_0_tau = 2,
               b_0_tau = 1,
               gamma = 0.5,
               phi = 10)

iters <- 5000

burnin <- 0

sim <- driver(parameters, priors, data, iters)

beta_post <- sim$beta[, (burnin+1):iters]
alpha_post <- sim$alpha[,(burnin+1):iters]
sigma_2_post <- sim$sigma_2[,(burnin+1):iters]
tau_2_post <- sim$tau_2[,(burnin+1):iters]
f_post <- sim$f[,(burnin+1):iters]
g_post <- sim$g[,(burnin+1):iters]

apply(beta_post, 1, mean)
apply(beta_post, 1, sd)

mean(sigma_2_post)
sd(sigma_2_post)

mean(tau_2_post)
sd(tau_2_post)

mean(alpha_post)
sd(alpha_post)


# Posterior Plots ----------------------------------------------------------------
# Posterior lambda for both sources
posterior_lambda <- matrix(NA, nrow = nrow(X_grid), ncol = (iters-burnin))

for(m in 1:(iters-burnin)){
  beta_m <- beta_post[,m]
  f_m <- f_post[, m]
  g_m <- g_post[, m]   # add source 2
  
  # log intensity = beta0 + beta1 * covariate + f + g
  log_lambda_m <- beta_m[1] + beta_m[2]*covariate + f_m + g_m
  
  posterior_lambda[, m] <- exp(log_lambda_m)
}

# Posterior mean intensity
lambda_mean <- rowMeans(posterior_lambda, na.rm = TRUE)

# Reshape to grid
lambda_mean_mat <- matrix(lambda_mean, 
                          nrow = grid_res, 
                          ncol = grid_res, 
                          byrow = FALSE)

# Plotting
par(mfrow = c(2,2))

image(x_seq, y_seq, log(lambda),
      main = "True Intensity",
      col = terrain.colors(50))

image.plot(x_seq, y_seq, log(lambda_mean_mat),
           main = "Posterior Mean",
           col = terrain.colors(50))

diff_mat <- log(lambda) - log(lambda_mean_mat)

image.plot(x_seq, y_seq, diff_mat,
           main = "Difference between True and Posterior",
           col = terrain.colors(50))

par(mfrow = c(1,1))

# Trace Plots -------------------------------------------------------------------

plot(sim$beta[1,], type = "l", main = "Beta Trace Plot")
plot(sim$f[1,], type = "l", main = "f trace plot")
plot(sim$g[1,], type = "l", main = "g trace plot")
plot(sim$sigma_2[1,], type = "l", main = "sigma_2 trace plot")
plot(sim$tau_2[1,], type = "l", main = "tau_2 trace plot")
plot(sim$alpha[1,], type = "l", main = "alpha trace plot")

# Testing -------------------------------------------------------------------------
loglike(parameters, data)
update_betas(parameters, priors, data)
update_f(parameters, priors, data)
update_g(parameters, priors, data)
update_sigma_2(parameters, priors, data)
update_alpha(parameters, priors, data)
update_tau_2(parameters, priors, data)
