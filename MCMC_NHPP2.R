# Packages
library(spatstat)
library(fields)
library(FastGP)

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
lambda_im <- im(lambda, xcol = x_seq, yrow = y_seq)

nhpp_sim <- rpoispp(lambda_im)

# Discretize using spatstat

nhpp_discretize <- pixellate(nhpp_sim, eps = 1)

par(mfrow = c(1,2))
image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
  points(nhpp_sim)
plot(nhpp_discretize)
par(mfrow = c(1,1))

# Source 1  --------------------------------------------------------------------

set.seed(111)
sigma_1 <- 0.01
nrow <- length(nhpp_discretize$yrow)
ncol <- length(nhpp_discretize$xcol)

mu_1 <- -0.5 * sigma_1^2
f_1 <- rnorm(nrow * ncol, mean = mu_1, sd = sigma_1)
exp_f_1 <- exp(f_1)

lambda_1 <- nhpp_discretize * exp_f_1

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

sigma_2 <- 1
alpha <- 0.75

f_2 <- rnorm(nrow * ncol, mean = alpha, sd = sigma_2)
exp_f_2 <- exp(f_2)

lambda_2 <- nhpp_discretize * exp_f_2

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

#Creating Dataframes -----------------------------------------------------------
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")
X_grid$covariate <- covariate

X_1 <- as.data.frame(nhpp_1)
nn_X_1 <- nncross(nhpp_1, cov_field_ppp)
X_1$covariate <- cov_field_ppp$marks[nn_X_1$which]


X_2 <- as.data.frame(nhpp_2)
nn_X_2 <- nncross(nhpp_2, cov_field_ppp)
X_2$covariate <- cov_field_ppp$marks[nn_X_2$which]


quilt.plot(X_1$x, X_1$y, X_1$covariate)
quilt.plot(X_2$x, X_2$y, X_2$covariate)
quilt.plot(X_grid$x, X_grid$y, X_grid$covariate)

#MCMC ---------------------------------------------------------------------------

#Log-Likelihood Functions

loglike_joint <- function(parameters, data) {
  
  log_lambda_points1 <- parameters$beta[1] + parameters$beta[2] * data$X_1$covariate + data$f_1[nn_X_1$which]
  term1 <- sum(log_lambda_points1)
  
  log_lambda_grid1 <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + data$f_1
  lambda_grid1 <- exp(log_lambda_grid1)
  term2 <- sum(lambda_grid1 * data$cell_area)

  log_lambda_points2 <- parameters$beta[1] + parameters$beta[2] * data$X_2$covariate + data$f_2[nn_X_2$which]
  term3 <- sum(log_lambda_points2)
  
  log_lambda_grid2 <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + data$f_2
  lambda_grid2 <- exp(log_lambda_grid2)
  term4 <- sum(lambda_grid2 * data$cell_area)
  
  likelihood <- (term1 - term2) + (term3 - term4)
  return(likelihood)
}

# log prior for both betas and f

logprior <- function(parameters, priors) {
  
  prior_beta <- sum(dnorm(parameters$beta, mean = priors$beta_mean, sd = priors$beta_sd, log = TRUE))
  
  prior_f_1 <- sum(dnorm(parameters$f_1_current, mean = priors$f_1_mean, sd = priors$f_1_sd, log = TRUE))
  
  prior_f_2 <- sum(dnorm(parameters$f_2_current, mean = priors$f_2_mean, sd = priors$f_2_sd, log = TRUE))
  
  return(as.numeric(prior_beta + prior_f_1 + prior_f_2))
}

# Update betas

update_betas_joint <- function(parameters, priors, data){
  
  beta_top <- parameters$beta
  beta_bottom <- rnorm(length(beta_top), mean = beta_top, sd = priors$beta_sd)
  
  params_top <- list(beta = beta_top, f_1_current = parameters$f_1_current, f_1_current = parameters$f_2_current)
  params_bottom <- list(beta = beta_bottom, f_1_current = parameters$f_1_current, f_1_current = parameters$f_2_current)
  
  # Posteriors
  post_top <- loglike_joint(params_top, data) 
  post_bottom <- loglike_joint(params_bottom, data) 
  + logprior(beta_bottom, f_1, f_2, sigma_1, sigma_2) 
  
  # Metropolis Hastings Ratio
  log_acceptance_ratio <- post_bottom - post_top
  
  if(log(runif(1)) < log_acceptance_ratio) {
    params$beta <- beta_bottom
  }
  
  return(params)
}

# Update f_1

update_f_1 <- function(beta, f_1_current, f_2, X_1, X_2, X_grid, cell_area, mu, sigma_1){
  
  # Choosing ellipse (nu) from prior (f)
  nu <- mvrnorm(n = 1, mu = rep(0, length(f_1_current)), Sigma = sigma_1)
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike_joint(beta,f_1_current, f_2, X_1, X_2, X_grid, cell_area) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
  repeat {
    # Calculate f'
    f_1_prime <- f_1_current*cos(theta) + nu*sin(theta)
    
    # Shrinking bracket
    if(loglike_joint(beta, f_1_prime,f_2, X_1,  X_2, X_grid, cell_area) > log_y){
      return(f_1_prime)} else {
        if(theta < 0) {
          theta_min <- theta
        } else {
          if (theta > 0) {
            theta_max <- theta
          }
        }
      }
    theta <- runif(1, theta_min, theta_max)
  }
}

# Update f_2

update_f_2 <- function(beta, f_2_current, f_1, X_1, X_2, X_grid, cell_area, mu, sigma_2){
  
  # Choosing ellipse (nu) from prior (f)
  nu <- mvrnorm(n = 1, mu = rep(0, length(f_2_current)), Sigma = sigma_2)
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike_joint(beta,f_2_current, f_1, X_1, X_2, X_grid, cell_area) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
  repeat {
    # Calculate f'
    f_2_prime <- f_2_current*cos(theta) + nu*sin(theta)
    
    # Shrinking bracket
    if(loglike_joint(beta, f_2_prime,f_1, X_1,  X_2, X_grid, cell_area) > log_y){
      return(f_2_prime)} else {
        if(theta < 0) {
          theta_min <- theta
        } else {
          if (theta > 0) {
            theta_max <- theta
          }
        }
      }
    theta <- runif(1, theta_min, theta_max)
  }
}

# Wrapper function
driver_joint <- function(X_1, X_2, X_grid, f_1, f_2, 
                         sigma_1, sigma_2, params, cell_area, priors,  proposal_sd, iters){
  out=list()
  out$params=params
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(params$beta),ncol = iters)
  out$f_1   <- matrix(NA, nrow = length(f_1), ncol = iters)
  out$f_2   <- matrix(NA, nrow = length(f_2), ncol = iters)
  
  for(k in 1:iters){
    params=update_betas_joint(params, priors, proposal_sd, f_1, f_2, sigma_1, sigma_2,
                              X_1,X_2, X_grid, cell_area)
    out$beta[,k]=params$beta
    
    f_1 <- update_f_1(params$beta, f_1, f_2, X_1, X_2, X_grid, cell_area, mu, sigma_1)
    out$f_1[, k]  <- f_1

    f_2 <- update_f_2(params$beta, f_2, f_1, X_1, X_2, X_grid, cell_area, mu, sigma_2)
    out$f_2[, k]  <- f_2
    
  }
  return(out)
}

# Attempting---------------------------------------------------------------------

data <- list(X_grid = X_grid, 
             X_1 = X_1, 
             X_2 = X_2,
             f_1 = f_1,
             f_2 = f_2,
             cell_area  = (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res))

parameters <- list(beta = c(0,0),
                   f_1_current = rep(0, nrow(X_grid)), 
                   f_2_current = rep(0, nrow(X_grid)),
                   f_1_sigma = diag(0.01, nrow(X_grid)),
                   f_2_sigma = diag(1, nrow(X_grid))) 

priors <- list(beta_mean = c(0,0), 
               beta_sd = c(0.7,0.7), 
               f_1_mean = 0, 
               f_1_sd = 0.1,
               f_2_mean = 0,
               f_2_sd = 2)

iters <- 7500


sim <- driver(parameters, priors, data, iters)



# Create posterior lambda

posterior_lambda <- matrix(NA, nrow = nrow(X_grid), ncol = iters)

for(m in 1:iters){
  beta_m <- sim$beta[,m]
  
  log_lambda_m <- beta_m[1] + beta_m[2]*covariate
  
  posterior_lambda[, m] <- exp(log_lambda_m)
}


lambda_mean <- rowMeans(posterior_lambda)

lambda_mean_mat <- matrix(lambda_mean, 
                          nrow = grid_res, 
                          ncol = grid_res, 
                          byrow = FALSE)

posterior_lambda_im <- im(mat = t(lambda_mean_mat), 
                          xcol = x_seq, 
                          yrow = y_seq)


par(mfrow = c(1,2))
plot(log(lambda), main = "True Mean" )
plot(log(posterior_lambda_im), main = "Posterior Mean")
par(mfrow = c(1,1))


plot(log(lambda$v),log(posterior_lambda_im$v))
  abline(0,1)

plot(sim$beta[1,],type = "l")
plot(sim$beta[2,], type = "l")

# Testing each function

loglike_joint(parameters, data)

logprior(parameters, priors)
