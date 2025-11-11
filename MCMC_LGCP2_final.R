# Packages
library(spatstat)
library(fields)
library(FastGP)
library(MASS)
library(ggplot2)


# True LGCP --------------------------------------------------------------------
# Simulate Covariate
par(mfrow = c(1,1))
win <- owin(xrange = c(0, 10), yrange = c(0, 10))

grid_res <- 10

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)

coords <- as.matrix(expand.grid(y_seq, x_seq))
grid_coords <- expand.grid(x = y_seq, y = x_seq)

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

#Simulate Gaussian random field

sigma_2 <- 0.25

S_z <-  sigma_2 * exp(-dists/1)

z <- as.vector(rcpp_rmvnorm(1,S_z,mu))

z_mat <- matrix(z, nrow = length(x_seq), ncol = length(y_seq))

z_ppp <- ppp(x = grid_coords$x,
             y = grid_coords$y,
             window = win,
             marks = z)

image.plot(x_seq, y_seq, z_mat, col = terrain.colors(100), main = "Simulated Exponential Latent Gaussian Field")

# Simulate LGCP
b_0 <- 1
b_1 <- 3

lambda <- exp(b_0 + b_1*(cov_field) + z)
lambda_im <- im(lambda, xcol = x_seq, yrow = y_seq)

lgcp_sim <- rpoispp(lambda_im)

plot(lgcp_sim)

# Discretize using spatstat

lgcp_discretize <- pixellate(lgcp_sim, eps = 1)

par(mfrow = c(1,2))
image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated LGCP")
points(lgcp_sim)
plot(lgcp_discretize)
par(mfrow = c(1,1))

# Source 1  --------------------------------------------------------------------

set.seed(111)
nrow <- length(lgcp_discretize$yrow)
ncol <- length(lgcp_discretize$xcol)

x_min_subwindow1 <- 0
x_max_subwindow1 <- 5
y_min_subwindow1 <- 5
y_max_subwindow1 <- 10

x_min_subwindow2 <- 5
x_max_subwindow2 <- 10
y_min_subwindow2 <- 0
y_max_subwindow2 <- 5

sub_window1 <- owin(xrange = c(x_min_subwindow1, x_max_subwindow1), yrange = c(y_min_subwindow1, y_max_subwindow1))

sub_window2 <- owin(xrange = c(x_min_subwindow2, x_max_subwindow2), yrange = c(y_min_subwindow2, y_max_subwindow2))

lambda_1 <- lgcp_discretize

lgcp_1_sub1 <- rpoispp(lambda_1[sub_window1])
lgcp_1_sub2 <- rpoispp(lambda_1[sub_window2])

lgcp_1 <- superimpose(lgcp_1_sub1, lgcp_1_sub2, W = owin(c(0,10), c(0,10)))

par(mfrow=(c(1,2)))
plot(lgcp_discretize)
plot(lambda_1)
par(mfrow=(c(1,1)))

par(mfrow=(c(1,2)))
plot(lgcp_sim)
plot(lgcp_1)
par(mfrow=(c(1,1)))

plot(win, main = "f measurement error field")
points(lgcp_1, col = "red", pch = 16)
rect(x_min_subwindow1, y_min_subwindow1, x_max_subwindow1, y_max_subwindow1, border = "blue", lwd = 2)
rect(x_min_subwindow2, y_min_subwindow2, x_max_subwindow2, y_max_subwindow2, border = "blue", lwd = 2)

# Source 2 ----------------------------------------------------------------------

tau_2 <- 0.4
S_g <- tau_2 * exp(-dists/1.5)
alpha <- -0.2

#g <- rnorm(nrow * ncol, alpha, tau_2)
g <- as.vector(rcpp_rmvnorm(1,S_g,alpha))

exp_g <- exp(g)

lambda_2 <- lgcp_discretize * exp_g

lgcp_2 <- rpoispp(lambda_2)

par(mfrow=(c(1,2)))
plot(lgcp_discretize)
plot(lambda_2)
par(mfrow=(c(1,1)))

par(mfrow=(c(1,2)))
plot(lgcp_sim)
plot(lgcp_2)
par(mfrow=(c(1,1)))

par(mfrow = (c(1,3)))
plot(lgcp_sim)
plot(win, main ="lgcp_1")
points(lgcp_1)
plot(lgcp_2)
par(mfrow=(c(1,1)))

# Data frame creation ----------------------------------------------------------

## X grid
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")
X_grid$covariate <- covariate

## Source 1 (small var)
X_1 <- as.data.frame(lgcp_1)
nn_X_1 <- nncross(lgcp_1, cov_field_ppp)
X_1$covariate <- cov_field_ppp$marks[nn_X_1$which]

# making a mask for X-grid to be used with source 1
inside_sub1 <- with(X_grid,
                    x >= x_min_subwindow1 & x <= x_max_subwindow1 &
                      y >= y_min_subwindow1 & y <= y_max_subwindow1)

inside_sub2 <- with(X_grid,
                    x >= x_min_subwindow2 & x <= x_max_subwindow2 &
                      y >= y_min_subwindow2 & y <= y_max_subwindow2)

X_grid$mask_source1 <- inside_sub1 | inside_sub2

## Source 2 (large var)
X_2 <- as.data.frame(lgcp_2)
nn_X_2 <- nncross(lgcp_2, cov_field_ppp)
X_2$covariate <- cov_field_ppp$marks[nn_X_2$which]


par(mfrow = (c(1,2)))
quilt.plot(X_1$x, X_1$y, X_1$covariate)
plot(lgcp_1)
par(mfrow = (c(1,1)))


par(mfrow = (c(1,2)))
quilt.plot(X_2$x, X_2$y, X_2$covariate)
plot(lgcp_2)
par(mfrow = (c(1,1)))

quilt.plot(X_grid$x, X_grid$y, X_grid$covariate)


# MCMC --------------------------------------------------------------------------

## log likelihood function
loglike <- function(parameters, data) {
  
  idx_mask1 <- which(X_grid$mask_source1)

  log_lambda_points1 <- parameters$beta[1] + parameters$beta[2] * data$X_1$covariate + parameters$z[data$nn_index_1]
  term1 <- sum(log_lambda_points1)

  log_lambda_grid1_full <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$z
  
  lambda_grid1_masked <- exp(log_lambda_grid1_full[idx_mask1])
  term2 <- sum(lambda_grid1_masked * data$cell_area)

  log_lambda_points2 <- parameters$beta[1] + parameters$beta[2] * data$X_2$covariate + parameters$g[data$nn_index_2] + parameters$z[data$nn_index_2]
  term3 <- sum(log_lambda_points2)

  log_lambda_grid2 <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$g + parameters$z
  lambda_grid2 <- exp(log_lambda_grid2)
  term4 <- sum(lambda_grid2 * data$cell_area)

  likelihood <- (term1 - term2) + (term3 - term4)
  return(likelihood)
}


## Updating slope estimates using Metropolis Hastings MCMC
update_betas<- function(parameters, priors, data){

  beta_cand <- rnorm(length(parameters$beta), mean = parameters$beta, sd = priors$beta_prop_sd)

  params_top <- list(beta = beta_cand, g = parameters$g, z = parameters$z)
  params_bottom <- list(beta = parameters$beta, g = parameters$g, z = parameters$z)

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

# Updating g - Source 2 measurement error
update_g <- function(parameters, priors, data){

  # Choosing ellipse (nu) from prior (g)
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$g)), Sigma = parameters$tau_2 * exp(-data$dists/1.5)))

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

update_z <- function(parameters, priors, data){

  # Choosing ellipse (nu) from prior (g)
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$z)), Sigma = parameters$sigma_2 * exp(-data$dists/1)))

  # Log likelihood threshold (finding log(y))

  u <- runif(1, min = 0, max = 1)

  log_y <- loglike(parameters, data) + log(u)

  # Draw an initial proposal for theta

  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta

  repeat {
    # Calculate z'
    z_prime <- as.vector(parameters$z*cos(theta) + nu*sin(theta))

    params_prime <- parameters
    params_prime$z <- z_prime

    # Shrinking bracket
    if(loglike(params_prime, data) > log_y){
      parameters$z <- z_prime
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

# Updating LGCP Gaussian Latent field variance - sigma_2

update_sigma_2 <- function(parameters, priors, data){
  n <- length(parameters$z)

  alpha_post <- priors$a_0_sigma + n/2
  beta_post  <- priors$b_0_sigma  + 0.5 * sum((parameters$z - priors$z_mean)^2)

  # Draw samples from Inverse-Gamma
  parameters$sigma_2 <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)

  return(parameters)

}

# Updating source 2 variance using Gibbs Sampling - tau_2

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
  mu_post <- (( sum(parameters$g) / parameters$tau_2 )) / denom
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
  out$g=matrix(NA, nrow = length(parameters$g), ncol = iters)
  out$z=matrix(NA, nrow = length(parameters$z), ncol = iters)
  out$sigma_2=matrix(NA, nrow = 1, ncol = iters)
  out$tau_2=matrix(NA, nrow = 1, ncol = iters)
  out$alpha=matrix(NA, nrow = 1, ncol = iters)

  for(k in 1:iters){
    parameters <- update_betas(parameters, priors, data)
    out$beta[,k]=parameters$beta

    parameters <- update_g(parameters, priors, data)
    out$g[,k] <- parameters$g

    parameters <- update_z(parameters, priors, data)
    out$z[,k] <- parameters$z

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
             nn_index_2 = nn_X_2$which,
             win = win,
             grid_res = 20,
             cell_size = cell_size,
             x_seq = x_seq,
             y_seq = y_seq,
             coords = as.matrix(expand.grid(y_seq, x_seq)),
             dists = as.matrix(dist(coords)))

parameters <- list(beta = c(0,0),
                   g = g,
                   z = z,
                   sigma_2 = sigma_2,
                   alpha = alpha,
                   tau_2 = tau_2)

priors <- list(beta_mean = c(0,0),
               beta_sd = c(10,10),
               beta_prop_sd = c(0.075, 0.075),
               z_mean = 0,
               a_0_sigma = 2,
               b_0_sigma = 1,
               a_0_tau = 2,
               b_0_tau = 1,
               phi = 10
)

iters <- 10000

burnin <- 3000

sim <- driver(parameters, priors, data, iters) # took 20 min to run with res = 20

beta_post <- sim$beta[, (burnin+1):iters]
alpha_post <- sim$alpha[,(burnin+1):iters]
sigma_2_post <- sim$sigma_2[,(burnin+1):iters]
tau_2_post <- sim$tau_2[,(burnin+1):iters]
g_post <- sim$g[,(burnin+1):iters]
z_post <- sim$z[,(burnin+1):iters]

apply(beta_post, 1, mean)
apply(beta_post, 1, sd)

mean(sigma_2_post)
sd(sigma_2_post)

mean(tau_2_post)
sd(tau_2_post)

mean(alpha_post)
sd(alpha_post)

mean(g_post)
sd(g_post)

mean(z_post)
sd(z_post)

# Posterior Plots ----------------------------------------------------------------
# Posterior lambda for both sources
posterior_lambda <- matrix(NA, nrow = nrow(X_grid), ncol = (iters-burnin))

for(m in 1:(iters-burnin)){
  beta_m <- beta_post[,m]
  
  log_lambda_m <- beta_m[1] + beta_m[2]*covariate + z_post[,m]
  
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

image.plot(x_seq, y_seq, log(t(lambda_mean_mat)),
           main = "Posterior Mean",
           col = terrain.colors(50))

diff_mat <- log(lambda) - log(t(lambda_mean_mat))

image.plot(x_seq, y_seq, diff_mat,
           main = "Difference between True and Posterior",
           col = terrain.colors(50))

par(mfrow = c(1,1))

# Trace Plots -------------------------------------------------------------------

plot(sim$beta[1,], type = "l", main = "Beta 1 Trace Plot")
plot(sim$beta[2,], type = "l", main = "Beta 2 Trace Plot")
plot(sim$z[1,], type = "l", main = "z trace plot")
plot(sim$g[1,], type = "l", main = "g trace plot")
plot(sim$sigma_2[1,], type = "l", main = "sigma_2 trace plot")
plot(sim$tau_2[1,], type = "l", main = "tau_2 trace plot")
plot(sim$alpha[1,], type = "l", main = "alpha trace plot")

# Running multiple simulations to show confidence intervals ---------------------

n_sims <- 100

n_sims_df <- data.frame()

for (i in 1:n_sims){
  cat("Running simulation", i, "of", n_sims, "\n")
  
  #Simulate Covariate
  win <- owin(xrange = c(0, 10), yrange = c(0, 10))
  
  grid_res <- 10
  
  cell_size <- diff(win$xrange) / grid_res
  
  x_seq <- seq(win$xrange[1] + cell_size/2,
               win$xrange[2] - cell_size/2,
               by = cell_size)
  y_seq <- seq(win$yrange[1] + cell_size/2,
               win$yrange[2] - cell_size/2,
               by = cell_size)
  
  coords <- as.matrix(expand.grid(y_seq, x_seq))
  grid_coords <- expand.grid(x = y_seq, y = x_seq)
  
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
  
  #Simulate Gaussian random field
  
  sigma_2 <- 0.25
  
  S_z <-  sigma_2 * exp(-dists/1)
  
  z <- as.vector(rcpp_rmvnorm(1,S_z,mu))
  
  z_mat <- matrix(z, nrow = length(x_seq), ncol = length(y_seq))
  
  z_ppp <- ppp(x = grid_coords$x,
               y = grid_coords$y,
               window = win,
               marks = z)
  
  # Simulate LGCP
  b_0 <- 1
  b_1 <- 3
  
  lambda <- exp(b_0 + b_1*(cov_field) + z)
  lambda_im <- im(lambda, xcol = x_seq, yrow = y_seq)
  
  lgcp_sim <- rpoispp(lambda_im)
  
  plot(lgcp_sim)
  
  # Discretize using spatstat
  
  lgcp_discretize <- pixellate(lgcp_sim, eps = 1)
  
  # Source 1  -----------
  nrow <- length(lgcp_discretize$yrow)
  ncol <- length(lgcp_discretize$xcol)
  
  x_min_subwindow1 <- 0
  x_max_subwindow1 <- 5
  y_min_subwindow1 <- 5
  y_max_subwindow1 <- 10
  
  x_min_subwindow2 <- 5
  x_max_subwindow2 <- 10
  y_min_subwindow2 <- 0
  y_max_subwindow2 <- 5
  
  sub_window1 <- owin(xrange = c(x_min_subwindow1, x_max_subwindow1), yrange = c(y_min_subwindow1, y_max_subwindow1))
  
  sub_window2 <- owin(xrange = c(x_min_subwindow2, x_max_subwindow2), yrange = c(y_min_subwindow2, y_max_subwindow2))
  
  lambda_1 <- lgcp_discretize
  
  lgcp_1_sub1 <- rpoispp(lambda_1[sub_window1])
  lgcp_1_sub2 <- rpoispp(lambda_1[sub_window2])
  
  lgcp_1 <- superimpose(lgcp_1_sub1, lgcp_1_sub2, W = owin(c(0,10), c(0,10)))
  
  # Source 2 --------
  
  tau_2 <- 0.4
  S_g <- tau_2 * exp(-dists/1.5)
  alpha <- -0.2
  
  #g <- rnorm(nrow * ncol, alpha, tau_2)
  g <- as.vector(rcpp_rmvnorm(1,S_g,alpha))
  
  exp_g <- exp(g)
  
  lambda_2 <- lgcp_discretize * exp_g
  
  lgcp_2 <- rpoispp(lambda_2)
  
  # Data frame creation -----
  
  ## X grid
  X_grid <- as.data.frame(coords)
  colnames(X_grid) <- c("x", "y")
  X_grid$covariate <- covariate
  
  ## Source 1 (small var)
  X_1 <- as.data.frame(lgcp_1)
  nn_X_1 <- nncross(lgcp_1, cov_field_ppp)
  X_1$covariate <- cov_field_ppp$marks[nn_X_1$which]
  
  # making a mask for X-grid to be used with source 1
  inside_sub1 <- with(X_grid,
                      x >= x_min_subwindow1 & x <= x_max_subwindow1 &
                        y >= y_min_subwindow1 & y <= y_max_subwindow1)
  
  inside_sub2 <- with(X_grid,
                      x >= x_min_subwindow2 & x <= x_max_subwindow2 &
                        y >= y_min_subwindow2 & y <= y_max_subwindow2)
  
  X_grid$mask_source1 <- inside_sub1 | inside_sub2
  
  ## Source 2 (large var)
  X_2 <- as.data.frame(lgcp_2)
  nn_X_2 <- nncross(lgcp_2, cov_field_ppp)
  X_2$covariate <- cov_field_ppp$marks[nn_X_2$which]
  
  data <- list(X_grid = X_grid,
               X_1 = X_1,
               X_2 = X_2,
               cell_area  = (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res),
               nn_index_1 = nn_X_1$which,
               nn_index_2 = nn_X_2$which,
               win = win,
               grid_res = 10,
               cell_size = cell_size,
               x_seq = x_seq,
               y_seq = y_seq,
               coords = as.matrix(expand.grid(y_seq, x_seq)),
               dists_z = as.matrix(dist(coords)))
  
  parameters <- list(beta = c(0,0),
                     g = g,
                     z = z,
                     sigma_2 = sigma_2,
                     alpha = alpha,
                     tau_2 = tau_2)
  
  priors <- list(beta_mean = c(0,0),
                 beta_sd = c(10,10),
                 beta_prop_sd = c(0.075, 0.075),
                 z_mean = 0,
                 a_0_sigma = 2,
                 b_0_sigma = 1,
                 a_0_tau = 2,
                 b_0_tau = 1,
                 phi = 10
  )
  
  iters <- 10000
  
  burnin <- 1000
  
  # Run simulation ----------
  
  sim <- driver(parameters, priors, data, iters)
  
  beta_post <- sim$beta[, (burnin+1):iters]
  alpha_post <- sim$alpha[, (burnin+1):iters]
  sigma_2_post <- sim$sigma_2[, (burnin+1):iters]
  tau_2_post <- sim$tau_2[, (burnin+1):iters]
  
  # Store just summary stats
   n_sims_df <- rbind(n_sims_df, data.frame(
    sim = i,
    
    # beta parameters
    beta0_mean = mean(beta_post[1,]), beta0_sd = sd(beta_post[1,]),
    beta1_mean = mean(beta_post[2,]), beta1_sd = sd(beta_post[2,]),
    
    # hyperparameters
    sigma2_mean = mean(sigma_2_post), sigma2_sd = sd(sigma_2_post),
    tau2_mean = mean(tau_2_post), tau2_sd = sd(tau_2_post),
    alpha_mean = mean(alpha_post), alpha_sd = sd(alpha_post)
    
  ))
}


# 95% Confidence Intervals for Estimates ----------------------------------------

#Beta 
n_sims_df$beta0_lower <- n_sims_df$beta0_mean - 1.96*n_sims_df$beta0_sd
n_sims_df$beta0_upper <- n_sims_df$beta0_mean + 1.96*n_sims_df$beta0_sd
n_sims_df$covered_beta0 <- (n_sims_df$beta0_lower <= 1 &
                              n_sims_df$beta0_upper >= 1)

n_sims_df$beta1_lower <- n_sims_df$beta1_mean - 1.96*n_sims_df$beta1_sd
n_sims_df$beta1_upper <- n_sims_df$beta1_mean + 1.96*n_sims_df$beta1_sd
n_sims_df$covered_beta1 <- (n_sims_df$beta1_lower <= 3 &
                              n_sims_df$beta1_upper >= 3)

beta0_covered <- mean(n_sims_df$beta0_lower <= 1 &
                        n_sims_df$beta0_upper >= 1)

beta1_covered <- mean(n_sims_df$beta1_lower <= 3 &
                        n_sims_df$beta1_upper >= 3)
beta0_covered
beta1_covered

#sigma_2
n_sims_df$sigma2_lower <- n_sims_df$sigma2_mean - 1.96*n_sims_df$sigma2_sd
n_sims_df$sigma2_upper <- n_sims_df$sigma2_mean + 1.96*n_sims_df$sigma2_sd
n_sims_df$covered_sigma2 <- (n_sims_df$sigma2_lower <= 0.25 &
                              n_sims_df$sigma2_upper >=0.25)

covered_sigma2 <- mean(n_sims_df$sigma2_lower <= 0.25 &
                     n_sims_df$sigma2_upper >=0.25)
covered_sigma2

#tau_2
n_sims_df$tau2_lower <- n_sims_df$tau2_mean - 1.96*n_sims_df$tau2_sd
n_sims_df$tau2_upper <- n_sims_df$tau2_mean + 1.96*n_sims_df$tau2_sd
n_sims_df$covered_tau2 <- (n_sims_df$tau2_lower <= 0.4 &
                               n_sims_df$tau2_upper >=0.4)

covered_tau2 <- mean(n_sims_df$tau2_lower <= 0.4 &
                   n_sims_df$tau2_upper >=0.4)
covered_tau2

#alpha
n_sims_df$alpha_lower <- n_sims_df$alpha_mean - 1.96*n_sims_df$alpha_sd
n_sims_df$alpha_upper <- n_sims_df$alpha_mean + 1.96*n_sims_df$alpha_sd
n_sims_df$covered_alpha <- (n_sims_df$alpha_lower <= -0.2 &
                             n_sims_df$alpha_upper >= -0.2)

covered_alpha <- mean(n_sims_df$alpha_lower <= -0.2 &
                    n_sims_df$alpha_upper >= -0.2)
covered_alpha

#Confidence Interval Plots -------------------------------------------------------
pdf("credible_intervals.pdf", width = 8, height = 6)  
#Betas
ggplot(n_sims_df, aes(x = sim, y = beta0_mean, 
                      color = covered_beta0)) +
  geom_pointrange(aes(ymin = beta0_lower, ymax = beta0_upper),
                  size = 1.1) +
  geom_hline(aes(yintercept = 1),
             linetype = "dashed", color = "black", size = 0.7) +
  scale_color_manual(values = c("red","green")) +
  labs(
    x = "Parameter",
    y = "Posterior Estimate",
    title = "95% Credible Intervals for" ~ beta[0],
    color = "Contains truth?"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Proportion covered: ", round(beta0_covered, 2)),
    hjust = 1.1, vjust = -1.5,
    size = 5
  ) +
  coord_cartesian(clip = "off")

ggplot(n_sims_df, aes(x = sim, y = beta1_mean, 
                      color = covered_beta1)) +
  geom_pointrange(aes(ymin = beta1_lower, ymax = beta1_upper),
                  size = 1.1) +
  geom_hline(aes(yintercept = 3),
             linetype = "dashed", color = "black", size = 0.7) +
  scale_color_manual(values = c("red","green")) +
  labs(
    x = "Parameter",
    y = "Posterior Estimate",
    title = "95% Credible Intervals for" ~ beta[1],
    color = "Contains truth?"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Proportion covered: ", round(beta1_covered, 2)),
    hjust = 1.1, vjust = -1.5,
    size = 5
  ) +
  coord_cartesian(clip = "off")

#Sigma2
ggplot(n_sims_df, aes(x = sim, y = sigma2_mean, 
                      color = covered_sigma2)) +
  geom_pointrange(aes(ymin = sigma2_lower, ymax = sigma2_upper),
                  size = 1.1) +
  geom_hline(aes(yintercept = 0.25),
             linetype = "dashed", color = "black", size = 0.7) +
  scale_color_manual(values = c("red","green")) +
  labs(
    x = "Parameter",
    y = "Posterior Estimate",
    title = "95% Credible Intervals for" ~ sigma^2,
    color = "Contains truth?"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Proportion covered: ", round(covered_sigma2, 2)),
    hjust = 1.1, vjust = -1.5,
    size = 5
  ) +
  coord_cartesian(clip = "off")

#Tau2
ggplot(n_sims_df, aes(x = sim, y = tau2_mean, 
                      color = covered_tau2)) +
  geom_pointrange(aes(ymin = tau2_lower, ymax = tau2_upper),
                  size = 1.1) +
  geom_hline(aes(yintercept = 0.4),
             linetype = "dashed", color = "black", size = 0.7) +
  scale_color_manual(values = c("green","red")) +
  labs(
    x = "Parameter",
    y = "Posterior Estimate",
    title = "95% Credible Intervals for" ~ tau^2,
    color = "Contains truth?"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Proportion covered: ", round(covered_tau2, 2)),
    hjust = 1.1, vjust = -1.5,
    size = 5
  ) +
  coord_cartesian(clip = "off")

#Alpha
ggplot(n_sims_df, aes(x = sim, y = alpha_mean, 
                      color = covered_alpha)) +
  geom_pointrange(aes(ymin = alpha_lower, ymax = alpha_upper),
                  size = 1.1) +
  geom_hline(aes(yintercept = -0.2),
             linetype = "dashed", color = "black", size = 0.7) +
  scale_color_manual(values = c("green","red")) +
  labs(
    x = "Parameter",
    y = "Posterior Estimate",
    title = "95% Credible Intervals for" ~ alpha,
    color = "Contains truth?"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Proportion covered: ", round(covered_alpha, 2)),
    hjust = 1.1, vjust = -1.5,
    size = 5
  ) +
  coord_cartesian(clip = "off")

dev.off()