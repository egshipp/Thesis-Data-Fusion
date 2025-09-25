# Packages
library(spatstat)
library(fields)
library(FastGP)
library(MASS)

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

cov_field <- matrix(covariate, nrow = grid_res, ncol = grid_res)
cov_field_ppp <- ppp(x = grid_coords$x,
                     y = grid_coords$y,
                     window = win,
                     marks = covariate)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated Exponential Covariate")

# Simulate NHPP

b_0 <- 1
b_1 <- 3

lambda <- exp(b_0 + b_1*(cov_field))

lambda_im <- im(t(lambda), xcol = x_seq, yrow = y_seq) # fixed transposition issue by transposing lambda

nhpp_sim <- rpoispp(lambda_im) # has to be an image to be passed through rpoispp

plot(nhpp_sim)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
  points(nhpp_sim)
  
# Discretizing
  
nhpp_discretize <- pixellate(nhpp_sim, eps = 1)

par(mfrow = c(1,2))
  image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
    points(nhpp_sim)
  plot(nhpp_discretize)
par(mfrow = c(1,1))

# Creating f 
f_Sigma <- 0.01
nrow <- length(nhpp_discretize$yrow)
ncol <- length(nhpp_discretize$xcol)

f_mu <- -0.5 * f_Sigma^2
f <- rnorm(nrow * ncol, mean = f_mu, sd = f_Sigma)
exp_f <- exp(f)

## MCMC 
  
# Creating X and X_grid data frames
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")
X_grid$covariate <- covariate

X <- as.data.frame(nhpp_sim)
nn_X <- nncross(nhpp_sim, cov_field_ppp)
X$covariate <- cov_field_ppp$marks[nn_X$which]

#Log-Likelihood Function

loglike <- function(parameters, data) {
  
  log_lambda_points <- parameters$beta[1] + parameters$beta[2] * data$X$covariate + parameters$f[data$nn_index]
  term1 <- sum(log_lambda_points)
  
  log_lambda_grid <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$f
  lambda_grid <- exp(log_lambda_grid)
  term2 <- sum(lambda_grid * data$cell_area)
  
  likelihood <- sum(term1 - term2)
  return(likelihood)
}

#Update betas

update_betas <- function(parameters, priors, data){
 
  #beta_top <- parameters$beta
  beta_cand <- rnorm(length(parameters$beta), mean = parameters$beta, sd = priors$beta_prop_sd) 
  
  params_top <- list(beta = beta_cand, f = parameters$f)
  params_bottom <- list(beta = parameters$beta, f = parameters$f)
  
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

#Update f

update_f <- function(parameters, priors, data){
  
  # Choosing ellipse (nu) from prior (f)
 nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$f)), Sigma = diag(parameters$Sigma_2, length(parameters$f))))
  
  # nu <- rnorm(length(parameters$f), 0, sqrt(0.2))
  
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

# Update source 2 variance

update_sigma_2 <- function(parameters, priors, data){
  n <- length(parameters$f)
  
  alpha_post <- priors$a_0 + n/2
  beta_post  <- priors$b_0  + 0.5 * sum((f - priors$f_mean)^2)
  
  # Draw samples from Inverse-Gamma
  parameters$sigma_2 <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)
  
  return(parameters)
  
}

# Wrapper function
driver <- function(parameters, priors, data, iters){
  out=list()
  out$params=parameters
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(parameters$beta),ncol = iters)
  out$f=matrix(NA, nrow = length(parameters$f), ncol = iters)
  out$sigma_2=matrix(NA, nrow = 1, ncol = iters)
  
  for(k in 1:iters){
    parameters <- update_betas(parameters, priors, data)
    out$beta[,k]=parameters$beta 
    
    parameters <- update_f(parameters, priors, data)
    out$f[,k] <- parameters$f
    
    parameters <- update_sigma_2(parameters, priors, data)
    out$sigma_2[,k] <- parameters$sigma_2
  }
  return(out)
}

# Attempting

data <- list(X_grid = X_grid, 
             X = X, 
             cell_area  = (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res),
             nn_index = nn_X$which)

parameters <- list(beta = c(0,0),
                   f = f,
                   Sigma_2 = 0.2) 

priors <- list(beta_mean = c(0,0), 
               beta_sd = c(10,10), 
               beta_prop_sd = c(0.1, 0.1),
               f_mean = 0, 
               a_0 = 0.5,
               b_0 = 0.5)

iters <- 7500

sim <- driver(parameters, priors, data, iters)

burnin <- 1000

beta_post <- sim$beta[, (burnin+1):iters]
f_post <- sim$f[, (burnin+1):iters]
sigma_2_post <- sim$sigma_2[,(burnin+1):iters]

apply(beta_post, 1, mean)
apply(beta_post, 1, sd)

mean(sigma_2_post)
sd(sigma_2_post)

# Calculating posterior lambda

posterior_lambda <- matrix(NA, nrow = nrow(X_grid), ncol = (iters-burnin))

for(m in 1:(iters-burnin)){
  beta_m <- beta_post[,m]
  f_m <- f_post[, m]
  
  log_lambda_m <- beta_m[1] + beta_m[2]*covariate + f_m
  
  posterior_lambda[, m] <- exp(log_lambda_m)
}

lambda_mean <- rowMeans(posterior_lambda, na.rm = TRUE)

lambda_mean_mat <- matrix(lambda_mean, 
                          nrow = grid_res, 
                          ncol = grid_res, 
                          byrow = FALSE)


par(mfrow = c(2,2))

zlim <- range(log(lambda), log(posterior_lambda))

image(x_seq, y_seq, log(lambda),
      main = "True Intensity",
      zlim = zlim,
      col = terrain.colors(50))

image.plot(x_seq, y_seq, log(lambda_mean_mat),
      main = "Posterior Mean",
      zlim = zlim,
      col = terrain.colors(50))


diff_mat <- log(lambda) - log(lambda_mean_mat)

image.plot(x_seq, y_seq, diff_mat,
           main = "Difference between True and Posterior",
           col = terrain.colors(50))

par(mfrow = c(1,1))

# Trace plots

par(mfrow = c(1,2))
plot(beta_post[1,],type = "l", main = "Beta 1 Trace Plot")
plot(beta_post[2,], type = "l", main = "Beta 2 Trace Plot")
par(mfrow = c(1,1))


plot(sim$sigma_2[1,],type = "l", main = "Sigma_2 Trace Plot")


