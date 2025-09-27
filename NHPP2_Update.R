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

nhpp_discretize <- pixellate(nhpp_sim, eps = 1)

par(mfrow = c(1,2))
image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
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

# MCMC --------------------------------------------------------------------------


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
                   tau_2 = 1) 

priors <- list(beta_mean = c(0,0), 
               beta_sd = c(10,10), 
               beta_prop_sd = c(0.1, 0.1),
               f_mean = 0, 
               a_0 = 0.5,
               b_0 = 0.5)

iters <- 7500

sim <- driver(parameters, priors, data, iters)

