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
