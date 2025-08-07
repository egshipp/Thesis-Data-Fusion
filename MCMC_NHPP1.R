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

lambda <- b_0 + b_1*(cov_im)

nhpp_sim <- rpoispp(lambda)

plot(nhpp_sim)

image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated NHPP")
  points(nhpp_sim)
