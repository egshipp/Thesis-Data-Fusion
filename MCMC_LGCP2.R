# Packages
library(spatstat)
library(fields)
library(FastGP)
library(MASS)


# True LGCP --------------------------------------------------------------------
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

z <- as.vector(rcpp_rmvnorm(1,S,mu))

# Simulate NHPP LGCP
b_0 <- 1
b_1 <- 3

lambda <- exp(b_0 + b_1*(cov_field) + z)
lambda_im <- im(lambda, xcol = x_seq, yrow = y_seq)

lgcp_sim <- rpoispp(lambda_im)

# Discretize using spatstat

lgcp_discretize <- pixellate(lgcp_sim, eps = 1)

par(mfrow = c(1,2))
image.plot(x_seq, y_seq, cov_field, col = terrain.colors(100), main = "Simulated LGCP")
points(lgcp_sim)
plot(lgcp_discretize)
par(mfrow = c(1,1))



# Source 1  --------------------------------------------------------------------

set.seed(111)
sigma_2 <- 0.1
nrow <- length(lgcp_discretize$yrow)
ncol <- length(lgcp_discretize$xcol)

mu_1 <- -0.5 * sigma_2^2
f <- rnorm(nrow * ncol, mean = 0, sd = 0)
exp_f <- exp(f)

lambda_1 <- lgcp_discretize * exp_f
lgcp_1 <- rpoispp(lambda_1)

par(mfrow=(c(1,2)))
plot(lgcp_discretize)
plot(lambda_1)
par(mfrow=(c(1,1)))

par(mfrow=(c(1,2)))
plot(lgcp_sim)
plot(lgcp_1)
par(mfrow=(c(1,1)))

# Source 2 ----------------------------------------------------------------------

tau_2 <- 0.4
alpha <- -0.5

g <- rnorm(nrow * ncol, mean = alpha, sd = tau_2)
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
plot(lgcp_1)
plot(lgcp_2)
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
