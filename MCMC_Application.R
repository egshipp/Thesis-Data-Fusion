# Packages
library(spatstat)
library(fields)
library(FastGP)
library(MASS)
library(maps)
library(sf)

# Massachusetts map --------------------------------------------------------------
mass_map <- maps::map("state", "massachusetts", plot = FALSE, fill = TRUE)

mass_poly <- sf::st_as_sf(mass_map) |> 
  st_transform(4326) 

# Data ---------------------------------------------------------------------------
automated_locations <- readRDS("automated_locations.rds")
manual_locations <- readRDS("manual_locations.rds")

auto_loc_mat <- as.data.frame(automated_locations)
  colnames(auto_loc_mat) = c("x","y")
  
manual_loc_mat <- as.data.frame(manual_locations) # Over entire window
  colnames(manual_loc_mat) = c("x","y")

x_all <- range(c(auto_loc_mat$x, manual_loc_mat$x))
y_all <- range(c(auto_loc_mat$y, manual_loc_mat$y))

par(mfrow = c(2,2))  
plot(auto_loc_mat$x, auto_loc_mat$y, main = "Automated Locations", 
     xlim = x_all, ylim = y_all)
  map("state", "massachusetts", add = TRUE)
plot(manual_loc_mat$x, manual_loc_mat$y, main = "Manual Locations", 
     xlim = x_all, ylim = y_all)
  map("state", "massachusetts", add = TRUE)
par(mfrow = c(1,1))

# Creating Window ---------------------------------------------------------------
x_bound <- c(-70.52832 - 0.05, -70.01345 + 0.05)
y_bound <- c(41.77413 - 0.05,  42.16180 + 0.05)

win <- owin(
  xrange = x_bound,
  yrange = y_bound
)

grid_res <- 20

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             length.out = grid_res)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             length.out = grid_res)

coords <- as.matrix(expand.grid(x_seq, y_seq))
grid_coords <- expand.grid(x = x_seq, y = y_seq)

grid_ppp <- ppp(x = grid_coords$x, y = grid_coords$y, window = win)

plot(grid_ppp, xlim = x_bound, ylim = y_bound, axes = TRUE)
  abline(v = c(-70.52832, -70.01345), col = "red")
  abline(h = c(41.77413,  42.16180), col = "red")

#Converting Data Matrices to Shape files ----------------------------------------
auto_sf <- st_as_sf(auto_loc_mat, coords = c("x", "y"), crs = 4326)
manual_sf <- st_as_sf(manual_loc_mat, coords = c("x", "y"), crs = 4326)
grid_sf <- st_as_sf(grid_coords, coords = c("x", "y"), crs = 4326)
  grid_sf$id <- seq_len(nrow(grid_sf))

auto_water <- auto_sf |>
  st_filter(mass_poly, .predicate = st_disjoint)
auto_water_xy <- as.data.frame(st_coordinates(auto_water))

manual_water <- manual_sf |>
  st_filter(mass_poly, .predicate = st_disjoint)
manual_water_xy <- as.data.frame(st_coordinates(manual_water))

grid_water <- grid_sf |>
  st_filter(mass_poly, .predicate = st_disjoint)
grid_water_xy <- as.data.frame(st_coordinates(grid_water))

grid_land <- grid_sf |> st_filter(mass_poly, .predicate = st_intersects)
grid_land_xy <- as.data.frame(st_coordinates(grid_land))
                               
par(mfrow = c(2,2))  
plot(auto_water_xy$X, auto_water_xy$Y, main = "Automated Locations", 
     xlim = x_all, ylim = y_all)
map("state", "massachusetts", add = TRUE)
plot(manual_water_xy$X, manual_water_xy$Y, main = "Manual Locations", 
     xlim = x_all, ylim = y_all)
map("state", "massachusetts", add = TRUE)
plot(grid_water_xy$X, grid_water_xy$Y, main = "Grid Locations", 
     xlim = x_all, ylim = y_all)
map("state", "massachusetts", add = TRUE)
par(mfrow = c(1,1))

# Dataframe creation for MCMC ---------------------------------------------------

## Grid
X_grid <- as.data.frame(grid_water_xy)
  colnames(X_grid) <- c("x", "y")
  
X_1 <- manual_water_xy
nn_X_1 <- nncross(manual_water_xy, grid_ppp)

## Source 2 (Automated Locations)

X_2 <- auto_water_xy
nn_X_2 <- nncross(auto_water_xy, grid_ppp)

par(mfrow = c(2,2))
plot(X_grid$x, X_grid$y, xlim = x_bound, ylim = y_bound, main = "Grid")
plot(X_1, xlim = x_bound, ylim = y_bound, main = "Source 1 Point Process")
plot(X_2, xlim = x_bound, ylim = y_bound, main = "Source 2 Point Process")
par(mfrow = c(1,1))

# MCMC --------------------------------------------------------------------------

## log likelihood function
loglike <- function(parameters, data) {
  
  log_lambda_points1 <- parameters$beta[1] + parameters$z[data$nn_index_1]
  term1 <- sum(log_lambda_points1)
  
  log_lambda_grid1_full <- parameters$beta[1] + parameters$z
  
  lambda_grid1_masked <- exp(log_lambda_grid1_full)
  term2 <- sum(lambda_grid1_masked * data$cell_area)
  
  log_lambda_points2 <- parameters$beta[1] + parameters$g[data$nn_index_2] + parameters$z[data$nn_index_2]
  term3 <- sum(log_lambda_points2)
  
  log_lambda_grid2 <- parameters$beta[1] + parameters$g + parameters$z
  lambda_grid2 <- exp(log_lambda_grid2)
  term4 <- sum(lambda_grid2 * data$cell_area)
  
  likelihood <- (term1 - term2) + (term3 - term4)
  return(likelihood)
}


## Updating slope estimates using Metropolis Hastings MCMC
update_betas<- function(parameters, priors, data){
  
  beta_cand <- rnorm(1, mean = parameters$beta, sd = priors$beta_prop_sd) 
  
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
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$g)), Sigma = parameters$tau_2 * exp(-data$dists/0.1)))
  
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
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$z)), Sigma = parameters$sigma_2 * exp(-data$dists/0.1)))
  
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

# Parameters, priors, and data ---------------------------------------------------

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
             coords = as.matrix(expand.grid(x_seq, y_seq)),
             dists = as.matrix(dist(coords)))


parameters <- list(
             beta = 0,       
             sigma_2 = 1,          
             alpha = 0,            
             tau_2 = 1,            
             g = rep(0, nrow(data$coords)),  
             z = rep(0, nrow(data$coords))   
)
priors <- list(beta_mean = 0, 
               beta_sd = 10, 
               beta_prop_sd = 0.075,
               z_mean = 0,
               a_0_sigma = 2,
               b_0_sigma = 1,
               a_0_tau = 2,
               b_0_tau = 1,
               phi = 0.1,
               z_var = 0.2)

iters <- 10000
burnin <- 3000

# Run driver -----------------------------------------------------------------------
sim <- driver(parameters, priors, data, iters)

save(sim, file = "application_sim_res20.RData")

beta_post <- sim$beta[, (burnin+1):iters]
alpha_post <- sim$alpha[,(burnin+1):iters]
sigma_2_post <- sim$sigma_2[,(burnin+1):iters]
tau_2_post <- sim$tau_2[,(burnin+1):iters]
g_post <- sim$g[,(burnin+1):iters]
g_post_water <- g_post[grid_water$id,]
z_post <- sim$z[,(burnin+1):iters]
z_post_water <- z_post[grid_water$id,]

mean(beta_post)
sd(beta_post)

mean(sigma_2_post)
sd(sigma_2_post)

mean(tau_2_post)
sd(tau_2_post)

mean(alpha_post)
sd(alpha_post)

mean(g_post_water)
sd(g_post_water)

mean(z_post_water)
sd(z_post_water)


# Posterior Plots ----------------------------------------------------------------
# Posterior intensity matrix
posterior_lambda <- matrix(NA, nrow = nrow(grid_coords), ncol = (iters - burnin))

for (m in 1:(iters - burnin)) {
  beta_m <- beta_post[m]
  log_lambda_m <- beta_m + z_post[, m] + g_post[, m]
  posterior_lambda[, m] <- exp(log_lambda_m)
}

# Posterior mean intensity
lambda_mean <- rowMeans(posterior_lambda)

lambda_mean[grid_land$id] <- 1e-3

# Reshape to spatial grid
lambda_mean_mat <- matrix(lambda_mean,
                          nrow = length(x_seq),
                          ncol = length(y_seq),
                          byrow = TRUE)

# Plot posterior mean intensity (on log scale)
par(mfrow = c(2,2))
image.plot(x_seq, y_seq, log(t(lambda_mean_mat)),
           main = "Posterior Mean Intensity (log scale)",
           xlab = "x", ylab = "y",
           col = terrain.colors(100), 
           xlim = x_bound, ylim = y_bound)
  map("state", "massachusetts", add = TRUE)
plot(auto_water_xy$X, auto_water_xy$Y, main = "Automated Locations", 
     xlim = x_bound, ylim = y_bound)
  map("state", "massachusetts", add = TRUE)
plot(manual_water_xy$X, manual_water_xy$Y, main = "Manual Locations", 
     xlim = x_bound, ylim = y_bound)
  map("state", "massachusetts", add = TRUE)

par(mfrow = c(1,1))

# Posterior g map ---------------------------------------------------------------
posterior_g <- matrix(NA, nrow = nrow(grid_coords), ncol = (iters - burnin))

for (m in 1:(iters - burnin)) {
  log_g_m <-  g_post[, m]
  posterior_g[, m] <- exp(log_g_m)
}

# Posterior mean intensity
g_mean <- rowMeans(posterior_g)

g_mean[grid_land$id] <- 1e-3

# Reshape to spatial grid
g_mean_mat <- matrix(g_mean,
                          nrow = length(x_seq),
                          ncol = length(y_seq),
                          byrow = TRUE)
zmin <- min(c(log(g_mean_mat), log(lambda_mean_mat), na.rm = TRUE))
zmax <- max(c(log(g_mean_mat), log(lambda_mean_mat), na.rm = TRUE))

# Plot posterior mean intensity (on log scale)
par(mfrow = c(2,2))
image.plot(x_seq, y_seq, log(t(g_mean_mat)),
           main = "Posterior g intensity (linking function)",
           xlab = "x", ylab = "y",
           col = terrain.colors(100), 
           zlim = c(zmin, zmax),
           xlim = x_bound, ylim = y_bound)

image.plot(x_seq, y_seq, log(t(lambda_mean_mat)),
           main = "Posterior Mean Intensity",
           xlab = "x", ylab = "y",
           col = terrain.colors(100), 
           zlim = c(zmin, zmax),
           xlim = x_bound, ylim = y_bound)
par(mfrow = c(1,1))

# Trace Plots -------------------------------------------------------------------
par(mfrow = c(1,1))
plot(sim$beta[1,], type = "l", main = "Beta 1 Trace Plot")
plot(sim$z[1,], type = "l", main = "z trace plot")
plot(sim$g[1,], type = "l", main = "g trace plot")
plot(sim$sigma_2[1,], type = "l", main = "sigma_2 trace plot")
plot(sim$tau_2[1,], type = "l", main = "tau_2 trace plot")
plot(sim$alpha[1,], type = "l", main = "alpha trace plot")

