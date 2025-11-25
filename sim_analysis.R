# Packages
library(spatstat)
library(fields)
library(FastGP)
library(MASS)
library(ggplot2)
library(coda)

# Individual sim analysis: ------------------------------------------------------
load("sim_full.RData")
load("sim_source1.RData")
load("sim_source2.RData")
load("sim_data.RData")
load("sim_parameters.RData")
load("sim_priors.Rdata")
load("sim_covariate.Rdata")
load("sim_lambda.RData")

win <- owin(xrange = c(0, 10), yrange = c(0, 10))

grid_res <- 10

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)
grid_coords <- expand.grid(x = y_seq, y = x_seq)

burnin <- 3000
iters <- 10000
## Mean, standard deviation, and 95% credible intervals for parameters

# Full sim
beta_post <- sim$beta[, (burnin+1):iters]
alpha_post <- sim$alpha[,(burnin+1):iters]
sigma_2_post <- sim$sigma_2[,(burnin+1):iters]
tau_2_post <- sim$tau_2[,(burnin+1):iters]
g_post <- sim$g[,(burnin+1):iters]
z_post <- sim$z[,(burnin+1):iters]

apply(beta_post, 1, mean)
apply(beta_post, 1, sd)
quantile(beta_post[1,], c(0.025, 0.975))
quantile(beta_post[2,], c(0.025, 0.975))

mean(sigma_2_post)
sd(sigma_2_post)
quantile(sigma_2_post, c(.025,0.975))

mean(tau_2_post)
sd(tau_2_post)
quantile(tau_2_post, c(.025,0.975))

mean(alpha_post)
sd(alpha_post)
quantile(alpha_post, c(.025,0.975))

mean(g_post)
sd(g_post)
quantile(g_post, c(.025,0.975))

mean(z_post)
sd(z_post)
quantile(z_post, c(.025,0.975))

# Source 1 Sim
beta_post1 <- sim_source1$beta[, (burnin+1):iters]
alpha_post1 <- sim_source1$alpha[,(burnin+1):iters]
sigma_2_post1 <- sim_source1$sigma_2[,(burnin+1):iters]
tau_2_post1 <- sim_source1$tau_2[,(burnin+1):iters]
g_post1 <- sim_source1$g[,(burnin+1):iters]
z_post1 <- sim_source1$z[,(burnin+1):iters]

apply(beta_post1, 1, mean)
apply(beta_post1, 1, sd)
quantile(beta_post1[1,], c(0.025, 0.975))
quantile(beta_post1[2,], c(0.025, 0.975))

mean(sigma_2_post1)
sd(sigma_2_post1)
quantile(sigma_2_post1, c(.025,0.975))

mean(z_post1)
sd(z_post1)
quantile(z_post1, c(.025,0.975))


# Source 2 Sim
beta_post2 <- sim_source2$beta[, (burnin+1):iters]
alpha_post2 <- sim_source2$alpha[,(burnin+1):iters]
sigma_2_post2 <- sim_source2$sigma_2[,(burnin+1):iters]
tau_2_post2 <- sim_source2$tau_2[,(burnin+1):iters]
g_post2 <- sim_source2$g[,(burnin+1):iters]
z_post2 <- sim_source2$z[,(burnin+1):iters]

apply(beta_post2, 1, mean)
apply(beta_post2, 1, sd)
quantile(beta_post2[1,], c(0.025, 0.975))
quantile(beta_post2[2,], c(0.025, 0.975))

mean(sigma_2_post2)
sd(sigma_2_post2)
quantile(sigma_2_post2, c(.025,0.975))

mean(z_post2)
sd(z_post2)
quantile(z_post2, c(.025,0.975))

## Posterior Plots

# Posterior lambda for both sources
posterior_lambda <- matrix(NA, nrow = 100, ncol = (iters-burnin))
    
    for(m in 1:(iters-burnin)){
      beta_m <- beta_post[,m]
      
      log_lambda_m <- beta_m[1] + beta_m[2]*covariate + z_post[,m]
      
      posterior_lambda[, m] <- exp(log_lambda_m)
    }
    
    # Posterior mean intensity
    lambda_mean <- rowMeans(posterior_lambda, na.rm = TRUE)
    
# Just source 1 
posterior_lambda_source1 <- matrix(NA, nrow = nrow(data$X_grid), ncol = (iters-burnin))
    
    for(m in 1:(iters-burnin)){
      beta_m1 <- beta_post1[,m]
      
      log_lambda_m1 <- beta_m1[1] + beta_m1[2]*covariate + z_post1[,m]
      
      posterior_lambda_source1[, m] <- exp(log_lambda_m1)
    }
    
    # Posterior mean intensity
    lambda_mean1 <- rowMeans(posterior_lambda_source1, na.rm = TRUE)
    
    # Reshape to grid
    lambda_mean_mat1 <- matrix(lambda_mean1, 
                               nrow = data$grid_res, 
                               ncol = data$grid_res, 
                               byrow = FALSE)
    
# Just source 2
posterior_lambda_source2 <- matrix(NA, nrow = nrow(data$X_grid), ncol = (iters-burnin))
    
    for(m in 1:(iters-burnin)){
      beta_m2 <- beta_post2[,m]
      
      log_lambda_m2 <- beta_m2[1] + beta_m2[2]*covariate + z_post2[,m]
      
      posterior_lambda_source2[, m] <- exp(log_lambda_m2)
    }
    
    # Posterior mean intensity
    lambda_mean2 <- rowMeans(posterior_lambda_source2, na.rm = TRUE)
    
    # Reshape to grid
    lambda_mean_mat2 <- matrix(lambda_mean2, 
                               nrow = data$grid_res, 
                               ncol = data$grid_res, 
                               byrow = FALSE)

# Plotting
par(mfrow = c(2,2))
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda), 
           nx = 10,
           ny = 10,
           col = terrain.colors(100),
           main = "True Intensity")

quilt.plot(grid_coords$y,
           grid_coords$x,
           log(lambda_mean),
           nx = 10, 
           ny = 10, 
           col = terrain.colors(100),
           main = "Posterior Mean Fused")

quilt.plot(grid_coords$y,
           grid_coords$x,
           log(lambda_mean1),
           nx = 10, 
           ny = 10, 
           col = terrain.colors(100),
           main = "Posterior Mean Source 1 Only")

quilt.plot(grid_coords$y,
           grid_coords$x,
           log(lambda_mean2),
           nx = 10, 
           ny = 10, 
           col = terrain.colors(100),
           main = "Posterior Mean Source 2 Only")
par(mfrow = c(1,1))

# Posterior g map ---------------------------------------------------------------
posterior_g <- matrix(NA, nrow = nrow(data$X_grid), ncol = (iters - burnin))
    
  for (m in 1:(iters - burnin)) {
    log_g_m <-  g_post[, m]
    posterior_g[, m] <- exp(log_g_m)
  }
    
  # Posterior mean intensity
  g_mean <- rowMeans(posterior_g)

  # Plot posterior mean intensity (on log scale)
  par(mfrow = c(2,2))
    quilt.plot(grid_coords$y, 
               grid_coords$x, 
               log(g_mean),
               nx = 10,
               ny = 10, 
               col = topo.colors(100))
    
    quilt.plot(grid_coords$y, 
               grid_coords$x, 
               log(lambda_mean),
               nx = 10,
               ny = 10, 
               col = terrain.colors(100),
               main = "Posterior Mean Intensity")
par(mfrow = c(1,1))
    
# Multiple simulations for credible intervals ------------------------------------
load("n_sims_df1.RData")
load("n_sims_df2.RData")
#Beta 

beta0_covered <- (n_sims_df$beta0_lower <= 1 &
                    n_sims_df$beta0_upper >= 1)

covered_beta0 <- mean(n_sims_df$beta0_lower <= 1 &
                        n_sims_df$beta0_upper >= 1)

beta1_covered <- (n_sims_df$beta1_lower <= 3 &
                    n_sims_df$beta1_upper >= 3)

covered_beta1 <- mean(n_sims_df$beta1_lower <= 3 &
                        n_sims_df$beta1_upper >= 3)
covered_beta0
covered_beta1

#sigma_2

covered_sigma2 <- mean(n_sims_df$sigma2_lower <= 0.25 &
                         n_sims_df$sigma2_upper >=0.25)
covered_sigma2

#tau_2
covered_tau2 <- mean(n_sims_df$tau2_lower <= 0.4 &
                       n_sims_df$tau2_upper >=0.4)
covered_tau2

#alpha

covered_alpha <- mean(n_sims_df$alpha_lower <= -0.2 &
                        n_sims_df$alpha_upper >= -0.2)
covered_alpha

#covered count

count_lower <- n_sims_df$expected_count - 1.96 * n_sims_df$expected_count_sd
count_upper <- n_sims_df$expected_count + 1.96 * n_sims_df$expected_count_sd
covered_count <- mean(count_lower <= n_sims_df$actual_count &
      count_upper >= n_sims_df$actual_count)

covered_count
