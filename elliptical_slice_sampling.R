# Creating log-likelihood function for multivariate normal

set.seed(111)
identity_matrix <- diag(10)
sigma_2 <- 0.5
S <- sigma_2 * identity_matrix
mu <- rep(0, 10)

f <- MASS::mvrnorm(n = 1, mu = mu, Sigma = S)

# instead of mvn do poisson

loglike <- function(beta, f, X, X_grid, cell_area) {
  log_lambda_points <- beta[1] + beta[2] * X$covariate + f
  term1 <- sum(log_lambda_points)
  
  log_lambda_grid <- beta[1] + beta[2] * X_grid$covariate + f
  lambda_grid <- exp(log_lambda_grid)
  term2 <- sum(lambda_grid * cell_area)
  
  likelihood <- sum(term1 - term2)
  return(likelihood)
}


# Creating the function -----------------------------------------------------------
ess_function <- function(beta, f_current, X, X_grid, cell_area, mu, Sigma){
  
  # Choosing ellipse (nu) from prior (f)
  nu <- mvrnorm(n = 1, mu = rep(0, length(f_current)), Sigma = Sigma)
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike(beta,f_current,X, X_grid, cell_area) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
repeat {
  # Calculate f'
  f_prime <- f_current*cos(theta) + nu*sin(theta)
  
  # Shrinking bracket
  if(loglik(beta, f_prime, X, X_grid, cell_area) > log_y){
    return(f_prime)} else {
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

n_iter <- 10000
samples <- matrix(NA, nrow = n_iter, ncol = length(f))
samples[1, ] <- f

for (i in 2:n_iter) {
  samples[i, ] <- ess_function(beta, samples[i-1, ], X, X_grid, cell_area, mu, Sigma)
}


plot(samples[,1], type = "l")
head(samples)

posterior_mean <- colMeans(samples)
posterior_mean

posterior_sd <- apply(samples, 2, sd)
posterior_sd

## Previous likelihood
loglike <- function(beta, X, X_grid, cell_area) {
  log_lambda_points <- beta[1] + beta[2] * X$covariate
  term1 <- sum(log_lambda_points)
  
  log_lambda_grid <- beta[1] + beta[2] * X_grid$covariate
  lambda_grid <- exp(log_lambda_grid)
  term2 <- sum(lambda_grid * cell_area)
  
  likelihood <- sum(term1 - term2)
  return(likelihood)
}