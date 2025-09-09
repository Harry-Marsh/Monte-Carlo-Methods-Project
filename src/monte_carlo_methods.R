# Monte Carlo Simulation Project
# Candidate number: 248605
# Part 1: Rejection Sampling and Importance Sampling
# Part 2: Markov Chain Simulation and MCMC
# Part 3: MCS in current research topics - Sonawane PDEs

##############################################################################

# Part 1
#Setup
set.seed(123)  # for reproducibility
alpha <- 2.5
k <- 1
mu <- 0

# Rejection Sampling Function
generate_custom_dist <- function(n, alpha, k, mu) {
  # n: number of samples to generate
  # alpha, k, mu: parameters of the target PDF
  # Target PDF: f(x) = (alpha/k) * ((k-mu+x)/k)^(-alpha-1), for x ≥ mu
  
  accepted <- numeric(0)
  c <- 10 * alpha            # c chosen so that f(x) ≤ c * g(x)
  g_density <- 1 / (10 * k)  # Uniform density ~ U[mu, mu+10k]
  
  while (length(accepted) < n) {
    Y <- runif(1, mu, mu + 10 * k)
    f_Y <- (alpha / k) * ((k - mu + Y)/k)^(-alpha - 1)
    U <- runif(1)
    
    # Accept with probability below
    if (U <= f_Y / (c * g_density)) {
      accepted <- c(accepted, Y)
    }
  }
  
  return(accepted)
}

# Generate samples
samples <- generate_custom_dist(n = 10000, alpha, k, mu)

# Plot Histogram with Theoretical PDF 
plot_hist_with_pdf <- function(samples, alpha, k, mu) {
  hist(samples, breaks = 50, probability = TRUE,
       main = "Histogram with Theoretical PDF",
       xlab = "x", col = "lightgray", border = "white")
  
  curve((alpha / k) * ((k - mu + x)/k)^(-alpha - 1),
        from = min(samples), to = max(samples),
        col = "red", lwd = 2, add = TRUE)
  
  legend("topright", legend = "Theoretical PDF", col = "red", lwd = 2)
}

plot_hist_with_pdf(samples, alpha, k, mu)

# Goodness-of-Fit (KS Test) 
cdf_theoretical <- function(x, alpha, k, mu) {
  ifelse(x < mu, 0, 1 - (k / (k - mu + x))^alpha)
}

ks_analysis <- function(sizes, alpha, k, mu) {
  results <- data.frame(SampleSize = sizes, D = NA, p_value = NA)
  
  for (i in seq_along(sizes)) {
    x <- generate_custom_dist(n = sizes[i], alpha, k, mu)
    ks <- suppressWarnings(ks.test(x, function(x) cdf_theoretical(x, alpha, k, mu)))
    results$D[i] <- ks$statistic
    results$p_value[i] <- ks$p.value
  }
  
  return(results)
}

ks_results <- ks_analysis(c(100, 500, 1000, 5000, 10000), alpha, k, mu)
print(ks_results)

# Importance Sampling for E[|X|]
importance_sampling_absX <- function(n, alpha, k, mu) {
  proposal <- runif(n, mu, mu + 10 * k)
  
  f_x <- (alpha / k) * ((k - mu + proposal)/k)^(-alpha - 1)
  g_x <- dunif(proposal, mu, mu + 10 * k)
  
  # For each sample size, generate data and perform a KS test
  weights <- f_x / g_x
  estimate <- sum(abs(proposal) * weights) / sum(weights)
  
  return(estimate)
}

importance_estimate <- importance_sampling_absX(n = 10000, alpha, k, mu)
print(paste("Estimated E[|X|] via importance sampling:", round(importance_estimate, 4)))

#################################################################################

#Part 2 


# Setup
set.seed(514)  
n_states <- 5

# Generate Transition Matrix
transition_matrix <- matrix(runif(n_states * n_states), nrow = n_states)
transition_matrix <- t(apply(transition_matrix, 1, function(x) x / sum(x)))


# Compute stationary distribution by solving (pi)P = pi
# Equivalent to finding the left eigenvector for eigenvalue 1

eigen_result <- eigen(t(transition_matrix))
stationary <- Re(eigen_result$vectors[,1])
stationary <- stationary / sum(stationary)

# Simulate Markov Chain 
simulate_chain <- function(P, n_steps, init_state = 1) {
  states <- numeric(n_steps)
  states[1] <- init_state
  for (t in 2:n_steps) {
    states[t] <- sample(1:nrow(P), size = 1, prob = P[states[t - 1], ])
  }
  return(states)
}

samples <- simulate_chain(transition_matrix, n_steps = 100000)

# Empirical Distribution 
empirical_dist <- table(samples) / length(samples)

# Compare Empirical vs Theoretical
print("Stationary distribution (theoretical):")
print(round(stationary, 4))

print("Empirical distribution (from simulation):")
print(round(empirical_dist, 4))

##############################################################################

# Part 3: 

rm(list = ls())
set.seed(514)

# 1D solver
mc1d <- function(N = 500, n = 31, method = "fixed") {
  # N: number of random walks per grid point
  # n: number of grid points (solution resolution)
  # method: "fixed" or "full" step size
  
  h <- 1 / (n - 1)
  xs <- seq(0, 1, length.out = n)
  u  <- numeric(n)  # store mean values for each x
  
  for (i in seq_len(n)) {
    x0 <- xs[i]
    hits <- numeric(N)
    
    for (k in seq_len(N)) {
      x <- x0
      while (TRUE) {
        x <- x + if (method == "fixed") sample(c(-h, h), 1) else runif(1, -h, h)
        if (x <= 0) { break }
        if (x >= 1) { hits[k] <- 1; break }
      }
    }
    u[i] <- mean(hits)
  }
  
  data.frame(x = xs, u = u)
}

# Parameters for 1D 
grid_pts   <- 31                    # Use 51 for full 
Ns_1d      <- c(200, 500, 1000)     # Full: c(1000, 5000, 10000)
methods_1d <- c("fixed", "full")

errors1d <- expand.grid(method = methods_1d, N = Ns_1d)
errors1d$error <- NA_real_

cat("\n--- 1D Error Table ---\n")
for (r in seq_len(nrow(errors1d))) {
  m <- errors1d$method[r]
  N <- errors1d$N[r]
  cat(sprintf("  %6s, N=%5d ... ", m, N))
  res <- mc1d(N = N, n = grid_pts, method = m)
  errors1d$error[r] <- max(abs(res$u - res$x))  # true solution is u(x)=x
  cat(sprintf("max error = %.5f\n", errors1d$error[r]))
}
print(errors1d)

# 1D log–log convergence plot
plot(log10(errors1d$N), log10(errors1d$error),
     col = as.numeric(factor(errors1d$method)),
     pch = 16, cex = 1.3,
     xlab = expression(log[10](N)), ylab = expression(log[10](E[N])),
     main = "1D Monte Carlo Convergence vs Sample Size"
     , las = 1)
legend("bottomleft", legend = methods_1d, col = 1:2, pch = 16, bty = "n")


# 2D solver [Laplace equation on [0,1]^2]
# This solves Laplace's equation ∇²u = 0 using random walks on a 2D grid
# Boundary condition: u(x=0,y)=0, u(x=1,y)=1, other sides reflect (Neumann)

mc2d <- function(N = 2000, nx = 15, ny = 15) {
  # N: number of walks per interior point
  # nx, ny: grid resolution (reduced: 15×15; full: 21×21)
  
  hx <- 1 / (nx - 1)
  hy <- 1 / (ny - 1)
  xs <- seq(0, 1, length.out = nx)
  ys <- seq(0, 1, length.out = ny)
  U  <- matrix(NA_real_, nx, ny)
  
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      x0 <- xs[i]; y0 <- ys[j]
      hits <- numeric(N)
      
      for (k in seq_len(N)) {
        x <- x0; y <- y0
        while (TRUE) {
          step <- sample(1:4, 1)
          if (step == 1) x <- x + hx
          if (step == 2) x <- x - hx
          if (step == 3) y <- y + hy
          if (step == 4) y <- y - hy
          
          if (y < 0) y <- -y
          if (y > 1) y <- 2 - y
          
          if (x <= 0) { break }
          if (x >= 1) { hits[k] <- 1; break }
        }
      }
      U[i, j] <- mean(hits)
    }
  }
  
  list(x = xs, y = ys, u = U)
}

# Parameters for 2D
nx <- ny <- 15                     # Full: nx = ny = 21
Ns_2d <- c(1000, 2000)             # Full: Ns_2d = c(3000, 10000)

errors2d <- data.frame(N = Ns_2d, error = NA_real_)

cat("\n--- 2D Error Table ---\n")
for (i in seq_along(Ns_2d)) {
  N <- Ns_2d[i]
  cat(sprintf("  N=%5d ... ", N))
  res2 <- mc2d(N = N, nx = nx, ny = ny)
  Utrue <- outer(res2$x, res2$y, function(x, y) x)  # u(x,y) = x
  errors2d$error[i] <- max(abs(res2$u - Utrue))
  cat(sprintf("max error = %.5f\n", errors2d$error[i]))
}
print(errors2d)

cat("\n Tables above and 1D plot displayed.\n")
