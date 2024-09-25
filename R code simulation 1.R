# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(MASS)

# Set seed for reproducibility
set.seed(123)

# Create an output directory
output_dir <- "Polya_Urn_Simulation_Outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Generalized Polya Urn Simulation Function
simulate_polya_urn <- function(K, N, delta_func, X0) {
  # K: Number of colors
  # N: Number of steps
  # delta_func: Function to calculate reinforcement
  # X0: Initial composition vector
  
  # Initialize variables
  X <- matrix(0, nrow = N+1, ncol = K)
  X[1, ] <- X0
  Tn <- sum(X0)
  Pn <- matrix(0, nrow = N+1, ncol = K)
  Pn[1, ] <- X0 / Tn
  
  # Store drawn colors
  Cn <- integer(N)
  
  for (n in 1:N) {
    # Current composition
    Xn <- X[n, ]
    Tn <- sum(Xn)
    Pn[n, ] <- Xn / Tn
    
    # Draw a ball based on current proportions
    Cn[n] <- sample(1:K, size = 1, prob = Xn / Tn)
    
    # Calculate reinforcement
    Delta <- delta_func(Cn[n], Xn)
    
    # Update composition
    X[n+1, ] <- Xn + Delta
  }
  
  # Remove the last row (since we simulated N steps)
  X <- X[-(N+1), ]
  Pn <- Pn[-(N+1), ]
  
  # Return results as a list
  return(list(
    X = X,
    Pn = Pn,
    Cn = Cn
  ))
}

# Define reinforcement functions
# Linear interaction model
delta_linear <- function(Cn, Xn) {
  K <- length(Xn)
  a <- matrix(1, nrow = K, ncol = K)  # Base reinforcement
  b <- matrix(0.1, nrow = K, ncol = K)  # Linear dependence
  
  Delta <- a[Cn, ] + b[Cn, ] * Xn
  return(Delta)
}

# Nonlinear interaction model
delta_nonlinear <- function(Cn, Xn) {
  K <- length(Xn)
  lambda <- matrix(0.05, nrow = K, ncol = K)
  alpha <- 0.5
  
  Delta <- lambda[Cn, ] * (Xn ^ alpha)
  return(Delta)
}

# Simulation parameters
K <- 3  # Number of colors
N <- 5000  # Number of steps
X0 <- c(10, 10, 10)  # Initial composition

# Simulate with linear reinforcement
simulation_linear <- simulate_polya_urn(K, N, delta_linear, X0)

# Simulate with nonlinear reinforcement
simulation_nonlinear <- simulate_polya_urn(K, N, delta_nonlinear, X0)

# Plotting function for proportions
plot_proportions <- function(Pn, title) {
  df <- data.frame(
    Step = 1:nrow(Pn),
    P1 = Pn[,1],
    P2 = Pn[,2],
    P3 = Pn[,3]
  )
  df_long <- gather(df, key = "Color", value = "Proportion", -Step)
  
  p <- ggplot(df_long, aes(x = Step, y = Proportion, color = Color)) +
    geom_line() +
    labs(title = title, x = "Step", y = "Proportion") +
    theme_minimal()
  
  return(p)
}

# Plot proportions over time
p_linear <- plot_proportions(simulation_linear$Pn, "Proportions Over Time (Linear Reinforcement)")
p_nonlinear <- plot_proportions(simulation_nonlinear$Pn, "Proportions Over Time (Nonlinear Reinforcement)")

# Save plots
ggsave(filename = file.path(output_dir, "Proportions_Linear.png"), plot = p_linear, width = 8, height = 5)
ggsave(filename = file.path(output_dir, "Proportions_Nonlinear.png"), plot = p_nonlinear, width = 8, height = 5)

# Display plots
grid.arrange(p_linear, p_nonlinear, ncol = 1)

# Analyze convergence
analyze_convergence <- function(Pn) {
  final_proportions <- Pn[nrow(Pn), ]
  return(final_proportions)
}

final_prop_linear <- analyze_convergence(simulation_linear$Pn)
final_prop_nonlinear <- analyze_convergence(simulation_nonlinear$Pn)

# Print final proportions
cat("Final Proportions (Linear Reinforcement):\n")
print(final_prop_linear)
cat("\nFinal Proportions (Nonlinear Reinforcement):\n")
print(final_prop_nonlinear)

# Validate CLT for scaled fluctuations
validate_clt <- function(simulation, title) {
  X <- simulation$X
  Pn <- simulation$Pn
  N <- nrow(X)
  
  # Total number of balls at each step
  Tn <- rowSums(X)
  
  # Compute scaled fluctuations
  Zn <- sqrt(Tn) * (Pn - matrix(Pn[N, ], nrow = N, ncol = K, byrow = TRUE))
  
  # Take last 1000 steps for analysis
  Zn_last <- Zn[(N-999):N, ]
  
  # Plot histograms of fluctuations
  df_Zn <- data.frame(
    Z1 = Zn_last[,1],
    Z2 = Zn_last[,2],
    Z3 = Zn_last[,3]
  )
  
  df_Zn_long <- gather(df_Zn, key = "Color", value = "Fluctuation")
  
  p_hist <- ggplot(df_Zn_long, aes(x = Fluctuation, fill = Color)) +
    geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
    facet_wrap(~ Color, scales = "free") +
    labs(title = paste("Histogram of Scaled Fluctuations -", title), x = "Fluctuation", y = "Frequency") +
    theme_minimal()
  
  # Q-Q plots
  p_qq <- ggplot(df_Zn_long, aes(sample = Fluctuation, color = Color)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ Color, scales = "free") +
    labs(title = paste("Q-Q Plot of Scaled Fluctuations -", title), x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Save plots
  ggsave(filename = file.path(output_dir, paste0("Histogram_Fluctuations_", gsub(" ", "_", title), ".png")), plot = p_hist, width = 8, height = 5)
  ggsave(filename = file.path(output_dir, paste0("QQPlot_Fluctuations_", gsub(" ", "_", title), ".png")), plot = p_qq, width = 8, height = 5)
  
  # Display plots
  grid.arrange(p_hist, p_qq, ncol = 1)
}

# Validate CLT for linear reinforcement
validate_clt(simulation_linear, "Linear Reinforcement")

# Validate CLT for nonlinear reinforcement
validate_clt(simulation_nonlinear, "Nonlinear Reinforcement")

# Application Example: Opinion Dynamics in Social Networks
# Simulate the spread of opinions using the generalized Pólya urn model

# Define reinforcement function for opinion dynamics
delta_opinion <- function(Cn, Xn) {
  K <- length(Xn)
  influence_matrix <- matrix(0.5, nrow = K, ncol = K)
  diag(influence_matrix) <- 1  # Strong self-reinforcement
  
  Delta <- influence_matrix[Cn, ] * Xn / sum(Xn)
  return(Delta)
}

# Initial composition: equal number of supporters for each opinion
X0_opinion <- c(100, 100, 100)

# Simulate opinion dynamics
simulation_opinion <- simulate_polya_urn(K, N, delta_opinion, X0_opinion)

# Plot proportions over time
p_opinion <- plot_proportions(simulation_opinion$Pn, "Opinion Dynamics Over Time")

# Save plot
ggsave(filename = file.path(output_dir, "Opinion_Dynamics.png"), plot = p_opinion, width = 8, height = 5)

# Display plot
print(p_opinion)

# Analyze final proportions
final_prop_opinion <- analyze_convergence(simulation_opinion$Pn)
cat("\nFinal Proportions (Opinion Dynamics):\n")
print(final_prop_opinion)

# Generate a table of final proportions for all simulations
final_proportions_table <- data.frame(
  Reinforcement = c("Linear", "Nonlinear", "Opinion Dynamics"),
  Color1 = c(final_prop_linear[1], final_prop_nonlinear[1], final_prop_opinion[1]),
  Color2 = c(final_prop_linear[2], final_prop_nonlinear[2], final_prop_opinion[2]),
  Color3 = c(final_prop_linear[3], final_prop_nonlinear[3], final_prop_opinion[3])
)

# Save table to CSV
write.csv(final_proportions_table, file.path(output_dir, "Final_Proportions.csv"), row.names = FALSE)

# Print table
print(final_proportions_table)

# Save all outputs in the output directory
# (All plots and tables are already saved)

# End of R code
