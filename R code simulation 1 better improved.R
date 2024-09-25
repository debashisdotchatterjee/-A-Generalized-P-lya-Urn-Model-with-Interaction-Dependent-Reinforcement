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

# Generalized Pólya Urn Simulation Function
simulate_polya_urn <- function(K, N, delta_func, X0, max_count = 1e6) {
  # Initialize variables
  X <- matrix(0, nrow = N+1, ncol = K)
  X[1, ] <- X0
  Pn <- matrix(0, nrow = N+1, ncol = K)
  Pn[1, ] <- X0 / sum(X0)
  
  # Store drawn colors
  Cn <- integer(N)
  
  for (n in 1:N) {
    # Current composition
    Xn <- X[n, ]
    Tn <- sum(Xn)
    
    # Ensure probabilities are valid
    if (any(is.na(Xn)) || any(Xn < 0) || Tn == 0) {
      warning(paste("Invalid composition at step", n))
      break
    }
    
    Pn[n, ] <- Xn / Tn
    
    # Draw a ball based on current proportions
    if (any(is.na(Pn[n, ])) || any(Pn[n, ] < 0)) {
      warning(paste("Invalid probabilities at step", n))
      break
    }
    
    Cn[n] <- sample(1:K, size = 1, prob = Pn[n, ])
    
    # Calculate reinforcement
    Delta <- delta_func(Cn[n], Xn)
    
    # Ensure Delta is non-negative and finite
    Delta[Delta < 0 | !is.finite(Delta)] <- 0
    
    # Update composition
    Xn_plus1 <- Xn + Delta
    
    # Prevent counts from becoming excessively large
    Xn_plus1[Xn_plus1 > max_count] <- max_count
    
    # Update the matrix
    X[n+1, ] <- Xn_plus1
  }
  
  # Determine the actual number of iterations completed
  actual_N <- min(n, N)
  
  # Trim the matrices to the actual number of iterations
  X <- X[1:(actual_N+1), , drop = FALSE]
  Pn <- Pn[1:actual_N, , drop = FALSE]
  Cn <- Cn[1:actual_N]
  
  # Return results as a list
  return(list(
    X = X,
    Pn = Pn,
    Cn = Cn
  ))
}

# Modified Linear interaction model with increased variability
delta_linear_modified <- function(Cn, Xn) {
  K <- length(Xn)
  # Asymmetric base reinforcement
  a <- matrix(c(1, 2, 1.5,
                1.5, 1, 2,
                2, 1.5, 1), nrow = K, ncol = K, byrow = TRUE)
  # Asymmetric linear dependence with increased variability
  b <- matrix(c(0.5, 0.3, 0.6,
                0.6, 0.5, 0.3,
                0.3, 0.6, 0.5), nrow = K, ncol = K, byrow = TRUE)
  
  # Add a small constant to Xn to prevent zero counts
  Xn_safe <- Xn + 1e-2
  
  Delta <- a[Cn, ] + b[Cn, ] * Xn_safe
  return(Delta)
}

# Simulation parameters
K <- 3  # Number of colors
N <- 2000  # Number of steps (reduced to prevent over-stabilization)
X0 <- c(50, 30, 20)  # Unequal initial composition

# Simulate with modified linear reinforcement
simulation_linear_mod <- simulate_polya_urn(K, N, delta_linear_modified, X0)

# For linear model
a_linear <- matrix(c(1, 2, 1.5,
                     1.5, 1, 2,
                     2, 1.5, 1), nrow = K, ncol = K, byrow = TRUE)
b_linear <- matrix(c(0.5, 0.3, 0.6,
                     0.6, 0.5, 0.3,
                     0.3, 0.6, 0.5), nrow = K, ncol = K, byrow = TRUE)

# Compute average reinforcement for each color
avg_reinforcement_linear <- rowMeans(a_linear + b_linear * mean(X0))

# Predicted limiting proportions (normalized)
predicted_proportions_linear <- avg_reinforcement_linear / sum(avg_reinforcement_linear)

# Plotting function for proportions with predicted proportions
plot_proportions_predicted <- function(Pn, predicted, title) {
  df <- data.frame(
    Step = 1:nrow(Pn),
    P1 = Pn[,1],
    P2 = Pn[,2],
    P3 = Pn[,3]
  )
  df_long <- gather(df, key = "Color", value = "Proportion", -Step)
  
  predicted_df <- data.frame(
    Color = c("P1", "P2", "P3"),
    Proportion = predicted
  )
  
  p <- ggplot(df_long, aes(x = Step, y = Proportion, color = Color)) +
    geom_line() +
    geom_hline(data = predicted_df, aes(yintercept = Proportion, color = Color), linetype = "dashed") +
    labs(title = title, x = "Step", y = "Proportion") +
    theme_minimal()
  
  return(p)
}

# Plot proportions over time with predicted proportions
p_linear_mod <- plot_proportions_predicted(simulation_linear_mod$Pn, predicted_proportions_linear, "Proportions Over Time (Modified Linear Reinforcement)")

# Save plot
ggsave(filename = file.path(output_dir, "Proportions_Linear_Modified.png"), plot = p_linear_mod, width = 8, height = 5)

# Display plot
print(p_linear_mod)

# Analyze convergence
analyze_convergence <- function(Pn) {
  final_proportions <- Pn[nrow(Pn), ]
  return(final_proportions)
}

final_prop_linear_mod <- analyze_convergence(simulation_linear_mod$Pn)

# Print final proportions and predicted proportions
cat("Final Proportions (Modified Linear Reinforcement):\n")
print(round(final_prop_linear_mod, 4))
cat("Predicted Proportions (Approximate):\n")
print(round(predicted_proportions_linear, 4))

# Validate CLT for scaled fluctuations with theoretical variance
validate_clt_improved <- function(simulation, title) {
  X <- simulation$X
  Pn <- simulation$Pn
  N <- nrow(Pn)
  
  # Total number of balls at each step
  Tn <- rowSums(X[1:N, ])
  
  # Ensure that Pn is finite and has no NaN values
  if (any(!is.finite(Pn)) || any(is.na(Pn))) {
    warning("Non-finite or NA values detected in Pn.")
    return(NULL)
  }
  
  # Compute scaled fluctuations (ensure Pn[N, ] is valid)
  Pn_last <- Pn[N, ]
  if (any(!is.finite(Pn_last)) || any(is.na(Pn_last))) {
    warning("Non-finite or NA values detected in the last row of Pn.")
    return(NULL)
  }
  
  Zn <- sqrt(Tn) * (Pn - matrix(Pn_last, nrow = N, ncol = K, byrow = TRUE))
  
  # Remove rows with non-finite values in Zn
  valid_rows <- apply(Zn, 1, function(row) all(is.finite(row)))
  Zn <- Zn[valid_rows, , drop = FALSE]
  Tn <- Tn[valid_rows]
  Pn <- Pn[valid_rows, , drop = FALSE]
  N <- nrow(Pn)
  
  # Check if Zn has valid rows after removing non-finite values
  if (N == 0) {
    warning("No valid data remaining after removing non-finite values.")
    return(NULL)
  }
  
  # Take last 1000 steps for analysis
  n_steps <- min(1000, N)
  Zn_last <- Zn[(N - n_steps + 1):N, , drop = FALSE]
  
  # Check for zero variance in columns
  sd_Zn_last <- apply(Zn_last, 2, sd)
  zero_sd_cols <- which(sd_Zn_last == 0)
  
  if (length(zero_sd_cols) == ncol(Zn_last)) {
    warning("All columns in Zn_last have zero variance. Cannot perform CLT validation.")
    return(NULL)
  } else if (length(zero_sd_cols) > 0) {
    # Remove columns with zero variance
    Zn_last <- Zn_last[, -zero_sd_cols, drop = FALSE]
    # Update column names accordingly
    col_names <- paste0("Z", setdiff(1:K, zero_sd_cols))
  } else {
    col_names <- paste0("Z", 1:K)
  }
  
  # Check if any columns are left
  if (ncol(Zn_last) == 0) {
    warning("No columns left in Zn_last after removing zero variance columns. Cannot perform CLT validation.")
    return(NULL)
  }
  
  # Standardize fluctuations
  Zn_std <- scale(Zn_last)
  
  # Convert to data frame
  df_Zn <- as.data.frame(Zn_std)
  names(df_Zn) <- col_names
  
  # Gather data for plotting
  df_Zn_long <- gather(df_Zn, key = "Color", value = "Fluctuation")
  
  # Remove any remaining non-finite values
  df_Zn_long <- df_Zn_long[is.finite(df_Zn_long$Fluctuation), ]
  
  # Check if df_Zn_long is empty
  if (nrow(df_Zn_long) == 0) {
    warning("No valid data to plot after removing non-finite values.")
    return(NULL)
  }
  
  # Plot histograms of standardized fluctuations
  p_hist <- ggplot(df_Zn_long, aes(x = Fluctuation, fill = Color)) +
    geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
    facet_wrap(~ Color, scales = "fixed") +
    labs(title = paste("Histogram of Standardized Scaled Fluctuations -", title), x = "Fluctuation", y = "Frequency") +
    theme_minimal() +
    guides(fill = FALSE)
  
  # Q-Q plots
  p_qq <- ggplot(df_Zn_long, aes(sample = Fluctuation)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ Color, scales = "fixed") +
    labs(title = paste("Q-Q Plot of Standardized Scaled Fluctuations -", title), x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Save plots
  ggsave(filename = file.path(output_dir, paste0("Histogram_Fluctuations_", gsub(" ", "_", title), ".png")), plot = p_hist, width = 8, height = 5)
  ggsave(filename = file.path(output_dir, paste0("QQPlot_Fluctuations_", gsub(" ", "_", title), ".png")), plot = p_qq, width = 8, height = 5)
  
  # Display plots
  grid.arrange(p_hist, p_qq, ncol = 1)
}

# Validate CLT for modified linear reinforcement
validate_clt_improved(simulation_linear_mod, "Modified Linear Reinforcement")

# Additional Analysis: Cumulative counts over time
plot_cumulative_counts <- function(X, title) {
  df <- data.frame(
    Step = 1:nrow(X),
    Count1 = X[,1],
    Count2 = X[,2],
    Count3 = X[,3]
  )
  df_long <- gather(df, key = "Color", value = "Count", -Step)
  
  p <- ggplot(df_long, aes(x = Step, y = Count, color = Color)) +
    geom_line() +
    labs(title = title, x = "Step", y = "Cumulative Count") +
    theme_minimal()
  
  return(p)
}

p_cumulative_counts <- plot_cumulative_counts(simulation_linear_mod$X, "Cumulative Counts Over Time (Modified Linear Reinforcement)")

# Save plot
ggsave(filename = file.path(output_dir, "Cumulative_Counts_Linear_Modified.png"), plot = p_cumulative_counts, width = 8, height = 5)

# Display plot
print(p_cumulative_counts)

# Additional Analysis: Reinforcement over time
plot_reinforcement_over_time <- function(simulation, title) {
  Delta_history <- simulation$X[-1, ] - simulation$X[-nrow(simulation$X), ]
  df <- data.frame(
    Step = 1:nrow(Delta_history),
    Delta1 = Delta_history[,1],
    Delta2 = Delta_history[,2],
    Delta3 = Delta_history[,3]
  )
  df_long <- gather(df, key = "Color", value = "Reinforcement", -Step)
  
  p <- ggplot(df_long, aes(x = Step, y = Reinforcement, color = Color)) +
    geom_line(alpha = 0.5) +
    labs(title = title, x = "Step", y = "Reinforcement") +
    theme_minimal()
  
  return(p)
}

p_reinforcement <- plot_reinforcement_over_time(simulation_linear_mod, "Reinforcement Over Time (Modified Linear Reinforcement)")

# Save plot
ggsave(filename = file.path(output_dir, "Reinforcement_Linear_Modified.png"), plot = p_reinforcement, width = 8, height = 5)

# Display plot
print(p_reinforcement)

# Additional Analysis: Table of Final Counts and Proportions
final_counts <- simulation_linear_mod$X[nrow(simulation_linear_mod$X), ]
final_proportions <- simulation_linear_mod$Pn[nrow(simulation_linear_mod$Pn), ]
predicted_proportions <- predicted_proportions_linear

result_table <- data.frame(
  Color = c("Color 1", "Color 2", "Color 3"),
  Final_Counts = round(final_counts),
  Simulated_Proportions = round(final_proportions, 4),
  Predicted_Proportions = round(predicted_proportions, 4)
)

print("Final Counts and Proportions:")
print(result_table)

# Save table to CSV
write.csv(result_table, file = file.path(output_dir, "Final_Counts_Proportions.csv"), row.names = FALSE)

# End of R code
