# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(xtable)
library(stats)

# Create output directory
output_dir <- "Bitcoin_Trust_Analysis"
if (!dir.exists(output_dir)) dir.create(output_dir)

###############################################
# Step 1: Load and Prepare the Data
###############################################

# Update this path as needed
data_path <- "C:/Users/User/Desktop/Debashis 2024/Bitcoin Python/soc-sign-bitcoinalpha.csv"
soc.sign.bitcoinalpha <- read.csv(data_path, 
                                  header = FALSE, 
                                  col.names = c("SOURCE", "TARGET", "RATING", "TIME"))

data <- soc.sign.bitcoinalpha

# Convert TIME to human-readable date using the Unix epoch
data$DATE <- as.POSIXct(data$TIME, origin = "1970-01-01", tz = "UTC")

# Save the cleaned dataset with readable dates
write.csv(data, file = file.path(output_dir, "cleaned_dataset.csv"), row.names = FALSE)

###############################################
# Step 2: Exploratory Data Analysis
###############################################

# Basic Statistical Attributes
basic_stats <- data %>%
  summarize(
    Total_Ratings = n(),
    Unique_Sources = n_distinct(SOURCE),
    Unique_Targets = n_distinct(TARGET),
    Min_Rating = min(RATING),
    Max_Rating = max(RATING),
    Avg_Rating = mean(RATING),
    Median_Rating = median(RATING),
    Std_Rating = sd(RATING),
    Start_Date = min(DATE),
    End_Date = max(DATE)
  )

write.csv(basic_stats, file = file.path(output_dir, "basic_statistics.csv"), row.names = FALSE)

# Print Basic Statistics
print("Basic Statistics:")
print(basic_stats)

# Plot 1: Distribution of Ratings
plot1 <- ggplot(data, aes(x = RATING)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  ggtitle("Distribution of Ratings in the Bitcoin Alpha Network") +
  xlab("Rating") +
  ylab("Frequency") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "distribution_of_ratings.png"), plot1, width=8, height=6)

# Plot 2: Ratings Over Time (Jitter)
plot2 <- ggplot(data, aes(x = DATE, y = RATING)) +
  geom_jitter(alpha = 0.5, color = "blue", size=0.7) +
  ggtitle("Ratings Over Time") +
  xlab("Date") +
  ylab("Rating") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "ratings_over_time.png"), plot2, width=8, height=6)

###############################################
# Step 3: Implementing the Pólya Urn Trust Model
###############################################

# We consider a simplified model:
# Each node is a color. Initially, assign each node a small positive weight.
# At each rating event: SOURCE rates TARGET with RATING.
# Probability of this event (TARGET chosen) under the model is proportional to the target's weight.
# Update the target's weight: X_{new}(target) = X_{old}(target) + alpha + beta * RATING.

# Key challenges:
# - The dataset is large, we will use a subset for demonstration.
# - Parameter estimation by MLE involves computing the probability of each observed TARGET given the history.

# For demonstration, let's take the first 5000 events sorted by time.
data_sub <- data %>%
  arrange(TIME) %>%
  slice(1:5000)

# Extract unique nodes and map them to indices
all_nodes <- unique(c(data_sub$SOURCE, data_sub$TARGET))
node_index <- seq_along(all_nodes)
names(node_index) <- all_nodes

# Map SOURCE and TARGET to indices
data_sub$SOURCE_ID <- node_index[as.character(data_sub$SOURCE)]
data_sub$TARGET_ID <- node_index[as.character(data_sub$TARGET)]

N <- length(all_nodes) # number of unique nodes in the subset
events <- nrow(data_sub)

# Initialize weights
# To ensure no zero probabilities at start, give each node a small baseline weight.
initial_weight <- 1
weights <- rep(initial_weight, N)
total_weight <- sum(weights)

# Define the log-likelihood function given alpha and beta:
# LL(alpha, beta) = sum over all events of log( Probability(TARGET chosen at that event) )
# Probability = X_{n,TARGET} / T_n
#
# After each event:
# X_{n+1,TARGET} = X_{n,TARGET} + alpha + beta*RATING

log_likelihood_fn <- function(params, data_events, N, init_weights) {
  alpha <- params[1]
  beta <- params[2]
  
  # Enforce positivity constraint: alpha + beta*RATING > 0 for min rating = -10
  # One simple approach: if not positive, return large negative LL
  # Ensure alpha + beta*(-10) > 0 => alpha + (-10)*beta > 0 => alpha > 10*beta if beta>0
  # We'll just penalize if any increment is negative.
  if (alpha + beta * (-10) <= 0) {
    return(-1e15) # large negative to reject these parameters
  }
  
  w <- init_weights
  T_w <- sum(w)
  
  ll <- 0
  for (i in seq_len(nrow(data_events))) {
    target_id <- data_events$TARGET_ID[i]
    rating <- data_events$RATING[i]
    
    # Probability of choosing this target
    p_target <- w[target_id] / T_w
    if (p_target <= 0) return(-1e15)
    
    ll <- ll + log(p_target)
    
    # Update weights
    increment <- alpha + beta * rating
    # If increment <= 0, penalize
    if (increment <= 0) {
      return(-1e15)
    }
    w[target_id] <- w[target_id] + increment
    T_w <- T_w + increment
  }
  
  return(ll)
}

# To estimate alpha and beta, we optimize the log-likelihood.
# We know alpha > 0 and beta > 0.
# We'll start with some initial guesses and use `optim`.

init_params <- c(alpha = 1, beta = 0.1)  # initial guess
lower_bounds <- c(0, 0)
# We can use a box-constrained optimization (method="L-BFGS-B")

fit <- optim(
  par = init_params,
  fn = function(par) -log_likelihood_fn(par, data_sub, N, rep(initial_weight, N)),
  lower = lower_bounds,
  method = "L-BFGS-B",
  control = list(maxit = 100, fnscale = 1)
)

fit

fitted_alpha <- fit$par[1]
fitted_beta <- fit$par[2]

fitted_params <- data.frame(Parameter=c("alpha","beta"),
                            Estimate=c(fitted_alpha, fitted_beta))
write.csv(fitted_params, file.path(output_dir, "fitted_parameters.csv"), row.names=FALSE)

print("Fitted Parameters:")
print(fitted_params)

###############################################
# Step 4: Model Fit Evaluation
###############################################

# With the estimated alpha and beta, let's compute the log-likelihood and
# see how well the model explains the observed sequence.

final_ll <- log_likelihood_fn(c(fitted_alpha, fitted_beta), data_sub, N, rep(initial_weight, N))
print(paste("Final Log-Likelihood:", final_ll))

# AIC and BIC as simple model selection criteria
k <- 2 # number of parameters
AIC_val <- -2*final_ll + 2*k
BIC_val <- -2*final_ll + k*log(events)

model_fit_stats <- data.frame(
  LogLikelihood = final_ll,
  AIC = AIC_val,
  BIC = BIC_val
)
write.csv(model_fit_stats, file.path(output_dir, "model_fit_stats.csv"), row.names=FALSE)

print("Model Fit Statistics:")
print(model_fit_stats)

###############################################
# Step 5: Analyzing the Implied Dynamics
###############################################

# With fitted alpha and beta, we can simulate the trajectory of node weights and identify influential nodes.
# Let's rerun through the data and track weights over time.

w <- rep(initial_weight, N)
T_w <- sum(w)

time_series <- data.frame(Event=integer(),
                          Target=integer(),
                          Rating=integer(),
                          Chosen_Prob=numeric(),
                          stringsAsFactors=FALSE)

for (i in seq_len(nrow(data_sub))) {
  target_id <- data_sub$TARGET_ID[i]
  rating <- data_sub$RATING[i]
  
  p_target <- w[target_id] / T_w
  time_series <- rbind(time_series, 
                       data.frame(Event=i, Target=target_id, Rating=rating, Chosen_Prob=p_target))
  
  increment <- fitted_alpha + fitted_beta * rating
  w[target_id] <- w[target_id] + increment
  T_w <- T_w + increment
}

# Identify top nodes by final weight
final_weights <- data.frame(Node=all_nodes, Final_Weight=w)
top_nodes <- final_weights %>%
  arrange(desc(Final_Weight)) %>%
  slice(1:10)

write.csv(top_nodes, file.path(output_dir, "top_influential_nodes.csv"), row.names=FALSE)

print("Top 10 Influential Nodes by Final Weight:")
print(top_nodes)

# Plot the distribution of chosen probabilities over events
plot_prob <- ggplot(time_series, aes(x=Event, y=Chosen_Prob)) +
  geom_line(color="darkgreen") +
  ggtitle("Evolution of Chosen Probability of the Observed Targets") +
  xlab("Event Index") +
  ylab("Probability") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "chosen_probability_evolution.png"), plot_prob, width=8, height=6)

# Visualize the final weight distribution (influence) of nodes
plot_weights <- ggplot(final_weights, aes(x=Final_Weight)) +
  geom_histogram(binwidth = 5, fill="orange", color="black") +
  ggtitle("Final Weight (Influence) Distribution of Nodes") +
  xlab("Final Weight") +
  ylab("Frequency") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "final_weight_distribution.png"), plot_weights, width=8, height=6)

###############################################
# Step 6: Interpretation and Insights
###############################################
# - The fitted alpha and beta parameters provide insights into how baseline reinforcement (alpha) and sensitivity
#   to rating values (beta) influence the trust dynamics.
# - The final weight distribution can highlight the emergence of "influential" nodes that accumulate higher weights.
# - Future extensions could include more complex reinforcement functions, time-varying parameters, or hierarchical modeling.

# Save summary report as a PDF (if desired, using rmarkdown or knitr for a publication-ready report).

print("Analysis complete. Check the output directory for plots and tables.")
