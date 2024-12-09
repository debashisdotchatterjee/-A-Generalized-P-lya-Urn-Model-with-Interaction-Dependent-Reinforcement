# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create output directory
output_dir <- "Bitcoin_Trust_Analysis"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Load dataset (update the file path accordingly)
soc.sign.bitcoinalpha <- read.csv("C:/Users/User/Desktop/Debashis 2024/Bitcoin Python/soc-sign-bitcoinalpha.csv", 
                                  header = FALSE, col.names = c("SOURCE", "TARGET", "RATING", "TIME"))
data <- soc.sign.bitcoinalpha

# Convert TIME to human-readable date using the Unix epoch
data$DATE <- as.POSIXct(data$TIME, origin = "1970-01-01", tz = "UTC")

# Save the cleaned dataset with readable dates
write.csv(data, file = file.path(output_dir, "cleaned_dataset.csv"), row.names = FALSE)

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
    Time_Span = paste(min(DATE), "to", max(DATE))
  )

# Save Basic Statistics
write.csv(basic_stats, file = file.path(output_dir, "basic_statistics.csv"), row.names = FALSE)

# Print Basic Statistics
print("Basic Statistics:")
print(basic_stats)

# Plot 1: Distribution of Ratings
plot1 <- ggplot(data, aes(x = RATING)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  ggtitle("Distribution of Ratings") +
  xlab("Rating") +
  ylab("Frequency") +
  theme_minimal()
ggsave(file.path(output_dir, "distribution_of_ratings.png"), plot1)

# Plot 2: Time Series of Ratings
plot2 <- ggplot(data, aes(x = DATE, y = RATING)) +
  geom_jitter(alpha = 0.5, color = "blue") +
  ggtitle("Ratings Over Time") +
  xlab("Date") +
  ylab("Rating") +
  theme_minimal()
ggsave(file.path(output_dir, "ratings_over_time.png"), plot2)

# Plot 3: Heatmap of Source-Target Interactions
interaction_matrix <- table(data$SOURCE, data$TARGET)
interaction_matrix_df <- as.data.frame(as.table(interaction_matrix))

plot3 <- ggplot(interaction_matrix_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  xlab("Source Nodes") +
  ylab("Target Nodes") +
  ggtitle("Heatmap of Source-Target Interactions")
ggsave(file.path(output_dir, "source_target_heatmap.png"), plot3)

# Plot 4: Top 10 Sources by Number of Ratings
top_sources <- data %>%
  group_by(SOURCE) %>%
  summarize(Ratings_Count = n()) %>%
  arrange(desc(Ratings_Count)) %>%
  head(10)

plot4 <- ggplot(top_sources, aes(x = reorder(SOURCE, Ratings_Count), y = Ratings_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  ggtitle("Top 10 Sources by Number of Ratings") +
  xlab("Source Nodes") +
  ylab("Number of Ratings") +
  theme_minimal()
ggsave(file.path(output_dir, "top_sources.png"), plot4)

# Plot 5: Top 10 Targets by Number of Ratings
top_targets <- data %>%
  group_by(TARGET) %>%
  summarize(Ratings_Count = n()) %>%
  arrange(desc(Ratings_Count)) %>%
  head(10)

plot5 <- ggplot(top_targets, aes(x = reorder(TARGET, Ratings_Count), y = Ratings_Count)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  ggtitle("Top 10 Targets by Number of Ratings") +
  xlab("Target Nodes") +
  ylab("Number of Ratings") +
  theme_minimal()
ggsave(file.path(output_dir, "top_targets.png"), plot5)

# Save tables for Top Sources and Top Targets
write.csv(top_sources, file = file.path(output_dir, "top_sources.csv"), row.names = FALSE)
write.csv(top_targets, file = file.path(output_dir, "top_targets.csv"), row.names = FALSE)

# Final Output Confirmation
print(paste("All plots and tables saved in the directory:", output_dir))

%%%%%%%%%%%%%%%%%%%%%

# Load necessary libraries
library(dplyr)
library(knitr)
library(kableExtra)

# Print and save Basic Statistics table
basic_stats_table <- basic_stats %>%
  as.data.frame() %>%
  rename(
    `Total Ratings` = Total_Ratings,
    `Unique Sources` = Unique_Sources,
    `Unique Targets` = Unique_Targets,
    `Minimum Rating` = Min_Rating,
    `Maximum Rating` = Max_Rating,
    `Average Rating` = Avg_Rating,
    `Median Rating` = Median_Rating,
    `Standard Deviation` = Std_Rating,
    `Time Span` = Time_Span
  )

# Save the table as a LaTeX file
sink(file = file.path(output_dir, "basic_statistics_table.tex"))
kable(basic_stats_table, format = "latex", booktabs = TRUE, caption = "Summary Statistics of the Bitcoin Alpha Trust Network Dataset") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
sink()

# Print the table in the console for quick reference
print("Basic Statistics Table:")
print(basic_stats_table)

# Print and save Top Sources table
top_sources_table <- top_sources %>%
  rename(`Source Node` = SOURCE, `Number of Ratings` = Ratings_Count)

sink(file = file.path(output_dir, "top_sources_table.tex"))
kable(top_sources_table, format = "latex", booktabs = TRUE, caption = "Top 10 Source Nodes by Number of Ratings") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
sink()

# Print the table in the console for quick reference
print("Top Sources Table:")
print(top_sources_table)

# Print and save Top Targets table
top_targets_table <- top_targets %>%
  rename(`Target Node` = TARGET, `Number of Ratings` = Ratings_Count)

sink(file = file.path(output_dir, "top_targets_table.tex"))
kable(top_targets_table, format = "latex", booktabs = TRUE, caption = "Top 10 Target Nodes by Number of Ratings") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
sink()

# Print the table in the console for quick reference
print("Top Targets Table:")
print(top_targets_table)

