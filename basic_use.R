basic_summary <- function(data,na.rm=T) {
  # Calculate summary statistics
  mean_val <- mean(data)
  median_val <- median(data)
  sd_val <- sd(data)
  se_val <- sd(data)/length(data)
  min_val <- min(data)
  max_val <- max(data)
  quartiles <- quantile(data)

  # Create a summary data frame
  summary_df <- data.frame(
    Mean = mean_val,
    Median = median_val,
    SD = sd_val,
    SE = se_val,
    Min = min_val,
    Max = max_val,
    Q1 = quartiles[2],
    Q3 = quartiles[4],
    N=length(data)
  )

  return(summary_df)
}



# Function to generate basic summary statistics for a data frame
basic_summary_df<- function(df) {
  # Filter only numeric columns
  numeric_cols <- sapply(df, is.numeric)
  numeric_df <- df[, numeric_cols]

  # Remove NAs from the numeric data frame
  numeric_df <- na.omit(numeric_df)

  # Calculate summary statistics
  summary_stats <- summary(numeric_df)

  return(summary_stats)
}
