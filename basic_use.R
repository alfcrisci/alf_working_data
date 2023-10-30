bas_summary <- function(data,na.rm=T) {
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

ext_summary <- function(data) {
  # Calculate summary statistics
  mean_val <- mean(data, na.rm = TRUE)
  median_val <- median(data, na.rm = TRUE)
  sd_val <- sd(data, na.rm = TRUE)
  se_val <- sd(data)/length(data)
  min_val <- min(data, na.rm = TRUE)
  max_val <- max(data, na.rm = TRUE)
  var_val <- var(data, na.rm = TRUE)
  range_val <- max_val - min_val
  missing_count <- sum(is.na(data))
  quartiles <- quantile(data)

  # Create a summary data frame
  summary_df <- data.frame(
    Mean = mean_val,
    Median = median_val,
    SD = sd_val,
    SE = se_val,
    Variance = var_val,
    Min = min_val,
    Max = max_val,
    Range = range_val,
    Missing_Values = missing_count,
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

ex_summary_df<- function(df) {
  # Filter only numeric columns
  numeric_cols <- sapply(df, is.numeric)
  numeric_df <- df[, numeric_cols]

  
  # Calculate summary statistics
  summary_stats <- do.call("rbind",apply(numeric_df,2,ext_summary))

  return(summary_stats)
}
