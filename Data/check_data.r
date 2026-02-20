rm(list = ls())
library(lattice)
new_files <- list.files("new")
old_files <- list.files("old")
new_data <- lapply(paste0('new/', new_files), read.csv)
new_data <- lapply(new_data, function(df){
  df$time <- df$doy + df$hour/24
  return(df)
})
old_data <- lapply(paste0('old/', old_files), read.csv)

i <- 3
# Extract column names
old_cols <- colnames(old_data[[i]])
new_cols <- colnames(new_data[[i]])

# Find columns in old_data but not in new_data
missing_in_new <- setdiff(old_cols, new_cols)

# Find columns in new_data but not in old_data
missing_in_old <- setdiff(new_cols, old_cols)

# Find columns in both (intersection)
common_cols <- intersect(old_cols, new_cols)

# Print results
cat("Columns in old_data but NOT in new_data:\n")
print(missing_in_new)

cat("\nColumns in new_data but NOT in old_data:\n")
print(missing_in_old)

# Initialize a list to store differences
differences <- list()

# Compare each common column
for (col in common_cols) {
  old_val <- old_data[[i]][[col]]
  new_val <- new_data[[i]][[col]]
  
  # Count differences
  num_diff <- sum(old_val != new_val, na.rm = TRUE)
  
  # Only print if there are differences
  if (num_diff > 1e-6) {
    cat("Column:", col, "\n")
    cat("  Number of differing values:", num_diff, 
        "out of", length(old_val), 
        sprintf("(%.2f%%)\n", 100 * num_diff / length(old_val)))
    
    # Show first few differences
    diff_idx <- which(old_val != new_val)
    show_n <- min(5, length(diff_idx))
    
    cat("  First", show_n, "differences (row index: old -> new):\n")
    for (j in 1:show_n) {
      idx <- diff_idx[j]
      cat("    Row", idx, ":", old_val[idx], "->", new_val[idx], "    ",
          100*(new_val[idx]-old_val[idx])/old_val[idx] , '% Change', "\n")
    }
    cat("\n")
  }
}


selected_cols <- c("day_tracker", "dusk_kick", "light", "solar")

# Get time variable
time_var <- new_data[[i]]$time

# Set up plotting area (2x2 grid)
par(mfrow = c(2, 2))
par(mar = c(4, 4, 2, 1))

# Calculate and plot percentage difference for each selected column
for (col in selected_cols) {
  old_val <- old_data[[i]][[col]]
  new_val <- new_data[[i]][[col]]
  
  # Calculate percentage difference: (new - old) / old * 100
  pct_diff <- ((new_val - old_val) / old_val) * 100
  
  plot(time_var, pct_diff, 
       type = "l",
       main = col,
       xlab = "Time",
       ylab = "% Difference",
       col = "blue",
       lwd = 2)
  grid()
}

# Reset plotting parameters
par(mfrow = c(1, 1))
