# Load necessary libraries (if uncommented)
#library(tidyverse)
#library(magrittr)

# Create new columns of locus tags, including those present in more than one copy.
# This will replace the original column in the dataset and allow insertion of 
# blanks when a locus tag does not appear more than once for certain samples

# Reduced dataset for testing
#df_1 <- read.delim("~/Desktop/tagsReduced.txt", header=FALSE)

# Load the full dataset
df_1 <- read.delim("~/Desktop/combinedWithOutHeaders.txt", header=FALSE)

# Initialize an empty list to store results of duplicated values
results_list <- list()

# Loop through each even column (containing locus tags) in the data frame
for (i in seq(2, ncol(df_1), by = 2)) {
  
  # Extract values from the current column and remove empty strings
  values <- df_1[, i]
  values <- values[values != " "]
  
  # Identify duplicated values in the column
  duplicated_values  <- values[duplicated(values) | duplicated(values, fromLast = TRUE)]
  
  # If there are duplicated values, create a data frame and add it to the results list
  if (length(duplicated_values ) > 0) {
    results <- data.frame(Value = duplicated_values , ColumnID = i, stringsAsFactors = FALSE)
    results_list[[length(results_list) + 1]] <- results
  }
}

# Combine all results into a single dataframe and remove empty values
final_results <- do.call(rbind, results_list)
final_results <- final_results[-which(final_results$Value == ""), ]

# Extract values from the first column (assumed to be the reference column)
col1_values <- df_1$V1
col1_values <- as.data.frame(col1_values)

# Count the number of occurrences of each locus tag per sample (column)
results_by_cols <- final_results %>%
  group_by(ColumnID, Value) %>%
  summarize(count = n(), .groups = 'drop') 

# Determine the maximum count of each locus tag across all columns
max_count_per_tag <- results_by_cols %>% 
  group_by(Value) %>%
  summarize(max_count = max(count), .groups = 'drop')

# Calculate the number of times each locus tag should be repeated
repeat_values <- max_count_per_tag %>%
  mutate(repeat_count = max_count - 1) # Subtract 1 from each max_count

# Sanity check
# For example, if 'locus_tag=lpt_01135' appears 7 times in column 74,
# then 'repeat_count' should be 6, meaning 6 additional copies should be added.
max(results_by_cols$count) # Get the maximum count
results_by_cols %>% filter(count == 7) # Find the column with the maximum count
repeat_values %>% filter(Value == 'locus_tag=lpt_01135')  # Check repeat count for specific tag

# Convert repeat counts to numeric and remove any NA values
repeat_values$repeat_count <- as.numeric(as.character(repeat_values$repeat_count))
repeat_values <- na.omit(repeat_values)

# Create a vector of repeated values based on the repeat counts
repeated_vector <- unlist(mapply(function(tag, times) rep(tag, times), repeat_values$Value, repeat_values$repeat_count))
repeated_df <- as.data.frame(repeated_vector)

# Combine the original locus tags with the repeated locus tags and sort
all_with_dups <- c(col1_values$col1_values, repeated_df$repeated_vector)
all_with_dups <- as.data.frame(all_with_dups)
all_sorted <- all_with_dups[order(all_with_dups$all_with_dups), ]
all_sorted <- as.data.frame(all_sorted)

# Add row names as an 'id' column to the sorted list of tags
all_sorted = cbind("id"=rownames(all_sorted),all_sorted)

# Add row names as an 'id' column to the original dataset
original_with_id  <- cbind("id"=rownames(df_1), df_1)

# Merge the datasets on the 'id' column, keeping all rows from both datasets
merged_data <- merge(all_sorted, original_with_id , all=T)

# Remove the 'id' and original column 'V1' from the merged dataset
merged_data <- subset(merged_data, select = -c(id, V1))

# Rename the column 'all_sorted' to 'V1' to standardize column names
names(merged_data)[names(merged_data) == "all_sorted"] <- "V1"

# Sort the data frame by the 'V1' column
sorted_data <- merged_data[order(merged_data$V1), ]

# Create a copy of the sorted data frame for updates
test_out_updated <- sorted_data

# Define the number of column pairs in the dataset
# Each pair consists of two columns, so calculate the number of pairs
num_pairs <- (ncol(sorted_data) - 1) / 2

# Iterate over each column pair
for (i in seq_len(num_pairs)) {
  
  # Generate column names for the current pair
  col2_name <- paste0("V", 2 * i)
  col3_name <- paste0("V", 2 * i + 1)

  # Initialize new columns with empty spaces
  new_col2 <- rep(" ", nrow(sorted_data))
  new_col3 <- rep(" ", nrow(sorted_data))
  
  # Initialize an index to track the position in the current pair of columns
  index <- 1
  
  # Iterate over each row of the data frame
  for (j in seq_len(nrow(sorted_data))) {
    
    # Check if the current index is within the bounds of the column length
    if (index <= length(sorted_data[[col2_name]])) {
    
        # Check if the values in 'V1' and the current pair of columns match
      if (!is.na(sorted_data$V1[j]) && !is.na(sorted_data[[col2_name]][index]) &&
          sorted_data$V1[j] == sorted_data[[col2_name]][index]) {
        
        # If values match, update the new columns with the corresponding values
        new_col2[j] <- sorted_data[[col2_name]][index]
        new_col3[j] <- sorted_data[[col3_name]][index]
        
        # Move to the next index in the column
        index <- index + 1
      } else {
        # If values do not match or are NA, leave the new columns as empty spaces
        new_col2[j] <- " "
        new_col3[j] <- " "
      }
    } else {
      # If index exceeds the column length, leave the new columns as empty spaces
      new_col2[j] <- " "
      new_col3[j] <- " "
    }
  }
  
  # Add the new columns to the updated data frame
  test_out_updated[[col2_name]] <- new_col2
  test_out_updated[[col3_name]] <- new_col3
}

# At this point, 'test_out_updated' contains updated columns based on the 
# pair-wise comparisons. This is the file to save 
output_file_path <- "~/Desktop/tagsUpdated.csv"
write.csv(test_out_updated, file = output_file_path, row.names = FALSE)
