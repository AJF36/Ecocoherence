###################################
#     Combine Results
###################################
# Script to combine the results from all 25 processes
rm(list= ls())
library(readr)
library(dplyr)
library(this.path)

this.dir = this.dir()
setwd(this.dir)

# Initialize empty dataframe for combined results
combined_df <- data.frame()

# Combine all result files
total_jobs <- 25
for (job_id in 1:total_jobs) {
  filename <- paste0("intersection_union_df_part_", job_id, ".tsv")
  
  # Check if file exists
  if (file.exists(filename)) {
    part_df <- read_tsv(filename,col_select = c("IDs","Familys","n_intersection","n_union"))
    combined_df <- rbind(combined_df, part_df)
    
    print(paste("Added results from job", job_id))
  } else {
    warning(paste("File not found:", filename))
  }
  
}

# Save the combined results
write_tsv(combined_df, "intersection_union_df_combined.tsv")
print(paste("Combined", nrow(combined_df), "rows of results from all jobs."))