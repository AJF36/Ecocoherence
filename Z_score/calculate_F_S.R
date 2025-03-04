########################################n#####################################
# This script will make the Z score of the observed vaulues for the families #
##############################################################################

# Clear workspace
rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(this.path)
# Define substrates
work_dir <- this.dir()
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")

# Main loop
for (substrate in substrates) {
  # Set working directory
  # path <- paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = "")
  path <- file.path(work_dir, "..","modules_table")
  path <- normalizePath(path)
  setwd(path)
  
  # Load data
  data <- read_tsv(paste("group_analysis_", substrate, ".tsv", sep = ""))
  df <- data
  
  # Select numeric columns
  numeric_columns <- df %>% select(matches("^\\d"))
  
  # Normalize numeric columns by n_members
  df <- df %>%
    mutate(across(
      .cols = matches("^\\d"),
      .fns = ~ ifelse(n_members != 0, .x / n_members, NA),
      .names = "{.col}" # Keeps original column names
    ))
  
  # Re-select normalized numeric columns
  numeric_columns <- df %>% select(matches("^\\d"))
  
  # Calculate Shannon entropy (S) and effective number of states (X)
  df$S <- -rowSums(numeric_columns * log(numeric_columns), na.rm = TRUE)
  df$X <- exp(df$S)
  # Save results
  setwd(work_dir)
  write_tsv(df,paste(substrate,"observed_S_X_family.tsv",sep="_"))
  # write_tsv(df, file.path(path, "observed_S_X_family.tsv"))
}
