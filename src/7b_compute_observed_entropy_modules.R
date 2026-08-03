########################################n#####################################
# This script will make the Z score of the observed vaulues for the families #
##############################################################################

# Clear workspace
rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(this.path)
library(purrr)
# Define substrates
script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "7b_compute_observed_entropy_modules"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")
### Create the lists that store the results
# results_S <- list()
results_X <- list()

# Main loop
for (substrate in substrates) {
  # Set working directory
  # path <- paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = "")
  path <- file.path(script_dir, "..","results","6a_build_family_module_table")
  path <- normalizePath(path)
  setwd(path)
  
  # Load data
  data <- read_tsv(paste("group_analysis_", substrate, ".tsv", sep = ""))
  df <- data
  rownames(df) <- data$Family
  length(data$Family)
  
  # Select numeric columns
  numeric_columns <- df %>% select(matches("^\\d"))
  
  # Normalize numeric columns by n_members
  df <- df %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ .x / sum(.x, na.rm = TRUE),
    .names = "{.col}"
  ))
  
  # Re-select normalized numeric columns
  numeric_columns <- df %>% select(matches("^\\d")) ### Select the modules from the colums
  # Calculate Shannon entropy (S) and effective number of states (X)
  S <- -colSums(numeric_columns * log(numeric_columns), na.rm = TRUE) 
  X <- exp(S) 
  #Store the results
  # results_S[[substrate]] <- S
  results_X[[substrate]] <- X

  # Save results
  # setwd(work_dir)
  # write_tsv(df,paste(substrate,"observed_S_X_family.tsv",sep="_"))
  # write_tsv(df, file.path(path, "observed_S_X_family.tsv"))
}

str(results_X)

### Create the df for each modules
final_df <- imap_dfr(results_X, function(mod_list, substrate) {
  tibble(
    substrate = substrate,
    module = names(mod_list),
    value = unname(mod_list)
  )
})

setwd(results_dir)

write_tsv(final_df,"entropy_by_modules.tsv")