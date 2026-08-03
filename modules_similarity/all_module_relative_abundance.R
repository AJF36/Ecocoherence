###############################################################################################
#This cript creates a tsv file with the information of the relative abundance of all          #
# the modules of the different substrates in order to plot it along the script heatmap_modules#
###############################################################################################



rm(list = ls())
library(readr)
library(tidyr)
library(dplyr)
library(this.path)

work_dir<-this.dir()

substrates <- c("Agarose","AgaroseAlginate","AgaroseChitosan","AgaroseCarrageenan","Carrageenan","Chitin","Alginate")

# Generate the first vector
vector1 <- as.vector(outer(seq(0, 84, by = 12), c("_A", "_B", "_C"), paste0))

# Generate the second vector
vector2 <- as.vector(outer(seq(108, 156, by = 24), c("_A", "_B", "_C"), paste0))
vector3 <- c("204_A","204_B","204_C")

column_names <- c(vector1,vector2,vector3)
numeric_part <- as.numeric(gsub("_[A-Z]", "", column_names))  # Remove letters to isolate numbers

column_names <- column_names[order(numeric_part)]
column_names <- c("module",column_names)

numeric_groups <- gsub("_[A-Z]$", "", column_names[-1])
grouped_columns <- unique(numeric_groups)

results_df <-data.frame(matrix(ncol = length(c(1,grouped_columns)) ))
colnames(results_df) <- c("module",grouped_columns)


for (substrate in substrates){
  
  # substrate <- "Agarose"
  ###FOrmatting of the table
  # substrate <- "AgaroseAlginate"
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))
  setwd(work_dir)
  abundance_modules <- read_tsv(paste("otu_table_",substrate,"_byModulesSize4.tsv",sep = ""))
  # colnames(abundance_modules) <- column_names
  colnames(abundance_modules)[1] <- "module"
  # abundance_modules <- separate(abundance_modules,"204_B", into = c("204_B","204_C"),sep = "\t")
  abundance_modules[-1] <- lapply(abundance_modules[-1], as.numeric)
  

  
  
  
  # numeric_groups <- gsub("_[A-Z]$", "", colnames(abundance_modules)[-1])
  numeric_groups <- gsub(".*_([0-9]+)_[A-Z]$", "\\1", colnames(abundance_modules)[-1])
  # numero <- gsub(".*_([0-9]+)$", "\\1", string)
  
  numeric_groups
  grouped_columns <- unique(numeric_groups)
  
  colnames(abundance_modules) <- c("module",gsub(".*_([0-9]+)_([A-Z])$", "\\1_\\2", colnames(abundance_modules)[-1]))
  # Step 1: Gather all time-replicate columns into a long format
  long_format <- abundance_modules %>%
    pivot_longer(
      cols = -module,  # Exclude "module" column
      names_to = c("Time", "Replicate"),  # Split column names into Time and Replicate
      names_sep = "_",  # Separator is "_"
      values_to = "Abundance"  # Name of the new column for values
    )
  
  # Step 2: Calculate the median of each group (Time) for each module
  median_df <- long_format %>%
    group_by(module, Time) %>%
    summarise(MedianAbundance = median(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Step 3: Pivot back to wide format if needed
  wide_format <- median_df %>%
    pivot_wider(
      names_from = Time,  # Create columns based on Time
      values_from = MedianAbundance  # Fill these columns with median values
    )
  
  
  # Assuming `abundance_modules` has the format: first column is `module`, rest are abundances.
  relative_abundance <- wide_format %>%
    mutate(across(-module, ~ .x / sum(.x, na.rm = TRUE)))  # Calculate relative abundance
  colSums(relative_abundance[,-1])
  # View the result
  relative_abundance$module <- gsub("mod_","",relative_abundance$module)
  relative_abundance$module <- paste(relative_abundance$module,substrate,sep= "_")
  
  results_df <- rbind(relative_abundance,results_df)
  
  
}
colSums(drop_na(results_df[,-1]))


###Now for carrageenan
setwd(work_dir)
carrageenan_data <- read_tsv("otu_table_Carrageenan_byModulesSize4.tsv")

# Generate the first vector
vector1 <- as.vector(outer(seq(0, 84, by = 12), c("_A", "_B", "_C"), paste0))

# Generate the second vector
vector2 <- as.vector(outer(seq(108, 156, by = 24), c("_A", "_B", "_C"), paste0))
vector3 <- c("204_A","204_B")

column_names <- c(vector1,vector2,vector3)
numeric_part <- as.numeric(gsub("_[A-Z]", "", column_names))  # Remove letters to isolate numbers
column_names <- column_names[order(numeric_part)]
column_names <- c("module",column_names)
column_names <- column_names[column_names != "24_C"]


colnames(carrageenan_data) <- column_names
# carrageenan_data <- separate(carrageenan_data,"204_B", into = c("204_B","204_C"),sep = "\t")
carrageenan_data[-1] <- lapply(carrageenan_data[-1], as.numeric)





numeric_groups <- gsub("_[A-Z]$", "", colnames(carrageenan_data)[-1])

grouped_columns <- unique(numeric_groups)

# Step 1: Gather all time-replicate columns into a long format
long_format_carrageenan <- carrageenan_data %>%
  pivot_longer(
    cols = -module,  # Exclude "module" column
    names_to = c("Time", "Replicate"),  # Split column names into Time and Replicate
    names_sep = "_",  # Separator is "_"
    values_to = "Abundance"  # Name of the new column for values
  )

# Step 2: Calculate the median of each group (Time) for each module
median_df_carrageenan <- long_format_carrageenan %>%
  group_by(module, Time) %>%
  summarise(MedianAbundance = median(Abundance, na.rm = TRUE), .groups = "drop")

# Step 3: Pivot back to wide format if needed
wide_format_carrageenan <- median_df_carrageenan %>%
  pivot_wider(
    names_from = Time,  # Create columns based on Time
    values_from = MedianAbundance  # Fill these columns with median values
  )


# Assuming `abundance_modules` has the format: first column is `module`, rest are abundances.
relative_abundance_carrageenan <- wide_format_carrageenan %>%
  mutate(across(-module, ~ .x / sum(.x, na.rm = TRUE)))  # Calculate relative abundance
colSums(relative_abundance[,-1])
# View the result
relative_abundance_carrageenan$module <- gsub("mod_","",relative_abundance_carrageenan$module)
relative_abundance_carrageenan$module <- paste(relative_abundance_carrageenan$module,"Carrageenan",sep= "_")

results_df <- rbind(relative_abundance_carrageenan,results_df)

results_df <- drop_na(results_df)
setwd(work_dir)
write_tsv(results_df,"all_modules_relativeabundance.tsv")

