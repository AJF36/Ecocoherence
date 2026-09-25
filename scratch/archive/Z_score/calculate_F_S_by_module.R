########################################n#####################################
# This script will make the Z score of the observed vaulues for the modules  #
##############################################################################

# Clear workspace
rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(this.path)
library(purrr)
# Define substrates
work_dir <- this.dir()
# substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")
substrate <- "Alginate"

###Load the tax table
# fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
# tax.pseq = tax_table(as.matrix(taxonomy))

### functionink table 
subs_path <- file.path(work_dir,"..","functionink",substrate)
  subs_path <- normalizePath(subs_path)
  setwd(subs_path)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  ### First we need the complete functionink file
  dir_path <- "functionink_tmp/"
  file <- list.files(
    path = dir_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv$"),
    full.names = TRUE
  )

functionink <- as.data.frame(read_tsv(file,skip = 9,col_names = F))
colnames(functionink) <- c("ESV","mod")

### Add the family info to the functionink df
functionink$family <- taxonomy$Family[match(functionink$ESV,rownames(taxonomy))]

### Store the total number of families
total_f <- length(unique(functionink$family))

list_n_families_for_each_module <- list()
### Now store the number of families in each of the modules
x <- 1
for (module in unique(functionink$mod)) {
  
  functionink_mod <- filter(functionink, mod == module)
  families_module <- length(unique(functionink_mod$family))
  list_n_families_for_each_module[[x]] <- families_module
  names(list_n_families_for_each_module)[x] <- module
  x <- x + 1
}
total_f

p_modules <- purrr::map(list_n_families_for_each_module, \(x) (-x/total_f * log(x/total_f)))



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
  # write_tsv(df,paste(substrate,"observed_S_X_family.tsv",sep="_"))
  # write_tsv(df, file.path(path, "observed_S_X_family.tsv"))
}
