###################################
#     name snippe
###################################
# Author: Adrian Jimenez Fernandez
# Copyright (c)  Adrian Jimenez Fernandez,  `r paste(format(Sys.Date(), "%Y"))`
# Web:  
#
# Date: `r paste(Sys.Date())`
# Script Name: name_of_script.R 
# Script Description:
#
#
# Notes:
#
#
# library(this.path)
# 
# 
# # SET WORKING DIRECTORY -----------------------------
# # this.dir = find_rstudio_root_file()
# this.dir = this.dir()
# # #this.dir=strsplit(rstudioapi::getActiveDocumentContext()\$path, "/src/")[[1]][1] # don't edit, just comment it if problems...
# # dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
# # dirInput=paste(this.dir,"/INPUT/",sep="") # Dir of input data CHANGE INPUT BY THE CORRECT DIRECTORY
# # dirOutput=paste(this.dir,"/OUTPUT/",sep="") # Dir of output data CHANGE OUTPUT BY THE CORRECT DIRECTORY
# # INSTALL PACKAGES & LOAD LIBRARIES -----------------
# cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
# packages <- c("readr", "tidyr", "dplyr") # list of packages to load
# n_packages <- length(packages) # count how many packages are required
# new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
# 
# # install missing packages
# if(length(new.pkg)){
#   install.packages(new.pkg)
# }
# 
# # load all requried libraries
# for(n in 1:n_packages){
#   cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
#   lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
#   eval(parse(text = lib_load)) # evaluate the string to load the library
# }
# 
# setwd(this.dir)
# ###Load the ID of the aligned ESV
# aligned_genomes <- read_tsv("matched_ESV_id_0.97.tsv",col_names = F)
# aligned_genomes <- select(aligned_genomes,c("X9","X10"))
# colnames(aligned_genomes) <- c("ASV","ID")
# aligned_genomes <- filter(aligned_genomes, ID != "*")
# 
# unique_ids <- unique(aligned_genomes$ID)
# aligned_genomes_unique <- data.frame("ID" = unique_ids)
# # aligned_genomes_unique$ID <- gsub("^.*?(GCF_|GCA_)", "\\1", aligned_genomes_unique$ID)
# # aligned_genomes_unique$ID <- gsub("\\..*$", "", aligned_genomes_unique$ID)
# # aligned_genomes_unique$ID <- gsub("_$", "", aligned_genomes_unique$ID)
# 
# ##Load metadata table
# 
# metadata_table <- read_tsv("bac120_metadata_r220.tsv")
# colnames(metadata_table)
# metadata_table_f <- filter(metadata_table,accession %in% aligned_genomes_unique$ID) %>%
#   select(c("accession","checkm2_completeness","checkm2_contamination","gc_percentage","genome_size","coding_density","protein_count"))
# 
# setwd(this.dir())
# # write_tsv(metadata_table_f,"genomes_metadata_filtered.tsv")
# 
# ##Load the kegg db and the table with the id and the kos
# table_ko <- read_tsv("list_id_with_ko.tsv")
# ko_db <- read_tsv("kegg_db.tsv")
# 
# ###Modify the id so they match 
# metadata_table_f$accession <- gsub("^.*?(GCF_|GCA_)", "\\1", metadata_table_f$accession)
# metadata_table_f$accession <- gsub("\\..*$", "", metadata_table_f$accession)
# metadata_table_f$accession <- gsub("_$", "", metadata_table_f$accession)
# colnames(metadata_table_f)[1] = "ID"
# ###join the two tables
# 
# full_table <- inner_join(table_ko, metadata_table_f, by = "ID")
# 
# length(full_table$Family)
# ##Filter familys of interest
# # familys_of_interest <- c("Helicobacteraceae","Nannocystaceae","Campylobacteraceae","Porphyromonadaceae","Colwelliaceae","Alteromonadaceae","Rhodobacteraceae"
# #                          ,"Flavobacteriaceae","Oceanospirillaceae","Syntrophaceae","Anaerolineaceae","Desulfobacteraceae","Psychromonadaceae")
# # coherency_families_of_interest <- c("early_c","late_c","early_c","late_c","incoherent","incoherent","incoherent","incoherent","incoherent","late_c",
# #                                     "late_c","incoherent","early_c")
# # df_for_color <- data.frame("family" = familys_of_interest,"coherency" = coherency_families_of_interest)
# 
# # length(unique(full_table$Family))
# # # full_table_f <- filter(full_table,Family %in% familys_of_interest)  
# # length(unique(full_table$Family)) # 8  of the  350  Familys compound 1/4 of the ID
# 
# 
# ###Generate the pair of ids
# # pairs_ID <- expand_grid("ID1" =full_table_f$ID,"ID2" =full_table_f$ID)
# pairs_ID <- expand_grid("ID1" =full_table$ID,"ID2" =full_table$ID)
# 
# pairs_ID <- filter(pairs_ID,ID1 != ID2)
# 
# ###Vector for filtering the metabolic KOs
# metabolic_KOs_df <- filter(ko_db,KO_family_L1 == "Metabolism")
# metabolic_KOs_vector <- metabolic_KOs_df$KO
# ###Loop for calculating union and intersection
# intersection_union_df <- data.frame("IDs" = character(),"Familys" = character(),"union" = character(),
#                                     "intersection" = character(),"n_union" = numeric(),
#                                     "n_intersection" = numeric())
# 
# 
# ###Loof for making the df
# 
# for (i in 1:length(pairs_ID$ID1)) {
#   print(i)
#   ID_1 <- pairs_ID$ID1[i]
#   ID_2 <- pairs_ID$ID2[i]
# 
#   fam_1 <- full_table$Family[match(ID_1, full_table$ID)]
#   fam_2 <- full_table$Family[match(ID_2, full_table$ID)]
# 
#   if (!is.na(fam_1) && !is.na(fam_2)) {
#     KO_ID_1 <- full_table$KO_list[match(ID_1, full_table$ID)]
#     KO_ID_1 <- strsplit(KO_ID_1, split = ";")[[1]]
# 
#     ##filter the metabolic KOs
#     KO_ID_1_m <- KO_ID_1[KO_ID_1 %in% metabolic_KOs_vector]
#     KO_ID_2 <- full_table$KO_list[match(ID_2, full_table$ID)]
#     KO_ID_2 <- strsplit(KO_ID_2, split = ";")[[1]]
#     KO_ID_2_m <- KO_ID_2[KO_ID_2 %in% metabolic_KOs_vector]
#     ##Calculate the union and the intersection
#     union <- paste(unique(c(KO_ID_1_m,KO_ID_2_m)),collapse = ";")
#     intersection <- paste(intersect(KO_ID_1_m, KO_ID_2_m),collapse = ";")
#     n_union <- length(unique(c(KO_ID_1_m,KO_ID_2_m)))
#     n_intersection <- length(intersect(KO_ID_1_m, KO_ID_2_m))
#     temp_df <- data.frame(
#       "IDs" = paste(min(ID_1, ID_2), max(ID_1, ID_2), sep = ";"),
#       "Familys" = paste(fam_1, fam_2, sep = "_"),
#       "union" = union,
#       "intersection" = intersection,
#       "n_union" = n_union,
#       "n_intersection" = n_intersection
#     )
# 
#     intersection_union_df <- rbind(intersection_union_df, temp_df)
#   }
# }
# 
# 
# setwd(this.dir())
# write_tsv(intersection_union_df,"intersection_union_df.tsv")


###################################
#     Parallelized KO Comparison
###################################
# Author: Based on Adrian Jimenez Fernandez's script
# Script Description: Split KO comparison calculations across 25 processes

library(this.path)
library(readr)
library(tidyr)
library(dplyr)

# Get job ID from command line argument (1-25)
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
total_jobs <- 25

# SET WORKING DIRECTORY -----------------------------
this.dir = this.dir()
setwd(this.dir)

# Load required data
# 1. Load genome alignments
aligned_genomes <- read_tsv("matched_ESV_id_0.97.tsv", col_names = FALSE)
aligned_genomes <- select(aligned_genomes, c("X9", "X10"))
colnames(aligned_genomes) <- c("ASV", "ID")
aligned_genomes <- filter(aligned_genomes, ID != "*")

unique_ids <- unique(aligned_genomes$ID)
aligned_genomes_unique <- data.frame("ID" = unique_ids)

# 2. Load metadata
metadata_table <- read_tsv("bac120_metadata_r220.tsv")
metadata_table_f <- filter(metadata_table, accession %in% aligned_genomes_unique$ID) %>%
  select(c("accession", "checkm2_completeness", "checkm2_contamination", "gc_percentage", 
           "genome_size", "coding_density", "protein_count"))

# 3. Load KEGG database
table_ko <- read_tsv("list_id_with_ko.tsv")
ko_db <- read_tsv("kegg_db.tsv")

# Clean ID format
metadata_table_f$accession <- gsub("^.*?(GCF_|GCA_)", "\\1", metadata_table_f$accession)
metadata_table_f$accession <- gsub("\\..*$", "", metadata_table_f$accession)
metadata_table_f$accession <- gsub("_$", "", metadata_table_f$accession)
colnames(metadata_table_f)[1] = "ID"

# Join tables
full_table <- inner_join(table_ko, metadata_table_f, by = "ID")

# Generate the pairs of IDs
pairs_ID <- expand_grid("ID1" = full_table$ID, "ID2" = full_table$ID)
pairs_ID <- filter(pairs_ID, ID1 != ID2)

# Divide work across processes
total_pairs <- nrow(pairs_ID)
pairs_per_job <- ceiling(total_pairs / total_jobs)
print(pairs_per_job)
# Calculate start and end indices for this job
start_idx <- (job_id - 1) * pairs_per_job + 1
print(start_idx)
end_idx <- min(job_id * pairs_per_job, total_pairs)

# Subset the pairs for this job
job_pairs <- pairs_ID[start_idx:end_idx, ]

# Filter for metabolic KOs
metabolic_KOs_df <- filter(ko_db, KO_family_L1 == "Metabolism")
metabolic_KOs_vector <- metabolic_KOs_df$KO

# Initialize results dataframe for this job
intersection_union_df <- data.frame(
  "IDs" = character(),
  "Familys" = character(),
  "union" = character(),
  "intersection" = character(),
  "n_union" = numeric(),
  "n_intersection" = numeric()
)

# Process pairs assigned to this job
for (i in 1:nrow(job_pairs)) {
  if (i %% 100 == 0) {
    print(paste("Job", job_id, "processing pair", i, "of", nrow(job_pairs)))
  }
  
  ID_1 <- job_pairs$ID1[i]
  ID_2 <- job_pairs$ID2[i]
  
  fam_1 <- full_table$Family[match(ID_1, full_table$ID)]
  fam_2 <- full_table$Family[match(ID_2, full_table$ID)]
  
  if (!is.na(fam_1) && !is.na(fam_2)) {
    KO_ID_1 <- full_table$KO_list[match(ID_1, full_table$ID)]
    KO_ID_1 <- strsplit(KO_ID_1, split = ";")[[1]]
    KO_ID_1_m <- KO_ID_1[KO_ID_1 %in% metabolic_KOs_vector]
    
    KO_ID_2 <- full_table$KO_list[match(ID_2, full_table$ID)]
    KO_ID_2 <- strsplit(KO_ID_2, split = ";")[[1]]
    KO_ID_2_m <- KO_ID_2[KO_ID_2 %in% metabolic_KOs_vector]
    
    union_kos <- paste(unique(c(KO_ID_1_m, KO_ID_2_m)), collapse = ";")
    intersection_kos <- paste(intersect(KO_ID_1_m, KO_ID_2_m), collapse = ";")
    n_union <- length(unique(c(KO_ID_1_m, KO_ID_2_m)))
    n_intersection <- length(intersect(KO_ID_1_m, KO_ID_2_m))
    
    temp_df <- data.frame(
      "IDs" = paste(min(ID_1, ID_2), max(ID_1, ID_2), sep = ";"),
      "Familys" = paste(fam_1, fam_2, sep = "_"),
      "union" = union_kos,
      "intersection" = intersection_kos,
      "n_union" = n_union,
      "n_intersection" = n_intersection
    )
    
    intersection_union_df <- rbind(intersection_union_df, temp_df)
  }
}

# Save results for this job
output_filename <- paste0("intersection_union_df_part_", job_id, ".tsv")
write_tsv(intersection_union_df, output_filename)

print(paste("Job", job_id, "completed successfully. Results saved to", output_filename))
