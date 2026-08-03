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
# Initialize results list for this job
intersection_union_list <- list()
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
    
    intersection_union_list[[length(intersection_union_list) + 1]] <- temp_df
  }
}
# Convert list to dataframe at the end
intersection_union_df <- do.call(rbind, intersection_union_list)
# Save results for this job
output_filename <- paste0("intersection_union_df_part_", job_id, ".tsv")
write_tsv(intersection_union_df, output_filename)
print(paste("Job", job_id, "completed successfully. Results saved to", output_filename))