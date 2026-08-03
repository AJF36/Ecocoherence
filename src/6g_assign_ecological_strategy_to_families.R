##############################################################################################################
# This script assigns an ecological strategy to each family, similarly to how it’s done with functionink modules
##############################################################################################################

# --- Environment setup ---
rm(list = ls())
library(tidyverse)
library(phyloseq)
library(this.path)
library(gplots)
library(BiotypeR)
library(pheatmap)
library(dendextend)
library(stringr)


script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "6g_assign_ecological_strategy_to_families"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "6g_assign_ecological_strategy_to_families"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
setwd(script_dir)
# --- Load data ---
fileOTU <- file.path(script_dir, "..", "data", "marine_particles_source_data", "count_table.ESV.4R.csv")
fileTaxonomy <- file.path(script_dir, "..", "data", "marine_particles_source_data", "sequence_table.ESV.fasta_RDPclassified.txt")

otu.in <- read.csv(normalizePath(fileOTU, mustWork = FALSE), row.names = 1)
taxa.in <- read.table(normalizePath(fileTaxonomy, mustWork = FALSE), sep = ";")

colnames(taxa.in) <- c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum",
                       "Class","sig_Class","Order","sig_Order","Family","sig_Family",
                       "Genus","sig_Genus")

taxonomy <- taxa.in %>%
  select(taxa_id, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  column_to_rownames("taxa_id") %>%
  select(-Kingdom, -Genus)

# --- Select substrate ---
#
substrates <- c("Agarose","Alginate","AgaroseAlginate","Carrageenan","AgaroseCarrageenan","Chitin","AgaroseChitosan")


### df to store the information about the strategies of the families
strategies_df <- data.frame("substrate" = character() , "family" = character(), "strategy" = character() ) ## df to store the results


for (substrate in substrates){
  print(substrate)
  setwd(script_dir)
  # --- Filter abundance table for substrate ---
otu.in.substrate <- otu.in %>%
  select(matches(paste0("beads.*_", substrate, "_")))

# --- Filter only ESVs present in functionink modules ---
dir_path <- file.path(script_dir,"..","results","5_filter_modules_by_size_and_plot_abundance",substrate)

file_fun <- list.files(
  path = dir_path,
  pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
  full.names = TRUE
)
functionink <- read.table(file_fun, sep = "\t", header = TRUE)
functionink$ESV <- rownames(functionink)
colnames(functionink) <- c("mod", "ESV")

otu.in.substrate.f <- otu.in.substrate[rownames(otu.in.substrate) %in% functionink$ESV, ]

# --- Compute average per timepoint (mean of replicates) ---

times <- unique(str_extract(names(otu.in.substrate.f), "(?<=Beads_)[A-Za-z]+_([0-9]+)") %>%
                  str_extract("[0-9]+"))

df_mean_by_time <- sapply(times, function(t) {
  cols <- grep(paste0("Beads_", substrate, "_", t, "_[A-Z]$"), names(otu.in.substrate.f), value = TRUE)
  rowMeans(otu.in.substrate.f[, cols], na.rm = TRUE)
}) %>% as.data.frame()

colnames(df_mean_by_time) <- paste0(substrate, "_", times)

# --- Normalize abundances per sample (columns sum to 1) ---
otu.rel <- sweep(df_mean_by_time, 2, colSums(df_mean_by_time), "/")

# --- Add family info ---
otu.rel$family <- taxonomy$Family[match(rownames(otu.rel), rownames(taxonomy))]

# --- Aggregate by family (mean abundance of all ESVs in each family) ---
substrate_df_filtered_grouped <- otu.rel %>%
  group_by(family) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()
colSums(substrate_df_filtered_grouped[,-1])
# --- Remove families with zero abundance ---
# substrate_df_filtered_grouped <- substrate_df_filtered_grouped %>%
#   filter(rowSums(across(where(is.numeric))) > 0)

# --- Prepare matrix for distance calculations ---
families <- substrate_df_filtered_grouped$family
substrate_df_filtered_grouped <- select(substrate_df_filtered_grouped, -family)

# --- Normalize rows (each family sums to 1, required for JSD) ---
substrate_df_filtered_grouped <- sweep(substrate_df_filtered_grouped, 1,
                                       rowSums(substrate_df_filtered_grouped), "/")

rownames(substrate_df_filtered_grouped) <- families
# --- Ensure columns are ordered by time ---
substrate_df_filtered_grouped <- substrate_df_filtered_grouped[, order(as.numeric(gsub("\\D", "", colnames(substrate_df_filtered_grouped))))]

# --- Create ideal ecological strategy vectors ---
Nsamp <- ncol(substrate_df_filtered_grouped)
profiles <- matrix(0, nrow = 5, ncol = Nsamp)

profiles[1, 1:1] <- 1                # attachment (early only)
profiles[2, 2:5] <- 1/4              # early stage
profiles[3, 6:7] <- 1/2              # mid
profiles[4, 8:10] <- 1/3             # late
profiles[5, ] <- 1/Nsamp             # all

profiles_df <- as.data.frame(profiles)
colnames(profiles_df) <- colnames(substrate_df_filtered_grouped)
profiles_df$family <- c("attachment","selection","transition","facilitation","generalist")

# --- Combine families + ideal vectors ---
combined_df <- rbind(
  profiles_df %>% select(-family),
  substrate_df_filtered_grouped
)

rownames(combined_df) <- c(profiles_df$family,
                           rownames(substrate_df_filtered_grouped))

# --- Compute Jensen-Shannon distance and cluster ---
dist_matrix <- t(as.matrix(combined_df))
dist_JSD <- dist.JSD(dist_matrix)

  
setwd(figures_dir)

  
# heatmap visualization
pdf(paste("heatmap_ecological_Strategy_families_",substrate,".pdf", sep = ""), 15,15)
heatmap.2(as.matrix(as.dist(dist_JSD)), trace = "none", cexCol = 0.1,cexRow = 0.5, dendrogram = "column",
density.info = "none", keysize = 0.9, hclustfun = function(x) hclust(x, method = "average"))
dev.off()

 dendogram_substrate <-   as.dendrogram(hclust(dist_JSD, method = "average"))

clusters <- cutree(dendogram_substrate, k = 1)

  

n <- 1
while (length(unique(clusters[c("attachment", "selection", "transition", "facilitation", "generalist")])) != 5 ) {
  n <- n + 1
 clusters <- cutree(dendogram_substrate,n)  
}
clusters
### Code for the dendogram
pdf(paste("dendogram_strategy_families_",substrate,".pdf",sep = ""),15,22)
dendogram_substrate %>% 
  dendextend::set("branches_k_color", k = n) |> 
  plot()
rect.dendrogram(dendogram_substrate,k = n , border = 8) 
dev.off()

###Store the info
  n ### n is the number of clusters
  
  for (i in 1:n){
    
  test<-names(clusters[clusters == i])  ### Method to select the different clusters
  test_strategy <- test[grepl("\\b(attachment|selection|transition|facilitation|generalist)\\b",test)]

    if (i <= 5 ){
      cluster_strategies_df <- data.frame("substrate" = substrate, "family" = test, "strategy" = test_strategy)
      strategies_df <- rbind(strategies_df,cluster_strategies_df)

    } else if (i > 5){
      cluster_no_strategy_df <- data.frame("substrate" = substrate, "family" = test , "strategy" = "no_strategy")
      strategies_df <- rbind(strategies_df, cluster_no_strategy_df)
    }
  }
}

strategies_df<- strategies_df %>%
    filter(!(family %in% c("attachment", "selection", "transition", "facilitation", "generalist")))
setwd(results_dir)
write_tsv(strategies_df,"ecological_strategies_families.tsv")
print("Finished!")
