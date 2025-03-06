########################################################################################################
#This script is used to calculate the significance of the in the interactions between distinct families
########################################################################################################
rm(list = ls())
library(readr)
library(tidyr)
library(dplyr)
library(gplots)
library(this.path)

work_dir <- this.dir()

fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))

substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")



permute_taxonomy <- function(taxonomy) {
  taxonomy_permuted <- taxonomy %>%
    mutate(across(everything(), ~ sample(.)))
  return(taxonomy_permuted)
}


links_df <- data_frame(Sp1 = character(), Sp2 = character())

for (substrate in substrates){
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))
  # network_df <- read_tsv(paste("interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep = ""))
  net_path <- file.path(work_dir,"..","filter_network_by_threshold")
  net_path <- normalizePath(net_path)
  setwd(net_path)
  network_df <- read_tsv(paste("interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep = ""))
  colnames(network_df) <- c("Sp1","Sp2","Cor","Type")
  ###Since we are going to study if some families goes together with a statistical significance, we only take the positive correlations
  network_df_positive <- network_df %>% filter(Cor > 0)
  
  # network_df_positive$Sp1 <- taxonomy$Family[match(network_df_positive$Sp1,rownames(taxonomy))]
  # network_df_positive$Sp2 <- taxonomy$Family[match(network_df_positive$Sp2,rownames(taxonomy))]
  
  network_df_positive <-  select(network_df_positive, c("Sp1","Sp2"))
  
  # links_df <- rbind(links_df,network_df_positive)




###Filter the taxonomy by ASV
ASV_to_filter <- unique(c(network_df_positive$Sp1,network_df_positive$Sp2))
taxonomy$ASV <- row.names(taxonomy)
filtered_taxonomy <- filter(taxonomy, ASV %in% ASV_to_filter)

###Change ASV for farmilies
network_df_positive$Sp1 <- taxonomy$Family[match(network_df_positive$Sp1,rownames(taxonomy))]
network_df_positive$Sp2 <- taxonomy$Family[match(network_df_positive$Sp2,rownames(taxonomy))]

###Now lets make the matrix of the counts
network_df_positive$Pair <- apply(network_df_positive, 1, function(x) paste(sort(x), collapse = "-"))

number_network_df_positive <-as.data.frame(table(network_df_positive$Pair))

families <- unique(c(network_df_positive$Sp1,network_df_positive$Sp2))


##Calculate the randomized counts
results_list <- list()
##Filter the taxnomy list so it onnly have the families that appears in the pairs

# filtered_taxonomy <- filter(taxonomy, Family %in% families) ###FIltering by family, better for ASV

for (iteration in 1:100) {
  iteration_df <-  data_frame(Sp1 = character(), Sp2 = character())
  # for (substrate in substrates){
    permutated_taxonomy <- permute_taxonomy(filtered_taxonomy)
    # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))
    # network_df <- read_tsv(paste("interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep = ""))
    net_path <- file.path(work_dir,"..","filter_network_by_threshold")
    net_path <- normalizePath(net_path)
    setwd(net_path)
    network_df_r <- read_tsv(paste("interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep = ""))
    colnames(network_df_r) <- c("Sp1","Sp2","Cor","Type")
    
    ###Since we are going to study if some families goes together with a statistical significance, we only take the positive correlations
    network_df_positive_r <- network_df_r %>% filter(Cor > 0)
    
    network_df_positive_r$Sp1 <- permutated_taxonomy$Family[match(network_df_positive_r$Sp1,rownames(permutated_taxonomy))]
    network_df_positive_r$Sp2 <- permutated_taxonomy$Family[match(network_df_positive_r$Sp2,rownames(permutated_taxonomy))]
    
    network_df_positive_r <-  select(network_df_positive_r, c("Sp1","Sp2"))
    
    iteration_df <- rbind(iteration_df,network_df_positive_r)
   
  
  iteration_df$Pair <- apply(iteration_df, 1, function(x) paste(sort(x), collapse = "-"))
  
  number_links_iteration_df <-as.data.frame(table(iteration_df$Pair))
  
  results_list[[paste0("perm_", iteration)]] <- number_links_iteration_df
  
  print(iteration)

}



# Combine all randomized iterations into a single dataframe
results_df <- bind_rows(results_list, .id = "Iteration")

# Extract the observed pairs and their frequencies
# pairs <- number_links_df$Var1
pairs <- number_network_df_positive$Var1
# observed_freqs <- number_links_df$Freq
observed_freqs <- number_network_df_positive$Freq
# Initialize a vector for pseudop-values
pvalues <- numeric(length(pairs))

# Calculate pseudop-values
for (i in seq_along(as.character(pairs))) {
  pair <- as.character(pairs[i])
  
  # Get randomized frequencies for the current pair
  randomized_freqs <- results_df$Freq[results_df$Var1 == pair]
  
  # Get observed frequency for the current pair
  observed_value <- observed_freqs[i]
  
  # Calculate the pseudop-value
  pvalue <- (sum(randomized_freqs >= observed_value) + 1) / (length(randomized_freqs) + 1)
  
  # Store the p-value
  pvalues[i] <- pvalue
}

# Create a dataframe with pairs and their p-values
pvalues_df <- data.frame(pair = pairs, pvalue = pvalues)
pvalues_df
families_bar_plot <- c("Alteromonadaceae","Alteromonadales","Campylobacteraceae","Cocleimonas","Colwelliaceae","Cryomorphaceae","Cytophagaceae","Desulfobacteraceae","Desulfobulbaceae","Flammeovirgaceae","Flavobacteraceae","Helicobacteraceae","Oceanospirillaceae","Prolixibacteraceae","Psychromonadaceae","Rhodobacteraceae","Saprospiraceae","Anaerolineaceae","Candidatus Carsonella","Nannocystaceae","Planctomycetaceae","Porphyromonadaceae","Verrumicrobiaceae")

length(families_bar_plot)
pvalues_df_separated<- separate(pvalues_df,col = pair, into = c("Sp1","Sp2"), sep = "-") %>%
  filter(Sp1 %in% families_bar_plot & Sp2 %in% families_bar_plot)



# Create a unique list of families
families <- unique(c(pvalues_df_separated$Sp1, pvalues_df_separated$Sp2))

# Initialize a square matrix with row/column names as families
similarity_matrix <- matrix(NA, nrow = length(families), ncol = length(families), 
                            dimnames = list(families, families))

# Populate the matrix with similarity values
for (i in 1:nrow(pvalues_df_separated)) {
  row <- pvalues_df_separated$Sp1[i]
  col <- pvalues_df_separated$Sp2[i]
  similarity_matrix[row, col] <- pvalues_df_separated$pvalue[i]
  similarity_matrix[col, row] <- pvalues_df_separated$pvalue[i]  # Make the matrix symmetric
}


# Replace NA with 0 (if needed for visualization)
similarity_matrix[is.na(similarity_matrix)] <- 1

# Print the matrix
print(similarity_matrix)

# pvalues_df_separated_wide <- pivot_wider(pvalues_df_separated,names_from = c("Sp1") , values_from = pvalue )

# pvalues_df_separated_wide[is.na(pvalues_df_separated_wide)] <- 1

# heatmap_matrix <- as.matrix(pvalues_df_separated_wide[,-1])
# rownames(heatmap_matrix) <- colnames(heatmap_matrix)

# setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc/figures/Heatmaps_families_correlated")
fig_dir <- file.path(work_dir,"..","figures")
fig_dir <- normalizePath(fig_dir)
setwd(fig_dir)
pdf(paste("Heatmap_all_modules_networtk_pair_families_signifcance.pdf",substrate,sep = "_"),12,12)
heatmap.2(
  log10(similarity_matrix) + 1*10^-8, 
  trace = "none",
  margins = c(12,12),
  cellnote = round(similarity_matrix,2),
  notecol = "black",
  key = T,
  key.title = "Log10 of the p_value"
  
  # scale = "row"
)
dev.off()
}
# 
# # Filter significant pairs (e.g., pvalue <= 0.05 for significance)
# significant_pairs <- filter(pvalues_df, pvalue <= 0.01)
# 
# # Output
# print("Number of significant pairs:")
# print(nrow(significant_pairs))
# 
# # Inspect significant pairs
# View(significant_pairs)
# length(significant_pairs$pair)
# length(pvalues_df$pair)
# 
# 
# # incoherent_families vector 
# incoherent_families <- c("Alteromonadaceae", "Alteromonadales_incertae_sedis", "Rhodobacteraceae","Flavobacteriaceae")
# 
# # Check if any family from incoherent_families is in the pair column
# significant_pairs_filtered <- significant_pairs %>%
#   filter(!apply(significant_pairs, 1, function(row) {
#     any(sapply(incoherent_families, function(fam) grepl(fam, row['pair'])))
#   }))
# 
# # View the filtered result
# View(significant_pairs_filtered)
# 
# 
# setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc/group_analysis")
# write_tsv(significant_pairs_filtered,"significant_paired_families.tsv")











##############################################################################################################3

# pvalues <- c()
# for (pair in pairs) {
#   pseudo_pvalue <- results_df$Freq[results_df$Var1 == pair]
#   sorted_pseudo_pvalue <- sort(pseudo_pvalue_test)
#   
#   observed_value <- number_links_df$Freq[number_links_df$Var1 == pair] 
#   pvalue <- sum(observed_value >= sorted_pseudo_pvalue) / length(sorted_pseudo_pvalue)
#   pvalues <- c(pvalues,pvalue)
# }
# 
# 
# pvalues_df <- data.frame(pair = pairs, pvalue = pvalues)
# 
# sum(pvalues_df$pvalue >= 0.99)
# length(pvalues_df$pvalue)
# pvalue
# 
# pvalues_df_filtered <- filter(pvalues_df, pvalue >= 0.99)
# View(pvalues_df_filtered)

####################################################################3

# results_df <- data.frame(Var1 = character(), Freq = numeric())
# 
# for (i in 1:length(results_list)) {
#   results_df <- rbind(results_list[[i]],results_df)
#   
# }
# results_df
# pairs <-  unique(number_links_df$Var1)
# 
# mean_sd_df <- data_frame(pair = character(),mean = numeric(),sd = numeric())
# z = 1
# for (pair in pairs) {
#   results_df_filtered <- results_df %>%
#     filter(Var1 == pair)
#   
#   mean_family <- mean(results_df_filtered$Freq)
#   sd_family <- sd(results_df_filtered$Freq)
#   
#   provisional_df <- data.frame(pair = pair, mean = mean_family,sd = sd_family)
#   mean_sd_df <- rbind(mean_sd_df, provisional_df)
#   z <- z + 1
#   print(z)
#   
# }
# 
# mean_sd_df_filtered <- filter(mean_sd_df, pair %in% number_links_df$Var1)
# 
# Z_df <- data.frame(pair= pairs, Z = (mean_sd_df_filtered$mean - number_links_df$Freq)/mean_sd_df_filtered$sd)
# length(Z_df$Z)
# View(Z_df)
# 
# library(pheatmap)
# any(!is.finite(Z_df$Z))
# sum(anyNA(Z_df$pair))
# which(!is.finite(Z_df$Z), arr.ind = TRUE)
# 
# Z_df$Z[which(!is.finite(Z_df$Z), arr.ind = TRUE)] <- 0
# Z_df <- filter(Z_df, abs(Z) > 2.5)
# sum(abs(Z_df$Z) > 2.5)


#### QUIZA HACER EL ANALISIS POR SUTRATO, Y NO EN GENERAL, ASI SE VERIA COMO LAS FAMILIAS SON COHERENTES EN TODOS LOS SUSTRATOS
####HABRIA QUE RANDOMIZAR DE UNA MANERA LA CUAL SE PUEDAN IGUALAR LOS PARES (FILTRAR LA TAXONOMIA DE MANERA QUE SOLO ESTEN LOS OTUS QUE VAMOS A USAR EN EL ANALISIS??)