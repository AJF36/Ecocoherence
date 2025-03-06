######################################################
#Script for comparing the modules between substrates 
#####################################################


rm(list = ls())
# Load necessary libraries
library(tidyr)
library(phyloseq)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(dendextend)
library(gplots)
library(gridExtra)
library(RColorBrewer)
library(this.path)

work_dir <- this.dir()
# Define substrates
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")
# substrate <- "Agarose"
# Read the all the families of all the modules

families <- c()
for (substrate in substrates) {
  path_file <- file.path(work_dir,"..","modules_table")
  path_file <- normalizePath(path_file)
  # families_table <- read_tsv(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, "/group_analysis_", substrate, ".tsv", sep = ""))
  families_table <- read_tsv(file.path(path_file,paste("group_analysis_", substrate, ".tsv", sep = "")))
  families_module <- families_table$Family
  families <- c(families,families_module)
                             
}
families <- unique(families_module)
# incoherent_families <- c("Alteromonadaceae", "Alteromonadales","Alteromonadales_incertae_sedis", "Rhodobacteraceae","Flavobacteriaceae","Pseudomonadales","Desulfobulbaceae","Candidatus Carsonella","Desulfobacteraceae","Hyphomonadaceae")
# families <- families[!(families %in% incoherent_families)] #### Removing the most incoherent families

all_modules_df <- data.frame(Family = families)
# families_table <- read_tsv(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", "Alginate", "/group_analysis_", "Alginate", ".tsv", sep = ""))
# all_modules_df <- data.frame("Family" = families_table$Family)

###Filter the incoherent families




# Process each substrate
for (substrate in substrates) {
  ## Load the table
  # group_analysis_table <- read_tsv(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, "/group_analysis_", substrate, ".tsv", sep = ""))
  path_file <- file.path(work_dir,"..","modules_table")
  path_file <- normalizePath(path_file)
  # families_table <- read_tsv(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, "/group_analysis_", substrate, ".tsv", sep = ""))
  group_analysis_table <- read_tsv(file.path(path_file,paste("group_analysis_", substrate, ".tsv", sep = "")))
  
  ## Select the columns with the modules
  numeric_columns <- names(group_analysis_table)[!is.na(as.numeric(names(group_analysis_table)))]
  df_numeric <- group_analysis_table[, numeric_columns, drop = FALSE]
  df_numeric<- df_numeric[, colSums(df_numeric) >= 15]
  ## Calculate the total members of each module
  total_members_modules <- colSums(df_numeric)
  total_members_modules <- total_members_modules[total_members_modules >= 15]
  

  
  ## Make the columns proportional
  df_proportional <- t(t(df_numeric) / total_members_modules)
  df_proportional <- as_data_frame(df_proportional)
  
  ## Add the family column
  df_proportional$Family <- group_analysis_table$Family
  
  # df_proportional <- filter(df_proportional, !(Family %in% incoherent_families))
  
  ## Rename module columns to include substrate name
  colnames(df_proportional)[-ncol(df_proportional)] <- paste0(colnames(df_proportional)[-ncol(df_proportional)], "_", substrate)
  
  ## Join with the final dataframe
  # all_modules_df <- inner_join(all_modules_df, df_proportional, by = "Family")
  # all_modules_df <- inner_join(all_modules_df, df_proportional, by = "Family") ####without the incoherent families
  all_modules_df <- full_join(all_modules_df, df_proportional, by = "Family") #### for all the families
}


# Replace NA with 0 in the data
all_modules_df[is.na(all_modules_df)] <- 0

# Extract the numeric matrix from the dataframe
data_matrix <- as.matrix(all_modules_df[,-1])

# Compute the distance matrix for rows (families) and columns (modules)
dist_rows <- dist(data_matrix)  # Distance for families (rows)
dist_cols <- dist(t(data_matrix))  # Distance for modules (columns)

# Perform hierarchical clustering for rows and columns
hclust_rows <- hclust(dist_rows, method = "average")  # For families
hclust_cols <- hclust(dist_cols, method = "average")  # For modules

hclust_av = function(x){hclust(x, method = "average")}

#Define custom color palette with more bins and higher contrast
custom_palette <- colorRampPalette(c("blue", "white", "red"))(200)  # 100 bins

# Optional: Adjust breaks for more noticeable transitions
# breaks <- seq(min(data_matrix), max(data_matrix), length.out = 201)



margins = c(12, 8)
fig_path <- file.path(work_dir,"..","figures")
fig_path  <- normalizePath(fig_path)
# setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc/figures")
setwd(fig_path)
pdf("heatmap_modules_families_all_modules_no_incoherent2.pdf", width = 10, height = 8)
#
# # png("heatmap_high_res_modules.png", width = 4000, height = 2000, res = 300)
heatmap.2(
  x = data_matrix,
  Rowv = T, #as.dendrogram(hclust_rows),  # Enable clustering for rows
  Colv = T, #as.dendrogram(hclust_cols),  # Enable clustering for columns
  col = custom_palette,  # Custom color palette
  hclustfun = hclust_av,
  # main = "Heatmap with Clustering",  # Title
  trace = "none",  # Remove trace lines inside the heatmap
  dendrogram = "both",  # Display dendrograms for rows and columns
  labRow = all_modules_df$Family,  # Row labels (Family names)
  labCol = colnames(all_modules_df)[-1],  # Column labels (Module names)
  cexRow = 0.4,  # Font size for row labels
  cexCol = 0.2,  # Font size for column labels
  key = TRUE,  # Show color key
  keysize = 1.2,  # Size of the color key
  key.title = "Relative abundance",  # Title for the legend
  key.xlab = "Z-score of the modules",  # X-axis label for the key
  margins = margins,  # Margins for row and column labels
  scale = "col"  # Scale rows before plotting
)

heatmap.2(
  x = data_matrix,
  Rowv = T, #as.dendrogram(hclust_rows),  # Enable clustering for rows
  Colv = T, #as.dendrogram(hclust_cols),  # Enable clustering for columns
  col = custom_palette,  # Custom color palette
  hclustfun = hclust_av,  
  # main = "Heatmap with Clustering",  # Title
  trace = "none",  # Remove trace lines inside the heatmap
  dendrogram = "both",  # Display dendrograms for rows and columns
  labRow = all_modules_df$Family,  # Row labels (Family names)
  labCol = colnames(all_modules_df)[-1],  # Column labels (Module names)
  cexRow = 0.4,  # Font size for row labels
  cexCol = 0.2,  # Font size for column labels
  key = TRUE,  # Show color key
  keysize = 1.2,  # Size of the color key
  key.title = "Relative abundance",  # Title for the legend
  key.xlab = "Z-score of the families",  # X-axis label for the key
  margins = margins,  # Margins for row and column labels
  scale = "row"  # Scale rows before plotting
)
dev.off()
pdf("heatmap_modules_families_all_modules_no_scale.pdf", width = 10, height = 8)
#
# # png("heatmap_high_res_modules.png", width = 4000, height = 2000, res = 300)
heatmap.2(
  x = -log10(data_matrix + 1e-08),
  Rowv = T, #as.dendrogram(hclust_rows),  # Enable clustering for rows
  Colv = T, #as.dendrogram(hclust_cols),  # Enable clustering for columns
  #col = custom_palette,  # Custom color palette
  hclustfun = hclust_av,
  # main = "Heatmap with Clustering",  # Title
  trace = "none",  # Remove trace lines inside the heatmap
  dendrogram = "both",  # Display dendrograms for rows and columns
  labRow = all_modules_df$Family,  # Row labels (Family names)
  labCol = colnames(all_modules_df)[-1],  # Column labels (Module names)
  cexRow = 0.4,  # Font size for row labels
  cexCol = 0.5,  # Font size for column labels
  key = TRUE,  # Show color key
  keysize = 1.2,  # Size of the color key
  key.title = "",  # Title for the legend
  key.xlab = "-log10(Relative abundances)",  # X-axis label for the key
  margins = margins,  # Margins for row and column labels
  #scale = "col"  # Scale rows before plotting
)
dev.off()

head(data_matrix)
colSums(data_matrix)
colSums(log10(data_matrix+1e-08))
# stop()
# # Compute the distance matrix between modules (transpose because we compare columns)
dist_matrix <- dist(t(as.matrix(all_modules_df[,-1])))
# dist_matrix_families <- dist(as.matrix(all_modules_df[,-1]))
# cluster_families <- hclust(dist_matrix_families,method = "average")
# # Define module names
# modules_names <- colnames(all_modules_df)[-1]
# 
# # Heatmap visualization


# Clustered dendrogram
par(mar = c(10, 5, 5, 5))  # Adjust margins
# par(cex = 0.9)            # Adjust label size
clustered_matrix <- hclust(dist_matrix, method = "average")
# plot(clustered_matrix, main = "Clustered Dendrogram")


dendogram <- as.dendrogram(clustered_matrix)
dendogram  %>% set("labels_cex" , 0.9) %>% set("labels_col",k=12) %>% plot
dendogram %>% rect.dendrogram(k=12, border = 8,lty = 5, lwd = 2 )
plot(dendogram)

###Get the module names from the cluster
# Cut the dendrogram into `k` clusters
k <- 12 # Number of clusters
clusters <- cutree(clustered_matrix, k = k)

# Create a list of labels for each cluster
cluster_labels <- split(names(clusters), clusters)

# Print the labels of each cluster
print(cluster_labels)

early_modules <- cluster_labels[[2]]
early_modules
df_early_modules <- select(all_modules_df, c("Family",early_modules))

colSums(df_early_modules[,-1])

late_modules <- cluster_labels[[1]]

df_late_modules <- select(all_modules_df, c("Family",late_modules))

colSums(df_early_modules[,-1])


####Make the barplots of the early and the late modules

tax_table <- tax_table(as.matrix(df_early_modules["Family"]))

# count_table_early <- otu_table(df_early_modules[,-1],taxa_are_rows = TRUE)
# phylo_obj <- merge_phyloseq(tax_table,count_table_early)
# threshold <- 0.03
# phylo_obj <- prune_taxa(apply(otu_table(phylo_obj), 1, max) >= threshold, phylo_obj)
# phylo_obj = transform_sample_counts(phylo_obj, function(x) 100 * x/sum(x))
# early_abundace_plot <- plot_bar(phylo_obj , fill = "Family") +
#   theme(
#     legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys
#     legend.text = element_text(size = 8),  # Adjust font size of legend text
#     legend.title = element_text(size = 10)  # Adjust font size of legend title
#   )
# 
# ###Lets make the second layer of the plot to see the abundance over time of the modules
# 
# count_table_late <- otu_table(df_late_modules[,-1],taxa_are_rows = TRUE)
# phylo_obj <- merge_phyloseq(tax_table,count_table_late)
# threshold <- 0.03
# phylo_obj <- prune_taxa(apply(otu_table(phylo_obj), 1, max) >= threshold, phylo_obj)
# phylo_obj = transform_sample_counts(phylo_obj, function(x) 100 * x/sum(x))
# late_abundace_plot <- plot_bar(phylo_obj , fill = "Family") +
#   theme(
#     legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys
#     legend.text = element_text(size = 8),  # Adjust font size of legend text
#     legend.title = element_text(size = 10)  # Adjust font size of legend title
#   )
# 


# Prepare early data
count_table_early <- otu_table(df_early_modules[,-1], taxa_are_rows = TRUE)
phylo_obj_early <- merge_phyloseq(tax_table, count_table_early)
threshold <- 0.03
phylo_obj_early <- prune_taxa(apply(otu_table(phylo_obj_early), 1, max) >= threshold, phylo_obj_early)
phylo_obj_early <- transform_sample_counts(phylo_obj_early, function(x) 100 * x / sum(x))

# Prepare late data
count_table_late <- otu_table(df_late_modules[,-1], taxa_are_rows = TRUE)
phylo_obj_late <- merge_phyloseq(tax_table, count_table_late)
phylo_obj_late <- prune_taxa(apply(otu_table(phylo_obj_late), 1, max) >= threshold, phylo_obj_late)
phylo_obj_late <- transform_sample_counts(phylo_obj_late, function(x) 100 * x / sum(x))

# Get the common levels of "Family" across both datasets
common_families <- union(tax_table(phylo_obj_early), tax_table(phylo_obj_late))
tax_table(phylo_obj_early)
# Create a color palette for all families
common_families

manual_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", 
  "#BCBD22", "#17BECF", "#FFB6C1", "#87CEFA", "#D3D3D3", "#FFD700", "#32CD32", "#8A2BE2", 
  "#FF1493", "#C71585", "#3CB371", "#FFD700", "#8B008B", "#FF6347", "#B22222", "#A52A2A", 
  "#5F9EA0", "#D2691E", "#CD5C5C", "#B8860B", "#2F4F4F", "#98FB98", "#AFEEEE", "#7FFF00","red","blue","green","violet"
)
family_colors <- setNames(manual_palette, common_families)

# Create the early plot
early_abundace_plot <- plot_bar(phylo_obj_early, fill = "Family") +
  scale_fill_manual(values = family_colors) +  # Apply the common color palette
  theme(
    legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys
    legend.text = element_text(size = 12),  # Adjust font size of legend text
    legend.title = element_text(size = 10)  # Adjust font size of legend title
  ) +
  xlab("Time (h)") +  # Example axis label
  ylab("Relative Abundance") +
  theme(
    legend.text = element_text(size = 12),        # Adjust font size of legend text
    legend.key.size = unit(1, "cm"),          # Adjust size of legend keys
    legend.title = element_text(size = 8)      # Adjust font size of legend title
  )

# Create the late plot
late_abundace_plot <- plot_bar(phylo_obj_late, fill = "Family") +
  scale_fill_manual(values = family_colors) +  # Apply the common color palette
  theme(
    legend.key.size = unit(1, "cm"),  # Adjust size of legend keys
    legend.text = element_text(size = 8),  # Adjust font size of legend text
    legend.title = element_text(size = 8)  # Adjust font size of legend title
  ) +
  xlab("Time (h)") +  # Example axis label
  ylab("Relative Abundance") 

###Load the abundance table of all the modules and make the plot
setwd(work_dir)
abundance_modules_df <- read_tsv("all_modules_relativeabundance.tsv")
early_modules
abundace_early_modules_df<- abundance_modules_df %>%
  filter(module %in% early_modules) %>% 
  pivot_longer(-module,names_to = "Time",values_to = "rel_abundance")

early_plot <-ggplot(abundace_early_modules_df,aes(x = as.numeric(Time), y =rel_abundance, colour = module)) + 
  geom_line(size = 2) + 
  theme_bw() +
  ylab("Relative abundance") +
  xlab("Time (h)")

abundace_late_modules_df<- abundance_modules_df %>%
  filter(module %in% late_modules) %>% 
  pivot_longer(-module,names_to = "Time",values_to = "rel_abundance")

discrete_palette <- c(
  "red", "blue", "green", "purple", "orange", "yellow", "pink", "brown", 
  "cyan", "magenta", "black", "gray", "violet", "lightgreen"
)


late_plot <- ggplot(abundace_late_modules_df, aes(x = as.numeric(Time), y = rel_abundance, colour = module)) + 
  geom_line(size = 2) + 
  theme_bw() +
  ylab("Relative abundance") +
  scale_colour_manual(values = discrete_palette) +
  xlab("Time (h)") +
  geom_point(data = subset(abundace_late_modules_df, module %in% c("1_Carrageenan","2_AgaroseAlginate","18_AgaroseCarrageenan","1_AgaroseChitosan","2_Alginate","9_Chitin","3_Agarose")),
             aes(x = as.numeric(Time), y = rel_abundance),  # Use rel_abundance for y
             colour = "red", size = 3)


file_path_figures <- file.path(work_dir,"..","figures")
file_path_figures <- normalizePath(file_path_figures)
setwd(file_path_figures)
# setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc/figures")
# pdf("Relative_abundance_composition_modules.pdf",width = 15,height = 9)
png("Relative_abundance_composition_modules_all_modules_2.png", 3000,2000)
grid.arrange(early_plot,late_plot,early_abundace_plot,late_abundace_plot,nrow = 2)
dev.off()


# 
# abundances_modules_df <- as.data.frame(matrix(ncol = length(times) + 1, nrow = 0))
# colnames(abundances_modules_df) <- c("module",as.character(times))
# abundances_modules_df
# 
# 
# 

# for (substrate in substrates){
#   setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))
#   abundance_modules_df <- read_tsv(paste("otu_table_",substrate,"_byModulesSize4.tsv", sep = ""))
#   abundance_modules_df
#   abundance_modules_df$module <- gsub("mod_",paste(substrate,"_",sep=""),abundance_modules_df$Beads_Agarose_0_C)
#   abundance_modules_df$module <-  gsub("^(\\w+)_(\\d+)$", "\\2_\\1", abundance_modules_df$module)
#   abundace_modules_df_percent<- sweep(abundance_modules_df[,-1], 2, colSums(abundance_modules_df[,-1]), FUN = "/") * 100
#   
#   class(abundance_modules_df)
#   
# 
#   #Extract time points from column names
#   time_points <- gsub(".*_(\\d+)_.*", "\\1", colnames(abundance_modules_df)
# 
#   # Create an empty data frame for results
#   unique_times <- unique(time_points)
#   result_df <- data.frame(matrix(nrow = nrow(df), ncol = length(unique_times)))
#   colnames(result_df) <- unique_times
# 
#   # Calculate mean for each time point and fill the result data frame
#   for (time in unique_times) {
#     replicate_columns <- grep(paste0("_", time, "_"), colnames(df))
#     result_df[[time]] <- rowMeans(df[, replicate_columns])
#   }
# 
#   # View the resulting data frame
#   print(result_df)
# 
#   # abundance_vector <- abundance_modules_df$module
#   # abundance_vector <- filter(abundance_modules_df,module %in% early_modules)
#   # abundance_vector <- as.vector(abundance_vector[[1:length(abundance_vector[[1]])]])
#   # abundance_vector
#   # class(abundance_vector)
# 
# }
# # 
# early_modules
# 
# 

# ####Make the barplots of the modules
# count_table_late <- otu_table(df_late_modules[,-1],taxa_are_rows = TRUE)
# phylo_obj <- merge_phyloseq(tax_table,count_table_late)
# plot_bar(phylo_obj , fill = "Family")
# 
# early_modules
# 
# 
# 
# 
# 
# combined_modules <- inner_join(df_early_modules, df_late_modules , by = "Family")
# 
# # Create Tax Table
# tax_table <- tax_table(as.matrix(combined_modules["Family"]))
# tax_table
# combined_modules <- combined_modules[,-1]
# # Create Count Table
# count_table <- otu_table(as.matrix(combined_modules), taxa_are_rows = TRUE)
# #Create the sample_metadata
# sample_data <- data.frame(module = colnames(count_table))
# 
# sample_data$Time_point <- c(
#   rep("Early", length(colnames(df_early_modules[,-1]))),
#   rep("Late", length(colnames(df_late_modules[,-1])))
# )
# sample_data <-sample_data(sample_data)
# 
# sample_names(sample_data) <- sample_data$module
# # rownames(count_table) <- tax_table
# 
# # Create Phyloseq Object
# phylo_obj <- phyloseq(tax_table, count_table,sample_data)
# # Prune taxa with less than 5% in all samples
# threshold <- 0.03
# phylo_obj <- prune_taxa(apply(otu_table(phylo_obj), 1, max) >= threshold, phylo_obj)
# phylo_obj = transform_sample_counts(phylo_obj, function(x) 100 * x/sum(x))
# 
# 
# # Plot Bar Chart
# manual_palette <- c(
#   "#4E79A7", # Blue
#   "#F28E2B", # Orange
#   "#E15759", # Red
#   "#76B7B2", # Teal
#   "#59A14F", # Green
#   "#EDC948", # Yellow
#   "#B07AA1", # Purple
#   "#FF9DA7", # Pink
#   "#9C755F", # Brown
#   "#BAB0AC", # Gray
#   "#86BCB6", # Light Teal
#   "#F1CE63", # Mustard
#   "#D37295", # Rose
#   "#FABFD2", # Light Pink
#   "#A8786E", # Mauve
#   "#8CD17D", # Lime
#   "#B6992D", # Olive
#   "#79706E", # Charcoal
#   "#FFBE7D", # Peach
#   "#6C4F30", # Coffee
#   "#AF7AA1", # Lavender
#   "#FF7A5C", # Coral
#   "#91B9D0", # Light Blue
#   "#FFA07A", # Salmon
#   "#C39BD3", # Orchid
#   "#A9DFBF", # Mint
#   "#E06D06", # Dark Orange
#   "#AD494A", # Brick Red
#   "#8B5E3B", # Cinnamon
#   "#729ECE", # Steel Blue
#   "#FFCB99", # Light Apricot
#   "#A1D490", # Sage Green
#   "#9EDAE5", # Powder Blue
#   "#C3C3C3", # Light Gray
#   "#D4A6C8", # Light Purple
#   "#6C6F30", # Olive Green
#   "#F8A5A5", # Soft Pink
#   "#8E0152", # Burgundy
#   "#276419"  # Deep Green
# )
# 
# plot_bar(phylo_obj, fill = "Family") +
#   scale_fill_manual(values = manual_palette) +
#   facet_wrap(~ Time_point)   # Facet by TimePoint





