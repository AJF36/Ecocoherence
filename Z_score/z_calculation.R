rm(list = ls())
library(readr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gplots)
library(ape)
library(data.tree)
library(dendextend)
library(this.path)

work_dir <- this.dir()
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")


###Load the randomized data
setwd(work_dir)
random_data <- read_tsv("X_mean_dv_randomized.tsv")


###almacen de resultados
results_df <- data.frame(Substrate = character(),
                         family = character(),
                         z = numeric())
for (x in substrates) {
   ###Load the observed values
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",x,sep=""))
  # observed_data <- read_tsv("observed_S_X_family.tsv")
  observed_data <- read_tsv(paste(x,"observed_S_X_family.tsv",sep = "_"))
  colnames(observed_data)[colnames(observed_data) == "Family"] <- "family"
  ###Filter the data of the substrates
  random_data_filtered <-  random_data %>% filter(substrate == x)
  
  ###Join the tables by family
  merged_tables <- inner_join(random_data_filtered, observed_data , by = "family")
  
  ### Calculate Z
  merged_tables <- merged_tables %>% mutate(Z = (merged_tables$X - merged_tables$mean_X)/merged_tables$sd_X)
  
  ###store the results 
  results_df <- rbind(results_df, data.frame(substrate = x,
                                                       family = merged_tables$family,
                                                       Z = merged_tables$Z))
}
View(results_df)
sum(results_df$Z > 0)
length(results_df$Z)
sum(results_df$Z < -2,5)

setwd(work_dir)
write_tsv(results_df, "z_score_families.tsv")
###Les try the heatmap

heatmap_data <- results_df %>%
  pivot_wider(names_from = substrate, values_from = Z)

# Convert to a matrix for visualization in heatmap.2
matrix_data <- as.matrix(heatmap_data[,-1])  # Remove family column
rownames(matrix_data) <- heatmap_data$family
matrix_data

heatmap_long <- results_df %>%
  mutate(Z = ifelse(is.na(Z), 0, Z))  # Replace NAs if needed

ggplot(heatmap_long, aes(x = substrate, y = family, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Z-scores",
       x = "Substrate",
       y = "Family",
       fill = "Z") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

matrix_data
View(matrix_data)
matrix_data[is.na(matrix_data)] <- 0
row_clust <- hclust(dist(matrix_data))
plot(row_clust)
cellnote_matrix <- ifelse(abs(matrix_data) > 2.5, "*","")



fig_dir <- file.path(work_dir,"..","figures")
fig_dir <- normalizePath(fig_dir)
setwd(fig_dir)

png("heatmap_coherent_families_inferno_test.png",1500,1500)
heatmap.2(matrix_data,
           col = viridis::inferno(50),         # Escala de colores
           trace = "none",             # Sin líneas de traza
           na.color = "white",
           key = T , 
           density.info = "none" , 
           key.par = list(cex = 1.1, font = 2) ,  
           key.title = "z score value", 
           key.xlab = NA, 
           margins = c(16, 20),          # Márgenes reducidos
           cexRow = 1.8,               # Tamaño del texto de las filas
           cexCol = 2,               # Tamaño del texto de las columnas
           lwid = c(3, 13),             # Relación del ancho de columnas
           lhei = c(1, 15),             # Relación del alto de filas
           srtCol = 45,
           cellnote = cellnote_matrix,
           notecol = "yellow",
           notecex = 3,
          hclustfun = function(x) hclust(x, method="complete"))                # R
dev.off()
# ?png
# ?heatmap.2



# heatmap_data_mod <- subset(heatmap_data,select= -family)

# df <- data.frame(lapply(heatmap_data_mod, function(x) {
#   case_when(
#     x >= 1 ~ -1,
#     x <= -1 ~ +1,
#     TRUE ~ 0
#   )
# }))
# rownames(df) <- heatmap_data$family


# names(rowSums(df)[rowSums(df) != 0])
# df
# plot(hclust(dist(df)))
# df_f <- filter(df,rownames(df) %in% names(rowSums(df)[rowSums(df) != 0]))
# plot(hclust(dist(df_f,method = "euclidian"),method = "average"))
# plot(hclust(dist(rowSums(df)[rowSums(df) != 0])))
# Handle NAs if needed
# df[is.na(df)] <- 0.5
# df <- data.frame(lapply(heatmap_data_mod, function(x) {
#   ifelse(x >= 2.5, 1, ifelse(x <= -2.5, -1, x))
# }))
# df[df != 1 & df != -1] <- 0
# df[is.na(df)] <- 2



# matrix_data_mod <- as.matrix(df_f)
# row_clust_mod <- hclust(dist(df_f))
# plot(row_clust_mod)

# getwd()
# png("heatmap_coherent_families_2.png",1500,1500)
# heatmap.2(matrix_data_mod,
#           col = bluered(100),         # Escala de colores
#           trace = "none",             # Sin líneas de traza
#           na.color = "white",         # Color para los valores NA
#           margins = c(12, 14),          # Márgenes reducidos
#           cexRow = 1.3,               # Tamaño del texto de las filas
#           cexCol = 1.5,               # Tamaño del texto de las columnas
#           lwid = c(1, 4),             # Relación del ancho de columnas
#           lhei = c(1, 4),             # Relación del alto de filas
#           srtCol = 45,
#           hclustfun = function(x) hclust(x, method="complete"))     
# dev.off()
# sapply


# # Create a custom color palette
# # For values -1, 0, 1
# my_palette <- colorRampPalette(c("red","white", "blue"))(3)
# # Add gray for value 2 (absent)
# my_colors <- c(my_palette, "gray50")

# # Create breaks to ensure proper color mapping
# my_breaks <- c(-1.5, -0.5, 0.5, 1.5, 2.5)

# png("heatmap_coherent_families_2.png", 1500, 1500)
# heatmap.2(matrix_data_mod,
#           col = my_colors,            # Custom color palette
#           breaks = my_breaks,         # Custom breaks for color mapping
#           trace = "none",             # No trace lines
#           key = TRUE,                 # Include color key
#           keysize = 1.5,              # Size of color key
#           density.info = "none",      # No density info
#           margins = c(12, 14),        # Margins
#           cexRow = 1.3,               # Row text size
#           cexCol = 1.5,               # Column text size
#           lwid = c(1, 4),             # Width ratio
#           lhei = c(1, 4),             # Height ratio
#           srtCol = 45,                # Column name rotation
#           key.title = "Legend",       # Legend title
#           key.xlab = "Value",         # X-axis label for legend
#           key.ylab = "Count",         # Y-axis label for legend
#           # Custom labels for the legend
#           key.xtickfun = function() {
#             list(at = c(-1, 0, 1, 2), 
#                  labels = c("Low (-1)", "Medium (0)", "High (1)", "Absent (2)"))
#           })
# dev.off()
# ####Comparation between clustering and taxonomy
# ??heatmap.2

##Load the tax table
# Cargar la tabla de taxonomía
# fileTaxonomy <- "/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
# taxa.in <- read.table(fileTaxonomy, sep = ";")
# colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
#                        "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
# taxonomy <- subset(taxa.in, select = c("taxa_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
# rownames(taxonomy) <- taxonomy$taxa_id
# taxonomy <- subset(taxonomy, select = -c(taxa_id))
# taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))
# 
# 
# ###Take only the families present in the clustering
# clustering_families <- rownames(matrix_data)
# taxonomy_filtered <- filter(taxonomy,Family %in% clustering_families )
# str(taxonomy_filtered)
# str(taxonomy)
# 
# families<- unique(taxonomy_filtered$Family)
# taxons <- colnames(taxonomy_filtered)
# taxons <- taxons[-4]
# taxons
# orders <-character()
# phylums <- character()
# classes <-character()
# i <- 1
# for (family in families) {
#   orders[i] <- unique(taxonomy$Order[taxonomy$Family == family])
#   classes[i] <- unique(taxonomy$Class[taxonomy$Family == family])
#   phylums[i] <- unique(taxonomy$Phylum[taxonomy$Family == family])
#   i <- i + 1
# }
# length(orders) == length(families)
# 
# tax_df <- data.frame(phylums,classes,orders,families)
# tax_df
# tax_matrix <- as.matrix(tax_df)
# tax_df$jerarquia <- apply(tax_matrix, 1, paste, collapse = "/")
# taxonomy_tree <- as.Node(data.frame(pathString = paste("Root", tax_df$jerarquia, sep = "/")))
# taxonomy_tree
# phylo_tree <- as.phylo(taxonomy_tree)
# plot(phylo_tree, show.node.label = TRUE , cex = 0.7)
# taxonomy_dendogram <- as.dendrogram(taxonomy_tree)
# cluster_dendogram <- as.dendrogram(row_clust)
# 
# list_dends <- dendlist(taxonomy_dendogram,cluster_dendogram)
# x11(width = 12, height = 8)
# 
# tax_df$orders
# 
# # Create a named vector of colors for each class
# unique_classes <- tax_df$orders
# 
# class_colors <- setNames(rainbow(length(unique_classes)), unique_classes)
# class_colors
# # Map families to their corresponding class colors
# family_colors <- class_colors[match(taxonomy$Order, names(class_colors))]
# names(family_colors) <- taxonomy$Family
# family_colors
# 
# # Get labels from the first dendrogram (taxonomy_dendogram)
# dendrogram_labels <- labels(taxonomy_dendogram)
# 
# # Assign colors based on the families in the taxonomy data
# link_colors <- family_colors[match(dendrogram_labels, names(family_colors))]
# 
# 
# 
# tanglegram(list_dends,
#            lab.cex = 1,                              # Adjust label size
#            margin_inner = 10,                          # Increase inner margins
#            color_lines = link_colors,                  # Apply class-based colors
#            highlight_distinct_edges = FALSE,           # Optional
#            common_subtrees_color_lines = TRUE)         # Optional
