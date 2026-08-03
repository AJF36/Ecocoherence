rm(list = ls())
### Clustering of the misosouop communties
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(dendextend)
library(yaml)
work_d <- "/home/ajf/Desktop/CNB/ecocoherence/metabolism/misosoup/results_misosoup"
setwd(work_d)

communties_df <-list.files(path = ".",pattern = ".tsv")

list_df <- list()
communties_df
for (df in communties_df){
  dummy_df <- read_tsv(df)
  dummy_df <- as.data.frame(dummy_df)
  class(dummy_df)
  list_df[[length(list_df ) + 1]] <- dummy_df
}

all_df <- bind_rows(list_df)


#Viena figure
#################
# 
# # 1. Suponiendo que 'Community' tiene etiquetas como '1_mod_1', '47_mod_48', etc.
# # Primero, guardamos esa columna como CommunityID
# all_df <- all_df %>%
#   mutate(CommunityID = Community)
# 
# # 2. Creamos una nueva columna con la categoría de comunidad: Early, Late, All
# all_df <- all_df %>%
#   mutate(Type = case_when(
#     grepl("_mod_48", CommunityID) ~ "Early",
#     grepl("_mod_1", CommunityID) ~ "Late",
#     grepl("_all", CommunityID) ~ "All",
#     TRUE ~ "Other"
#   ))
# 
# # 3. Ahora contamos cuántos miembros tiene cada comunidad individual
# community_sizes <- all_df %>%
#   group_by(Type, CommunityID) %>%
#   summarise(n_members = n(), .groups = "drop")
# 
# # 4. Contamos cuántas comunidades hay de cada tamaño, para cada tipo
# plot_df <- community_sizes %>%
#   group_by(Type, n_members) %>%
#   summarise(n_communities = n(), .groups = "drop")
# 
# # 5. Finalmente, hacemos el gráfico
# ggplot(plot_df, aes(x = n_members, y = n_communities, fill = Type)) +
#   geom_col(position = "dodge") +
#   labs(title = "",
#        x = "Number of Members in Community",
#        y = "Number of Communities") +
#   scale_fill_manual(values = c("Early" = "blue", "Late" = "red", "All" = "orange")) +
#   theme_minimal() +
#   theme(
#     axis.title.x = element_text(size = 16),
#     axis.title.y = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 16),
#     legend.text = element_text(size = 14),
#     plot.title = element_text(size = 18, face = "bold")
#   ) +
#   guides(fill = guide_legend(title = "Modules"))
# 
# 
# ggplot(plot_df, aes(x = n_members, y = n_communities, color = Type, group = Type)) +
#   geom_col(aes(fill = Type), position = "dodge", width = 0.8, alpha = 0.6) +
#   geom_smooth(stat = "identity", size = 1.2) +
#   geom_point(size = 2) +
#   labs(title = "Community Size Distribution by Module",
#        x = "Number of Members in Community",
#        y = "Number of Communities") +
#   scale_fill_manual(values = c("Early" = "blue", "Late" = "red", "All" = "orange")) +
#   scale_color_manual(values = c("Early" = "blue", "Late" = "red", "All" = "orange")) +
#   theme_minimal() +
#   theme(
#     axis.title.x = element_text(size = 16),
#     axis.title.y = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 16),
#     legend.text = element_text(size = 14),
#     plot.title = element_text(size = 18, face = "bold")
#   ) +
#   guides(fill = guide_legend(title = "Modules"), color = guide_legend(title = "Modules"))
# 
#                                 
# ###figure viena
# plot_df <- all_df %>%
#   group_by(Community, n_members) %>%
#   summarise(n_communities = n(), .groups = "drop")
# 
# # Step 1: Clean up community type labels
# all_df$Community <- gsub("\\d+_mod_48", "Early", all_df$Community)
# all_df$Community <- gsub("\\d+_mod_1", "Late", all_df$Community)
# all_df$Community <- gsub("\\d+_all", "All", all_df$Community)
# 
# # Step 2: Create community IDs (optional but useful)
# # Let's assume you originally had a column that identified each unique community (e.g. "1_mod_1")
# # If not, recreate that now:
# # all_df$CommunityID <- paste0(all_df$OriginalCommunityID)  # Replace with actual column if needed
# 
# # Step 3: Count members per community
# community_sizes <- all_df %>%
#   group_by(Community, Models) %>%
#   summarise(n_members = n(), .groups = "drop")
# 
# # Step 4: Count how many communities of each type have n members
# plot_df <- community_sizes %>%
#   group_by(Community, n_members) %>%
#   summarise(n_communities = n(), .groups = "drop")
# 
# # Step 5: Plot
# ggplot(plot_df, aes(x = n_members, y = n_communities, fill = Community)) +
#   geom_col(position = "dodge") +
#   labs(title = "",
#        x = "Number of members in the community",
#        y = "Number of communities") +
#   scale_fill_manual(values = c("Early" = "blue", "Late" = "red", "All" = "orange")) +
#   theme_minimal() +
#   theme(
#     axis.title.x = element_text(size = 16),
#     axis.title.y = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 16),
#     legend.text = element_text(size = 14),
#     plot.title = element_text(size = 18, face = "bold")
#   )
# 
# 
# # Plot it
# ggplot(plot_df, aes(x = n_members, y = n_members, fill = Community)) +
#   geom_col(position = "dodge") +
#   labs(title = "Number of communities by size and type",
#        x = "Number of members in community",
#        y = "Number of communities") +
#   scale_fill_brewer(palette = "Set2") +
#   theme_minimal()
###Extract all the models (rows of the clustering)

models <- unique(all_df$Models)

#### Extract all the communities

communities <- unique(all_df$Community)

### Create the matrix

clustering_matrix <- matrix(0,nrow = length(models),ncol = length(communities) )
rownames(clustering_matrix) <- models 
colnames(clustering_matrix) <- communities


grouped_df <- all_df %>%
  group_by(Models) %>% 
  summarize(
    communities = paste(Community,collapse = ";")
  )
rownames(grouped_df) <- grouped_df$Models

### Fill the matrix
# model <- "GCA_002869565"
for (model in models){
  model_communities <- grouped_df[model,]
  model_communities <- grouped_df$communities[match(model,grouped_df$Models)]
  model_communities <- strsplit(model_communities,split = ";")
  model_communities <- model_communities[[1]]
  
  for (communty in model_communities){
  clustering_matrix[model,communty] <- clustering_matrix[model,communty] +1
  }
    
}

View(clustering_matrix)
distance_jaccard <- vegdist(t(clustering_matrix),method = "jaccard")
clustered_communities <- hclust(distance_jaccard)

###Filter the communities that are equal between sets
dist.matrix <- as.matrix(distance_jaccard)
col.sums.dist.matrix <- colSums(as.matrix(distance_jaccard))
index.matrix <- col.sums.dist.matrix != round(col.sums.dist.matrix)

dist.matrix <- dist.matrix[index.matrix,index.matrix]
plot(hclust(as.dist(dist.matrix)))
clustered_communities<- hclust(as.dist(dist.matrix))


###Check the number of clusters
dendo <- as.dendrogram(clustered_communities)
dendo %>% set("branches_k_color" , k = 6) %>%
  plot()

clusters <- cutree(clustered_communities, k = 6)
table(clusters)
#### Namew of the communities that are interesting for the analysis
clusters.names <- names(clusters)
clusters.names
clusters

setwd("/home/ajf/Desktop/CNB/ecocoherence/misosoupR/minimal_communities")
files <- c("misosoup_Agarose.yaml","misosoup_Agarose_mod_48.yaml")

for (cluster in unique(clusters)){
  setwd("/home/ajf/Desktop/CNB/ecocoherence/misosoupR/minimal_communities")
  
  cluster.names <- clusters.names[clusters == cluster]
  cluster.names
  list.communities.all <- list()
  list.communities.48 <- list()
  for (file in files){
    if (file == "misosoup_Agarose_mod_48.yaml"){
      pattern <- "mod_48"
      clusters.names.filter <- cluster.names[grepl(pattern, cluster.names)]
      clusters.names.filter <- gsub("(\\d+).*","\\1",clusters.names.filter)
      clusters.names.filter <- as.numeric(clusters.names.filter)
      yaml.file <- yaml.load_file(file)
      communities.yaml <- yaml.file[[1]]$min
      list.communities.48 <- communities.yaml[clusters.names.filter]
    }
    if (file == "misosoup_Agarose.yaml"){
      pattern <- "all"
      clusters.names.filter <- cluster.names[grepl(pattern, cluster.names)]
      clusters.names.filter <- gsub("(\\d+).*","\\1",clusters.names.filter)
      clusters.names.filter <- as.numeric(clusters.names.filter)
      yaml.file <- yaml.load_file(file)
      communities.yaml <- yaml.file[[1]]$min
      list.communities.all <- communities.yaml[clusters.names.filter]
    }
  }
  length(list.communities.all)
  length(list.communities.48)
  cluster
  full.list <- c(list.communities.48,list.communities.all)
  full.list.min <- list("min" = full.list)
  full.list.gal <- list("gal" = full.list.min)
  length(full.list)
  yaml_string <- as.yaml(full.list.gal)
  setwd("/home/ajf/Desktop/CNB/ecocoherence/misosoupR/minimal_communities/cluster_communities")
  write_lines(yaml_string,paste("cluster_",cluster,sep=""))
  
}


#### I delete cluster 1 and 5 because they are irrelevant for the ananylisis (2 clusters of the same)

test.list <- list("min" = list.communities.all)
test.list2 <- list("gal" = test.list)
