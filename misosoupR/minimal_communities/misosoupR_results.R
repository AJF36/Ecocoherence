rm(list = ls())

library(readr)
library(tidyr)
library(dplyr)
library(this.path)
library(gt)
library(pheatmap)
library(vegan)

work.dir <- this.dir()
setwd(work.dir)


df.list <- list()

file.list <- list.files(pattern = "yaml.tsv")

z = 1

for (file in file.list){
  
  file.df <-read_tsv(file)
  file.df$modules <- gsub(".tsv","",file)
  df.list[[z]] <- file.df
  names(df.list)[z] <- file
  z <- z + 1
}

### Create the df variables
list2env(df.list,envir = .GlobalEnv) ##this function extracts the list elements to the env
glob.df <- bind_rows(df.list)

grouped.df <- glob.df %>%
  group_by(modules) %>%
  summarize("max_N_sp" = max(N_sp),"mean_N_sp"=mean(N_sp),"median_N_sp"=median(N_sp),"max_N_eff_sp" = max(N_eff_sp),"mean_N_eff_sp" = mean(N_eff_sp),"median_N_eff_sp"=median(N_eff_sp),
            "max_N_reac" = max(N_react),"max_eff_reac" = max(N_eff_react),
            "max_eff_reac_pos" = max(N_eff_react_pos),"max_eff_reac_neg" = max(N_eff_react_neg))
grouped.df[1,1] <- "All"
grouped.df[2,1] <- "Early"
class(grouped.df)

late.df <- data.frame("modules" = "Late","max_N_sp" = 1,"mean_N_sp" = 1,"median_N_sp"= 1,"max_N_eff_sp" = 1 ,"mean_N_eff_sp" = 1,"median_N_eff_sp"= 1,
            "max_N_reac" = 1,"max_eff_reac" = 1,
            "max_eff_reac_pos" = 1,"max_eff_reac_neg" = 1)

final.df <- rbind(grouped.df,late.df)
final.df <- mutate(final.df, "n_communities" = c(nrow(stats_df_misosoup_Agarose.yaml.tsv),
nrow(stats_df_misosoup_Agarose_mod_48.yaml.tsv),6), .after = modules)


gt.table <- gt(as.data.frame(final.df))
gtsave(gt.table,"stats_summary.html")


### Now the analysis of the networks

##Add the family info

##table of aligments (Model-ESV connection)
# 1. Load genome alignments
aligned_genomes <- read_tsv("/home/ajf/Desktop/CNB/ecocoherence/metabolism/genome_aligment/matched_ESV_id_0.97.tsv", col_names = FALSE)
aligned_genomes <- select(aligned_genomes, c("X9", "X10"))
colnames(aligned_genomes) <- c("ASV", "ID")
aligned_genomes <- filter(aligned_genomes, ID != "*")
unique_ids <- unique(aligned_genomes$ID)
aligned_genomes_unique <- data.frame("ID" = unique_ids)
aligned_genomes$ID <- gsub("[A-Z]{2}_([A-Z]{3}_[0-9]{9})(\\.\\d+)?", "\\1", aligned_genomes$ID)
# 2 taxonomy table

fileTaxonomy <- file.path(work.dir, "..","..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
taxonomy$otu <- taxonomy$taxa_id
taxonomy <- subset(taxonomy, select = c('otu', 'Family', 'Genus', 'Order', 'Class'))



### Creation of net
net.all <- read_tsv("misosoup_Agarose/network-count_gal_min_sols_aggr.csv")


net.all$nodeB <- gsub("_model","",net.all$nodeB)
any(aligned_genomes$ID %in% net.all$nodeB)

net.all$ASV <- aligned_genomes$ASV[match(net.all$nodeB,aligned_genomes$ID)]
net.all$fam <- taxonomy$Family[match(net.all$ASV,taxonomy$otu)]



###Lets try making a heatmap

head(net.all)

net.all$nodeB <- paste(net.all$nodeB,net.all$fam)

net.all.wider <- pivot_wider(net.all,names_from = c(nodeB),values_from = weight,id_cols = nodeA)
# rownames(net.all.wider) <- net.all.wider$nodeB
# net.all.wider <- subset(net.all.wider,select = -c(nodeB))
net.all.wider[is.na(net.all.wider)] <- 0



heatmap.mat <- as.matrix(net.all.wider[,-1])

test.all <- ifelse(heatmap.mat > 0,1,0)
dist.jaccard.all <- vegdist(t(test.all),method = "jaccard")
pheatmap(test.all, fontsize = 6,fontsize_col = 6)
plot(hclust(dist.jaccard.all))


dist.all <- dist(t(heatmap.mat),method = "euclidian")
dendo.all <- hclust(dist.all)
plot(dendo.all)
# rownames(heatmap.mat) <- net.all.wider$nodeA
pheatmap(heatmap.mat,fontsize = 9)


############## For early modules
### Creation of net
net.early <- read_tsv("misosoup_Agarose_mod_48/network-count_gal_min_sols_aggr.csv")


net.early$nodeB <- gsub("_model","",net.early$nodeB)
any(aligned_genomes$ID %in% net.early$nodeB)

net.early$ASV <- aligned_genomes$ASV[match(net.early$nodeB,aligned_genomes$ID)]
net.early$fam <- taxonomy$Family[match(net.early$ASV,taxonomy$otu)]



###Lets try making a heatmap

head(net.early)

net.early$nodeB <- paste(net.early$nodeB,net.early$fam)

net.early.wider <- pivot_wider(net.early,names_from = c(nodeB),values_from = weight,id_cols = nodeA)
# rownames(net.early.wider) <- net.early.wider$nodeB
# net.early.wider <- subset(net.early.wider,select = -c(nodeB))
net.early.wider[is.na(net.early.wider)] <- 0



heatmap.mat.early <- as.matrix(net.early.wider[,-1])
rownames(heatmap.mat.early) <- net.early.wider$nodeA

test <- ifelse(heatmap.mat.early > 0,1,0)
dist.jaccard.early <- vegdist(t(test),method = "jaccard")
plot(hclust(dist.jaccard.early))

dist.early <- dist(t(heatmap.mat.early),method = "euclidian")
dendo.early <- hclust(dist.early)
plot(dendo.early)
pheatmap(heatmap.mat.early,fontsize = 6)
# rownames(heatmap.mat) <- net.all.wider$nodeA
# pheatmap(heatmap.mat,fontsize = 9)



