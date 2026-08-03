#########################################################
# This script make the cord plots of all the substrates #
#########################################################
rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(circlize)
library(this.path)

work_dir <- this.dir()

# fileTaxonomy <- "/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
fileTaxonomy <- file.path(work_dir,"..","data","marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy)
taxa.in <- read.table(fileTaxonomy, sep = ";")
colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
                       "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
taxonomy <- subset(taxa.in, select = c("taxa_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
rownames(taxonomy) <- taxonomy$taxa_id
taxonomy <- subset(taxonomy, select = -c(taxa_id))
taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))


substrates <- c("Agarose","Carrageenan","Alginate","AgaroseCarrageenan","AgaroseChitosan","AgaroseAlginate","Chitin")
eco_strategy_df <- data.frame("Module" = character(), "Strategy" = character())
for (substrate in substrates){
  setwd(this.dir())
  subs_df <- read_tsv(paste("modules_classified_internal",substrate,".tsv",sep  = ""))
  #format the module column as in the df 
  subs_df$Module <- gsub("mod_","",subs_df$Module)
  subs_df$Module <- paste(subs_df$Module,"_",substrate,sep="")
  #bind the dataframes
  eco_strategy_df <- rbind(eco_strategy_df,subs_df)
}


strategy_color_palette <- c(
  facilitation = "red",
  selection = "blue",
  transition = "yellow",
  generalist = "green",
  not_defined = "gray"
)


for (substrate in substrates) {
  functionink_path <- file.path(work_dir,"..","functionink")
  functionink_path <- normalizePath(functionink_path)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp",sep = ""))
  setwd(paste(functionink_path,substrate,"functionink_tmp",sep="/"))
  
  # List all matching files based on the pattern
  files <- list.files(
    path = ".",  # Directory to search in
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE  # Get the full file path
  )
  
  
  
  data <- read.table(files, sep = " ",header = TRUE)
  
  data2 <- separate(data, guild, sep = "\t" , into = c("ASV","Guild"))
  
  
  data2$ASV <- taxonomy$Family[match(data2$ASV,rownames(taxonomy))] 
  data2$Guild <- gsub("mod_","",data2$Guild)
  data2$Guild <- paste(data2$Guild,substrate,sep="_")
  ##adding to data2  the ecological strategy info
  data2$eco_strategy <- eco_strategy_df$Strategy[match(data2$Guild,eco_strategy_df$Module)]
  data2 <- select(data2,-Guild)
  data2$ASV <- gsub("_incertae_sedis","",data2$ASV)
  data2$ASV <- gsub("_incertae sedis","",data2$ASV)
  freq_table <- as.data.frame(table(data2))
  
  colnames(freq_table) <- c("from","to","value")
  freq_table_filtered <- filter(freq_table,value >=6) ### Keep modules with at least 6 members
  
  fig_path <- file.path(work_dir,"..","figures")
  fig_path <- normalizePath(fig_path)
  setwd(fig_path)
  dir.create("chord_plots",showWarnings = F)
  setwd("./chord_plots")
  png(paste("circular_plot_",substrate,".png",sep = ""),1000,1000)
  # Increase margins to provide more room for labels
  # par(oma = c(20, 20, 20, 20) + 0.1) # Bottom, left, top, right margins
  # Set text size for sector labels (e.g., "A", "B", "X", etc.)
  circos.par(track.height = 0.6, cell.padding = c(1,1,1,1), gap.degree = 1)
  chordDiagram(freq_table_filtered, annotationTrack = "grid",scale = T, grid.col = strategy_color_palette)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2.7] , CELL_META$sector.index, 
                facing = "inside", niceFacing = TRUE, adj = c(0.5, 0,5),cex = 1 , font = 2)
  }, bg.border = NA) # here set bg.border to NA is important
  
  # # Adjust text size for sector labels
  # circos.track(track.index = 1, panel.fun = function(x, y) {
  #   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
  #               facing = "clockwise", niceFacing = TRUE, cex = 2)
  # }
  
  
  
  legend("bottomleft",pch = 1, legend = "Number of member of each family that goes to each module",cex = 1)
  title(substrate , cex = 2)
  
  # Clear circos parameters after plotting
  circos.clear()
  dev.off()
  
}


