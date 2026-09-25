#########################################################
# This script make the cord plots of all the substrates #
#########################################################
rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(circlize)
library(this.path)

script_dir <- this.dir()
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "7d_plot_chord_diagrams_family_module"), mustWork = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
producer_dir <- normalizePath(file.path(script_dir, "..", "results", "6i_classify_modules_by_ecological_strategy_combined"), mustWork = FALSE)

# fileTaxonomy <- "/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
fileTaxonomy <- file.path(script_dir,"..","data","marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy)
taxa.in <- read.table(fileTaxonomy, sep = ";")
colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
                       "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
taxonomy <- subset(taxa.in, select = c("taxa_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
rownames(taxonomy) <- taxonomy$taxa_id
taxonomy <- subset(taxonomy, select = -c(taxa_id))
taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))





eco_strategy_df <- data.frame("Module" = character(), "Strategy" = character())

substrates <- c("Agarose","Carrageenan","Alginate","AgaroseCarrageenan","AgaroseChitosan","AgaroseAlginate","Chitin")
for (substrate in substrates){
  setwd(producer_dir)
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


substrates <- c("Agarose","Carrageenan","Alginate","AgaroseCarrageenan","AgaroseChitosan","AgaroseAlginate","Chitin")

for (substrate in substrates) {
  functionink_path <- file.path(script_dir,"..","results","5_filter_modules_by_size_and_plot_abundance")
  functionink_path <- normalizePath(functionink_path)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp",sep = ""))
  setwd(paste(functionink_path,substrate,sep="/"))
  
  # List all matching files based on the pattern
  files <- list.files(
    path = ".",  # Directory to search in
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE  # Get the full file path
  )
  
  
  
  data <- read.table(files, sep = " ",header = TRUE)
  
  data_test<- separate(data, guild, sep = "\t" , into = c("ASV","Guild"))
  
  
  data_test$ASV <- taxonomy$Family[match(data_test$ASV,rownames(taxonomy))] 
  data_test$Guild <- gsub("mod_","",data_test$Guild)
  data_test$Guild <- paste(data_test$Guild,substrate,sep="_")
  # ##adding to data2  the ecological strategy info
  # data2$eco_strategy <- eco_strategy_df$Strategy[match(data2$Guild,eco_strategy_df$Module)]
  # data2 <- select(data2,-Guild)
  data_test$ASV <- gsub("_incertae_sedis","",data_test$ASV)
  data_test$ASV <- gsub("_incertae sedis","",data_test$ASV)
  freq_table <- as.data.frame(table(data_test))
  
  colnames(freq_table) <- c("from","to","value")
  freq_table_filtered <- filter(freq_table,value >0) ### Keep modules with at least 6 members
  
  ### Crear el vector de coloración

  vec_for_color <- c(1:length(unique(freq_table_filtered$to)))
  names(vec_for_color) <- unique(freq_table_filtered$to)
  
  ### Now i need to change the value of each of the elements of the vector for the colour of the strategy
  vec_for_color

  for (mod in names(vec_for_color)) {
    
    mod_strategy <- eco_strategy_df$Strategy[match(mod,eco_strategy_df$Module)]
    vec_for_color[mod] <- strategy_color_palette[mod_strategy]

  }

  filtered_vec_for_color <- na.omit(vec_for_color)  ## Drop the modules that does not have ecological strategy assigned 
                                                    ## due to its size, we can also use it to filter the df for the 
                                                    ## cord plot
  filtered_vec_for_color

  freq_table_filtered <- filter(freq_table_filtered, to %in% names(filtered_vec_for_color))

  
  ### What families should we drop for the figure. Lets say at least 5 members for the moment
  suma_df <- freq_table_filtered %>%
    group_by(from) %>%
    summarise(total_members = sum(value)) %>%
    filter(total_members >= 5)


  ## Now we filter the families of suma_df

  freq_table_filtered <- filter(freq_table_filtered, from %in% suma_df$from)

   ## Create a vector of colour for the paths
  path_colour_df <- freq_table_filtered 

  path_colour_df$color <- filtered_vec_for_color[match(path_colour_df$to,names(filtered_vec_for_color))]



  setwd(figures_dir)

  # png(paste("circular_plot_modified_",substrate,".png",sep = ""),1000,1000)
  # #Increase margins to provide more room for labels
  # par(oma = c(20, 20, 20, 20) + 0.1) # Bottom, left, top, right margins
  # #Set text size for sector labels (e.g., "A", "B", "X", etc.)
  # circos.par(track.height = 0.6, cell.padding = c(1,1,1,1), gap.degree = 1)

png(paste("circular_plot_modified_", substrate, ".png", sep = ""),
    width = 2000, height = 2000, res = 300)  # High resolution, square canvas

# Remove extra margins and set equal aspect ratio
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = TRUE)

# Adjust circlize parameters to make better use of space
circos.clear()
circos.par(     # rotate if needed
       # reduce gaps between sectors
           track.margin = c(0, 0),
           cell.padding = c(0, 0, 0, 0),
           canvas.xlim = c(-1.4, 1.5),  # enlarge plot area
           canvas.ylim = c(-1.4, 1.5))  # enlarge plot area



  chordDiagram(freq_table_filtered, annotationTrack = "grid",scale = F, grid.col = filtered_vec_for_color,
   col = path_colour_df$color, big.gap = 30)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ycenter + 7.5, CELL_META$sector.index , 
                facing = "reverse.clockwise", niceFacing = T, adj = c(0.2, 0,5),cex = 0.7, font = 2)
  }, bg.border = NA) # here set bg.border to NA is important
  
  # # Adjust text size for sector labels
  # circos.track(track.index = 1, panel.fun = function(x, y) {
  #   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
  #               facing = "clockwise", niceFacing = TRUE, cex = 2)
  # }
  

  # legend("bottomleft",pch = 1, legend = "Number of member of each family that goes to each module",cex = 1)
  legend("topright", legend = names(strategy_color_palette), fill = strategy_color_palette,
       title = "Ecological strategies", inset = c(-0.00001, 0.0), xpd = TRUE , cex = 1 , text.font = 2 ,text.width = 0.4,
      bty = "n")
  # title(, cex = 3)
  
  # Clear circos parameters after plotting
  circos.clear()
  dev.off()
  
}

