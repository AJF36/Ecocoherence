#####This script filter all the ESV of the modules of one substrates and
#####filter the ids for moving the genomes to a folder for the substrate
#####create the metabolic models with the 

rm(list = ls())
library(readr)
library(tidyr)
library(dplyr)
library(this.path)
work_dir <- this.dir()
substrates <- c("Alginate","Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")


###Load the id-ESV table
matched_esv_file <- file.path(work_dir,"genome_aligment","matched_ESV_id_0.97.tsv")

# ESV_ID_table <- read_tsv("/home/ajf/Desktop/CNB/ecocoherence_metabolism/genome_aligment/matched_ESV_id_0.97.tsv",col_names = F)
ESV_ID_table <- read_tsv(matched_esv_file,col_names = F)
##Load the abundance table to know how much abundance do we explain 
count_table_file <- file.path(work_dir,"..","data","marine_particles_source_data")
count_table_file <- normalizePath(count_table_file)
# count_table <- read.table(file = "/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv" , sep = ",")
count_table <- read.table(file = file.path(count_table_file,"count_table.ESV.4R.csv") , sep = "," , header = T)
###Process it so it only haves ESV with matched entry and it has the correct format
ESV_ID_table <- ESV_ID_table %>%
  select(c("X9","X10")) %>% 
  filter(X10 != "*")
colnames(ESV_ID_table) <- c("ESV","ID_entry")
ESV_ID_table$ID_entry <- gsub(".*([A-Z]{3}_\\d+).*", "\\1", ESV_ID_table$ID_entry)


for (substrate in substrates) {
  ###Load the functionk table
  ###Read the module table
  subs_path <- file.path(work_dir,"..","functionink",substrate)
  subs_path <- normalizePath(subs_path)
  setwd(subs_path)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  dir_path <- "functionink_tmp/"
  file <- list.files(
    path = dir_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE
  )
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))###Filter the ASV of the modules
  # # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp",sep = ""))
  # file_name <- list.files(
  #   path = paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp",sep = ""),  # Set to the directory containing the files
  #   pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
  #   full.names = TRUE
  # )
  network_functionink_df <- read_tsv(file)
  network_functionink_df <- separate(network_functionink_df,col = guild,into = c("ASV","module"), sep = "\t")
  
  ###Filter by the ones that are has an access ID
  ESV_ID_table_f_substrate <- filter(ESV_ID_table,ESV %in% network_functionink_df$ASV)
 
  ###Continue creating the metabolic models
  # setwd("/home/ajf/Desktop/CNB/ecocoherence_metabolism/carveme_smetana")
  setwd(work_dir)
  system(paste("mkdir","carveme_smetana",sep = " "))
  carveme_dir <- file.path(work_dir,"carveme_smetana")
  setwd(carveme_dir)
  system(paste("mkdir ",substrate,sep=""))
  for (id in ESV_ID_table_f_substrate$ID_entry){
    # setwd("/home/ajf/Desktop/globdb_database/globdb_r220_protein_faa")
    glod_db_dir <- normalizePath(file.path(work_dir,"..","data","globdb_r220_protein_faa"))
    setwd(glod_db_dir)
    print(id)
    ###Copy the genomes in the substrate directory
    # system(paste("ls | grep ",id," | xargs cp -t ../../CNB/ecocoherence_metabolism/carveme_smetana/",substrate,sep = ""))
    system(paste("ls | grep ",id," | xargs cp -t ","../../metabolism/carveme_smetana/",substrate,sep = ""))
    
    
  }
  carve_sme_dir <- file.path(work_dir,"carveme_smetana")
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_metabolism/carveme_smetana/",substrate,sep=""))
  setwd(paste(carve_sme_dir,substrate,sep="/"))
  
  system("gunzip *gz ") ##Descomprimimos
  system("ls > number_id.tsv")
  ########################################## How much abundance do we lose
  proportion_module_explained <- data.frame(module =character(), abundance_explained = numeric())
  number_id <- read_tsv("number_id.tsv",col_names = F)
  number_id$X1 <- gsub(".faa","",number_id$X1)
  ESV_in_db <- ESV_ID_table_f_substrate$ESV[match(number_id$X1,ESV_ID_table_f_substrate$ID_entry)]
  count_table_substrate <- count_table[,grepl(paste("Beads_",substrate,"_",sep = ""),colnames(count_table))] ##filter the abundance_table
  count_table_substrate$ASV <- rownames(count_table_substrate)
  # guild= "mod_2" ###DELETE THIS
  for (guild in unique(network_functionink_df$module)){
    network_functionink_df_module <- filter(network_functionink_df, module == guild)
    count_table_substrate_filtered <- filter(count_table_substrate,ASV %in% network_functionink_df_module$ASV)
    length(colnames(count_table_substrate_filtered))
    abundance_relative <- sweep(count_table_substrate_filtered[-length(colnames(count_table_substrate_filtered))], 2, colSums(count_table_substrate_filtered[,-length(colnames(count_table_substrate_filtered))]), FUN = "/")
    colSums(abundance_relative)
    
    abundance_relative$ESV <- rownames(abundance_relative)
    abundance_relative_f <- filter(abundance_relative,ESV %in% ESV_in_db)
    abundance_relative_f <-  abundance_relative_f[, colSums(is.na(abundance_relative_f)) == 0]
    mean(colSums(abundance_relative_f[,-length(colnames(abundance_relative_f))]))
    

    temporal_df <- data_frame(module = guild, abundance_explained = mean(colSums(abundance_relative_f[,-length(colnames(abundance_relative_f))])))
    proportion_module_explained<- rbind(temporal_df,proportion_module_explained) ###Nas por quee???
    
  }
  
  write_tsv(proportion_module_explained, "proportion_abundance_modules_explained.tsv")
  
  

}





