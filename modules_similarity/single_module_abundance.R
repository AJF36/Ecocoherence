rm( list = ls())

library(readr)
library(tidyr)
library(this.path)

work_dir <- this.dir()
substrates <- c("Agarose","AgaroseAlginate","Alginate","Chitin","AgaroseChitosan","AgaroseCarrageenan","Carrageenan")

for (substrate in substrates){
  ##load the otu table
  otu_table_dir <- file.path(work_dir,"..","network_creation_sparcc",substrate)
  otu_table_dir <- normalizePath(otu_table_dir)
  setwd(otu_table_dir)
  
  otu_table <- read_tsv(paste("otu_table",substrate, sep = "_"))
  
  funink_table_dir <- file.path(work_dir,"..","functionink",substrate)
  funink_table_dir <- normalizePath(funink_table_dir)
  setwd(funink_table_dir)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  dir_path <- "functionink_tmp/"
  file <- list.files(
    path = dir_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE
  )
  funink_table <- read.table(file , sep = "\t" , header = T)
  funink_table$ESV <- rownames(funink_table)
  colnames(funink_table) <- c("guild","ESV")
  
  
  ###Filter the otu table by the ASV in the funinkl table
  colnames(otu_table)[1] <- "ESV"
  otu_table_f <- filter(otu_table, ESV %in% funink_table$ESV)
  
  ###Add the module info to the otu table
  otu_table_f$Module <- funink_table$guild[match(otu_table_f$ESV,funink_table$ESV)]
  
  otu_table_grouped <- otu_table_f %>%
    group_by(Module) %>% 
    summarise(across(where(is.numeric),sum,na.rm = T))
  
  setwd(work_dir)
  write_tsv(otu_table_grouped,paste("otu_table_",substrate,"_byModulesSize4.tsv",sep = ""))

}