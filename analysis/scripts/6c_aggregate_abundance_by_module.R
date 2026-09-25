rm( list = ls())

library(readr)
library(tidyr)
library(this.path)

script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "6c_aggregate_abundance_by_module"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
substrates <- c("Agarose","AgaroseAlginate","Alginate","Chitin","AgaroseChitosan","AgaroseCarrageenan","Carrageenan")

for (substrate in substrates){
  ##load the otu table
  otu_table_dir <- file.path(script_dir,"..","results","1_filter_asv_table_and_build_sparcc_network",substrate)
  otu_table_dir <- normalizePath(otu_table_dir)
  setwd(otu_table_dir)

  otu_table <- read_tsv(paste("otu_table",substrate, sep = "_"))

  funink_table_dir <- file.path(script_dir,"..","results","5_filter_modules_by_size_and_plot_abundance",substrate)
  funink_table_dir <- normalizePath(funink_table_dir)
  setwd(funink_table_dir)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  dir_path <- "."
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
  
  setwd(results_dir)
  write_tsv(otu_table_grouped,paste("otu_table_",substrate,"_byModulesSize4.tsv",sep = ""))

}