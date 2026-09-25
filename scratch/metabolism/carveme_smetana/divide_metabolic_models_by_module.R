rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(this.path)


work_dir <- this.dir()
# substrates <- c("Alginate","Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")
substrates <- "AgaroseChitosan"
number_modules_id <- c()
# matched_esv <- read_tsv("/home/ajf/Desktop/CNB/ecocoherence_metabolism/genome_aligment/matched_ESV_id_0.97.tsv",col_names = F)
matched_esv_file <- file.path(work_dir,"..","genome_aligment","matched_ESV_id_0.97.tsv")
matched_esv_file <- normalizePath(matched_esv_file)
matched_esv <- read_tsv(matched_esv_file,col_names = F)
matched_esv <- select(matched_esv, c("X9","X10"))
colnames(matched_esv) <- c("ASV","ID")
matched_esv <- filter(matched_esv,ID != "*")
matched_esv$ID <- gsub(".*([A-Z]{3}_\\d+).*", "\\1", matched_esv$ID)
for (substrate in substrates) {
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))###Filter the ASV of the modules
  # # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp",sep = ""))
  # file_name <- list.files(
  #   path = paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp",sep = ""),  # Set to the directory containing the files
  #   pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
  #   full.names = TRUE
  # )
  # setwd(this.dir())
  
  print(this.dir())
  subs_path <- file.path(this.dir(),"..","..","functionink",substrate)
  subs_path <- normalizePath(subs_path)
  setwd(subs_path)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  dir_path <- "functionink_tmp/"
  file <- list.files(
    path = dir_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE
  )
  network_functionink_df <- read_tsv(file)
  network_functionink_df <- separate(network_functionink_df,col = guild,into = c("ASV","module"), sep = "\t")
  network_functionink_df_f <- filter(network_functionink_df,ASV %in% matched_esv$ASV)
  
  asv_module_id_df <- inner_join(network_functionink_df_f,matched_esv, by = "ASV")
  
  list_modules <- unique(asv_module_id_df$module)
  
#   for (guild in list_modules) {
#     setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_metabolism/carveme_smetana/",substrate,"/metabolic_models",sep = ""))
#     system(paste("mkdir ",guild,sep = ""))
#     module_df <- filter(asv_module_id_df, module == guild)
#     for (id in module_df$ID) {
#       system(paste("cp ", id,"* ./",guild,sep=""))
#     }
#     
#   }
# }
  

  library(dplyr)
  
  # Definir ruta base
  # base_path <- "/home/ajf/Desktop/CNB/ecocoherence_metabolism/carveme_smetana/"
  base_path <- file.path(this.dir())
  
  # Iterar sobre módulos
  for (guild in list_modules) {
    
    # Definir directorios
    work_dir <- file.path(base_path, substrate, "metabolic_models")
    module_dir <- file.path(work_dir, guild)
  
    # Crear directorio si no existe
    dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Filtrar datos del módulo
    module_df <- filter(asv_module_id_df, module == guild)
    
    # Iterar sobre cada ID en el módulo
    for (id in module_df$ID) {
      file.path(this.dir(),substrate,"metabolic_models")
      # Buscar archivos que comiencen con la ID
      list.files(path = file.path(this.dir(),substrate,"metabolic_models"), pattern = paste0("^", id, ".xml"), full.names = TRUE)
      files_to_copy <- list.files(path = file.path(this.dir(),substrate,"metabolic_models"), pattern = paste0("^", id, "_model.xml"), full.names = TRUE)
      # Copiar archivos al directorio del módulo
      file.copy(files_to_copy, module_dir, overwrite = TRUE)
    }
  
    }
  
  for (guild in list_modules){
    work_dir <- file.path(base_path, substrate, "metabolic_models")
    module_dir <- file.path(work_dir, guild)
    module_dir
    length(list.files(module_dir,pattern = "*.xml"))
     x <- length(list.files(module_dir,pattern = "*.xml"))
     number_modules_id <- c(x,number_modules_id)
  }
}  

hist(number_modules_id)
View(table(number_modules_id))
number_modules_id_f <- number_modules_id[number_modules_id > 30]
table(number_modules_id_f)

hist(number_modules_id_f)