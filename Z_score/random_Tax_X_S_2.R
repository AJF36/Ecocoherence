###########################################################################################################
#This script makes the table with the number of members of each family in each module and then calculates #
# 50 times and then calculates the mean S (shannon entropy) and X (exp^S)                                 #
###########################################################################################################



rm(list = ls())
library(readr)
library(phyloseq)
library(dplyr)
library(grid)
library(gridExtra)
library(this.path)
work_dir <- this.dir()
substrates <- c("Agarose")

# Cargar la tabla de taxonomía
# fileTaxonomy <- "/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
# taxa.in <- read.table(fileTaxonomy, sep = ";")
# colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
#                        "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
# taxonomy <- subset(taxa.in, select = c("taxa_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
# rownames(taxonomy) <- taxonomy$taxa_id
# taxonomy <- subset(taxonomy, select = -c(taxa_id))
# taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))
fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))

# Función para permutar todas las columnas de la tabla de taxonomía
permute_taxonomy <- function(taxonomy) {
  taxonomy_permuted <- taxonomy %>%
    mutate(across(everything(), ~ sample(.)))
  return(taxonomy_permuted)
}

# Crear una lista vacía para almacenar los resultados
results_list <- list()

# Bucle para realizar el proceso 50 veces
for (iteration in 1:50) {
  # Crear una sublista para almacenar los resultados de cada sustrato en esta permutación
  perm_results <- list()
  # Iterar sobre cada sustrato
  for (substrate in substrates) {
    substrate <- "Agarose"
    # Cargar la tabla de FunctionInk
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
    # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
    # dir_path <- "functionink_tmp/"
    # file <- list.files(
    #   path = dir_path,
    #   pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    #   full.names = TRUE
    # )
    functionink <- read.table(file, sep = "\t", header = TRUE)
    functionink$ESV <- rownames(functionink)
    colnames(functionink) <- c("guild", "ESV")
    
    # Filtrar la tabla de taxonomía para que contenga solo los OTUs presentes en FunctionInk
    tax.pseq.families <- as.data.frame(subset(taxonomy, rownames(taxonomy) %in% functionink$ESV, select = "Family"))
    tax.pseq.families$ESV <- rownames(tax.pseq.families)
    tax.pseq.families
    # Permutar la tabla de taxonomía filtrada
    taxonomy_permuted <- permute_taxonomy(tax.pseq.families)
    taxonomy_permuted

    # Fusionar FunctionInk con la tabla de familias permutada
    merged_table <- merge(functionink, taxonomy_permuted, by = "ESV")
    merged_table$guild <- as.numeric(gsub("mod_", "", merged_table$guild))
    merged_table <- arrange(merged_table, guild)
    merged_table$Family <- gsub("_incertae_sedis|_Incertae Sedis XI", "", merged_table$Family)
    
    families <- unique(merged_table$Family)
    heatmap_groups <- unique(merged_table$guild) # Grupos de guilds
    
    # Crear la matriz guild_matrix
    guild_matrix <- matrix(0, nrow = length(families), ncol = length(heatmap_groups))
    rownames(guild_matrix) <- families
    colnames(guild_matrix) <- heatmap_groups
    
    # Llenar la matriz con el conteo de OTUs por familia y guild
    for (fam in families) {
      for (group in heatmap_groups) {
        guild_matrix[fam, as.character(group)] <- sum(merged_table$Family == fam & merged_table$guild == group)
      }
    }
    
    # Crear DataFrame con el conteo y clasificaciones de estrategias
    guild_df <- as.data.frame(guild_matrix)
    guild_df$n_members <- rowSums(guild_df)
    
    # Filtrar y ordenar la tabla
    guild_df <- guild_df %>%
      arrange(desc(n_members)) %>%
      filter(n_members >= 4)
    guild_df$Family <- row.names(guild_df)
    guild_df
    # Guardar el resultado en la sublista de la permutación actual
    perm_results[[substrate]] <- guild_df
  }
  
  # Almacenar la sublista de la permutación actual en la lista principal de resultados
  results_list[[paste0("perm_", iteration)]] <- perm_results
}

# Ahora `results_list` contiene los DataFrames de cada permutación y sustrato.
results_list[[1]]

# Ahora procederemos a calcular el estadístico X para cada sustrato y cada familia como antes.

# Crear un dataframe vacío para almacenar los resultados
results_summary <- data.frame(substrate = character(),
                              family = character(),
                              mean_X = numeric(),
                              sd_X = numeric(),
                              stringsAsFactors = FALSE)

# Iterar sobre cada sustrato en perm_results (suponiendo que perm_results tiene las tablas de cada sustrato)
for (i in 1:length(results_list)) {
  
  # Tomar el dataframe correspondiente al sustrato (ahora es un dataframe de 50 tablas por sustrato)
  substrate_df_list <- results_list[[i]]
  z = 1
  # Para cada tabla dentro de este sustrato (todas las 50 tablas)
  for (j in 1:length(substrate_df_list)) {
    # Extraer el dataframe para la j-ésima tabla dentro de este sustrato
    substrate_df <- substrate_df_list[[j]]
    
    
    
    substrate_df <- substrate_df %>%
      mutate(across(
        .cols = matches("^\\d"),
        .fns = ~ ifelse(n_members != 0, .x / n_members, NA),
        .names = "{.col}" # Keeps original column names
      ))
    
  
    
    
    # Seleccionar las columnas numéricas (equivalente a las columnas con módulos)
    numeric_columns <- substrate_df %>% select(matches("^\\d"))
    
    # Asegurarnos de que las columnas no tengan valores nulos en 'n_members'
    # substrate_df[names(numeric_columns)] <- lapply(names(numeric_columns), function(col) {
    #   ifelse(substrate_df$n_members != 0, substrate_df[[col]] / substrate_df$n_members, NA)
    # })
    
    
    
    
    
    
    # Calcular el valor de S para esta tabla
    substrate_df$S <- -rowSums(numeric_columns * log(numeric_columns), na.rm = TRUE)
    # Calcular el valor de X para esta tabla
    substrate_df$X <- exp(substrate_df$S)
    
    # Agrupar por familia
    family_groups <- unique(substrate_df$Family)
    
    # Para cada familia, calcular X
    for (family in family_groups) {
      family_data <- substrate_df[substrate_df$Family == family, ]
    family_data
      # Agregar el valor de X de esta familia en esta tabla
      for (x_value in family_data$X) {
        results_summary <- rbind(results_summary, data.frame(substrate = substrates[z],
                                                             family = family,
                                                             X = x_value))
        
      
      }
    
    }
    z <- z + 1
  }
}

# Ahora calculamos la media y desviación típica de X por familia y por sustrato
final_summary <- results_summary %>%
  group_by(substrate, family) %>%
  summarise(mean_X = mean(X, na.rm = TRUE),
            sd_X = sd(X, na.rm = TRUE))

# setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc/group_analysis")
setwd(work_dir)
# write_tsv(final_summary,"X_mean_dv_randomized.tsv")
write_tsv(final_summary,"X_mean_dv_randomized.tsv")

# Ver los resultados
View(final_summary)




