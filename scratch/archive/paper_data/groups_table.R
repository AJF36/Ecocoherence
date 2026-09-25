#This script make the table with the number of members of each family in each module and each ecological strategy
#For this it takes a tax table and the functionink table which the heatmap code of alberto produces and the module -ecological strategy of the file modules_classified
#from the cript modules_to_ecologicalstrategy.R
rm(list=ls())
library(readr)
library(phyloseq)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(gt)
library(webshot)
library(this.path)

setwd("/home/ajf/Desktop/CNB/ecocoherence/paper_data")

tax.pseq <- readRDS("taxonomy.rds")
tax.pseq.families <- subset(tax.pseq, select = "Family")
tax.pseq.families <- as.data.frame(tax.pseq.families)
tax.pseq.families$ESV <- rownames(tax.pseq.families)


  
# functionink <- read.table("Total_Max_Density_ S1841.tsv", sep = "\t",header = TRUE)
functionink <- read_tsv("Internal_Max_Density_ S1841.tsv" , skip = 9, col_names = F)
# functionink$ESV <- rownames(functionink)
colnames(functionink) <- c("ESV","guild")
any(functionink$ESV %in% tax.pseq.families$ESV)
  
  
  #Table with the OTU family and group information
  merged_table <- merge(functionink, tax.pseq.families, by = "ESV")
  
  merged_table$guild <- as.numeric(gsub("mod_", "", merged_table$guild))
  
  # Verificar que la columna guild ahora solo contiene números
  str(merged_table)
  # Mostrar las primeras filas para verificar el cambio
  merged_table <- arrange(merged_table,guild)
  
  # View(merged_table)
  merged_table$Family <- gsub("_incertae_sedis","",merged_table$Family)
  merged_table$Family <- gsub("_Incertae Sedis XI","",merged_table$Family)
  merged_table <- drop_na(merged_table)
  # ####This chunck try to automatize the selection of the modules
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep= ""))
  # modules_classified <- read_tsv(paste("modules_classified_",substrate,".tsv",sep=""))
  # modules_classified$Module <-as.numeric(gsub("mod_", "", modules_classified$Module))
  
  
  
  #Now the table with the information of the heatmap
  ##Uncomment for doing only with the heatmap groups
  # heatmap_groups <- modules_classified$Module
  heatmap_groups <- unique(merged_table$guild)
  families <- unique(merged_table$Family)
  num_members <- vector()
  for (i in 1:length(families)) {
    print(families[i])
    members <- sum(merged_table$Family == families[i])
    sum(merged_table$Family == families[i])
    merged_table$Family
    members
    num_members[i] <- members
  }
  
  
  
  #Now im going to create a matrix with the counts of each guild for each family
  
  guild_matrix <- matrix(0,nrow= length(families), ncol= length(heatmap_groups))
  heatmap_groups
  rownames(guild_matrix) <- families
  colnames(guild_matrix) <- heatmap_groups
  
  # Bucle para contar los OTUs por cada combinación de familia y guild
  for (i in 1:length(families)) {
    for (j in 1:length(heatmap_groups)) {
      # Sumar los OTUs que pertenecen a la familia y al guild específicos
      guild_matrix[i, j] <- sum(merged_table$Family == families[i] & 
                                  merged_table$guild == heatmap_groups[j])
    }
  }
  
  
  # Mostrar los resultados
  
  
  guild_df <- as.data.frame(guild_matrix)
  guild_df$n_members <- num_members
  numeric_columns <- guild_df %>% select(matches("^\\d"))
  
  # #classify the guilds
  sums_df <- guild_df
  # sums_df$generalist <- 0
  # sums_df$facilitation <- 0
  # sums_df$selection <- 0
  # sums_df$transition <- 0
  # sums_df$attachment <- 0
  # 
  # generalist <- modules_classified$Module[modules_classified$Strategy == "generalist"]
  # facilitation <- modules_classified$Module[modules_classified$Strategy == "facilitation"]
  # selection <- modules_classified$Module[modules_classified$Strategy == "selection"]
  # attachment <-  modules_classified$Module[modules_classified$Strategy == "attachment"]
  # transition <- modules_classified$Module[modules_classified$Strategy == "transition"]
  
  
  # View(sums_df)
  # for (i in 1:nrow(sums_df)) {
  #   # Para cada fila, sumar los valores de las columnas de guilds que pertenecen a cada categoría
  #   sums_df$generalist[i] <- sum(sums_df[i, as.character(generalist)], na.rm = TRUE)
  #   sums_df$facilitation[i] <- sum(sums_df[i, as.character(facilitation)], na.rm = TRUE)
  #   sums_df$selection[i] <- sum(sums_df[i, as.character(selection)], na.rm = TRUE)
  #   sums_df$transition[i] <- sum(sums_df[i, as.character(transition)], na.rm = TRUE)
  #   sums_df$attachment[i] <- sum(sums_df[i, as.character(attachment)], na.rm = TRUE)
  # }
  
  
  sums_df_ordered <- arrange(sums_df, desc(n_members) )
  sums_df_ordered <- sums_df_ordered %>% filter(n_members >= 4)
  sums_df_ordered$Family <- row.names(sums_df_ordered)
  # View(sums_df_ordered)
  setwd(work_dir)
  write_tsv(sums_df_ordered,paste("group_analysis_",substrate,".tsv", sep = ""))
}


familys_of_interest <- c("Helicobacteraceae","Nannocystaceae","Campylobacteraceae","Porphyromonadaceae","Colwelliaceae","Alteromonadaceae","Rhodobacteraceae"
                         ,"Flavobacteriaceae","Oceanospirillaceae","Syntrophaceae","Anaerolineaceae","Desulfobacteraceae","Psychromonadaceae")
gt_table_df <- subset(sums_df_ordered[colnames(numeric_columns)])
row_sums <- rowSums(gt_table_df)
gt_table_df <- sweep(gt_table_df,1,row_sums,"/")
gt_table_df <- sweep(gt_table_df,1,100,"*")
gt_table_df <- round(gt_table_df,2)

gt_table_df$Family <- rownames(gt_table_df)
gt_table_df <- filter(gt_table_df,Family %in% familys_of_interest)
gt_table_df_f <- gt_table_df[,!colSums(gt_table_df[,-length(colnames(gt_table_df))]) == 0]


gt_table <- gt(gt_table_df_f)
gt_table <- 
  gt_table %>%
  tab_spanner(label = "Modules",
              columns = colnames(numeric_columns)[colnames(numeric_columns) %in% colnames(gt_table_df_f)]
  ) %>% tab_header(
    title = md("**Percentage of members of each family in the different modules**"),
    subtitle = substrate
  ) %>% data_color(columns = colnames(gt_table_df_f)[-length(colnames(gt_table_df_f))],
                   fn = function(x) ifelse(x >= 70 , "blue","white"))



setwd(this.dir())
gtsave(gt_table,path = this.dir(),filename = paste(substrate,"modules_table.html",sep = "_"))

# Agregar un título