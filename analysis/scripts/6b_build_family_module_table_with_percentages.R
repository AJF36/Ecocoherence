#This script make the table with the number of members of each family in each module and each ecological strategy
#For this it takes a tax table and the functionink table which the heatmap code of alberto produces and the module -ecological strategy of the file modules_classified
#from the cript modules_to_ecologicalstrategy.R
#03/10/2025 Modified to know the % of members of a certain family that belong to modules >= 4
rm(list=ls())
library(readr)
library(phyloseq)
library(dplyr)
library(grid)
library(gridExtra)
library(gt)
# library(webshot)
library(this.path)
library(gt)
library(gtExtras)


script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "6b_build_family_module_table_with_percentages"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "6b_build_family_module_table_with_percentages"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")

###Load the tax table
# fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
fileTaxonomy <- file.path(script_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
tax.pseq = tax_table(as.matrix(taxonomy))
# fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
# taxa.in=read.table(fileTaxonomy,sep=";") 
# colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
#                     "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
# taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
# rownames(taxonomy)=taxonomy$taxa_id
# taxonomy=subset(taxonomy,select=-c(taxa_id))
# tax.pseq = tax_table(as.matrix(taxonomy))
tax.pseq.families <- subset(tax.pseq, select = "Family")
tax.pseq.families <- as.data.frame(tax.pseq.families)
tax.pseq.families$ESV <- rownames(tax.pseq.families)

list_ratios_families <- list() ## To store the ratio families keeped/families eliminated of the analysis

for (substrate in substrates) {
  #Functionink table
  
  raw_partition_path <- file.path(script_dir,"..","results","4a_detect_modules_functionink",substrate,"functionink_tmp")
  raw_partition_path <- normalizePath(raw_partition_path)
  guild_filtered_path <- file.path(script_dir,"..","results","5_filter_modules_by_size_and_plot_abundance",substrate)
  guild_filtered_path <- normalizePath(guild_filtered_path)
  ### First we need the complete functionink file
  file <- list.files(
    path = raw_partition_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv$"),
    full.names = TRUE
  )

  file_functionink_filtered <- list.files(
    path = guild_filtered_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE
  )
  file_functionink_filtered
  
  # functionink <- read.table(file, sep = "\t",header = TRUE)
  functionink <- as.data.frame(read_tsv(file,skip = 9,col_names = F))

  ###Extract the modules with at least 4 members to filter later the df
  functionink_filtered <- read.table(file_functionink_filtered, sep = "\t",header = TRUE)
  big_modules <- unique(functionink_filtered$guild)
  
  # functionink$ESV <- rownames(functionink)
  colnames(functionink) <- c("ESV","guild")
  rownames(functionink) <- functionink$ESV

  
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
  sums_df_ordered_filtered <- sums_df_ordered %>% filter(n_members >= 4) ### Filter families that have at least 4 members

  sums_df_ordered$Family <- row.names(sums_df_ordered)


###Extract the n_members in modules that have at least 4 members
  sums_df_ordered_big_modules <- sums_df_ordered_filtered %>% 
    select(gsub("mod_","",big_modules))
  sums_df_ordered$Family <- row.names(sums_df_ordered)
  n_members_in_big_modules <- rowSums(sums_df_ordered_big_modules)

###Add the column of the %
sums_df_ordered_filtered$Family <- rownames(sums_df_ordered_filtered)
sums_df_ordered_filtered$n_members_big_modules <- n_members_in_big_modules ### Number of members in big_modules
sums_df_ordered_filtered <- sums_df_ordered_filtered %>%
  mutate("percentage_members_in_big_modules" = round((n_members_in_big_modules/n_members)* 100, 2)) %>%
  relocate(n_members, .before = colnames(sums_df_ordered_filtered)[1]) %>%
  relocate(percentage_members_in_big_modules, .after = n_members) %>%
  relocate(Family, .before = n_members ) %>%
  relocate(n_members_big_modules, .before = percentage_members_in_big_modules)

### This is a change needed for the code for the z_score calculation code to work
colnames(sums_df_ordered_filtered)[1:4] <- c("Family","total_members","n_members","percentage_members")
colnames(sums_df_ordered_filtered)
big_modules_subed <- gsub("mod_","",big_modules)
big_modules_df_final <- select(sums_df_ordered_filtered,c("Family","total_members","n_members","percentage_members",big_modules_subed))  ### Filter modules
big_modules_df_final <- filter(big_modules_df_final, n_members >= 4)
  
### Compute the difference of families keeped and eliminated
  ratio_families_substrate <- nrow(big_modules_df_final)/nrow(guild_df)
  list_ratios_families[[substrate]] <- ratio_families_substrate * 100
  # View(sums_df_ordered)
  setwd(results_dir)
  write_tsv(big_modules_df_final,paste("group_analysis_",substrate,".tsv", sep = ""))
}

### Code for making the table figure
# familys_of_interest <- c("Helicobacteraceae","Nannocystaceae","Campylobacteraceae","Porphyromonadaceae","Colwelliaceae","Alteromonadaceae","Rhodobacteraceae"
#                          ,"Flavobacteriaceae","Oceanospirillaceae","Syntrophaceae","Anaerolineaceae","Desulfobacteraceae","Psychromonadaceae")
# gt_table_df <- subset(sums_df_ordered[colnames(numeric_columns)])
# row_sums <- rowSums(gt_table_df)
# gt_table_df <- sweep(gt_table_df,1,row_sums,"/")
# gt_table_df <- sweep(gt_table_df,1,100,"*")
# gt_table_df <- round(gt_table_df,2)

# gt_table_df$Family <- rownames(gt_table_df)
# gt_table_df <- filter(gt_table_df,Family %in% familys_of_interest)
# gt_table_df_f <- gt_table_df[,!colSums(gt_table_df[,-length(colnames(gt_table_df))]) == 0]


# gt_table <- gt(gt_table_df_f)
# gt_table <- 
#   gt_table %>%
#   tab_spanner(label = "Modules",
#               columns = colnames(numeric_columns)[colnames(numeric_columns) %in% colnames(gt_table_df_f)]
#               ) %>% tab_header(
#                 title = md("**Percentage of members of each family in the different modules**"),
#                 subtitle = substrate
#               ) %>% data_color(columns = colnames(gt_table_df_f)[-length(colnames(gt_table_df_f))],
#                                fn = function(x) ifelse(x >= 70 , "blue","white"))
  


# setwd(this.dir())
# gtsave(gt_table,path = this.dir(),filename = paste(substrate,"modules_table.html",sep = "_"))
setwd(figures_dir)
ratio_df <- as.data.frame(list_ratios_families)
rownames(ratio_df) <- "ratio_families_used"
library(webshot2)
gt_table <-gt(ratio_df)
gt_table <- gt(ratio_df) %>%
  # Make table and container expand fully
  tab_options(
    table.width = pct(100),
    container.width = pct(100),
    container.overflow.x = F,
    table.align = "center",
    data_row.padding = px(3)
  )

gtsave(gt_table,"ratio_df.html")


df_to_calculate_percentages <- select(big_modules_df_final, matches("^\\d"))
View(df_to_calculate_percentages)


df_normalized <- df_to_calculate_percentages %>%
  mutate(across(
    matches("^\\d+$"),  # columnas cuyos nombres son solo números
    ~ .x / rowSums(across(matches("^\\d+$"))),  # divide cada valor entre la suma de su fila
    .names = "{.col}"   # mantiene los mismos nombres
  ))

max_value <- apply(df_normalized, 1, max)



df_for_figure_table <- select(big_modules_df_final, c("Family","total_members","n_members","percentage_members"))
df_for_figure_table$max_value <- round(max_value * 100,2)
# row.names(df_for_figure_table) <- rownames(big_modules_df_final)

color_palette <- c("#CC79A7", "#009E73")
module_table_figure <-gt(df_for_figure_table) %>%
  tab_header(title = md("**Representation of the families in the modules**") ,
  subtitle = substrate) %>%
  cols_label(Family = "Family", total_members = "Total members" , n_members = "Number of members in big modules", percentage_members = "Percentage of members in big modules",max_value = "max percentage in one module") %>%
  gt_color_rows(
    columns = c(percentage_members,max_value), 
    domain = c(0, 100),
    palette = color_palette
  ) %>%
  tab_footnote(
    footnote = "Modules with at least 4 members",
    locations = cells_column_labels(columns = n_members)
  )

module_table_figure

gtsave_extra(module_table_figure,"module_table_figure.png", zoom = 1, vwidth = 1200)
