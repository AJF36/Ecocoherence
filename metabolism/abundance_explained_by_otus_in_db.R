rm(list= ls())

library(readr)
library(dplyr)
library(tidyr)



#Load OTU count table
count_table <- read.table(file = "/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv" , sep = ",")

#Transform to relabundance
# Convertir las abundancias absolutas a relativas
abundance_relative <- count_table
abundance_relative <- sweep(count_table, 2, colSums(count_table), FUN = "/")
abundance_relative$OTU <- rownames(abundance_relative)

##Now load the OTUs to filter

otus_to_filter <- read.table("/home/ajf/Desktop/CNB/ecocoherence_metabolism/genome_aligment/aligned_ssu.tsv",sep = "\t")

otus_to_filter <- otus_to_filter$V1
otus_to_filter

abundance_relative_filtered <- filter(abundance_relative, OTU %in% otus_to_filter)
length(colnames(abundance_relative_filtered))
mean(colSums(abundance_relative_filtered[,-348])) ###% explained of abundance explained in different samples
colSums(abundance_relative_filtered[,-348])
abundance_relative_filtered_b <- abundance_relative_filtered[,!grepl("Seawater",names(abundance_relative_filtered))]
abundance_relative_filtered_b_c <- abundance_relative_filtered_b[,!grepl("Chitosan",names(abundance_relative_filtered_b))]
length(abundance_relative_filtered_b_c)
mean(colSums(abundance_relative_filtered_b_c[,-216])) #### the mean across al the substrates is 0.49 without seawater and chitosan
##How about the different modules






####

#####
modules_algilante <- read.table("/home/ajf/Desktop/CNB/ecocoherence_sparcc/Alginate/functionink_tmp/Partition-NL_Average_StopStep-1046_interactions_filtered_p0.01_threshold_Alginate_.tsv_guildGT4.txt",sep= " ",skip = 1)
modules_algilante_separated <- separate(modules_algilante,col = V1,into = c("ESV","guild"), sep = "\t")
modules_algilante_separated_filtered <- filter(modules_algilante_separated, ESV %in% otus_to_filter)

length(modules_algilante_separated_filtered$ESV)/length(modules_algilante_separated$ESV)
?sweep


###What is the % of the abundace of the modules explained

rownames(count_table) %in% modules_algilante_separated$ESV
abundance_relative_module <- count_table[rownames(count_table) %in% modules_algilante_separated$ESV,]
abundance_relative_module <- sweep(abundance_relative_module,2,colSums(abundance_relative_module),"/")
colSums(abundance_relative_module)

abundace_relative_module_filtered<- filter(abundance_relative_module, rownames(abundance_relative_module) %in% otus_to_filter)
mean(colSums(abundace_relative_module_filtered))

abundace_relative_module_filtered_b <- abundace_relative_module_filtered[,!grepl("Seawater",colnames(abundace_relative_module_filtered))]
abundace_relative_module_filtered_b_c <- abundace_relative_module_filtered_b[,!grepl("Seawater",colnames(abundace_relative_module_filtered_b))]
# mean(colSums(abundace_relative_module_filtered_b_c))
colSums(abundace_relative_module_filtered_b_c)