rm(list = ls())

library(readr)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(KEGGREST)
library(this.path)

work_dir <- this.dir()
metadata_filepath <- file.path(work_dir,"genome_aligment","matched_ESV_id_0.97.tsv")
metadata_filepath <- normalizePath(metadata_filepath)
# metadata_genomes <- read_tsv("/home/ajf/Desktop/CNB/ecocoherence_metabolism/genome_aligment/matched_ESV_id_0.97.tsv",col_names = F)
metadata_genomes <- read_tsv(metadata_filepath,col_names = F)
metadata_genomes <- select(metadata_genomes,c("X9","X10"))
colnames(metadata_genomes) <- c("ASV","ID")
metadata_genomes <- filter(metadata_genomes, ID != "*")


#tax table for the families, we will use the same tax tamble as other analisys instead the taxonomy of the metadata
# I try with the manual upload
fileTaxonomy<- file.path(work_dir,"..","data","marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(tax_file_path)
# fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))


###Adding the family column
metadata_genomes$family <- taxonomy$Family[match(metadata_genomes$ASV,rownames(taxonomy))]
taxonomy$Family[match(metadata_genomes$ASV,rownames(taxonomy))]

##Make a df only with the unique ID (just for efficiency)
length(metadata_genomes$ID)
length(unique(metadata_genomes$ID))
unique_ids <- unique(metadata_genomes$ID)
metadata_genomes_unique <- data.frame("ID" = unique_ids)
metadata_genomes_unique$Family <- NA
metadata_genomes_unique$Family <- metadata_genomes$family[match(metadata_genomes_unique$ID,metadata_genomes$ID)]
metadata_genomes_unique$ID <- gsub("^.*?(GCF_|GCA_)", "\\1", metadata_genomes_unique$ID)
metadata_genomes_unique$ID <- gsub("\\..*$", "", metadata_genomes_unique$ID)
metadata_genomes_unique$ID <- gsub("_$", "", metadata_genomes_unique$ID)


###Now we need to add the ko column
##Lets test
metadata_genomes_unique$KO_list <- NA
metadata_genomes_df <- as.data.frame(metadata_genomes_unique)

rownames(metadata_genomes_df) <- metadata_genomes_unique$ID




##Uncomment to copy the kegg files to the directory of work
# for (accesion in metadata_genomes_df$ID){
#   print(accesion)
#   setwd("/home/ajf/Desktop/globdb_database/globdb_r220_gff_kegg")
#   system(paste("cp ",paste(accesion,"_kegg.gff.gz",sep="")," /home/ajf/Desktop/CNB/ecocoherence_metabolism/ko_analysis/kegg_files",sep = ""))
# }



###Add the ko from the sff3 file of globdb
for (accesion in metadata_genomes_df$ID){
  print(accesion)
  kegg_test <- read.table(paste("/home/ajf/Desktop/CNB/ecocoherence_metabolism/ko_analysis/kegg_files/",accesion,"_kegg.gff",sep=""),sep = "\t",header = F)
  kegg_test2 <-separate(kegg_test,col = V9, sep = ";", into = c("ID","name","db_xref_and_KO_fam","product"))
  kegg_test2$name <- gsub("Name=","",kegg_test2$name)
  kegg_test2_no_na <- drop_na(kegg_test2)
  ##Now create the string with the KO names and add it to the metadata_genomes_table
  list_ko <-paste(kegg_test2_no_na$name,collapse = ";")
  metadata_genomes_df[accesion,"KO_list"] <- list_ko

}

setwd("/home/ajf/Desktop/CNB/ecocoherence_metabolism/ko_analysis")
write_tsv(metadata_genomes_df,"list_id_with_ko.tsv")


##Take only the unique ko to make the db
list_ko <- paste(metadata_genomes_df$KO_list,collapse = ";")
list_ko_separated <- strsplit(list_ko,split = ";")[[1]]
class(list_ko_separated)

length(unique(list_ko_separated))
unique_ko <- unique(list_ko_separated)
unique_ko[1]
###Get the function of the BRITE KEGG database

# kegg_db <- data.frame("KO" = unique_ko, "KO_family" = NA)
# i= 1
# for (ko in unique_ko){
#   print(ko)
#   kegg_object <- keggGet(ko)
#   kegg_object[[1]]$BRITE
#   kegg_db$KO_family[i] <- kegg_object[[1]]$BRITE[2]
#   i = i +1 
# }

# ##Create the db of the ko using keggrest
# kegg_db <- data.frame("KO" = unique_ko, "KO_family" = NA, stringsAsFactors = FALSE)
# i = 1
# 
# for (i in 1:length(unique_ko)) {
#   ko <- unique_ko[i]
#   print(i)
#   
#   # Intentar obtener la información con tryCatch
#   kegg_object <- tryCatch({
#     keggGet(ko)
#   }, error = function(e) {
#     # Si ocurre un error, devolver NULL
#     return(NULL)
#   })
#   
#   # Verificar si keggGet fue exitoso
#   if (!is.null(kegg_object)) {
#     # Si la URL fue encontrada, extraer la familia KO
#     brite_info <- kegg_object[[1]]$BRITE
#     if (length(brite_info) >= 2) {
#       kegg_db$KO_family[i] <- brite_info[2]  # Usamos el segundo valor de BRITE
#     } else {
#       kegg_db$KO_family[i] <- NA  # Si no hay suficiente información, dejamos NA
#     }
#   } else {
#     # Si no se encuentra el KO, asignar NA
#     kegg_db$KO_family[i] <- NA
#   }
# }

# # write the db of the kos
# setwd("/home/ajf/Desktop/CNB/ecocoherence_metabolism/ko_analysis/")
# write.table(kegg_db,"kegg_db.tsv",sep = "\t" ,row.names = F)


# ###Lets start comparing the genome length of the distinct families
# grouped_data <- metadata_genomes_df %>%
#   group_by(family) %>%
#   summarize("mean_genome_size" = mean(genome_size))
# 
# ggplot(grouped_data,aes(family, mean_genome_size, )) +
#   geom_point() +
#   theme(axis.text = element_text(size = 5)) 
# 
# 
# ###Compare the genomes
# 
# pair_wise_combination_ESV <- expand_grid("ESV_1" =metadata_genomes_df$ESV_list,"ESV_2" = metadata_genomes_df$ESV_list)
# pair_wise_combination_ESV$overlap <- NA
# pair_wise_combination_ESV$family
# 
# similarity_jaccard <- function(vec1, vec2) {
#   length(intersect(vec1, vec2)) / length(unique(c(vec1, vec2)))
# }

# 
# for (n in 1:nrow(pair_wise_combination_ESV)) {
#   print(n)
#   OTU_1<-pair_wise_combination_ESV$ESV_1[n]
#   OTU_2<-pair_wise_combination_ESV$ESV_2[n]
#   family1 <-metadata_genomes_df$family[metadata_genomes_df$ESV_list == OTU_1]
#   family2 <-metadata_genomes_df$family[metadata_genomes_df$ESV_list == OTU_2]
#   KO_OTU_1 <- metadata_genomes_df$KO_list[metadata_genomes_df$ESV_list == OTU_1]
#   KO_OTU_2 <- metadata_genomes_df$KO_list[metadata_genomes_df$ESV_list == OTU_2]
#   KO_OTU_1 <- unlist(strsplit(KO_OTU_1,split = ";"))
#   KO_OTU_2 <- unlist(strsplit(KO_OTU_2,split = ";"))
#   
#   pair_wise_combination_ESV$overlap[n] <- similarity_jaccard(KO_OTU_1,KO_OTU_2)
#   
#   pair_wise_combination_ESV$family[n] <-ifelse(family1 == family2,family1,"different_family")
#   ###Habria que añadir la informacion de las familias
# }
# 
# 
# pair_wise_combination_ESV_gruoped_by_family <- pair_wise_combination_ESV %>%
#   filter(overlap != 1) %>%
#   group_by(family) 
#   # summarise("mean_overlap" = median(overlap))
# 
# 
# ggplot(pair_wise_combination_ESV_gruoped_by_family,aes(family,overlap)) +
#   geom_point()
