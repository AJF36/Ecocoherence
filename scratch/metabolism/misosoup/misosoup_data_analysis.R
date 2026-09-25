rm(list = ls())

# install.packages("yaml")

library(readr)
library(tidyr)
library(dplyr)
library(yaml)
library(ggplot2)
library(pheatmap)
setwd("/home/ajf/Desktop/CNB/ecocoherence/metabolism/misosoup/mod_48")

my_file <- yaml.load_file("misosoup_Agarose_mod_48.yaml")

my_file$gal
str(my_file)

names(my_file[[1]][[1]][[1]])
length(names(my_file[[1]][[1]][[1]]))


number_of_communities <- length(my_file[[1]][[1]])


###Construction of the df of the communities
# names(my_file[[1]][[1]][[52]])[grep("Growth",names(my_file[[1]][[1]][[52]]))]
df <- bind_cols(lapply(my_file[[1]][[1]][[1]], as.data.frame))
community_df <- data.frame("Community" = numeric(),"Models" = character())


for (i in 1:number_of_communities){
  print(i)
  print(names(my_file[[1]][[1]][[i]])[grep("Growth",names(my_file[[1]][[1]][[i]]))])
  dummy_df <- data.frame("Community" = i,"Models" = (names(my_file[[1]][[1]][[i]])[grep("Growth",names(my_file[[1]][[1]][[i]]))]))
  community_df <- rbind(dummy_df,community_df)
}


###Format the id

community_df$Models <- gsub("Growth_","",community_df$Models)
community_df$Models <- gsub("_model","",community_df$Models)


####Add the  module info

setwd("/home/ajf/Desktop/CNB/ecocoherence/functionink/Agarose/functionink_tmp")
###functionink df
functionink_df <- read.table("Partition-NL_Average_StopStep-640_interactions_filtered_p0.01_threshold_Agarose_.tsv_guildGT4.txt")
functionink_df$ESV <- rownames(functionink_df)

setwd("/home/ajf/Desktop/CNB/ecocoherence/metabolism/genome_aligment")
###esv-model df
matched_esv <- read_tsv("matched_ESV_id_0.97.tsv",col_names = F)
matched_esv <- subset(matched_esv,select = c(X9,X10))
matched_esv <- filter(matched_esv,X10 != "*")
matched_esv$X10 <- gsub("[A-Z]{2}_([A-Z]{3}_[0-9]+).[0-9]","\\1",matched_esv$X10)

###family_esv
setwd("/home/ajf/Desktop/CNB/ecocoherence/data/marine_particles_source_data")
fileTaxonomy <-"sequence_table.ESV.fasta_RDPclassified.txt"
# fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
taxonomy$otu <- taxonomy$taxa_id
taxonomy <- subset(taxonomy, select = c('otu', 'Family', 'Genus', 'Order', 'Class'))



community_df$ESV <- matched_esv$X9[match(community_df$Models,matched_esv$X10)]
community_df$Guild <- functionink_df$guild[match(community_df$ESV,functionink_df$ESV)]
community_df$Family <-  taxonomy$Family[match(community_df$ESV,taxonomy$otu)]
number_of_members_per_community <- table(community_df$Community)

community_df$n_members <- number_of_members_per_community[match(community_df$Community,names(number_of_members_per_community))]
# community_df$Community,names(number_of_members_per_community)
community_df <- as.data.frame(community_df)

sum(length(unique(community_df$Community)))

communitys_of_1 <- filter(community_df,n_members == 1)
# community_17_df <- filter(community_df,n_members == 17)
# View(table(community_17_df$Guild))
modules <- "mod_48"
modules
length(unique(community_df$Models))
community_df$n_members <- as.numeric(community_df$n_members)
community_df$Community <- paste(community_df$Community,modules ,sep = "_") 
setwd("/home/ajf/Desktop/CNB/ecocoherence/metabolism/misosoup/results_misosoup")
write_tsv(community_df,paste("communities_misosoup_df_",modules,".tsv",sep = ""))
paste(community_df,modules,".tsv",sep = "")
# pair_wise_familys <- expand_grid("F1" = community_df$Family,"F2" = community_df$Family)


######Analysis of the structure of the communities

number_of_members_community <- 2

community_structure_df <- filter(community_df, n_members == number_of_members_community)


matrix_communities_structure <- matrix(0,ncol = length(unique(community_structure_df$Community)),nrow = length(unique(community_structure_df$Family)))

rownames(matrix_communities_structure) <- unique(community_structure_df$Family)
colnames(matrix_communities_structure) <- unique(community_structure_df$Community)

for (community in unique(community_structure_df$Community)){
  community_bucle_df <- filter(community_structure_df, Community == community)
  familys <- community_bucle_df$Family
  
  for (family in familys){
    print(family)
    print(community)
    
    matrix_communities_structure[family,as.character(community)] <-  matrix_communities_structure[family,as.character(community)] + 1
  }
  
}

pheatmap(matrix_communities_structure)



####
#Analysis of the modules
modules_df <- data.frame("modules"= character())

for (C in unique(community_structure_df$Community)){
  loop_modules_df <- filter(community_structure_df, Community == C)
  modules <- sort(loop_modules_df$Guild)
  modules<- paste(modules, collapse = ";")
  dummy_df <- data.frame("modules" = modules)
  modules_df <- rbind(modules_df,dummy_df)

}
modules_sums <- as.data.frame(table(modules_df))



###### Heatmap for the number of times that 2 species appear in the same community
# number_of_members_community <- 4

community_df_f <- filter(community_df, n_members == number_of_members_community) %>% 
  group_by(Community) %>%
  summarize("members" = toString(Family)) 
  
community_df_f <- separate(community_df_f,col = members, into = paste0("F",1:number_of_members_community),sep = ", ")

family_matrix <- matrix(0,nrow = length(unique(community_df$Family)),ncol = length(unique(community_df$Family)))
rownames(family_matrix) <- unique(community_df$Family)
colnames(family_matrix) <- unique(community_df$Family)
any(family_matrix != 0)

for (i in 1:length(community_df_f$F1)){
  test <- community_df_f[i,-1]
  combinaciones <- as.data.frame(combn(test,2))
  
  for (z in 1:length(colnames(combinaciones))){
    print(z)
    combination <- combinaciones[,z]
    F1 <- combination[[1]]
    F2 <- combination[[2]]
    print(paste(F1,F2))
    family_matrix[F1,F2] <- family_matrix[F1,F2] + 1
    
    if (F1 != F2){
      family_matrix[F2,F1] <- family_matrix[F2,F1] + 1
    }

    sum(family_matrix == 2)

  }
 
  

}
par(mar = c(10, 6, 6, 8))  # c(inferior, izquierdo, superior, derecho)

max_family_matrix <-max(family_matrix) + 1

pheatmap(family_matrix,
         color = colorRampPalette(c("blue", "red"))(max_family_matrix),
         )




