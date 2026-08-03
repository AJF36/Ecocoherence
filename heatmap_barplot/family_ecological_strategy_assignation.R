
##############################################################################################################
#This script will asign a ecological strategy to each family in a similar way that is made with the modules  #
##############################################################################################################
rm(list = ls())
library(tidyverse)
library(phyloseq)
library(this.path)
library(gplots)
library(BiotypeR)
library(pheatmap)
work_dir <- this.dir()

###Load the abundance table
# fileOTU="/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv
fileOTU <- file.path(work_dir, "..", "data", "marine_particles_source_data", "count_table.ESV.4R.csv")
fileOTU <- normalizePath(fileOTU, mustWork = FALSE)

# --- Load OTUs
otu.in=read.csv(fileOTU)

### Charge the tax table

fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))


### Add the family info to the otu table
# head(otu.in)
# head(test)
# Nfam=ntaxa(particles.rar.mod.aggr) # your modules
# Nsamp=nsamples(particles.rar.mod.aggr)

substrate <- "Alginate"




### Filter by substrate
### Lets fix an arbitrari number of rowSums, ill change it later

otu.in.substrate <- otu.in %>%
  select(matches(paste("_",substrate,"_", sep = "")) & matches("beads")) 

## Load functionink to filter the ASV
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

    functionink <- read.table(file, sep = "\t", header = TRUE)
    functionink$ESV <- rownames(functionink)
    colnames(functionink) <- c("mod", "ESV")

otu.in.substrate.f <- filter(otu.in.substrate,rownames(otu.in.substrate) %in% functionink$ESV)

library(stringr)

#### Code to make the mean of the replicates
times <- unique(str_extract(names(otu.in.substrate.f), "(?<=Beads_)[A-Za-z]+_([0-9]+)") |> 
                  str_extract("[0-9]+"))


df_mean_by_time <- sapply(times, function(t) {
  # construimos el patrón dinámico usando la variable substrate
  cols <- grep(paste0("Beads_", substrate, "_", t, "_[A-Z]$"), names(otu.in.substrate.f), value = TRUE)
  rowMeans(otu.in.substrate.f[, cols])
}) %>% as.data.frame()

# renombramos las columnas con el tiempo (opcional)
colnames(df_mean_by_time) <- paste0(substrate, "_", times)
####


otu.in.substrate.f.relative <- sweep(df_mean_by_time,2,colSums(df_mean_by_time),"/")

otu.in.substrate.f.relative$family <- taxonomy$Family[match(rownames(otu.in.substrate.f.relative),rownames(taxonomy))]

substrate_df_filtered_grouped <- otu.in.substrate.f.relative %>%
  group_by(family) %>%
  summarise(across(everything(),mean))

colSums(substrate_df_filtered_grouped[,-1])
rowSums(substrate_df_filtered_grouped[,-1])
# substrate_df_filtered <- substrate_df %>%
#   filter(rowSums(across(where(is.numeric))) != 0)



rownames(substrate_df_filtered_grouped) <- substrate_df_filtered_grouped$family
# Nmod = 200 # your modules ### En este caso seran familias
Nmod <- nrow(substrate_df_filtered_grouped)
# Esto sera el numero de columnas de samples
Nsamp <- ncol(substrate_df_filtered_grouped) -1

# --- Create association vectors to determine the preferred stage
Nrep=1 # number replicates
Nstage=5 # Number artificial vectors we create to test associations
founder=1*Nrep  # Number of bins in each stage, 1/founder is the null probablity
if(substrate == "Carrageenan"){ # one sample missing, 24C
  early=(5*Nrep)-1 # same for other stages
}else{
  early=5*Nrep # same for other stages
}
mid=2*Nrep
late=4*Nrep
all=12*Nrep
1/early
# --- Now create the vectors,
# ..... first the null vectors
# profiles=matrix(0,nrow=(Nmod+Nstage),ncol=Nsamp) # Matrix of rows all the families + the artificial vectors 
### Voy a cambiarlo, y lo que voy a hacer es hacer un df con los vectores y unirlo por filas al df de la abundancia
profiles=matrix(0,nrow=(5),ncol=Nsamp) # Matrix of rows all the families + the artificial vectors 

profiles[1,(1:founder)]=1/founder        # And as columns samples
profiles[2,(founder+1):(founder+early)]=1/early 
profiles[3,(founder+early+1):(founder+early+mid)]=1/mid 
profiles[4,(founder+early+mid+1):(founder+early+mid+late)]=1/late 
profiles[5,]=1/Nsamp 


colnames(profiles) <- colnames(substrate_df_filtered_grouped)[-1]

# rownames(profiles)[6:nrow(profiles)] <- substrate_df_filtered_grouped$family
profiles_df <- as.data.frame(profiles)
dummy_df <- data.frame(family = c("attachment","selection","mid","late","all"))
profiles_df <- cbind(profiles_df,dummy_df)

substrate_df_filtered_grouped <- rbind(profiles_df,substrate_df_filtered_grouped)
substrate_df_filtered_grouped <- select(substrate_df_filtered_grouped,-family)

dist_matrix <- t(as.matrix(substrate_df_filtered_grouped))
dist_JSD <- dist.JSD(dist_matrix)
plot(hclust(dist_JSD, method = "average"))

# test_4 <- as.matrix(dist_JSD)
# heatmap.2(as.matrix(test_4),
# trace = "none",
# cexCol = 0.8)


