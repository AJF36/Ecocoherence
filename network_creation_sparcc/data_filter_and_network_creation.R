#Script para limpiar el dataset
rm(list = ls())
library(phyloseq)
library(readr)
library(dplyr)
library(this.path)

work_dir <- this.dir()
## Este codigo es de Alberto
# --- Load OTUs
fileOTU <- file.path(work_dir, "..", "data", "marine_particles_source_data", "count_table.ESV.4R.csv")
fileOTU <- normalizePath(fileOTU, mustWork = FALSE)
print(fileOTU)
otu.in=read.csv(fileOTU)
otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)
fileOTU
# --- Load samples metadata
# fileSample="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
fileSample = file.path(work_dir, "..", "data", "marine_particles_source_data", "samples_properties","samples_metadata_AltSubstrSyntax.4R.tsv")
fileSample  <- normalizePath(fileSample, mustWork = FALSE)
sample_metadata = import_qiime_sample_data(fileSample)
#Que 

# --- Load taxonomy and clean it 

# I try with the manual upload
# fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
tax.pseq = tax_table(as.matrix(taxonomy))

##final del codigo de Alberto

#Creo el objeto de phyloseq

phylo_obj <- merge_phyloseq(sample_metadata, tax.pseq, otu.pseq)

#First we remove seawater

phylo_obj_B <- prune_samples(sample_data(phylo_obj)$Media == "Beads", phylo_obj)


#--------------------------------------------------------------------
substrates <- unique(sample_data(phylo_obj_B)$Substrate) #Hago un vector con los sustratos

length <- length(substrates)
#Now starts the for loop
for (i in 1:2) {
  setwd(work_dir)
  substrate <- substrates[i]
  system2("mkdir",substrate)
  setwd(paste(work_dir,substrate,sep = "/"))
  
  #Now we take only 1 substatre
  
  phylo_obj_B_substrate <- prune_samples(sample_data(phylo_obj_B)$Substrate == substrate ,phylo_obj_B)
  
  #---------------------------------------------------------------------
  
  #Now we take the samples that have at least 5000 reads
  
  data <- prune_samples(sample_sums(phylo_obj_B_substrate) > 5000, phylo_obj_B_substrate)
  
  #We take the otus that dont have at least 100 reads
  
  data <- prune_taxa(taxa_sums(data) > 10,data)
  
  
  #Eliminar los otus que no esten en al menos 10 muestras
  #---------------------------------------------------------------------------------------------------------
  # Paso 1: Obtener la tabla de OTUs (abundancias)
  otu_table_data <- otu_table(data)
  
  # Paso 2: Contar cuántas muestras tienen presencia de cada OTU (valores > 0)
  otu_presence <- apply(otu_table_data, 1, function(x) sum(x > 0))
  
  # Paso 3: Filtrar los OTUs presentes en al menos 10 muestras
  otus_to_keep <- names(otu_presence[otu_presence >= 10])
  
  # Paso 4: Usar prune_taxa() para quedarte solo con los OTUs filtrados
  data2<- prune_taxa(otus_to_keep, data)
  #---------------------------------------------------------------------------------------------------------
  otu_df <- as.data.frame(otu_table(data2))
  
  setwd(paste(work_dir,substrate,sep = "/"))
  
  write.table(otu_df, file = paste("otu_table",substrate,sep = "_"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t" )
  #Arreglamos el formato biom para fast_spar
  system(paste("sed -i '1s/^/#OTU_ID\t/'",paste("otu_table",substrate,sep = "_"),sep = " "))
  
  #Ahora vamos con fast_spar------------------------------------------------------------------------------------
  #Calculo de las correlaciones y de las covariancias
  
  command <- paste("fastspar --threshold 0.3 --otu_table ",paste("otu_table",substrate,sep = "_")," --correlation median_correlation.tsv --covariance median_covariance.tsv",sep = "")
  system(command)
  
  #bootstraping
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc",substrate,sep = "/"))
  setwd(paste(work_dir,substrate,sep = "/"))
  
  system("mkdir bootstrap_counts")
  y <- paste("fastspar_bootstrap --otu_table ",paste("otu_table",substrate,sep = "_")," --number 1000"," --prefix bootstrap_counts/",substrate,sep = "")
  system(y)
  #infer correlations of the boostrap counts
  system("mkdir bootstrap_correlation")
  system("parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*")
  
  #Calculating the p-values

  z <- paste("fastspar_pvalues --otu_table ", paste("otu_table",substrate,sep = "_"),
             " --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_", substrate, "_",
             " --permutations 1000 --outfile pvalues.tsv", sep="")
  system(z)
  #fin fastsparcc--------------------------------------------------------
  
  #Este codigo es de diego, se utiliza para filtrar los p-valores significativos---------------------------------------
  setwd(paste(work_dir,substrate,sep = "/"))
  p_values <- read_tsv("pvalues.tsv")
  
  p_values <- as.data.frame(p_values)
  
  row.names(p_values) = colnames(p_values)[2:length(colnames(p_values))]
  p_values <- p_values %>% select(-c("#OTU ID"))
  
  
  cor <- read.table("median_correlation.tsv", header = FALSE, sep = "\t", row.names = 1)
  colnames(cor) = rownames(cor)
  
  df = data.frame()
  for (i in 1:(nrow(cor)-1))
    for (j in (i+1):ncol(cor))
    {
      if (p_values[i,j] <= 0.01)
      {
        df[nrow(df) +1 ,"Sp1"] = rownames(cor)[i]
        df[nrow(df),"Sp2"] = rownames(cor)[j]
        df[nrow(df),"Cor"] = cor[i,j]
      }
    }
  
  write.table(df, paste("interactions_filtered_0.01.",substrate,".tsv", sep = "_"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  
}
  #---------------------------------------------------------------------------------
  
  
  
  ###All the code below its before the implementation of the dual model
  
  
#   #Añadir la tabla de interacciones los tipos y filtrar las correlaciones menores a 0,7
#   setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc",substrate,sep = "/"))
#   df$type <- ifelse(df$Cor > 0, 1, 2)
#   str(df)
#   
#   df <- df %>% filter(abs(Cor) > 0.7 )
#   str(df)
#   
#   write_tsv(df,paste("interactions_filtered_0.01","0.7",substrate,".tsv", sep = "_"))
#   
#   
#   
#   #---------------------------------------------------------------------------------
#   #Por ultimo functionink
#   fileNet <- paste("interactions_filtered_0.01","0.7",substrate,".tsv", sep = "_")
#   pathNet <- paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate, sep ="")
#   pathRepo <- "/home/ajf/functionInk"
#   run_pipeline = function(fileNet,pathNet,pathRepo, # mandatory
#                           weighted=FALSE,directed=FALSE,
#                           types=FALSE,method="Average",mode="all"){
#     # ... set up the environment and source dependencies
#     setwd(pathRepo)
#     src.dir=paste("scripts","analysis_R",sep="/")
#     setwd(src.dir)
#     source("extractPartDensity.R")
#     setwd(pathRepo)
#     setwd(pathNet)
#     dir.create("functionink_tmp")
#     setwd("functionink_tmp")
#     fileNetPath=paste0("../",fileNet)
#     
#     # ... process arguments
#     if(weighted == FALSE){par_w = 0}else{par_w = 1}
#     if(directed == FALSE){par_d = 0}else{par_d = 1}
#     if(types == FALSE){par_t = 0}else{par_t = 1}
#     
#     # ... build commands basic run
#     # ...... Node similarity
#     script=paste(pathRepo,"NodeSimilarity.pl",sep="/")
#     options=paste("-w",par_w,"-d",par_d,"-t",par_t,"-f",fileNetPath)
#     comm_sim=paste(script,options)
#     # ..... Node linkage
#     fileSim=paste0("Nodes-Similarities_",fileNet)
#     script=paste(pathRepo,"NodeLinkage.pl",sep="/")
#     options=paste("-fn",fileNetPath,"-fs",fileSim,"-a",method)
#     comm_link_base=paste(script,options)
#     file.hist=paste0("HistCompact-NL_",method,"_NoStop_",fileNet) # expected history file output
#     
#     # --- Run the first analysis
#     system(comm_sim) # compute nodes' similarities
#     system(comm_link_base) # cluster nodes
#     hist.comp=read.table(file=file.hist,sep="\t",header = TRUE) # read history file
#     part_density=extractPartDensity(hist.comp) # extract partition densities
#     
#     # --- Run the extraction of the communities
#     if(mode != "none"){ # if the user wants to retrieve the communities
#       # .... determine the criteria to  be used
#       labels=c("total_dens_step","int_dens_step","ext_dens_step")
#       if(mode == "total"){
#         idx=1
#       }else if(mode == "internal"){
#         idx=2
#       }else if(mode == "external"){
#         idx=3
#       }else{ # mode "all", we identify the maximum among modes
#         idx=which.max(c(part_density$total_dens,
#                         part_density$int_dens,
#                         part_density$ext_dens))
#       }
#       # ... find the step of the peak and create a new command
#       value=part_density[labels[idx]]
#       comm_link_spec=paste(comm_link_base,"-s step -v",value)
#       # ... finally run
#       system(comm_link_spec) # cluster nodes and extract partition
#     }
#     return(part_density)
#   }
#   
#   run_pipeline(fileNet, pathNet, pathRepo, weighted = TRUE,types = TRUE)
# }
