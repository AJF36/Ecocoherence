##################################################
# phyloseq_2functionInk.R
##################################################
# In this script I aim to plot the abundances of the ESVs belonging to
# the different guilds identified with functionInk for each substrate.
#
# Zürich, October 2020
# Theoretical Biology, ETH
# alberto.pascual.garcia@gmail.com
###################################################

rm(list=ls())
library(phyloseq)
library(reshape2)
#library(RDPutils) # Notused
library(ggplot2)
library(ggpubr) # for ggboxplot
library(BiotypeR)
library(readr)
library(dendextend)
library(dplyr)
library(this.path)
#extra=0 # Fixing this parameter to 1 will generate additional plots and print file outputs


# Edit options -----------------
# --- Select substrate that will be analysed
# substrate ="Agarose"
# stop.step= 237
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin","Carrageenan")
min.module <- 4
results_df <- data.frame(substrate = character(), module = character(), strategy = character(), stringsAsFactors = FALSE)
for (substrate in substrates){
  script_dir<- this.dir()
  results_dir <- normalizePath(file.path(script_dir, "..", "results", "6i_classify_modules_by_ecological_strategy_combined"), mustWork = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.source=script_dir
  home_dir <-  Sys.getenv("HOME")
  dir.funcInk <- file.path(home_dir,"functionInk")
  dir.funcInk <- file.path(script_dir,"..","results","4a_detect_modules_functionink",substrate,"functionink_tmp")
  dir.funcInk <- normalizePath(dir.funcInk)
  # --- Set working pathways
  substrate
  # dir.source="/home/ajf/Desktop/CNB/ecocoherence_sparcc"
  # dir.funcInk=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc",substrate,sep="/")
  # dirOut=paste(dir.funcInk,"figures",sep="/")
  dir.funcInk
  # --- Fix files

  # fileOTU="/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv"
  # fileSample="/home/ajf/Desktop/CNB/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
  fileOTU <- file.path(script_dir, "..", "data", "marine_particles_source_data", "count_table.ESV.4R.csv")
  fileOTU <- normalizePath(fileOTU, mustWork = FALSE)
  fileSample = file.path(script_dir, "..", "data", "marine_particles_source_data", "samples_properties","samples_metadata_AltSubstrSyntax.4R.tsv")
  fileSample  <- normalizePath(fileSample, mustWork = FALSE)
  # fileFun=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp/Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep="")
  
  ###############
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  subs_path <- file.path(script_dir,"..","results","5_filter_modules_by_size_and_plot_abundance",substrate)
  subs_path <- normalizePath(subs_path)
  setwd(subs_path)
  subs_path
  dir_path <- "."
  file <- list.files(
    path = dir_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE
  )
  fileFun <- file
  ###############
  # Primero, verifica que los archivos coincidan con el patrón esperado
  # fileFun <- list.files(path = dir_path, pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv$"), full.names = TRUE)
  
  # Imprime los archivos encontrados para depuración
  print(paste("Archivos encontrados para", substrate, ":"))
  print(fileFun)
  
  # Asegúrate de que fileFun no esté vacío
  if(length(fileFun) == 0) {
    stop(paste("No se encontraron archivos para el sustrato", substrate))
  }
  
  # Cargar el archivo de función si se encuentran archivos
  partition = read.table(fileFun, row.names = 1)
  colnames(partition) = "modules"
  partition$modules = paste("mod", partition$modules, sep = "_")
  
  # fileFun <- list.files(pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv$"),full.names = TRUE)
  # print(fileFun)
  file <- list.files(
    path = dir.funcInk,
    pattern = paste("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_",substrate,"_.tsv$",sep=""),
    full.names = F
  )
  # fileFun <- paste("Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep="")
  fileFun <- file
  ######## STOP EDITING
  
  # Load data ------------
  setwd(dir.source)
  
  # --- Load OTUs
  otu.in=read.csv(fileOTU)
  otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)
  
  # --- Load samples metadata
  sample_metadata = import_qiime_sample_data(fileSample)
  
  # --- Load partition of functionInk
  # ... we are going to use it as a taxonomy
  setwd(dir.funcInk)
  partition=read.table(fileFun,row.names = 1)
  colnames(partition)="modules"
  partition$modules=paste("mod",partition$modules,sep="_")
  
  
  #Until here seems good------------------------------------------------------------
  
  
  # Preprocessing (controls) ------------------
  # --- Control that there are no empty rows/columns 
  sum(taxa_sums(otu.pseq) == 0) # 2046 have no observation
  otu.pseq = prune_taxa(taxa_sums(otu.pseq) > 0, otu.pseq) # delete these guys
  any(taxa_sums(otu.pseq) == 0) # double check (FALSE)
  any(sample_sums(otu.pseq) == 0) # double check (FALSE)
  ntaxa(otu.pseq)
  
  # Process partition and create phyloseq object ----
  # ... we create a vector with only the most relevant modules
  partIds=unique(partition$modules)
  partIds
  part.count=vector(mode = "numeric",length = length(partIds))
  part.count
  names(part.count)=partIds 
  #This creates a vector with the counts of each module
  for(part in partition$modules){
    part.count[part]=part.count[part]+1
  }
  i=0
  for(part in partition$modules){
    i=i+1
    if(part.count[part] < min.module){
      partition$modules[i]="none"
    }
  }
  
  
  part.count
  fileOut=file.path(results_dir,paste(fileFun,"_guildGT",min.module,".txt",sep=""))
  fileOut
  partition.out=as.data.frame(partition[which(partition$modules != "none"),])
  rownames(partition.out)=rownames(partition)[which(partition$modules != "none")]
  colnames(partition.out)="guild"
  write.table(partition.out,file=fileOut,quote=FALSE,sep="\t")
  partition.out
  # ... then we map these modules to the ESVs
  partition.all=matrix("none",nrow=dim(otu.in)[1],ncol=1)
  rownames(partition.all)=rownames(otu.in)
  matched=match(rownames(partition.all),rownames(partition))
  partition.all[!is.na(matched),]=partition$modules[matched[!is.na(matched)]]
  colnames(partition.all)="Modules"
  
  # ... let us say that the partition is the taxonomy
  tax.pseq = tax_table(as.matrix(partition.all))
  
  # ... Finally, create a single phyloseq object
  particles=merge_phyloseq(otu.pseq,sample_metadata,tax.pseq) 
  idx=sort(as.numeric(levels(particles@sam_data$Time)),index.return=TRUE)$ix # reorder levels
  particles@sam_data$Time=factor(particles@sam_data$Time,
                                 levels(particles@sam_data$Time)[idx])
  nsamples(particles)
  particles@sam_data
  
  # Process the whole dataset ---------
  # ... Select a subset
  particles.nosea = subset_samples(particles, Media != "Seawater") # exclude seawater from the analysis
  particles.select = subset_samples(particles.nosea, Substrate == substrate) #"Alginate")
  nsamples(particles.select)
  
  # --- Before rarefying, we want to aggregate data within the same module
  # .....  to create a new count table for LV estimation
  # ..... The next steps are needed because there is a bug in the version of phyloseq I am using:
  # ..... https://github.com/joey711/phyloseq/issues/223
  fake.tax=as.matrix(particles.select@tax_table)
  fake.tax=cbind(fake.tax,fake.tax)
  fake.tax
  colnames(fake.tax)=c("Modules","modules")
  tax_table(particles.select) <- fake.tax
  particles.select.aggr=tax_glom(particles.select, taxrank = "Modules") # agglomerate, by default returns one ESV belonging to that class
  matched=match(rownames(particles.select.aggr@otu_table),
                rownames(particles.select@tax_table))# and we want the classes in row names, so look for the class for each ESV
  
  aggr.names=particles.select@tax_table[matched,1] # and use the classes as names for the otu table
  

  rownames(particles.select.aggr@otu_table)=aggr.names 
  particles.select.aggr@otu_table
  fileOut=file.path(results_dir,paste("otu_table_",substrate,"_byModulesSize",min.module,".tsv",sep=""))
  write.table(particles.select.aggr@otu_table,file=fileOut,quote = FALSE, sep="\t")
  #plot(particles.select.aggr@otu_table[3,],particles.select.aggr@otu_table[2,])
  
  # --- Now yes, we work with a rarefied dataset
  sampling=1000 # select the size here, 1018 is the minimum sampling sites, and the alpha diversity pattern is already there
  set.seed(30052018) # Today's date 30/05/2018. Stored for reproducibility
  particles.rar = rarefy_even_depth(particles.select, sample.size = sampling)
  nsamples(particles.rar)
  
  # ..... reorder the levels for time in particles.rar
  particles.rar@sam_data$Time=as.factor(particles.rar@sam_data$Time)
  idx=sort(as.numeric(levels(particles.rar@sam_data$Time)),index.return=TRUE)$ix # reorder levels
  particles.rar@sam_data$Time=factor(particles.rar@sam_data$Time,
                                     levels(particles.rar@sam_data$Time)[idx]) # reorder levels
  
  
  #
  # --- Create factors to separate replicates
  #
  
  sample_data(particles.rar)$RepTime <- as.factor(paste0(sample_data(particles.rar)$Time, sample_data(particles.rar)$Replica))
  
  # #particles.rar.byRepTime =merge_samples(particles.rar, "RepTime") # don't merge
  # ..... Now we want to sort this factor, but it is hard cause has numbers and letters, e.g. 204A
  # ..... so we create a vector with the order we want, and then we match it
  ii=sort(as.numeric(levels(particles.rar@sam_data$Time)),index.return=TRUE)$ix # Reorder time levels
  times=as.character(levels(particles.rar@sam_data$Time)[ii]) # ordered times are our reference
  order.fac=c()
  for(i in 1:length(times)){
    tmp.A=paste(times[i],"A",sep="") # Now create the strings per replica, e.g. 204A
    tmp.B=paste(times[i],"B",sep="") # 204B
    tmp.C=paste(times[i],"C",sep="") # 204C
    tmp.vec=c(tmp.A,tmp.B,tmp.C)  
    order.fac=c(order.fac,tmp.vec) # And the final vector grows every cycle
  }
  # ..... And we can match now to this vector
  matched=match(order.fac,levels(particles.rar@sam_data$RepTime))
  particles.rar@sam_data$RepTime=factor(particles.rar@sam_data$RepTime,
                                        levels(particles.rar@sam_data$RepTime)[matched]) # reorder level
  particles.rar.byRepTime=particles.rar
  
  # .... transform to proportions 
  particles.rar.freq= transform_sample_counts(particles.rar, function(x) 100 * x/sum(x))
  particles.rar.byRepTime.freq = transform_sample_counts(particles.rar, function(x) 100 * x/sum(x))
  
  # .... Perform some double checks
  # plot(colSums(particles.rar.byRepTime.freq@otu_table)) # all should be 100
  length(which(particles.rar.byRepTime.freq@tax_table == "none")) # 
  dim(particles.rar.byRepTime.freq@otu_table) # 13120 ESVs
  max(part.count) # the module with more ESVs has 66
  length(which(particles.rar.byRepTime.freq@otu_table[,1] == 0)) # how many are zero in a random sample
  idx.check=head(which(particles.rar.byRepTime.freq@otu_table[,1] == 0))
  particles.rar.byRepTime.freq@otu_table[idx.check,1] # check that is indeed the case, are these represented in the bar plot at all?
  # Associate preferences to the modules -----
  # --- Extract the table for ESVs in modules
  ESVs.modules.ids=which(particles.rar@tax_table != "none")
  ESVs.modules=rownames(particles.rar@otu_table)[ESVs.modules.ids]
  particles.rar.mod=prune_taxa(ESVs.modules, particles.rar)
  ntaxa(particles.rar.mod)
  nsamples(particles.rar.mod)
  
  # .... aggregate taxa within the same module
  # The next steps are needed because there is a bug in the version of phyloseq I am using:
  # https://github.com/joey711/phyloseq/issues/223
  fake.tax=as.matrix(particles.rar.mod@tax_table)
  fake.tax=cbind(fake.tax,fake.tax)
  fake.tax = subset(fake.tax, select = c("Modules","modules"))
  fake.tax
  colnames(fake.tax)=c("Modules","modules")
  tax_table(particles.rar.mod) <- fake.tax
  tax_table(particles.rar.mod)
  particles.rar.mod.aggr=tax_glom(particles.rar.mod, taxrank = "Modules")
  rank_names(particles.rar.mod)
  Nmod=ntaxa(particles.rar.mod.aggr) # your modules
  Nsamp=nsamples(particles.rar.mod.aggr)
  
  # --- Compute the proportions in each replicate/time point
  ESVs.modules.totals=rowSums(particles.rar.mod.aggr@otu_table)
  particles.rar.mod.aggr.prop=particles.rar.mod.aggr
  particles.rar.mod.aggr.prop@otu_table=particles.rar.mod.aggr@otu_table/ESVs.modules.totals
  # plot(rowSums(particles.rar.mod.aggr.prop@otu_table)) # double check, should sum up to one
  
  # --- Create association vectors to determine the preferred stage
  Nrep=3 # number replicates
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
  
  # --- Now create the vectors,
  # ..... first the null vectors
  profiles=matrix(0,nrow=(Nmod+Nstage),ncol=Nsamp)
  profiles[1,(1:founder)]=1/founder
  profiles[2,(founder+1):(founder+early)]=1/early
  profiles[3,(founder+early+1):(founder+early+mid)]=1/mid
  profiles[4,(founder+early+mid+1):(founder+early+mid+late)]=1/late
  profiles[5,]=1/Nsamp
  profiles
  # .... then add the real ones
  profiles[(Nstage+1):(Nmod+Nstage),1:Nsamp]=particles.rar.mod.aggr.prop@otu_table
  particles.rar.mod.aggr.prop@tax_table[,1]
  rownames(profiles)=c("attachment","selection","transition","facilitation","generalist",
                       particles.rar.mod.aggr.prop@tax_table[,1])
  
  
  
  dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    inMatrix <- t(inMatrix)
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
    
    for(i in 1:matrixColSize) {
      for(j in 1:matrixColSize) { 
        resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                               as.vector(inMatrix[,j]))
      }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix) 
  }
  
  
  
  
  # 
  # d <- dist.JSD(profiles)
  # d
  # number_of_clusters = 10
  # cluster <- hclust(d, method = "average")
  # dendrogram <- as.dendrogram(cluster)
  # 
  # dendrogram %>% set("labels_cex", 1) %>% set("labels_col", value = c(1,2,3,4,5,6,7,8,9,10,11), k=number_of_clusters) %>%set("branches_k_color",k =number_of_clusters) %>%
  #   plot(main = "Color labels \nper cluster")
  # plot(dendrogram)
  # abline(h = 0.53, lty = 3)
  # 
  # 
  # clusters <- cutree(dendrogram, k = 10)
  # clusters
  # 
  # 
  # 
  # attachment_modules <- c()
  # selection_modules <- c()
  # transition_modules <- c()
  # facilitation_modules <- c()
  # generalist_modules <- c()
  # 
  # 
  # attachment_modules <- clusters == 1
  # attachment_modules <- clusters[attachment_modules]
  # attachment_modules
  # 
  # selection_modules <- clusters == 2
  # selection_modules <- clusters[selection_modules]
  # selection_modules
  # transition_modules <- clusters == 3
  # transition_modules <- clusters[transition_modules]
  # 
  # facilitation_modules <- clusters == 4
  # facilitation_modules <- clusters[facilitation_modules]
  # facilitation_modules
  # 
  # generalist_modules <- clusters == 5
  # generalist_modules <- clusters[generalist_modules]
  # generalist_modules
  # 
  # all_modules <- c(attachment_modules, selection_modules, transition_modules, facilitation_modules, generalist_modules)
  # all_modules
  # modules_classified_df <- data.frame(names(all_modules),all_modules)
  # modules_classified_df <- modules_classified_df %>% select(all_modules)
  # colnames(modules_classified_df) <- "Module"
  # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 1,"attachment",modules_classified_df$Module)
  # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 2,"selection",modules_classified_df$Module)
  # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 3,"transition",modules_classified_df$Module)
  # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 4,"facilitation",modules_classified_df$Module)
  # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 5,"generalist",modules_classified_df$Module)
  # modules_classified_df$Strategy <- modules_classified_df$Module
  # modules_classified_df$Module <- row.names(modules_classified_df)
  # modules_classified_df <- modules_classified_df %>% 
  #   filter(!(Module %in% c("attachment", "selection", "transition", "facilitation", "generalist")))
  # modules_classified_df
  # class(modules_classified_df)
  # getwd()
  # write_tsv(modules_classified_df,"modules_classified.tsv")
  
  # Crear el data.frame a partir del vector
  
  
  # Calcula la distancia y realiza el clustering
  d <- dist.JSD(profiles)
  cluster <- hclust(d, method = "average")
  dendrogram <- as.dendrogram(cluster)
  
  # Define los nombres de los módulos
  module_names <- c("attachment", "selection", "transition", "facilitation", "generalist")
  
  # Función para verificar que cada módulo esté en un cluster único
  find_optimal_clusters <- function(dendrogram, module_names) {
    num_clusters <- length(module_names) # Inicia con al menos el número de módulos
    found_unique_clusters <- FALSE
    
    while (!found_unique_clusters) {
      clusters <- cutree(dendrogram, k = num_clusters)
      unique_modules <- unique(clusters[names(clusters) %in% module_names])
      
      # Verifica que cada módulo esté en un cluster diferente
      if (length(unique_modules) == length(module_names)) {
        found_unique_clusters <- TRUE
      } else {
        num_clusters <- num_clusters + 1
      }
    }
    
    return(num_clusters)
  }
  
  # Encuentra el número óptimo de clusters
  number_of_clusters <- find_optimal_clusters(dendrogram, module_names)
  
  # Aplica el número óptimo de clusters y asigna los módulos
  clusters <- cutree(dendrogram, k = number_of_clusters)
  
  # Visualización del dendrograma
  # dendrogram %>%
  #   set("labels_cex", 1) %>% 
  #   set("labels_col", value = 1:number_of_clusters, k = number_of_clusters) %>%
  #   set("branches_k_color", k = number_of_clusters) %>%
  #   plot(main = "Color labels \nper cluster")
  # abline(h = 0.53, lty = 3)
  
  # Asignación de módulos clasificados (similar al paso de clasificación previo)
  modules_classified_df <- data.frame(substrate = character(), module = character(), strategy = character(), stringsAsFactors = FALSE)
  for (i in seq_along(module_names)) {
    module <- module_names[i]
    module_clusters <- clusters == i
    cluster_modules <- clusters[module_clusters]
    
    if (length(cluster_modules) > 0) {
      temp_df <- data.frame(substrate = substrate,module = names(cluster_modules), strategy = module, stringsAsFactors = FALSE)
      modules_classified_df <- rbind(modules_classified_df, temp_df)
    }
  }
  getwd()
  # Exporta el dataframe a un archivo .tsv
  
  modules_classified_df <- modules_classified_df %>%
    filter(!(module %in% c("attachment", "selection", "transition", "facilitation", "generalist")))
  results_df <- rbind(results_df,modules_classified_df)


  setwd(results_dir)
  # write_tsv(modules_classified_df, paste("modules_classified_internal",substrate,".tsv",sep = ""))

}
setwd(results_dir)
write_tsv(results_df, "ecological_strategies_modules.tsv")


# substrate <- "Alginate"
# min.module=4 
# 
# # minimum size of a functionInk module to be plotted
# 
# # --- Set working pathways
# dir.source="/home/ajf/Desktop/CNB/ecocoherence_sparcc"
# dir.funcInk=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc",substrate,sep="/")
# dirOut=paste(dir.funcInk,"figures",sep="/")
# dirOut
# dir.funcInk
# # --- Fix files
# fileOTU="/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv"
# fileSample="/home/ajf/Desktop/CNB/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
# # fileFun=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp/Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep="")
# dir_path <- "functionink_tmp/"
# fileFun <- list.files(
#   path = dir_path,
#   pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv$"),
#   full.names = TRUE
# )
# fileFun
# ######## STOP EDITING
# 
# # Load data ------------
# setwd(dir.source)
# 
# # --- Load OTUs
# otu.in=read.csv(fileOTU)
# otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)
# 
# # --- Load samples metadata
# sample_metadata = import_qiime_sample_data(fileSample)
# 
# # --- Load partition of functionInk
# # ... we are going to use it as a taxonomy
# setwd(dir.funcInk)
# partition=read.table(fileFun,row.names = 1)
# colnames(partition)="modules"
# partition$modules=paste("mod",partition$modules,sep="_")
# 
# 
# #Until here seems good------------------------------------------------------------
# 
# 
# # Preprocessing (controls) ------------------
# # --- Control that there are no empty rows/columns 
# sum(taxa_sums(otu.pseq) == 0) # 2046 have no observation
# otu.pseq = prune_taxa(taxa_sums(otu.pseq) > 0, otu.pseq) # delete these guys
# any(taxa_sums(otu.pseq) == 0) # double check (FALSE)
# any(sample_sums(otu.pseq) == 0) # double check (FALSE)
# ntaxa(otu.pseq)
# 
# # Process partition and create phyloseq object ----
# # ... we create a vector with only the most relevant modules
# partIds=unique(partition$modules)
# partIds
# part.count=vector(mode = "numeric",length = length(partIds))
# part.count
# names(part.count)=partIds 
# #This creates a vector with the counts of each module
# for(part in partition$modules){
#   part.count[part]=part.count[part]+1
# }
# i=0
# for(part in partition$modules){
#   i=i+1
#   if(part.count[part] < min.module){
#     partition$modules[i]="none"
#   }
# }
# 
# 
# part.count
# fileOut=paste(fileFun,"_guildGT",min.module,".txt",sep="")
# fileOut
# partition.out=as.data.frame(partition[which(partition$modules != "none"),])
# rownames(partition.out)=rownames(partition)[which(partition$modules != "none")]
# colnames(partition.out)="guild"
# write.table(partition.out,file=fileOut,quote=FALSE,sep="\t")
# partition.out
# # ... then we map these modules to the ESVs
# partition.all=matrix("none",nrow=dim(otu.in)[1],ncol=1)
# rownames(partition.all)=rownames(otu.in)
# matched=match(rownames(partition.all),rownames(partition))
# partition.all[!is.na(matched),]=partition$modules[matched[!is.na(matched)]]
# colnames(partition.all)="Modules"
# partition.all
# # ... let us say that the partition is the taxonomy
# tax.pseq = tax_table(as.matrix(partition.all))
# 
# # ... Finally, create a single phyloseq object
# particles=merge_phyloseq(otu.pseq,sample_metadata,tax.pseq) 
# idx=sort(as.numeric(levels(particles@sam_data$Time)),index.return=TRUE)$ix # reorder levels
# particles@sam_data$Time=factor(particles@sam_data$Time,
#                                levels(particles@sam_data$Time)[idx])
# nsamples(particles)
# particles@sam_data
# 
# # Process the whole dataset ---------
# # ... Select a subset
# particles.nosea = subset_samples(particles, Media != "Seawater") # exclude seawater from the analysis
# particles.select = subset_samples(particles.nosea, Substrate == substrate) #"Alginate")
# nsamples(particles.select)
# 
# # --- Before rarefying, we want to aggregate data within the same module
# # .....  to create a new count table for LV estimation
# # ..... The next steps are needed because there is a bug in the version of phyloseq I am using:
# # ..... https://github.com/joey711/phyloseq/issues/223
# fake.tax=as.matrix(particles.select@tax_table)
# fake.tax=cbind(fake.tax,fake.tax)
# fake.tax
# colnames(fake.tax)=c("Modules","modules")
# tax_table(particles.select) <- fake.tax
# particles.select.aggr=tax_glom(particles.select, taxrank = "Modules") # agglomerate, by default returns one ESV belonging to that class
# matched=match(rownames(particles.select.aggr@otu_table),
#               rownames(particles.select@tax_table))# and we want the classes in row names, so look for the class for each ESV
# 
# aggr.names=particles.select@tax_table[matched,1] # and use the classes as names for the otu table
# 
# aggr.names
# rownames(particles.select.aggr@otu_table)=aggr.names 
# particles.select.aggr@otu_table
# fileOut=paste("otu_table_",substrate,"_byModulesSize",min.module,".tsv",sep="")
# write.table(particles.select.aggr@otu_table,file=fileOut,quote = FALSE, sep="\t")
# #plot(particles.select.aggr@otu_table[3,],particles.select.aggr@otu_table[2,])
# 
# # --- Now yes, we work with a rarefied dataset
# sampling=1000 # select the size here, 1018 is the minimum sampling sites, and the alpha diversity pattern is already there
# set.seed(30052018) # Today's date 30/05/2018. Stored for reproducibility
# particles.rar = rarefy_even_depth(particles.select, sample.size = sampling)
# nsamples(particles.rar)
# 
# # ..... reorder the levels for time in particles.rar
# particles.rar@sam_data$Time=as.factor(particles.rar@sam_data$Time)
# idx=sort(as.numeric(levels(particles.rar@sam_data$Time)),index.return=TRUE)$ix # reorder levels
# particles.rar@sam_data$Time=factor(particles.rar@sam_data$Time,
#                                    levels(particles.rar@sam_data$Time)[idx]) # reorder levels
# 
# 
# #
# # --- Create factors to separate replicates
# #
# 
# sample_data(particles.rar)$RepTime <- as.factor(paste0(sample_data(particles.rar)$Time, sample_data(particles.rar)$Replica))
# 
# # #particles.rar.byRepTime =merge_samples(particles.rar, "RepTime") # don't merge
# # ..... Now we want to sort this factor, but it is hard cause has numbers and letters, e.g. 204A
# # ..... so we create a vector with the order we want, and then we match it
# ii=sort(as.numeric(levels(particles.rar@sam_data$Time)),index.return=TRUE)$ix # Reorder time levels
# times=as.character(levels(particles.rar@sam_data$Time)[ii]) # ordered times are our reference
# order.fac=c()
# for(i in 1:length(times)){
#   tmp.A=paste(times[i],"A",sep="") # Now create the strings per replica, e.g. 204A
#   tmp.B=paste(times[i],"B",sep="") # 204B
#   tmp.C=paste(times[i],"C",sep="") # 204C
#   tmp.vec=c(tmp.A,tmp.B,tmp.C)  
#   order.fac=c(order.fac,tmp.vec) # And the final vector grows every cycle
# }
# # ..... And we can match now to this vector
# matched=match(order.fac,levels(particles.rar@sam_data$RepTime))
# particles.rar@sam_data$RepTime=factor(particles.rar@sam_data$RepTime,
#                                       levels(particles.rar@sam_data$RepTime)[matched]) # reorder level
# particles.rar.byRepTime=particles.rar
# 
# # .... transform to proportions 
# particles.rar.freq= transform_sample_counts(particles.rar, function(x) 100 * x/sum(x))
# particles.rar.byRepTime.freq = transform_sample_counts(particles.rar, function(x) 100 * x/sum(x))
# 
# # .... Perform some double checks
# plot(colSums(particles.rar.byRepTime.freq@otu_table)) # all should be 100
# length(which(particles.rar.byRepTime.freq@tax_table == "none")) # 
# dim(particles.rar.byRepTime.freq@otu_table) # 13120 ESVs
# max(part.count) # the module with more ESVs has 66
# length(which(particles.rar.byRepTime.freq@otu_table[,1] == 0)) # how many are zero in a random sample
# idx.check=head(which(particles.rar.byRepTime.freq@otu_table[,1] == 0))
# particles.rar.byRepTime.freq@otu_table[idx.check,1] # check that is indeed the case, are these represented in the bar plot at all?
# # Associate preferences to the modules -----
# # --- Extract the table for ESVs in modules
# ESVs.modules.ids=which(particles.rar@tax_table != "none")
# ESVs.modules=rownames(particles.rar@otu_table)[ESVs.modules.ids]
# particles.rar.mod=prune_taxa(ESVs.modules, particles.rar)
# ntaxa(particles.rar.mod)
# nsamples(particles.rar.mod)
# 
# # .... aggregate taxa within the same module
# # The next steps are needed because there is a bug in the version of phyloseq I am using:
# # https://github.com/joey711/phyloseq/issues/223
# fake.tax=as.matrix(particles.rar.mod@tax_table)
# fake.tax=cbind(fake.tax,fake.tax)
# fake.tax = subset(fake.tax, select = c("Modules","modules"))
# fake.tax
# colnames(fake.tax)=c("Modules","modules")
# View(fake.tax)
# tax_table(particles.rar.mod) <- fake.tax
# tax_table(particles.rar.mod)
# particles.rar.mod.aggr=tax_glom(particles.rar.mod, taxrank = "Modules")
# rank_names(particles.rar.mod)
# Nmod=ntaxa(particles.rar.mod.aggr) # your modules
# Nsamp=nsamples(particles.rar.mod.aggr)
# 
# # --- Compute the proportions in each replicate/time point
# ESVs.modules.totals=rowSums(particles.rar.mod.aggr@otu_table)
# particles.rar.mod.aggr.prop=particles.rar.mod.aggr
# particles.rar.mod.aggr.prop@otu_table=particles.rar.mod.aggr@otu_table/ESVs.modules.totals
# plot(rowSums(particles.rar.mod.aggr.prop@otu_table)) # double check, should sum up to one
# 
# # --- Create association vectors to determine the preferred stage
# Nrep=3 # number replicates
# Nstage=5 # Number artificial vectors we create to test associations
# founder=1*Nrep  # Number of bins in each stage, 1/founder is the null probablity
# if(substrate == "Carrageenan"){ # one sample missing, 24C
#   early=(5*Nrep)-1 # same for other stages
# }else{
#   early=5*Nrep # same for other stages
# }
# mid=2*Nrep
# late=4*Nrep
# all=12*Nrep
# 
# # --- Now create the vectors,
# # ..... first the null vectors
# profiles=matrix(0,nrow=(Nmod+Nstage),ncol=Nsamp)
# profiles[1,(1:founder)]=1/founder
# profiles[2,(founder+1):(founder+early)]=1/early
# profiles[3,(founder+early+1):(founder+early+mid)]=1/mid
# profiles[4,(founder+early+mid+1):(founder+early+mid+late)]=1/late
# profiles[5,]=1/Nsamp
# profiles
# # .... then add the real ones
# profiles[(Nstage+1):(Nmod+Nstage),1:Nsamp]=particles.rar.mod.aggr.prop@otu_table
# particles.rar.mod.aggr.prop@tax_table[,1]
# rownames(profiles)=c("attachment","selection","transition","facilitation","generalist",
#                      particles.rar.mod.aggr.prop@tax_table[,1])
# 
# View(profiles)
# 
# 
# dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
#   inMatrix <- t(inMatrix)
#   KLD <- function(x,y) sum(x *log(x/y))
#   JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
#   matrixColSize <- length(colnames(inMatrix))
#   matrixRowSize <- length(rownames(inMatrix))
#   colnames <- colnames(inMatrix)
#   resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
#   
#   inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
#   
#   for(i in 1:matrixColSize) {
#     for(j in 1:matrixColSize) { 
#       resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
#                              as.vector(inMatrix[,j]))
#     }
#   }
#   colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
#   as.dist(resultsMatrix)->resultsMatrix
#   attr(resultsMatrix, "method") <- "dist"
#   return(resultsMatrix) 
# }
# 
# 
# 
# 
# # 
# # d <- dist.JSD(profiles)
# # d
# # number_of_clusters = 10
# # cluster <- hclust(d, method = "average")
# # dendrogram <- as.dendrogram(cluster)
# # 
# # dendrogram %>% set("labels_cex", 1) %>% set("labels_col", value = c(1,2,3,4,5,6,7,8,9,10,11), k=number_of_clusters) %>%set("branches_k_color",k =number_of_clusters) %>%
# #   plot(main = "Color labels \nper cluster")
# # plot(dendrogram)
# # abline(h = 0.53, lty = 3)
# # 
# # 
# # clusters <- cutree(dendrogram, k = 10)
# # clusters
# # 
# # 
# # 
# # attachment_modules <- c()
# # selection_modules <- c()
# # transition_modules <- c()
# # facilitation_modules <- c()
# # generalist_modules <- c()
# # 
# # 
# # attachment_modules <- clusters == 1
# # attachment_modules <- clusters[attachment_modules]
# # attachment_modules
# # 
# # selection_modules <- clusters == 2
# # selection_modules <- clusters[selection_modules]
# # selection_modules
# # transition_modules <- clusters == 3
# # transition_modules <- clusters[transition_modules]
# # 
# # facilitation_modules <- clusters == 4
# # facilitation_modules <- clusters[facilitation_modules]
# # facilitation_modules
# # 
# # generalist_modules <- clusters == 5
# # generalist_modules <- clusters[generalist_modules]
# # generalist_modules
# # 
# # all_modules <- c(attachment_modules, selection_modules, transition_modules, facilitation_modules, generalist_modules)
# # all_modules
# # modules_classified_df <- data.frame(names(all_modules),all_modules)
# # modules_classified_df <- modules_classified_df %>% select(all_modules)
# # colnames(modules_classified_df) <- "Module"
# # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 1,"attachment",modules_classified_df$Module)
# # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 2,"selection",modules_classified_df$Module)
# # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 3,"transition",modules_classified_df$Module)
# # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 4,"facilitation",modules_classified_df$Module)
# # modules_classified_df$Module <- ifelse(modules_classified_df$Module == 5,"generalist",modules_classified_df$Module)
# # modules_classified_df$Strategy <- modules_classified_df$Module
# # modules_classified_df$Module <- row.names(modules_classified_df)
# # modules_classified_df <- modules_classified_df %>% 
# #   filter(!(Module %in% c("attachment", "selection", "transition", "facilitation", "generalist")))
# # modules_classified_df
# # class(modules_classified_df)
# # getwd()
# # write_tsv(modules_classified_df,"modules_classified.tsv")
# 
# # Crear el data.frame a partir del vector
# 
# 
# # Calcula la distancia y realiza el clustering
# d <- dist.JSD(profiles)
# cluster <- hclust(d, method = "average")
# dendrogram <- as.dendrogram(cluster)
# 
# # Define los nombres de los módulos
# module_names <- c("attachment", "selection", "transition", "facilitation", "generalist")
# 
# # Función para verificar que cada módulo esté en un cluster único
# find_optimal_clusters <- function(dendrogram, module_names) {
#   num_clusters <- length(module_names) # Inicia con al menos el número de módulos
#   found_unique_clusters <- FALSE
#   
#   while (!found_unique_clusters) {
#     clusters <- cutree(dendrogram, k = num_clusters)
#     unique_modules <- unique(clusters[names(clusters) %in% module_names])
#     
#     # Verifica que cada módulo esté en un cluster diferente
#     if (length(unique_modules) == length(module_names)) {
#       found_unique_clusters <- TRUE
#     } else {
#       num_clusters <- num_clusters + 1
#     }
#   }
#   
#   return(num_clusters)
# }
# 
# # Encuentra el número óptimo de clusters
# number_of_clusters <- find_optimal_clusters(dendrogram, module_names)
# 
# # Aplica el número óptimo de clusters y asigna los módulos
# clusters <- cutree(dendrogram, k = number_of_clusters)
# 
# # Visualización del dendrograma
# dendrogram %>%
#   set("labels_cex", 1) %>% 
#   set("labels_col", value = 1:number_of_clusters, k = number_of_clusters) %>%
#   set("branches_k_color", k = number_of_clusters) %>%
#   plot(main = "Color labels \nper cluster")
# abline(h = 0.53, lty = 3)
# 
# # Asignación de módulos clasificados (similar al paso de clasificación previo)
# modules_classified_df <- data.frame(Module = character(), Strategy = character(), stringsAsFactors = FALSE)
# for (i in seq_along(module_names)) {
#   module <- module_names[i]
#   module_clusters <- clusters == i
#   cluster_modules <- clusters[module_clusters]
#   
#   if (length(cluster_modules) > 0) {
#     temp_df <- data.frame(Module = names(cluster_modules), Strategy = module, stringsAsFactors = FALSE)
#     modules_classified_df <- rbind(modules_classified_df, temp_df)
#   }
#   # Exporta el dataframe a un archivo .tsv
#   getwd()
#   modules_classified_df
#   write_tsv(modules_classified_df, "modules_classified.tsv")
# }
# 
# 
