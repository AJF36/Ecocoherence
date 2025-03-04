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
library(this.path)
#extra=0 # Fixing this parameter to 1 will generate additional plots and print file outputs

# Edit options -----------------
# --- Select substrate that will be analysed
# substrate ="Agarose"
# stop.step= 266

min.module=4
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")

for (substrate in substrates){
# minimum size of a functionInk module to be plotted

# --- Set working pathways
# dir.source="/home/ajf/Desktop/CNB/ecocoherence_sparcc"
# dir.funcInk=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc",substrate,sep="/")
# dirOut=paste(dir.funcInk,"figures",sep="/")
# dirOut
# dir.funcInk
# # --- Fix files
# fileOTU="/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv"
# fileSample="/home/ajf/Desktop/CNB/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
# fileFun=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp/Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep="")
work_dir<- this.dir()


dir.source=work_dir

home_dir <-  Sys.getenv("HOME")
dir.funcInk <- file.path(home_dir,"functionInk")
# dir.funcInk=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc",substrate,sep="/")
dir.funcInk <- file.path(work_dir,"..","functionink",substrate,"functionink_tmp")
dir.funcInk <- normalizePath(dir.funcInk)
# dirOut=paste(dir.funcInk,"figures",sep="/")
dirOut <- file.path(work_dir,"..","figures",substrate)
dirOut <- normalizePath(dirOut)
dir.funcInk
# --- Fix files
# fileOTU="/home/ajf/Desktop/CNB/marine_particles_source_data/count_table.ESV.4R.csv
fileOTU <- file.path(work_dir, "..", "data", "marine_particles_source_data", "count_table.ESV.4R.csv")
fileOTU <- normalizePath(fileOTU, mustWork = FALSE)
# fileSample="/home/ajf/Desktop/CNB/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
fileSample = file.path(work_dir, "..", "data", "marine_particles_source_data", "samples_properties","samples_metadata_AltSubstrSyntax.4R.tsv")
fileSample  <- normalizePath(fileSample, mustWork = FALSE)
# fileFun=paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,"/functionink_tmp/Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_0.01_0.7_",substrate,"_.tsv",sep="")
# fileFun <- file.path(work_dir,"..","functionink",substrate,"functionink_tmp","functionink_tmp",
#                      paste("Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_0.01_threshold_",substrate,"_.tsv",sep=""))
# fileFun <- normalizePath(fileFun)

subs_path <- file.path(work_dir,"..","functionink",substrate)
subs_path <- normalizePath(subs_path)
setwd(subs_path)
# setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
dir_path <- "functionink_tmp/"
file <- list.files(
  path = dir_path,
  pattern = paste("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_",substrate,"_.tsv$",sep=""),
  full.names = F
)
# fileFun <- paste("Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep="")
fileFun <- file
# paste("Partition-NL_Average_StopStep-",stop.step,"_interactions_filtered_p0.01_threshold_",substrate,"_.tsv",sep="") == "Partition-NL_Average_StopStep-608_interactions_filtered_p0.01_threshold_Agarose_.tsv"
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
length(partition$modules)
any(partition$modules == "none")
length(part.count)
sum(part.count > 4)
part.count
fileOut=paste(fileFun,"_guildGT",min.module,".txt",sep="")
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
partition.all
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

aggr.names
rownames(particles.select.aggr@otu_table)=aggr.names 
particles.select.aggr@otu_table
fileOut=paste("otu_table_",substrate,"_byModulesSize",min.module,".tsv",sep="")
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
plot(colSums(particles.rar.byRepTime.freq@otu_table)) # all should be 100
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
View(fake.tax)
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
plot(rowSums(particles.rar.mod.aggr.prop@otu_table)) # double check, should sum up to one

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




particles.rar.mod.aggr.prop@tax_table
# --- Compare all against all profiles
#source("/home/apascual/APASCUAL/Research/Programs/libraries/R/dist.JSD.R")


#---------------------------------------
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








#----------------------------------------
Ncomp=(Nmod+Nstage)
dim(profiles)
dim(t(profiles))
JSD=dist.JSD(profiles)


JSD.mat=as.matrix(JSD)
colnames(JSD.mat)=rownames(profiles)
rownames(JSD.mat)=colnames(JSD.mat)

# Plots ----------------
command=paste("mkdir",dirOut)
system(command)
setwd(dirOut)
# .... Create a palette for the plots
library(RColorBrewer) #
library(colorRamps)
library(viridis)
# .... Colors for barplot
Ncol=length(unique(partition.all))
#col_vector = primary.colors(Ncol)
qual_col_pals = brewer.pal.info[(brewer.pal.info$category == 'qual'),] # Get qualitative palettes
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                            rownames(qual_col_pals)))
col_vector[8]="red2" #col_vector[60] # and then I assigned to the two guys with grey positions
col_vector[7]="gold" 
#col_vector[59] # And to one of the pink ones, too much pink around
# col_vector[13]="Orange"
col_vector[16]="Green"
# col_vector[17]="Red"
col_vector[Ncol-2]="cyan"
col_vector[Ncol-1]="violet"
col_vector[Ncol]="White" # should be none

# .... Colors for heatmap
col_heat=inferno(50,direction = -1)

# --- Select plot analysis
select.plot.data=particles.rar.byRepTime.freq    #particles.rar.byRepTime.freq

# .... plot by replicate and time
Analysis=substrate
byAxes="ReplicatesAndTime-V2" # Add a label here describing the sample factor used to plot (e.g. MediaSubAndTime)
byTaxonomy="Genus"
plotTitle=paste("BarPlot_diversity4-",Analysis,"_By",byAxes,"_By",byTaxonomy,".pdf",sep="")
pdf(file=plotTitle,width = 12,height = 11)
title = Analysis

#####This is the part I had to correct

seq1 <- seq(0,84,12)
seq2 <- seq(108,156,24)
seq3 <- c(204,204,204)
final_times <-sort( c(seq1,seq1,seq1,seq2,seq2,seq2,seq3))
### Uncomment only for carrageenan since has 1 time less
if (substrate == "Carrageenan"){
  final_times <-final_times[-7] 
} else {
  final_times <- final_times
}

# final_times <-final_times[-7]
select.plot.data@sam_data$Time = as.numeric(final_times)
test = paste(select.plot.data@sam_data$Time,"",select.plot.data@sam_data$Replica,"")


select.plot.data@sam_data$RepTime = factor(test,levels = test)
select.plot.data@sam_data
sum(sample_metadata$Time == 24)
#The time is set to NA so that makes the problem in the plot , im going to fix it

p = plot_bar(select.plot.data, "RepTime", fill = "Modules", title = title)+ # coord_flip() + labs(colour = "family")
  scale_fill_manual(values = col_vector)+ 
  theme(axis.text.y = element_text(size=22),axis.title = element_text(size=26), 
        axis.text.x = element_text(size=19,vjust=0.4),#angle=60),
        legend.text = element_text(size=16,face="italic"),
        legend.title = element_text(size=26),legend.title.align =0.4,
        strip.text = element_text(size=24),title = element_text(size=22))+
  scale_y_reverse()+
  scale_x_discrete(labels=c("0A"="/","0B"="0","0C"="\\",
                            "12A"="/","12B"="12","12C"="\\","24A"="/","24B"="24","24C"="\\",
                            "36A"="/","36B"="36","36C"="\\","48A"="/","48B"="48","48C"="\\",
                            "60A"="/","60B"="60","60C"="\\","72A"="/","72B"="72","72C"="\\",
                            "84A"="/","84B"="84","84C"="\\","108A"="/","108B"="108","108C"="\\",
                            "132A"="/","132B"="132","132C"="\\","156A"="/","156B"="156","156C"="\\",
                            "204A"="/","204B"="204","204C"="\\")) # final choice


# ... Now include facetting
p=p+ ylab("Percentage of Sequences")+xlab("Time (h)/Replicate")+labs("Modules"="Guilds") # single substrate
print(p)
dev.off()

#--- Compare all against all profiles
work_dir
source(paste(work_dir,"heatmap.2.mod.R" , sep = "/"))
plotTitle=paste("heatmap_modulesVsecologicalPrefs4","_",substrate,".pdf",sep="")
pdf(file=plotTitle,width = 12,height = 12)
heatmap.2.mod(JSD.mat,
              density.info="none",
              trace="none",
              #scale="column",
              #xlab="",ylab="Genes",
              cexRow = 0.8,cexCol=1.2,#cex.lab=6,
              margins=c(10,10),#mgps=c(12,0.2,0),xmgp=12,ymgp=9,F
              #key.xlab = "Z-score",
              #key.title = "Jensen-Shannon dist.",
              #margins=c(8,7),
              keysize=1, # For Z-scores
              #key.xlab = "log(counts)", #Metagenomics
              symkey = FALSE,
              key.xlab = "Jensen-Shannon dist.",
              key.title = "", #key.title = "Function Values",
              key.par = list(cex.axis=1.2,cex.lab=1.5),
             col=col_heat)
dev.off()
}