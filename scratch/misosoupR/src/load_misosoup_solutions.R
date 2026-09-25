# ****************************************
# load_misosoup_solutions.R
# ****************************************
# 
# 
# author = Alberto Pascual-García
# contact = apascualgarcia.github.io
# date = June 2021 
# description = This script loads a set  of models in yaml format describing 
#   communities coexisting in defined media. For each model, it first relabel the
#   fluxes as secretions or consumptions. Then it inspects all models and it
#   computes the number of times each flux appears across all comunities, and
#   the mean value of the flux (the mean among the models in which it appears)
# usage = tbf
# 
rm(list=ls())
library(yaml)
library(reshape)

######## START EDITING

# ... set vars for input solutions file (it has the form C_src.ext)
C_src="ac"
ext=".yaml"
strain="A1R12"

# ... set name of minimal medium file
fileMedium="medium_MBM_no_co2_hco3.yaml"

######## STOP EDITING

# Set up environment -------------
# --- Set directories
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/scripts/")[[1]][1] # don't edit, just comment it if problems...

dirSrc=paste(this.dir,"/scripts/",sep="") # Directory where the code is
dirModels=paste(this.dir,"/data/misosoup_results/cache/",strain,sep="") # directory for the models
dirOut=paste(this.dir,"/results/networks/",sep="") # directory for the output
dirMedium=paste(this.dir,"/data/",sep="") # directory for the minimal media

# --- Load some additional functions
setwd(dirSrc)
source("heatmap.2.mod.R")

# --- Load minimal medium
setwd(dirMedium)
met_medium = read_yaml(file = fileMedium)

# --- Load MisoSoup output
setwd(dirModels)
fileIn=paste(C_src,ext,sep="")
# ... This is a nested list with structure sols$focalStrain$C_source[[solutions]]
sols = read_yaml(file = fileIn)
# ** CREATE FIRST FUNCTION UP TO HERE, JUST LOADING SOLUTIONS ** 

# Start computation -------
Nsols=length(sols[[1]][[1]]) # count how many solutions 

for(i in 1:Nsols){
  sol_tmp=unlist(sols[[1]][[1]][[i]]) # convert first solution in a named vector
  #idx_medium=grep("_e$",names(sol_tmp)) # identify metabolites not produced by the bugs
  #met_medium=names(sol_tmp[idx_medium]) # are either input metabolites or final byproducts like h2o
  idx_ex=grep("_e_",names(sol_tmp)) # identify exchanged reactions
  sol_tmp=sol_tmp[idx_ex] # select only those
  idx_ex=c() # among reactions
  for(met in met_medium){ # exclude those involving metab present in the base medium
    idx_ex_tmp=grep(met,names(sol_tmp))
    idx_ex=c(idx_ex,idx_ex_tmp)
  }
  matched=match(seq(1,length(sol_tmp)),idx_ex) # take a list of indexes of sol_tmp and check which are not in idx_ex
  sol_tmp=sol_tmp[is.na(matched)] # we are looking for those not found 
  #names(sol_tmp) # double check
  idx_scr=which(sol_tmp > 0) # identify secretions
  idx_con=which(sol_tmp < 0) # identify consumption
  names(sol_tmp)[idx_scr]=paste(names(sol_tmp)[idx_scr],"_t_SCR",sep="") # rename secretion
  names(sol_tmp)[idx_con]=paste(names(sol_tmp)[idx_con],"_t_CON",sep="") # rename consumption
  if(i == 1){
    sol_df=as.data.frame(sol_tmp)
    colnames(sol_df)="flux"
    count_df=sol_df
    colnames(count_df)="count"
    count_df[count_df !=0]=1
  }else{
    # --- First identify existing and new solutions
    matched=match(names(sol_tmp),rownames(sol_df)) # identify those sols already in the df
    found=which(!is.na(matched)) # these are the positions where there is a sol in the df 
    # --- Add the fluxes
    sol_df[matched[found],]=sol_df[matched[found],]+sol_tmp[found] # add the fluxes
    sol_tmp_new=as.data.frame(sol_tmp[is.na(matched)]) # finally incorporate new entries
    colnames(sol_tmp_new)="flux"
    sol_df=rbind(sol_df,sol_tmp_new)
    # .... repeat for counts
    count_tmp=sign(sol_tmp)
    #count_tmp[count_tmp != 0]=1
    count_df[matched[found],]=count_df[matched[found],]+count_tmp[found] # add the fluxes
    count_tmp_new=as.data.frame(count_tmp[is.na(matched)]) # finally incorporate new entries
    colnames(count_tmp_new)="count"
    count_df=rbind(count_df,count_tmp_new)
  }
}
# .... compute means and rescale
sol_df$flux=sol_df$flux/Nsols
#quantile(sol_df$flux) # check the distribution
#sol_df$flux=log(abs(sol_df$flux)) # since we keep the id as con or scr we can remove the sign
#quantile(count_df$count)
# ** CREATE SECOND FUNCTION UP TO HERE, CREATES SOL_DF AND COL_DF ** 

# --- clean the names
names_tmp1=t(as.data.frame(strsplit(rownames(sol_df),split = "_e_")))
names_tmp2=t(as.data.frame(strsplit(names_tmp1[,2],split = "_i_t_")))
names_tmp3=t(as.data.frame(strsplit(names_tmp1[,1],split = "_EX_")))
names_final_undir=cbind(paste(names_tmp3[,2],names_tmp2[,2],sep="_"),names_tmp2[,1])
names_final_dir=cbind(names_tmp3[,2],names_tmp2[,1])

# --- create an undirected network (metabolites appear with label CON and SCR)
sol_df_undir=cbind(sol_df,names_final_undir)
colnames(sol_df_undir)=c("flux","metab","strain")
count_df_undir=cbind(count_df,names_final_undir)
colnames(count_df_undir)=c("count","metab","strain")

# --- create a directed network (metabolites appear without label, so direction is needed)
sol_df_dir=cbind(sol_df,names_final_dir,sign(sol_df$flux))
colnames(sol_df_dir)=c("flux","metab","strain","type")
count_df_dir=cbind(count_df,names_final_dir,sign(count_df$count))
colnames(count_df_dir)=c("count","metab","strain","type")


# --- Create a metadata table and write outputs
setwd(dirOut)
# .... first the undirected network
metabs=unique(count_df_undir$metab); strains=unique(count_df_undir$strain)
meta.metab=data.frame(metabs,"metab","null",stringsAsFactors = FALSE)
idx_con=grep("_CON",metabs)
idx_scr=grep("_SCR",metabs)
meta.metab[idx_con,3]="CON"
meta.metab[idx_scr,3]="SCR"
colnames(meta.metab)=c("node","entity","direction")
meta.strain=data.frame(strains,"strain","none")
colnames(meta.strain)=c("node","entity","direction")
metadata_undir=rbind(meta.metab,meta.strain)

fileOut=paste("Metadata_",strain,"_",C_src,"_undir.txt",sep="")
write.table(metadata_undir,file=fileOut,sep="\t",quote = FALSE,row.names = FALSE)

# ... then the directed network
metabs=unique(count_df_dir$metab); strains=unique(count_df_dir$strain)
meta.metab=data.frame(metabs,"metab",stringsAsFactors = FALSE)
colnames(meta.metab)=c("node","entity")
meta.strain=data.frame(strains,"strain")
colnames(meta.strain)=c("node","entity")
metadata_dir=rbind(meta.metab,meta.strain)

fileOut=paste("Metadata_",strain,"_",C_src,"_dir.txt",sep="")
write.table(metadata_dir,file=fileOut,sep="\t",quote = FALSE,row.names = FALSE)

# --- For the undirected network create links between the same metabolite consumed and secreted
metabs=unique(count_df_dir$metab);

k=0
for(met in metabs){
  idx_met=grep(met,metadata_undir$node)
  if(length(idx_met)==2){
    k=k+1
    metA_metB_tmp=data.frame(0.01,metadata_undir$node[idx_met[1]],
                               metadata_undir$node[idx_met[2]]) # a link with a small value to have little influence in functionink
    colnames(metA_metB_tmp)=c("value","nodeA","nodeB")
    if(k ==1){
      metA_metB=metA_metB_tmp
    }else{
      metA_metB=rbind(metA_metB,metA_metB_tmp)                         
    }
  }
}
colnames(sol_df_undir)=c("value","nodeA","nodeB")
sol_df_undir=rbind(sol_df_undir,metA_metB)
colnames(count_df_undir)=c("value","nodeA","nodeB")
count_df_undir=rbind(count_df_undir,metA_metB)

# --- Print the networks
fileOut=paste("Network_",strain,"_",C_src,"_undir_byFlux.txt",sep="")
write.table(sol_df_undir[,c(2,3,1)],file=fileOut,sep="\t",quote = FALSE,row.names = FALSE)
fileOut=paste("Network_",strain,"_",C_src,"_undir_byCount.txt",sep="")
write.table(count_df_undir[,c(2,3,1)],file=fileOut,sep="\t",quote = FALSE,row.names = FALSE)
fileOut=paste("Network_",strain,"_",C_src,"_dir_byFlux.txt",sep="")
write.table(sol_df_dir[,c(2,3,1,4)],file=fileOut,sep="\t",quote = FALSE,row.names = FALSE)
fileOut=paste("Network_",strain,"_",C_src,"_dir_byCount.txt",sep="")
write.table(count_df_dir[,c(2,3,1,4)],file=fileOut,sep="\t",quote = FALSE,row.names = FALSE)

# Plot -------
# --- Reshape to plot
sol_mat=reshape(sol_df_undir,idvar="nodeA",timevar="nodeB",direction="wide")
rownames(sol_mat)=sol_mat$nodeA
sol_mat=as.matrix(subset(sol_mat,select = -c(nodeA)))
sol_mat[is.na(sol_mat)]=0

# ..... same for counts
count_mat=reshape(count_df_undir,idvar="nodeA",timevar="nodeB",direction="wide")
rownames(count_mat)=count_mat$nodeA
count_mat=as.matrix(subset(count_mat,select = -c(nodeA)))
count_mat[is.na(count_mat)]=0


plotVec=c("flux","count")
for(plotType in plotVec){
  if(plotType == "flux"){
    mat_out=t(sol_mat)
  }else if(plotType == "count"){
    mat_out=t(count_mat)
  }
  fileOut=paste("heatmap_",strain,"_",C_src,"_",plotType,".pdf",sep="")
  pdf(fileOut,width=30,height=15)
  heatmap.2.mod(mat_out,
                trace = "none",
                main = C_src,
                cexRow = 1.5,cexCol = 1.5,
                margins = c(17,15),
                keysize = 0.5,
                xlab = "strains",ylab="metabolites")#,
  #col="bluered")
  dev.off()
}

