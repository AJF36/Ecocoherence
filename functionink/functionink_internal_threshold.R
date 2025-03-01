#################################################################################
#
#For this script to work functioink must be installed in the home directory
#
#################################################################################
rm(list=ls())
library(readr)
library(dplyr)
library(tidyr)
library(this.path)
library(fs)


# install.packages("fs")

substrates <- c("Agarose")

run_pipeline = function(fileNet,pathNet,pathRepo, # mandatory
                        weighted=FALSE,directed=FALSE,
                        types=FALSE,method="Average",mode="all"){
  # ... set up the environment and source dependencies
  setwd(pathRepo)
  src.dir=paste("scripts","analysis_R",sep="/")
  setwd(src.dir)
  source("extractPartDensity.R")
  setwd(pathRepo)
  setwd(pathNet)
  dir.create("functionink_tmp")
  setwd("functionink_tmp")
  fileNetPath=paste0("../",fileNet)
  
  # ... process arguments
  if(weighted == FALSE){par_w = 0}else{par_w = 1}
  if(directed == FALSE){par_d = 0}else{par_d = 1}
  if(types == FALSE){par_t = 0}else{par_t = 1}
  
  # ... build commands basic run
  # ...... Node similarity
  script=paste(pathRepo,"NodeSimilarity.pl",sep="/")
  options=paste("-w",par_w,"-d",par_d,"-t",par_t,"-f",fileNetPath)
  comm_sim=paste(script,options)
  # ..... Node linkage
  fileSim=paste0("Nodes-Similarities_",fileNet)
  script=paste(pathRepo,"NodeLinkage.pl",sep="/")
  options=paste("-fn",fileNetPath,"-fs",fileSim,"-a",method)
  comm_link_base=paste(script,options)
  file.hist=paste0("HistCompact-NL_",method,"_NoStop_",fileNet) # expected history file output
  
  # --- Run the first analysis
  system(comm_sim) # compute nodes' similarities
  system(comm_link_base) # cluster nodes
  hist.comp=read.table(file=file.hist,sep="\t",header = TRUE) # read history file
  part_density=extractPartDensity(hist.comp) # extract partition densities
  
  # --- Run the extraction of the communities
  if(mode != "none"){ # if the user wants to retrieve the communities
    # .... determine the criteria to  be used
    labels=c("total_dens_step","int_dens_step","ext_dens_step")
    if(mode == "total"){
      idx=1
    }else if(mode == "internal"){
      idx=2
    }else if(mode == "external"){
      idx=3
    }else{ # mode "all", we identify the maximum among modes
      idx=which.max(c(part_density$total_dens,
                      part_density$int_dens,
                      part_density$ext_dens))
    }
    # ... find the step of the peak and create a new command
    value=part_density[labels[idx]]
    comm_link_spec=paste(comm_link_base,"-s step -v",value)
    # ... finally run
    system(comm_link_spec) # cluster nodes and extract partition
  }
  return(part_density)
}

home_dir <-  Sys.getenv("HOME")
work_dir <- this.dir()
path_to_input <- file.path(work_dir,"..","filter_network_by_threshold")
path_to_input <- normalizePath(path_to_input)
setwd(work_dir)
for (substrate in substrates){
  system(paste("mkdir ",substrate))
  #Por ultimo functionink
  fileNet <- paste("interactions_filtered_p0.01","threshold",substrate,".tsv", sep = "_")
  # pathNet <- paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate, sep ="")
  pathNet <- path_to_input
  pathRepo <- file.path(home_dir,"functionInk")
  run_pipeline(fileNet = fileNet, pathRepo = pathRepo, pathNet = pathNet, weighted = TRUE, types = TRUE, mode = "internal")
  
  ###moving the directories
  destiny <- file.path(work_dir,substrate)
  origin <- file.path(work_dir,"..","filter_network_by_threshold","functionink_tmp")
  origin <- normalizePath(origin)
  file_move(origin,destiny)
}



for(substrate in substrates){
 
}
