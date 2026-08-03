rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
# #Ahora vamos con fast_spar------------------------------------------------------------------------------------
setwd("/home/ajf/Desktop/CNB/ecocoherence/paper_data")
# 
# ##Fixinf the format
#  system(paste("sed -i '1s/^/#OTU_ID\t/'","data_filtered_for_sparcc_paper_data.tsv",sep = " "))
#  #Calculo de las correlaciones y de las covariancias
#  command <- paste("fastspar --threshold 0.3 --otu_table ","data_filtered_for_sparcc_paper_data.tsv"," --correlation median_correlation.tsv --covariance median_covariance.tsv",sep = "")
# 
#  system(command)
# 
#  #bootstraping
#  system("mkdir bootstrap_counts")
#  y <- paste("fastspar_bootstrap --otu_table ","data_filtered_for_sparcc_paper_data.tsv"," --number 1000"," --prefix bootstrap_counts/",sep = "")
# 
#  system(y)
#  #infer correlations of the boostrap counts
#  system("mkdir bootstrap_correlation")
#  system("parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*")


##HAD to delete the cov files
setwd("/home/ajf/Desktop/CNB/ecocoherence/paper_data/bootstrap_correlation")
system("rm cov*")

setwd("/home/ajf/Desktop/CNB/ecocoherence/paper_data")

z <- paste("fastspar_pvalues --otu_table ", "data_filtered_for_sparcc_paper_data.tsv",
           " --correlation median_correlation.tsv --prefix bootstrap_correlation/",
           " --permutations 1000 --outfile pvalues.tsv", sep="")
system(z)
#fin fastsparcc--------------------------------------------------------

p_values <- read_tsv("pvalues.tsv")

p_values <- as.data.frame(p_values)

row.names(p_values) = colnames(p_values)[2:length(colnames(p_values))]
p_values <- p_values %>% select(-c("#OTU ID"))


cor <- read.table("median_correlation.tsv", header = FALSE, sep = "\t", row.names = 1)
colnames(cor) = rownames(cor)

df = data.frame()
for (i in 1:(nrow(cor)-1)){
  for (j in (i+1):ncol(cor))
  {
    if (p_values[i,j] <= 0.01)
    {
      df[nrow(df) +1 ,"Sp1"] = rownames(cor)[i]
      df[nrow(df),"Sp2"] = rownames(cor)[j]
      df[nrow(df),"Cor"] = cor[i,j]
    }
  }
}
# 
write.table(df,"interactions_filtered_0.01.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

###Formating correctly the file for functionink
functionink.input <- read_tsv("interactions_filtered_0.01.tsv")
colnames(functionink.input) <- c("#SpeciesA","SpeciesB","interaction")
functionink.input.type <- functionink.input %>%
  mutate("Type" = ifelse(interaction > 0 ,1,0))
# 
write.table(functionink.input.type,"interactions_filtered_0.01_functioink.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
write_tsv(functionink.input.type,"interactions_filtered_0.01_functioink2.tsv")
write.csv(functionink.input.type,
          "interactions_filtered_0.01_functioink.csv",
          row.names = FALSE)



## Functionink
  fileNet <- "interactions_filtered_0.01_functioink.tsv"
  pathNet <- "/home/ajf/Desktop/CNB/ecocoherence/paper_data/"
  pathRepo <- "/home/ajf/functionInk"
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

  run_pipeline(fileNet, pathNet, pathRepo, weighted = TRUE,types = TRUE)


# 
