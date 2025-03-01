rm(list = ls())

library(readr)
library(this.path)
library(dplyr)
library(tidyr)

substrates <- c("Agarose")

work_dir <- this.dir()
setwd(work_dir)
###Load the thresholds file
thres_file <- file.path(work_dir,"..","Threshold_election","minimun_correlations_threshold.tsv")
thres_file <- normalizePath(thres_file)

thres_data <- read_tsv(thres_file)
rownames(thres_data)  <- thres_data$substrates

for (substrate in substrates){
  net_file <- file.path(work_dir,"..","network_creation_sparcc",substrate,paste("interactions_filtered_0.01._",substrate,"_.tsv", sep = ""))
  net_file <- paste(normalizePath(net_file))
  substrate_net <- read_tsv(net_file)
  threshold <- thres_data$minimun_values[match(substrate,thres_data$substrates)]
  net_file_f <- filter(substrate_net,abs(Cor) >= threshold)
  net_file_f$type <- ifelse(net_file_f$Cor > 0,0,1)
  write_tsv(net_file_f,paste("interactions_filtered_p0.01","threshold",substrate,".tsv", sep = "_"))
}

