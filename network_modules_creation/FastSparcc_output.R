rm(list = ls())

library(dplyr)

# args = commandArgs(trailingOnly = TRUE)
# 
# out = args[1]

###you may have to make some changes to dhis code, specially in the rownames and colnames part
setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc/tara_oceans_data/OM-RGC_v2_taxonomic_profiles")
# pvalues <- read.table("pvalues.tsv", header = FALSE, sep = "\t", row.names = 1)
pvalues <- read_tsv("pvalues.tsv")
rownames(pvalues) = colnames(pvalues)[-1]
# cor <- read.table("median_correlation.tsv", header = FALSE, sep = "\t", row.names = 1)
cor <- read_tsv("median_correlation.tsv")
rownames(cor) = colnames(cor)[-1]

df = data.frame()
for (i in 1:(nrow(cor)-1))
  for (j in (i+1):ncol(cor))
  {
    if (pvalues[i,j] <= 0.05)
    {
      df[nrow(df) +1 ,"Sp1"] = rownames(cor)[i]
      df[nrow(df),"Sp2"] = rownames(cor)[j]
      df[nrow(df),"Cor"] = cor[i,j]
    }
  }

##Remove the columns that have NA as family because is causing a problem in functionink
df$Type <- ifelse(df$Cor > 0 ,0,1)
df_filtered <- filter(df, SpeciesB != "NA")
colnames(df_filtered) <- c("#SpeciesA","SpeciesB","interaction","type")
write.table(df_filtered, "tara_oceans_network.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
