rm(list = ls())

library(this.path)
library(readr)
library(dplyr)
library(tidyr)
library(phyloseq)

script_dir <- this.dir()
data_dir <- normalizePath(file.path(script_dir, "..", "..", "paper_data"), mustWork = FALSE)
results_dir <- normalizePath(file.path(script_dir, "..", "..", "results", "paper_data", "1_filter_asv_table"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setwd(data_dir)

### Load counts
file.counts <- "counts.rds"
counts.table <- readRDS(file.counts)


### Load taxonomy
file.tax <- "taxonomy.rds"
tax.table <- readRDS(file.tax)


### Load metadata
file.metadata <- "metadata.csv"
metadata.table <- read_csv(file.metadata)



#### To do list
# Filter samples
# check metadata
# Make the same filtering analisys as in the original


##Filter samples
# We only need the samples from the UE east coast

metadata.table.f <- filter(metadata.table, geo_region == "East Coast U.S")
samples.to.keep <- metadata.table.f$sample_name

#Filter the samples from the OTU table

counts.table.f <- counts.table[rownames(counts.table) %in% samples.to.keep, ]

## Check metadata

metadata.columns <- colnames(metadata.table.f)
test3 <- as.vector(metadata.table.f)
metadata.columns
metadata.list.unique <- list()
x <- 1
length(test3)
for (col in test3){
  col.unique <- unique(test3[[x]])
  metadata.list.unique[[x]] <- col.unique
  x <- x + 1
}

names(metadata.list.unique) <-metadata.columns

## Same filtering analysis as in the original

#make the phyloseq objects

tax.phylo <- tax_table(as.matrix(tax.table))
counts.phylo <- otu_table(t(counts.table.f),taxa_are_rows = T)
View(counts.phylo)
metadata.phylo <- sample_data(metadata.table.f)
rownames(metadata.phylo) <- metadata.phylo$sample_name

phylo.obj <- merge_phyloseq(tax.phylo,counts.phylo,metadata.phylo)

#Now we take the samples that have at least 5000 reads
phylo.obj
min(sample_sums(phylo.obj))
sum(sample_sums(phylo.obj) > 5000)

data <- prune_samples(sample_sums(phylo.obj) > 5000, phylo.obj)
data
min(sample_sums(data))

#We take the otus that dont have at least 100 reads
min(taxa_sums(data))

data <- prune_taxa(taxa_sums(data) > 100,data)
data
min(taxa_sums(data))

#Filter the ASV that are not at least in 10 samples
#---------------------------------------------------------------------------------------------------------
# Paso 1: Obtener la tabla de OTUs (abundancias)
otu_table_data <- otu_table(data)

# Paso 2: Contar cuántas muestras tienen presencia de cada OTU (valores > 0)
otu_presence <- apply(otu_table_data, 1, function(x) sum(x > 0))
?apply
# Paso 3: Filtrar los OTUs presentes en al menos 10 muestras
otus_to_keep <- names(otu_presence[otu_presence >= 10])

# Paso 4: Usar prune_taxa() para quedarte solo con los OTUs filtrados
data2<- prune_taxa(otus_to_keep, data)


otu_data_frame <- as.data.frame(otu_table(data2,taxa_are_rows = TRUE))

View(otu_data_frame)

median(taxa_sums(data2))

# write.table(otu_data_frame, file = "data_filtered_for_sparcc_paper_data.tsv", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t" )

### There is some OTUs that cause a problem in fastSpar, lets filter them

otus.to.filter <- c(180, 579,917,957,995,1002,1237,1308,1339,
                    1389,1721,1794,1809,1845,1936,1937,1938,1939,1940)
df_problematic_otus <- otu_table(otu_data_frame[otus.to.filter,],taxa_are_rows = T)
taxa_sums(df_problematic_otus)
otu_data_frame_f <- otu_data_frame[-otus.to.filter,]


write.table(otu_data_frame_f, file = file.path(results_dir, "data_filtered_for_sparcc_paper_data.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t" )
