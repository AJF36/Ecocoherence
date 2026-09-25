##################################################
# Untitled-1
##################################################
# Mia analysis for ecococoherence
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())


library(mia)
library(tidyverse)
library(this.path)
library(phyloseq)

source("../general_functions/modified_plotabundance.R")


script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "mia_analysis"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "mia_analysis"), mustWork = FALSE)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

## Este codigo es de Alberto
# --- Load OTUs
fileOTU <- file.path(script_dir, "..", "data", "marine_particles_source_data", "count_table.ESV.4R.csv")
fileOTU <- normalizePath(fileOTU, mustWork = FALSE)
print(fileOTU)
otu.in=read.csv(fileOTU)
otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)
fileOTU
# --- Load samples metadata
# fileSample="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
fileSample = file.path(script_dir, "..", "data", "marine_particles_source_data", "samples_properties","samples_metadata_AltSubstrSyntax.4R.tsv")
fileSample  <- normalizePath(fileSample, mustWork = FALSE)
sample_metadata = import_qiime_sample_data(fileSample)

# --- Load taxonomy and clean it

# I try with the manual upload
# fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
fileTaxonomy <- file.path(script_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
tax.pseq = tax_table(as.matrix(taxonomy))

##final del codigo de Alberto

phylo_obj <- merge_phyloseq(sample_metadata, tax.pseq, otu.pseq)
phylo_obj_B <- prune_samples(sample_data(phylo_obj)$Media == "Beads", phylo_obj)

# Conver to mia
mia_main_b = mia::convertFromPhyloseq(phylo_obj_B)

# Apply the same filter
samples_to_filter = colSums(assay(mia_main_b,"counts")) > 5000 # Samples with at least 5000 reads
otu_to_filter = rowSums(assay(mia_main_b,"counts")) > 10 # ASVs that have at least 10 counts
otu_to_filter
mia_main_f = mia_main_b[otu_to_filter,samples_to_filter]

# Now prune the asvs that are not present in at least 10 samples

mia_main_f = mia::transformAssay(mia_main_f,method = "pa")
pa_matrix = assay(mia_main_f,"pa") # Addd a precense/abscense matrix
otu_to_filter2 = rowSums(pa_matrix) > 10 # ASVs that are present in more than 10 samples

mia_main_final = mia_main_f[otu_to_filter2,]

tax_table = rowData(mia_main_final)

tax_table = tax_table |> 
  as.data.frame() |> 
  rownames_to_column("ASV")

rowData(mia_main_final) = tax_table

# Drop the blank trailing metadata columns from the QIIME import (empty
# header cells get auto-named "X"/"X.1"). miaViz::plotAbundance reserves
# the colname "X" for its own internal plotting data and errors if it's
# already present in colData.
colData(mia_main_final) = colData(mia_main_final)[setdiff(colnames(colData(mia_main_final)), c("X", "X.1"))]
colData(mia_main_final)[["Time"]] = as.numeric(colData(mia_main_final)[["Time"]])
# Now since every substrate has its own network/functionink clustering, split by substrate

saveRDS(mia_main_final,file ="analysis/results/mia_analysis/mia_main_final.RDS")

mia_substrate_list = mia::splitOn(mia_main_final,"Substrate")

functionink_dirs = fs::dir_ls("functionink/")

dirs = fs::dir_ls(functionink_dirs[[1]],recurse = T)

dirs[grepl("Partition-.*guildGT4.*",dirs,perl = T)]

load_functionink = function(path) {

  dirs = fs::dir_ls(path,recurse = T)
  funk = readr::read_tsv(dirs[grepl("Partition-.*guildGT4.*",dirs,perl = T)] , skip = 1 , col_names = c("ASV","module"))
  return(funk)
}

names(functionink_dirs) = basename(functionink_dirs) # Charge the names to map it 

funk_substrates = map(functionink_dirs,load_functionink)

imap(mia_substrate_list,function(mia,substrate) {
  glue::glue("Processing {substrate}")
  # Get their functionink table
  funk_substrate = funk_substrates[[substrate]]

  if (is.null(funk_substrate)) {
    warning(glue::glue("No functionink data for substrate '{substrate}', skipping"))
    return(NULL)
  }

  tax = rowData(mia)
  
  full_tax = dplyr::left_join(as.data.frame(tax),funk_substrate,by = "ASV")
  
  # print(full_tax)

  rowData(mia) = full_tax 
  
  plot = modified_plotabundance(mia,variable_to_aggregate = "module", order.col.by = "Time") +
    ggtitle(label = glue::glue("Module barplot: {substrate}"))
  ggsave(filename = glue::glue("analysis/figures/mia_analysis/barplot_modules_{substrate}.pdf"))
}
)

matched_asvs = readr::read_tsv("metabolism/genome_aligment/matched_ESV_id_0.97.tsv",col_names = F)
colnames(matched_asvs) <- c(
  "record_type", "cluster", "length", "pct_identity", "strand",
  "field6", "field7", "cigar", "ASV", "ID"
)

# Selected modules for misosoup analysis

selected_modules_misosoup = list("Alginate" = c("mod_4","mod_6"),
"Agarose" = c("mod_1","mod_48"), 
"Chitin" = c("mod_2","mod_9"))

# Move the metabolic models and run misosoup using all the models of each substrate
funk_substrates[["Alginate"]]
# Prepare the tables for the general_functions/move_metabolic_models.sh
purrr::iwalk(selected_modules_misosoup,function(modules,substrate) {
  print(glue::glue("Processing modules {modules} for {substrate}"))
  funk_substrate = funk_substrates[[substrate]]
  
  funk_modules = funk_substrate |> 
    filter(module %in% modules )
  
  # join with the matched ASVs
  joined_table = dplyr::inner_join(funk_modules,matched_asvs,by = "ASV")

  final_table = joined_table |> 
    filter(ID != "*") |> 
    select("ASV","ID")

  write_tsv(final_table,file = glue::glue("analysis/results/mia_analysis/metabolic_models_{substrate}.tsv"))

})

# Run general_scripts/move_metabolic_models.sh
substrates = c("Agarose","Alginate","Chitin")

walk(substrates, function(substrate) {

models_folder = glue::glue("metabolism/carveme_smetana/{substrate}/metabolic_models/")
community_file = glue::glue("analysis/results/mia_analysis/metabolic_models_{substrate}.tsv")
target_folder = glue::glue("analysis/results/mia_analysis/misosoup_analysis_{substrate}")

system2(command = "../general_scripts/move_metabolic_models.sh" , args = c(
  models_folder,
  community_file,
  target_folder
))

})


# Run general_scripts/run_misosoup.sh

run_misosoup_script = normalizePath(file.path("..", "general_scripts", "run_misosoup.sh"))
media_file = normalizePath(file.path("metabolism", "misosoup", "media.yaml"))

misosoup_media_dict = list("Agarose" = "gal", "Alginate" = "alg" , "Chitin" = "chitin")

walk(substrates , function(substrate) {

  models_folder = normalizePath(glue::glue("analysis/results/mia_analysis/misosoup_analysis_{substrate}"))

  misosoup_media_select = misosoup_media_dict[[substrate]]
  # move_metabolic_models.sh copies models in as "<ID>.xml" (e.g.
  # "RS_GCF_...xml"), not "*genomic.xml", so run_misosoup.sh's own
  # auto-discovery won't see them — pass the list explicitly instead.
  metabolic_models = fs::dir_ls(models_folder, glob = "*.xml")

  if (length(metabolic_models) == 0) {
    warning(glue::glue("No metabolic models found for substrate '{substrate}', skipping misosoup"))
    return(NULL)
  }

  old_wd = getwd()
  on.exit(setwd(old_wd), add = TRUE)
  # run_misosoup.sh must be run from inside the models folder: it names
  # its output/log from basename(getwd()), so cd into it instead of
  # passing it as an arg.
  setwd(models_folder)

  system2(command = run_misosoup_script , args = c(".", media_file,
   "--models", metabolic_models,
  "--media-select",misosoup_media_select))

})