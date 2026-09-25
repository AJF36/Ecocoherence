##################################################
# tara_oceans_analysis.R
##################################################
# Analysis of the tara oceans data from Salazar G., Paoli L., et al. (2019)
# to see if we find similar modules to the ones we have in our dataset
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())

library(mia)
library(tidyverse)
library(readxl)
library(this.path)

script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "tara_oceans_analysis"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "tara_oceans_analysis"), mustWork = FALSE)
data_dir <- normalizePath(file.path(script_dir, "..", "data", "tara_oceans"), mustWork = FALSE)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# mitags_tab_otu.tsv.gz has ~24000 columns (one per OTU); readr's default
# connection buffer (128KB) is too small to buffer a single line of it
Sys.setenv(VROOM_CONNECTION_SIZE = 5 * 1024 * 1024)

metadata_tara = suppressMessages(readxl::read_excel(file.path(data_dir, "Salazar_et_al_2019_Suppl_Info.xlsx"), sheet = "Table_W1"))
otu_raw = suppressMessages(read_tsv(file.path(data_dir, "OM-RGC_v2_taxonomic_profiles", "mitags_tab_otu.tsv.gz")))

# --- Reshape OTU table: samples-as-rows/OTUs-as-columns -> OTUs-as-rows/
# samples-as-columns, and split the taxonomy lineage out of each OTU's
# column header (format: "<accession> Domain;Phylum;Class;Order;Family;Genus[;more]")

sample_col = colnames(otu_raw)[1]
otu_raw = rename(otu_raw, PANGAEA_sample = !!sample_col)

col_headers = colnames(otu_raw)[-1]
lineage_str = str_trim(str_remove(col_headers, "^\\S+"))
# This is an NCBI-style variable-depth lineage (not the marine dataset's
# fixed-rank RDP taxonomy): Domain..Genus are assigned positionally from
# the first 6 ";"-separated tokens, an approximation that misassigns rank
# names when a lineage collapses to fewer/more levels. The untouched
# lineage string is kept in `Lineage` for reference.
lineage_split = str_split_fixed(lineage_str, ";", 20)[, 1:6]
colnames(lineage_split) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
lineage_split[lineage_split == ""] = NA

asv_ids = paste("ASV", seq_along(col_headers), sep = "_")

tax_tara = as.data.frame(lineage_split)
tax_tara$OTU_id = str_extract(col_headers, "^\\S+")
tax_tara$Lineage = lineage_str
rownames(tax_tara) = asv_ids

count_table = otu_raw |> select(-PANGAEA_sample) |> as.matrix() |> t()
rownames(count_table) = asv_ids
colnames(count_table) = otu_raw$PANGAEA_sample

# --- Simplify metadata labels ---
metadata_tara = metadata_tara |>
  # Table_W1 has one row per PANGAEA sample id per sequencing type
  # (MetaG/MetaT); the miTAG OTU table is metagenomic, so keep MetaG rows
  # only, otherwise the join key below is duplicated
  filter(`MetaG/MetaT` == "MetaG") |>
  rename(
    Sample = `PANGAEA sample id`,
    DateTime = `Date/Time`,
    Depth_m = `Depth, nominal`,
    Ocean_region = `OS region`
  ) |>
  as.data.frame()

rownames(metadata_tara) = metadata_tara$Sample

# Only samples with miTAG profiles are usable for the mia object
metadata_tara = metadata_tara[colnames(count_table), ]

# --- Create mia object ---

tse_tara = TreeSummarizedExperiment(
  assays = list(counts = as.matrix(count_table)),
  rowData = tax_tara,
  colData = metadata_tara
)

# --- Keep only samples that come from the Atlantic ocean ---

tse_tara_atlantic = tse_tara[, grepl("Atlantic", colData(tse_tara)$Ocean_region)]

# --- Apply the same prevalence/abundance filters used for the marine dataset ---
# (see src/1_filter_asv_table_and_build_sparcc_network.R): samples with >5000
# reads, ASVs with >10 total reads, ASVs present in >=10 samples

samples_to_keep = colSums(assay(tse_tara_atlantic, "counts")) > 5000
tse_f = tse_tara_atlantic[, samples_to_keep]

asv_to_keep = rowSums(assay(tse_f, "counts")) > 10
tse_f = tse_f[asv_to_keep, ]

tse_f = mia::transformAssay(tse_f, method = "pa")
pa_matrix = assay(tse_f, "pa")
asv_to_keep2 = rowSums(pa_matrix) >= 10
tse_atlantic_final = tse_f[asv_to_keep2, ]

glue::glue("Atlantic subset after filtering: {nrow(tse_atlantic_final)} ASVs x {ncol(tse_atlantic_final)} samples")

# --- Write the OTU table in BIOM-compatible format for fastspar ---

pipeline_dir = file.path(results_dir, "atlantic")
dir.create(pipeline_dir, recursive = TRUE, showWarnings = FALSE)

otu_out = as.data.frame(assay(tse_atlantic_final, "counts")) |>
  tibble::rownames_to_column(var = "#OTU_ID")
otu_table_file = file.path(pipeline_dir, "raw_counts_asv_f.tsv")
write_tsv(otu_out, file = otu_table_file)

# --- Run fastspar + functionInk (general_scripts/fastsparcc_functionink_pipeline.sh) ---

pipeline_script = normalizePath(file.path(script_dir, "..", "..", "general_scripts", "fastsparcc_functionink_pipeline.sh"), mustWork = FALSE)

old_wd = getwd()
setwd(pipeline_dir)

system2("bash", args = c(pipeline_script,
  "--input", "raw_counts_asv_f.tsv",
  "--threshold", "0.3",
  "--bootstraps", "1000",
  "--permutations", "1000",
  "--weighted", "TRUE",
  "--directed", "FALSE",
  "--types", "TRUE"
))

setwd(old_wd)

## How similar are the datasets in terms of taxonomic appeareance

# Load network
sparcc_net =readr::read_tsv( "results/tara_oceans_analysis/atlantic/functionink_input.tsv")

# Load funk partition

funk_part = readr::read_tsv("results/tara_oceans_analysis/atlantic/functionink_tmp/Partition-NL_Average_StopStep-1579_functionink_input.tsv",
skip = 9, col_names = c("ASV","cluster"))
funk_part
colnames(tax_tara)

tax_tara = tax_tara |>
  rownames_to_column("ASV")

funk_tax = inner_join(funk_part,tax_tara,by = "ASV")
