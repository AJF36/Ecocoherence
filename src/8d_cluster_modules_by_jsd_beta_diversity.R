##################################################
# 8d_cluster_modules_by_jsd_beta_diversity.R
##################################################
# Cluster every module, across all substrates at once, by the
# Jensen-Shannon Divergence between their family-composition profiles (a
# beta-diversity distance), using every guildGT4 module as-is (min size 4,
# no further filtering) - unlike the earlier Aitchison-distance version of
# this step, no module or family is dropped for being small/sparse.
#
# Built on the miaverse (TreeSummarizedExperiment + mia::getDissimilarity)
# instead of phyloseq: the family x module matrix (rows = families, columns
# = "<module>_<substrate>", one column per module from every substrate,
# values = count of that module's ESVs belonging to that family) is stored
# as a TreeSummarizedExperiment's counts assay, and mia's own JSD
# implementation computes the pairwise distance between columns (modules).
#
# Since all substrates' modules are clustered together in one distance
# matrix, the same cluster id means the same taxonomic-composition
# neighborhood regardless of substrate - that's what makes the resulting
# ids usable as comparable colors/groups across per-substrate plots.
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list = ls())

library(readr)
library(dplyr)
library(this.path)
library(mia)
library(dendextend)

script_dir <- this.dir()
repo_root  <- normalizePath(file.path(script_dir, ".."))

results_dir <- normalizePath(file.path(repo_root, "results", "8d_cluster_modules_by_jsd_beta_diversity"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(repo_root, "figures", "8d_cluster_modules_by_jsd_beta_diversity"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# NOTE: same caveat as coarse_grained_analysis.R - functionInk output under
# the new src/ pipeline (results/4a_detect_modules_functionink/) has only
# been rerun for Alginate so far, so this reads from the legacy
# pre-reorg functionink/ folder (complete for all 7 real substrates).
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan",
                "AgaroseChitosan", "Chitin", "Carrageenan")
k <- 14 # kept from the earlier Aitchison version for continuity - not re-derived here

fileTaxonomy <- normalizePath(file.path(repo_root, "data", "marine_particles_source_data", "sequence_table.ESV.fasta_RDPclassified.txt"))
taxa.in <- read.table(fileTaxonomy, sep = ";")
colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
                        "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
taxonomy <- taxa.in[, c("taxa_id", "Family")]
taxonomy$Family <- gsub("_incertae_sedis", "", taxonomy$Family)
taxonomy$Family <- gsub("_Incertae Sedis XI", "", taxonomy$Family)

# Build one family x module count matrix, pooling every substrate's guildGT4
# modules into the same set of columns (module numbers are only unique
# within a substrate, hence the "<module>_<substrate>" column naming).
module_family_counts <- list()
for (substrate in substrates) {
  dir.funcInk <- normalizePath(file.path(repo_root, "functionink", substrate, "functionink_tmp"))
  fileFun <- list.files(dir.funcInk, pattern = "Partition.+guildGT4\\.txt$", full.names = TRUE)

  partition <- read.table(fileFun, sep = "\t", header = TRUE)
  partition$ESV <- rownames(partition)
  colnames(partition) <- c("guild", "ESV")

  merged <- merge(partition, taxonomy, by.x = "ESV", by.y = "taxa_id")
  merged$guild <- gsub("mod_", "", merged$guild)

  counts <- table(merged$Family, merged$guild)
  colnames(counts) <- paste0(colnames(counts), "_", substrate)
  module_family_counts[[substrate]] <- as.data.frame.matrix(counts)
}

all_families <- unique(unlist(lapply(module_family_counts, rownames)))
combined <- matrix(0, nrow = length(all_families), ncol = 0, dimnames = list(all_families, character(0)))
for (substrate in substrates) {
  m <- module_family_counts[[substrate]]
  filled <- matrix(0, nrow = length(all_families), ncol = ncol(m), dimnames = list(all_families, colnames(m)))
  filled[rownames(m), ] <- as.matrix(m)
  combined <- cbind(combined, filled)
}

tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = list(counts = combined))

jsd_dist <- mia::getDissimilarity(tse, method = "jsd", assay.type = "counts")

clustered <- hclust(jsd_dist, method = "ward.D2")
clusters <- cutree(clustered, k = k)

# Dendrogram of every module (leaf label = "<module>_<substrate>"), branches
# and labels colored by the k = 14 cutree clusters, so cluster membership
# can be read directly off the tree instead of only from the results table.
dendrogram <- as.dendrogram(clustered) %>%
  dendextend::set("labels_cex", 1) %>%
  dendextend::set("branches_k_color", k = k) %>%
  dendextend::set("labels_col", k = k) %>%
  dendextend::set("branches_lwd", 2)

pdf(file.path(figures_dir, "module_jsd_dendrogram.pdf"), width = 40, height = 20, pointsize = 9)
par(mar = c(20, 2, 4, 2))
plot(dendrogram, main = "Module clustering by JSD beta-diversity (ward.D2, k = 14)")
dendrogram %>% rect.dendrogram(k = k, border = 1, lty = 4, lwd = 2, lower_rect = -0.05)
dev.off()

module_number <- sub("_.*$", "", names(clusters))
substrate <- sub("^[^_]*_", "", names(clusters))

results_df <- data.frame(
  substrate  = substrate,
  module     = paste0("mod_", module_number),
  cluster_id = as.integer(clusters),
  stringsAsFactors = FALSE
)

write_tsv(results_df, file.path(results_dir, "module_jsd_clusters.tsv"))
