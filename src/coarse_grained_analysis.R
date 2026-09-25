##################################################
# coarse_grained_analysis.R
##################################################
# Build the coarse grained (module-level) networks of each substrate, to
# check whether relations between modules are sustained across substrates.
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################


rm(list = ls())

library(this.path)
library(tidyverse)
library(NetCoMi)

script_dir <- this.dir()
repo_root  <- normalizePath(file.path(script_dir, ".."))

source(file.path(repo_root, "..", "Resultados_Ana_Cuenda", "src", "CoarseGrain.R"))

results_dir <- normalizePath(file.path(repo_root, "results", "coarse_grained_analysis"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(repo_root, "figures", "coarse_grained_analysis"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# NOTE: the reorganized src/ pipeline (results/4a_.../5_...) has only been
# rerun for Alginate so far, so every other substrate's functionInk +
# guildGT4 output still only exists in the legacy pre-reorg folders
# (network_creation_sparcc/, filter_network_by_threshold/, functionink/).
# Reading from there until the full pipeline is rerun for all substrates
# under src/ - switch these three path helpers to the results/<step>/
# folders once that happens.
substrates <- c("Agarose", "Alginate", "AgaroseAlginate", "AgaroseCarrageenan",
                "AgaroseChitosan", "Chitin", "Carrageenan")

count_path_for <- function(s) {
  file.path(repo_root, "network_creation_sparcc", s, paste0("otu_table_", s))
}

# The threshold-filtered network (not the raw step-1 output, which has no
# "type" column and is what made Cor_CoarseGrain() fail before: it expects
# a 4-column Node1/Node2/Interaction/Type table). This is also the network
# functionInk actually partitioned, so it's the one consistent with the
# guildGT4 module assignments below.
net_path_for <- function(s) {
  file.path(repo_root, "filter_network_by_threshold",
            paste0("interactions_filtered_p0.01_threshold_", s, "_.tsv"))
}

partition_path_for <- function(s) {
  fs::dir_ls(file.path(repo_root, "functionink", s), recurse = TRUE,
             regexp = "Partition.+guildGT4\\.txt$")
}

# Beta-diversity (JSD on family composition, ward.D2) clustering of modules,
# computed once across all substrates together by
# 8d_cluster_modules_by_jsd_beta_diversity.R (k = 14) - every guildGT4
# module is included, none dropped. Because that clustering pools every
# substrate's modules into one distance matrix, the same cluster id means
# the same taxonomic-composition neighborhood in every substrate's plot
# below - that's what makes the node colors comparable across networks.
jsd_clusters <- suppressMessages(readr::read_tsv(
  file.path(repo_root, "results", "8d_cluster_modules_by_jsd_beta_diversity", "module_jsd_clusters.tsv")
))


# One fixed id -> color table, shared by every substrate's plot below.
# NetCoMi's own nodeColor = "cluster" re-ranks and recolors cluster ids
# from scratch for each network, based only on which ids happen to be
# present in *that* network (e.g. cluster id 4 might be plotted red in one
# substrate and green in another, purely because of which other ids are or
# aren't present alongside it there) - that silently breaks the whole point
# of using one shared cross-substrate clustering. Building the palette once
# here and passing it via nodeColor = "colorVec" instead keeps a given
# cluster id's color identical across every plot.
cluster_ids <- sort(unique(jsd_clusters$cluster_id))
cluster_palette <- setNames(c("grey80", grDevices::rainbow(length(cluster_ids))),
                             c("0", as.character(cluster_ids)))

coarse_grained_nets <- purrr::map(substrates, function(s) {

  # otu table is ASVs x samples with a leading #OTU_ID column; CoarseGrain's
  # abun expects a numeric samples x ASVs matrix with no ID column.
  count_sub <- suppressMessages(readr::read_tsv(count_path_for(s)))
  count_sub <- tibble::column_to_rownames(count_sub, var = colnames(count_sub)[1])
  count_sub <- t(as.matrix(count_sub))

  net_sub <- suppressMessages(readr::read_tsv(net_path_for(s)))
  # step 3's own postprocessing prefixed every column name with "#"
  colnames(net_sub) <- sub("^#", "", colnames(net_sub))

  # skip = 1 to drop the "guild" header line only - the file has no other
  # preamble, so the previous skip = 9 silently dropped 8 real ASV/module
  # rows for every substrate.
  partition_sub <- suppressMessages(
    readr::read_tsv(partition_path_for(s), skip = 1, col_names = c("ASV", "module"))
  )

  Cor_CoarseGrain(cluster = partition_sub, cor = net_sub, abun = count_sub)
})
names(coarse_grained_nets) <- substrates

coarse_grained_nets[[1]]
# Save the module-level interaction tables, and build/plot a NetCoMi network
# for each substrate from the coarse-grained (module x module) matrix.

purrr::iwalk(coarse_grained_nets, function(cg_net, substrate) {

  readr::write_tsv(cg_net, file.path(results_dir, paste0("coarse_grained_network_", substrate, ".tsv")))

  modules <- union(cg_net$N1, cg_net$N2)
  mat <- matrix(0, length(modules), length(modules), dimnames = list(modules, modules))
  for (i in seq_len(nrow(cg_net))) {
    mat[cg_net$N1[i], cg_net$N2[i]] <- cg_net$Interaction[i]
    mat[cg_net$N2[i], cg_net$N1[i]] <- cg_net$Interaction[i]
  }

  # data is already a module-module association matrix, so skip NetCoMi's
  # own association estimation/sparsification (dataType = "association",
  # sparsMethod = "none") - it would otherwise try to treat mat as raw counts.
  netcomi_sub <- NetCoMi::netConstruct(data = mat, dataType = "association",
                                        sparsMethod = "none", weighted = TRUE)
  props_sub <- NetCoMi::netAnalyze(netcomi_sub)

  # Cross-substrate JSD beta-diversity cluster id per module, colored via
  # the fixed cluster_palette above so the same id always gets the same
  # color regardless of which other ids happen to be present in this plot.
  substrate_clusters <- jsd_clusters[jsd_clusters$substrate == substrate, ]
  module_cluster <- substrate_clusters$cluster_id[match(modules, substrate_clusters$module)]
  module_cluster[is.na(module_cluster)] <- 0
  node_colors <- setNames(cluster_palette[as.character(module_cluster)], modules)
  
  #Save the analysis of the net

  saveRDS(object = props_sub ,glue::glue("results/coarse_grained_analysis/coarse_grained_net_analysis_{substrate}.RDS"))
  
  pdf(file.path(figures_dir, paste0("coarse_grained_network_", substrate, ".pdf")))
  plot(props_sub,
    sameLayout = TRUE,
    repulsion = 1.5,
    rmSingles = "inboth",
    labelScale = FALSE,
    nodeSize = "betweenness",
    nodeSizeSpread = 2.5,
    nodeColor = "colorVec",
    colorVec = node_colors,
    nodeTransp = 0,
    sameColThresh = 2,
    hubBorderCol = "darkgray",
    cexNodes = 2,
    edgeTranspHigh = 20,
    title1 = paste("Coarse-grained network -", substrate),
    showTitle = TRUE,
    cexTitle = 2,
    mar = c(1, 4, 4, 4)
  )
  dev.off()
})
