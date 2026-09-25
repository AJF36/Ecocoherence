##################################################
# 8f_build_pooled_module_network.R
##################################################
# Pooled module-level network: every guildGT4 module from every substrate,
# drawn as its own node (no JSD clustering/aggregation, unlike
# 8e_build_module_cluster_meta_network.R) - one combined view of all 7
# substrates' coarse-grained networks side by side, to look at individual
# module-module connectivity directly instead of through the cluster
# summary. There is no data linking a module in one substrate to a module
# in another (correlations are computed within a substrate's own samples),
# so this is necessarily 7 disconnected subgraphs drawn together, not a
# network with genuine cross-substrate edges.
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

results_dir <- normalizePath(file.path(repo_root, "results", "8f_build_pooled_module_network"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(repo_root, "figures", "8f_build_pooled_module_network"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

substrates <- c("Agarose", "Alginate", "AgaroseAlginate", "AgaroseCarrageenan",
                "AgaroseChitosan", "Chitin", "Carrageenan")

coarse_grained_dir <- file.path(repo_root, "results", "coarse_grained_analysis")

# Pool every substrate's coarse-grained edges, renaming module ids to
# "<module>_<substrate>" so every node across all substrates is unique.
all_edges <- purrr::map_dfr(substrates, function(substrate) {
  read_tsv(file.path(coarse_grained_dir, paste0("coarse_grained_network_", substrate, ".tsv")), show_col_types = FALSE) %>%
    transmute(
      substrate = substrate,
      N1 = paste0(N1, "_", substrate),
      N2 = paste0(N2, "_", substrate),
      Interaction, Type
    )
})
write_tsv(all_edges, file.path(results_dir, "pooled_module_network_edges.tsv"))

nodes <- union(all_edges$N1, all_edges$N2)
mat <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
for (i in seq_len(nrow(all_edges))) {
  mat[all_edges$N1[i], all_edges$N2[i]] <- mat[all_edges$N2[i], all_edges$N1[i]] <- all_edges$Interaction[i]
}
diag(mat) <- 0

netcomi_pooled <- NetCoMi::netConstruct(data = mat, dataType = "association",
                                         sparsMethod = "none", weighted = TRUE)
props_pooled <- NetCoMi::netAnalyze(netcomi_pooled)

# Color each module node by its own ecological strategy label (not a
# cluster majority, since there's no clustering step here), coarsened the
# same way as 8e: attachment/selection -> early, transition -> mid,
# facilitation -> late, generalist -> generalist. Modules without a label
# (most of them - that file only covers 117/260 modules, see 8e) are grey.
ecological_strategies <- read_tsv(
  file.path(repo_root, "heatmap_barplot", "ecological_strategies_modules.tsv"),
  show_col_types = FALSE
)
strategy_bucket <- c(attachment = "early", selection = "early", transition = "mid",
                      facilitation = "late", generalist = "generalist")
strategy_colors <- c(early = "#2166AC", mid = "#F4A582", late = "#B2182B",
                      generalist = "#762A83", unknown = "grey80")

node_module <- sub("_[^_]+$", "", nodes)
node_substrate <- sub("^[^_]+_", "", nodes)
node_strategy <- ecological_strategies$strategy[match(paste(node_substrate, node_module), paste(ecological_strategies$substrate, ecological_strategies$module))]
node_bucket <- strategy_bucket[node_strategy]
node_bucket[is.na(node_bucket)] <- "unknown"
node_colors <- setNames(strategy_colors[node_bucket], nodes)

pdf(file.path(figures_dir, "pooled_module_network.pdf"), width = 34, height = 30, pointsize = 9)
plot(props_pooled,
  repulsion = 1,
  rmSingles = "all",
  labelScale = FALSE,
  nodeSize = "degree",
  nodeSizeSpread = 3,
  nodeColor = "colorVec",
  colorVec = node_colors,
  nodeTransp = 0,
  cexNodes = 1.5,
  cexLabels = 0.6,
  edgeTranspHigh = 20,
  title1 = "Pooled module network - all substrates, individual modules (no clustering)",
  showTitle = TRUE,
  cexTitle = 2,
  mar = c(1, 4, 6, 4)
)
legend(x = par("usr")[1], y = par("usr")[4], legend = names(strategy_colors), fill = strategy_colors,
       title = "Strategy", bty = "n", cex = 2, xpd = NA)
dev.off()
