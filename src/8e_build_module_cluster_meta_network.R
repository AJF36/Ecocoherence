##################################################
# 8e_build_module_cluster_meta_network.R
##################################################
# Cluster-level meta-network: aggregate the per-substrate coarse-grained
# (module-level) network edges from coarse_grained_analysis.R by the JSD
# beta-diversity cluster id of each endpoint module (8d), pooling across
# all 7 substrates. Each edge in the resulting network is a cluster pair
# (including self-pairs, i.e. two different modules from the same
# cluster), summarizing how modules of that composition type relate to
# modules of another composition type, wherever both appear together in a
# substrate's network - this is the direct test of whether a module-module
# relation is sustained across substrates rather than a one-off in a
# single substrate.
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

results_dir <- normalizePath(file.path(repo_root, "results", "8e_build_module_cluster_meta_network"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(repo_root, "figures", "8e_build_module_cluster_meta_network"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

substrates <- c("Agarose", "Alginate", "AgaroseAlginate", "AgaroseCarrageenan",
                "AgaroseChitosan", "Chitin", "Carrageenan")

jsd_clusters <- read_tsv(
  file.path(repo_root, "results", "8d_cluster_modules_by_jsd_beta_diversity", "module_jsd_clusters.tsv"),
  show_col_types = FALSE
)

# Majority ecological strategy per cluster, coarsened from the 5 named
# strategies (attachment/selection/transition/facilitation/generalist,
# from the legacy heatmap_barplot/ecological_strategies_modules.tsv) into
# early/mid/late/generalist. Coverage is partial (that file only labels
# modules that collapsed into one of those 5 reference clusters), so this
# is a lean based on however many labeled modules a cluster has, not a
# certainty - clusters with no labeled modules at all fall back to "unknown".
ecological_strategies <- read_tsv(
  file.path(repo_root, "heatmap_barplot", "ecological_strategies_modules.tsv"),
  show_col_types = FALSE
)
strategy_bucket <- c(attachment = "early", selection = "early", transition = "mid",
                      facilitation = "late", generalist = "generalist")

cluster_strategy <- jsd_clusters %>%
  left_join(ecological_strategies, by = c("substrate", "module")) %>%
  mutate(bucket = strategy_bucket[strategy]) %>%
  filter(!is.na(bucket)) %>%
  count(cluster_id, bucket) %>%
  group_by(cluster_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster_id, majority_strategy = bucket)

strategy_colors <- c(early = "#2166AC", mid = "#F4A582", late = "#B2182B",
                      generalist = "#762A83", unknown = "grey80")

coarse_grained_dir <- file.path(repo_root, "results", "coarse_grained_analysis")

all_edges <- purrr::map_dfr(substrates, function(substrate) {
  cg_net <- read_tsv(file.path(coarse_grained_dir, paste0("coarse_grained_network_", substrate, ".tsv")), show_col_types = FALSE)
  sub_clusters <- jsd_clusters[jsd_clusters$substrate == substrate, ]

  cg_net %>%
    mutate(
      substrate = substrate,
      cluster1 = sub_clusters$cluster_id[match(N1, sub_clusters$module)],
      cluster2 = sub_clusters$cluster_id[match(N2, sub_clusters$module)]
    ) %>%
    filter(!is.na(cluster1), !is.na(cluster2))
})

# Unordered cluster pair per module-module edge (self-pairs kept: two
# different modules from the same cluster, still an informative edge).
all_edges <- all_edges %>%
  mutate(cluster_A = pmin(cluster1, cluster2), cluster_B = pmax(cluster1, cluster2))

meta_edges <- all_edges %>%
  group_by(cluster_A, cluster_B) %>%
  summarise(
    mean_interaction = mean(Interaction),
    sd_interaction    = sd(Interaction),
    n_module_pairs    = n(),
    n_substrates      = n_distinct(substrate),
    substrates        = paste(sort(unique(substrate)), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(n_substrates), desc(abs(mean_interaction)))

write_tsv(meta_edges, file.path(results_dir, "module_cluster_meta_network.tsv"))

# Build/plot the meta-network: nodes = clusters, edge weight = mean
# interaction between their member modules, pooled over every substrate
# where both clusters co-occur.
clusters_present <- sort(unique(c(meta_edges$cluster_A, meta_edges$cluster_B)))
mat <- matrix(0, length(clusters_present), length(clusters_present),
               dimnames = list(clusters_present, clusters_present))
for (i in seq_len(nrow(meta_edges))) {
  a <- as.character(meta_edges$cluster_A[i])
  b <- as.character(meta_edges$cluster_B[i])
  mat[a, b] <- mat[b, a] <- meta_edges$mean_interaction[i]
}
diag(mat) <- 0

netcomi_meta <- NetCoMi::netConstruct(data = mat, dataType = "association",
                                       sparsMethod = "none", weighted = TRUE)
props_meta <- NetCoMi::netAnalyze(netcomi_meta)

# How many substrates actually support each cluster's edges at all (breadth
# of evidence) - reported in the results table rather than as node size,
# since NetCoMi's nodeSize options don't accept an arbitrary external vector.
support_by_cluster <- bind_rows(
  all_edges %>% transmute(cluster = cluster1, substrate),
  all_edges %>% transmute(cluster = cluster2, substrate)
) %>%
  group_by(cluster) %>%
  summarise(n_substrates_supporting = n_distinct(substrate), .groups = "drop") %>%
  arrange(cluster)
write_tsv(support_by_cluster, file.path(results_dir, "cluster_substrate_support.tsv"))

node_strategy <- cluster_strategy$majority_strategy[match(clusters_present, cluster_strategy$cluster_id)]
node_strategy[is.na(node_strategy)] <- "unknown"
node_colors <- setNames(strategy_colors[node_strategy], clusters_present)

pdf(file.path(figures_dir, "module_cluster_meta_network.pdf"), width = 9, height = 8)
plot(props_meta,
  repulsion = 1.2,
  rmSingles = "all",
  labelScale = FALSE,
  # betweenness is nearly uninformative here - the meta-network is almost
  # fully connected (105/105 possible cluster pairs), so betweenness
  # collapses toward zero for every node. Weighted degree (strength) varies
  # meaningfully instead.
  nodeSize = "degree",
  nodeSizeSpread = 3,
  nodeColor = "colorVec",
  colorVec = node_colors,
  nodeTransp = 0,
  cexNodes = 3,
  edgeTranspHigh = 0,
  title1 = "Module cluster meta-network (colored by majority ecological strategy)",
  showTitle = TRUE,
  cexTitle = 1.1,
  mar = c(1, 4, 6, 4)
)
# Keyword positions (e.g. "bottomright") get clipped by qgraph's own layout,
# so anchor the legend with explicit coordinates inside the confirmed
# post-plot par("usr") range instead, with xpd = NA to allow it to sit
# just past the node circles without being cut off.
legend(x = -1.4, y = 1.6, legend = names(strategy_colors), fill = strategy_colors,
       title = "Majority strategy", bty = "n", cex = 0.9, xpd = NA)
dev.off()
