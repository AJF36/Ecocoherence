##################################################
# centrality_vs_ngenes.R
##################################################
# Test if the centrality of the different modules is related to the number of genes
# related to the substrate degradation/synthesis
#
# For every substrate and every centrality measure of the coarse-grained
# (module-level) network, fit centrality ~ n_genes_norm, where n_genes_norm is
# the number of dbCAN genes per genome in the module (n_genes / n_ids). The
# model summary is written to results/ and annotated on the figures.
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())

library(tidyverse)
library(broom)

results_dir = "results/centrality_vs_ngenes"

figures_dir = "figures/centrality_vs_ngenes"

for (dir in c(results_dir, figures_dir)) {
  if(!(fs::dir_exists(dir))) {
    fs::dir_create(dir)
  }
}

# Load the analysis of the coarse-grained networks

coarse_grained_analysis_paths = fs::dir_ls("results/coarse_grained_analysis/",regexp = "coarse_grained_net_analysis_.+\\.RDS$")
names(coarse_grained_analysis_paths) = gsub(".+_analysis_([A-Za-z]+)\\.RDS","\\1",coarse_grained_analysis_paths)

coarse_grained_analysis_list = map(coarse_grained_analysis_paths, function(path) {
  readRDS(path)
})

# NetCoMi stores centralities for network 1 and network 2 of the comparison;
# only the "1" slots are filled here (a single network per substrate), the
# "2" slots are NULL.
centralities = c("degree1","between1","close1","eigenv1")

centrality_labels = c(degree1  = "Degree",
                      between1 = "Betweenness",
                      close1   = "Closeness",
                      eigenv1  = "Eigenvector")

# One long table per substrate: module x centrality measure
centralities_df_list = map(coarse_grained_analysis_list, function(coarse_net) {
  map(set_names(centralities), function(centrality) {
    enframe(coarse_net[["centralities"]][[centrality]], name = "module", value = "centrality")
  }) |>
    list_rbind(names_to = "centrality_measure")
})

# Load the results of the n_genes

n_genes_list = readRDS("results/dbcan_analysis/results_ngenes.RDS")

n_genes_summary_list = map(n_genes_list, "summary")

# Only the substrates present in both sources (the coarse-grained analysis also
# covers AgaroseCarrageenan, which has no dbCAN gene counts).
substrates = intersect(names(n_genes_summary_list), names(centralities_df_list))

# Now for each substrate join the dataframes.
# inner_join, not full_join: a few modules exist in the coarse-grained network
# but have no dbCAN gene summary (and vice versa). With full_join those rows
# enter the data as NA and are silently dropped by lm(), which makes the
# reported n disagree with the number of rows actually plotted.
join_df = map(set_names(substrates), function(substrate) {
  gen_sub = n_genes_summary_list[[substrate]]
  centralities_sub = centralities_df_list[[substrate]]
  inner_join(gen_sub, centralities_sub, by = "module")
}) |>
  list_rbind(names_to = "substrate") |>
  mutate(centrality_measure = factor(centrality_measure, levels = centralities,
                                     labels = centrality_labels[centralities]))

# Report the modules lost in the join, so the drop is explicit and not silent
walk(substrates, function(substrate) {
  cg_modules = unique(centralities_df_list[[substrate]]$module)
  gene_modules = n_genes_summary_list[[substrate]]$module
  dropped = union(setdiff(cg_modules, gene_modules), setdiff(gene_modules, cg_modules))
  if (length(dropped) > 0) {
    message(glue::glue("{substrate}: {length(dropped)} module(s) dropped in the join: ",
                       "{paste(dropped, collapse = ', ')}"))
  }
})

##################################################
# Linear models: centrality ~ n_genes_norm
##################################################
# centrality is the response and the gene density the predictor (we ask whether
# a module's gene content explains how central it is). n_genes_norm = genes per
# genome in the module, so module size does not drive the relation by itself.

fits = join_df |>
  nest(data = -c(substrate, centrality_measure)) |>
  mutate(fit = map(data, ~ lm(centrality ~ n_genes_norm, data = .x)))

lm_stats = fits |>
  mutate(coefs = map(fit, ~ tidy(.x) |> filter(term == "n_genes_norm")),
         gof   = map(fit, glance),
         # Spearman as a rank-based check: n_genes_norm is right-skewed and the
         # centralities are bounded, so a single module can dominate the OLS fit
         spearman_rho = map_dbl(data, ~ suppressWarnings(cor(.x$n_genes_norm, .x$centrality, method = "spearman")))) |>
  unnest(coefs) |>
  unnest(gof, names_sep = "_") |>
  transmute(substrate,
            centrality_measure,
            n = gof_nobs,
            slope = estimate,
            std_error = std.error,
            statistic,
            p_value = p.value,
            r_squared = gof_r.squared,
            adj_r_squared = gof_adj.r.squared,
            spearman_rho) |>
  # 6 substrates x 4 centralities = 24 tests
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

print(lm_stats, n = Inf)

write_tsv(lm_stats, file.path(results_dir, "centrality_vs_ngenes_lm_stats.tsv"))
write_tsv(join_df, file.path(results_dir, "centrality_vs_ngenes_data.tsv"))
saveRDS(list(data = join_df, lm_stats = lm_stats, fits = fits),
        file.path(results_dir, "centrality_vs_ngenes.RDS"))

##################################################
# Figures
##################################################

# Annotation placed in the top-left corner of each panel
annotation_df = lm_stats |>
  mutate(label = glue::glue("slope = {signif(slope, 3)}\n",
                            "R2 = {signif(r_squared, 3)}\n",
                            "p = {signif(p_value, 3)}\n",
                            "rho = {signif(spearman_rho, 2)}\n",
                            "n = {n}"))

plot_substrate = function(substrate_name) {
  dat = join_df |> filter(substrate == substrate_name)
  ann = annotation_df |> filter(substrate == substrate_name)

  ggplot(dat, aes(x = n_genes_norm, y = centrality)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.1, size = 3.2, lineheight = 0.95,
              inherit.aes = FALSE) +
    facet_wrap(~ centrality_measure, scales = "free_y") +
    # headroom so the annotation does not sit on top of the points
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    theme_classic(base_size = 14) +
    labs(title = substrate_name,
         x = "dbCAN genes per genome (n_genes / n_ids)",
         y = "Centrality")
}

walk(substrates, function(substrate_name) {
  p = plot_substrate(substrate_name)
  ggsave(file.path(figures_dir, glue::glue("{substrate_name}_centrality_vs_ngenes.pdf")),
         p, width = 9, height = 7)
})

# Combined overview: centrality measure (rows) x substrate (columns)
combined_plot = ggplot(join_df, aes(x = n_genes_norm, y = centrality)) +
  geom_point(alpha = 0.7, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
  geom_text(data = annotation_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.1, size = 2.2, lineheight = 0.95,
            inherit.aes = FALSE) +
  facet_grid(centrality_measure ~ substrate, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  theme_bw(base_size = 12) +
  labs(x = "dbCAN genes per genome (n_genes / n_ids)",
       y = "Centrality")

ggsave(file.path(figures_dir, "all_substrates_centrality_vs_ngenes.pdf"),
       combined_plot, width = 16, height = 10)
