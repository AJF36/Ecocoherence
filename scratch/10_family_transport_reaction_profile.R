# Builds a family-level transport-reaction repertoire from the CarveMe genome
# models (transport reactions = SBO:0000185, restricted to those touching the
# extracellular compartment "_e", i.e. import/secretion with the environment,
# excluding purely intracellular periplasm<->cytoplasm shuttling) and checks
# whether the number of distinct transport reactions per family relates to
# that family's ecological-coherence Z-score (results/8a_.../Z_score).

rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(this.path)

script_dir <- this.dir()
results_dir <- normalizePath(file.path(script_dir, "..", "results", "10_family_transport_reaction_profile"), mustWork = FALSE)
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "10_family_transport_reaction_profile"), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

metabolism_dir <- normalizePath(file.path(script_dir, "..", "metabolism"), mustWork = FALSE)

# --- ASV -> Family (RDP taxonomy, same source/cleaning as steps 1 & 6a) ----
fileTaxonomy <- normalizePath(file.path(script_dir, "..", "data", "marine_particles_source_data", "sequence_table.ESV.fasta_RDPclassified.txt"))
taxa.in <- read.table(fileTaxonomy, sep = ";")
colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
                        "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
asv_family <- taxa.in %>%
  transmute(ASV = taxa_id,
            Family = Family %>% gsub("_incertae_sedis", "", .) %>% gsub("_Incertae Sedis XI", "", .))

# --- ASV -> genome id (vsearch ESV<->GTDB genome match table) --------------
fileMatch <- normalizePath(file.path(metabolism_dir, "genome_aligment", "matched_ESV_id_0.97.tsv"))
matched_esv <- read_tsv(fileMatch, col_names = FALSE, show_col_types = FALSE) %>%
  transmute(ASV = X9, genome_raw = X10) %>%
  filter(genome_raw != "*") %>%
  mutate(genome_id = gsub(".*([A-Z]{3}_\\d+).*", "\\1", genome_raw))

# --- Family -> genome ids ----------------------------------------------
family_genomes <- asv_family %>%
  inner_join(matched_esv, by = "ASV") %>%
  distinct(Family, genome_id)

cat(sprintf("Families with >=1 matched genome: %d\n", n_distinct(family_genomes$Family)))
cat(sprintf("Distinct genomes involved: %d\n", n_distinct(family_genomes$genome_id)))

# --- Locate one model file per genome id --------------------------------
# Use the flat metabolic_models/ dirs (pre module-split, one file per genome)
# with metabolism/misosoup/ as a fallback for any genome not found there.
model_files <- list.files(metabolism_dir, pattern = "_model\\.xml$", recursive = TRUE, full.names = TRUE)
flat_models <- model_files[grepl("metabolic_models/[^/]+_model\\.xml$", model_files) |
                              grepl("/misosoup/[^/]+_model\\.xml$", model_files)]
model_lookup <- tibble(path = flat_models,
                        genome_id = str_extract(basename(flat_models), "^[A-Z]{3}_\\d+")) %>%
  distinct(genome_id, .keep_all = TRUE)

missing_genomes <- setdiff(unique(family_genomes$genome_id), model_lookup$genome_id)
cat(sprintf("Genomes missing a model file: %d of %d\n", length(missing_genomes), n_distinct(family_genomes$genome_id)))

genomes_to_parse <- model_lookup %>% filter(genome_id %in% family_genomes$genome_id)

# --- Parse each model: transport reactions crossing the extracellular space --
# Regex-based on the raw text instead of per-node xml2 XPath: CarveMe/libSBML
# output is regular enough for this, and it's >100x faster than walking each
# <reaction> node with xml_find_all() (which took 30+ min for 666 models and
# still hadn't finished).
extract_transport_reactions <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  blocks <- unlist(str_extract_all(txt, regex("<reaction\\b.*?</reaction>|<reaction\\b[^>]*/>", dotall = TRUE)))
  if (length(blocks) == 0) return(character(0))
  sbo <- str_match(blocks, 'sboTerm="([^"]+)"')[, 2]
  blocks <- blocks[!is.na(sbo) & sbo == "SBO:0000185"]
  if (length(blocks) == 0) return(character(0))
  ids <- str_match(blocks, '(?<![\\w])id="([^"]+)"')[, 2]
  crosses_e <- vapply(blocks, function(b) {
    sp <- str_match_all(b, 'species="([^"]+)"')[[1]][, 2]
    any(grepl("_e$", sp))
  }, logical(1))
  ids[crosses_e]
}

cat("Parsing genome models for transport reactions...\n")
genome_reactions <- pmap_dfr(genomes_to_parse, function(path, genome_id) {
  rxns <- tryCatch(extract_transport_reactions(path), error = function(e) character(0))
  if (length(rxns) == 0) return(tibble(genome_id = character(0), reaction_id = character(0)))
  tibble(genome_id = genome_id, reaction_id = rxns)
})

write_tsv(genome_reactions, file.path(results_dir, "genome_transport_reactions.tsv"))

# --- Family-level matrix: union of transport reactions across its genomes --
family_reactions <- family_genomes %>%
  inner_join(genome_reactions, by = "genome_id", relationship = "many-to-many") %>%
  distinct(Family, reaction_id)

family_matrix <- family_reactions %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = reaction_id, values_from = present, values_fill = 0L)
write_tsv(family_matrix, file.path(results_dir, "family_transport_reaction_matrix.tsv"))

family_profile <- family_genomes %>%
  group_by(Family) %>%
  summarise(n_genomes = n_distinct(genome_id), .groups = "drop") %>%
  left_join(
    family_reactions %>% group_by(Family) %>% summarise(n_unique_transport_reactions = n_distinct(reaction_id), .groups = "drop"),
    by = "Family"
  ) %>%
  mutate(n_unique_transport_reactions = replace_na(n_unique_transport_reactions, 0L))
write_tsv(family_profile, file.path(results_dir, "family_transport_reaction_profile.tsv"))

# --- Compare to coherence Z-scores (results/8a_compute_zscore_families) ----
# NB: the reorganized results/ tree only has this populated for substrates
# rerun since the src/ reorg; the legacy Z_score/ folder is complete for all
# 7 substrates, so prefer it when the new one is missing/incomplete.
z_new <- file.path(script_dir, "..", "results", "8a_compute_zscore_families", "z_score_families.tsv")
z_legacy <- file.path(script_dir, "..", "Z_score", "z_score_families.tsv")
z_file <- if (file.exists(z_new) && file.info(z_new)$size > 0) z_new else z_legacy
cat(sprintf("Using Z-score file: %s\n", z_file))
z_scores <- read_tsv(z_file, show_col_types = FALSE) %>%
  rename(Family = family)

merged <- inner_join(z_scores, family_profile, by = "Family")
write_tsv(merged, file.path(results_dir, "transport_reactions_vs_coherence.tsv"))

# Primary test: one point per family (mean Z across substrates) avoids
# pseudoreplication, since the transport profile doesn't vary by substrate.
family_mean_z <- merged %>%
  group_by(Family, n_unique_transport_reactions, n_genomes) %>%
  summarise(mean_Z = mean(Z), .groups = "drop")

cor_test_mean <- cor.test(family_mean_z$n_unique_transport_reactions, family_mean_z$mean_Z, method = "spearman")
cat("Spearman correlation, mean Z per family vs n_unique_transport_reactions:\n")
print(cor_test_mean)

p_mean <- ggplot(family_mean_z, aes(x = n_unique_transport_reactions, y = mean_Z)) +
  geom_point(aes(size = n_genomes), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  labs(x = "Unique transport reactions (family, union across genomes)",
       y = "Mean coherence Z-score across substrates (low = more coherent)",
       size = "# genomes",
       title = "Family transport-reaction repertoire vs ecological coherence",
       subtitle = sprintf("Spearman rho = %.2f, p = %.3g", cor_test_mean$estimate, cor_test_mean$p.value)) +
  theme_bw()
ggsave(file.path(figures_dir, "transport_reactions_vs_mean_coherence.pdf"), p_mean, width = 7, height = 5)

# Secondary view: per-substrate, in case the relationship is substrate-specific
p_facet <- ggplot(merged, aes(x = n_unique_transport_reactions, y = Z)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  facet_wrap(~substrate) +
  labs(x = "Unique transport reactions (family, union across genomes)",
       y = "Coherence Z-score (low = more coherent)",
       title = "Family transport-reaction repertoire vs coherence, per substrate") +
  theme_bw()
ggsave(file.path(figures_dir, "transport_reactions_vs_coherence_per_substrate.pdf"), p_facet, width = 10, height = 7)

per_substrate_cor <- merged %>%
  group_by(substrate) %>%
  summarise(
    n = n(),
    rho = suppressWarnings(cor.test(n_unique_transport_reactions, Z, method = "spearman")$estimate),
    p_value = suppressWarnings(cor.test(n_unique_transport_reactions, Z, method = "spearman")$p.value),
    .groups = "drop"
  )
write_tsv(per_substrate_cor, file.path(results_dir, "per_substrate_correlation.tsv"))
print(per_substrate_cor)
