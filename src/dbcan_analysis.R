##################################################
# dbcan_analysis.R
##################################################
# In this script I perform ...
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())

results_dir = "results/dbcan_analysis"

figures_dir = "figures/dbcan_analysis"

results_dir = c(results_dir,figures_dir)

for (dir in results_dir) {
  
  if(!(fs::dir_exists(dir))) {
      fs::dir_create(dir) 
  }
}

substrates = c("Alginate","Agarose","Chitin","AgaroseAlginate","AgaroseChitosan","Carrageenan")

combined_substrates = c("AgaroseAlginate","AgaroseChitosan")

dbcan_path = "metabolism/dbcan_all/"

# Load mapping ASV-genome
map_asv_genome = read_tsv("metabolism/genome_aligment/matched_ESV_id_0.97.tsv",col_names = F)

map_asv_genome = map_asv_genome |> 
  select(X9,X10) |> 
  filter(X10 != "*")
colnames(map_asv_genome) = c("ASV","ID")

clean_ids = gsub("(.+)\\.[0-9]+","\\1",map_asv_genome$ID)
clean_ids = gsub("GB_|RS_","",clean_ids)


dbcan_dirs = fs::dir_ls("metabolism/dbcan_all/",regexp = "[A-Z]{3}_[0-9]{9}")

all(basename(dbcan_dirs) %in% clean_ids)

map_asv_genome$ID = clean_ids

map_asv_genome_f = map_asv_genome |>
  filter(ID %in% basename(dbcan_dirs))

# Load ASV taxonomy (same RDP-classified table used in mia_analysis.R) and
# derive the family of each genome via the ASV that matched it
taxa_in = read.table("data/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt", sep=";")
colnames(taxa_in) = c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                       "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
asv_taxonomy = taxa_in |>
  select(ASV = taxa_id, family = Family)

genome_family = map_asv_genome_f |>
  left_join(asv_taxonomy, by = "ASV") |>
  select(ID, family) |>
  distinct(ID, .keep_all = TRUE)

list_df = purrr::map(substrates, function(subs) {
 # For each substrate, compute the number of genes related to the degradation of the module
  # Load the funk partition

  # Combined substrates (e.g. "AgaroseAlginate") don't exist as such in the
  # dbCAN Substrate column, so match genes annotated to any of the substrates
  # that make up the combination
  subs_terms = if (subs %in% combined_substrates) {
    tolower(stringr::str_extract_all(subs, "[A-Z][a-z]*")[[1]])
  } else {
    tolower(subs)
  }
  subs_pattern = paste(subs_terms, collapse = "|")

  funk_file = fs::dir_ls(fs::path("functionink",subs),recurse = T ,regexp = "Partition.+GT4.txt$")
  
  funk = read_tsv(file = funk_file, skip = 1,col_names = F)
  colnames(funk) = c("ASV","module")

  # Now join the two tables
  subs_table = inner_join(funk,map_asv_genome_f, by = "ASV")
  print(head(subs_table))

  modules = unique(subs_table[["module"]])

  module_results = map(modules, function(mod) {

    module_table = subs_table |>
      filter(module == mod)

    sub_id = module_table$ID
    unique_id = unique(sub_id)

    subs_df = map_df(unique_id, function(id) {
      # Search for the summary.tsv
      overview_file = fs::path(dbcan_path,id,"overview.tsv")
      suppressMessages(readr::read_tsv(overview_file)) |>
        mutate(ID = id)

    })

    # Number of substrate-related genes per genome, then mapped back onto every
    # member (ASV) of the module so ASVs sharing a genome each get its count
    colnames(subs_df)[2] = "EC"
    colnames(subs_df)[6] = "n_tools"
    
    genome_counts = subs_df |>
      filter(n_tools == 3) |> 
      mutate(is_hit = stringr::str_detect(Substrate, subs_pattern)) |>
      group_by(ID) |>
      summarise(n_genes = sum(is_hit), .groups = "drop")

    per_id_counts = tibble(ID = sub_id) |>
      left_join(genome_counts, by = "ID") |>
      mutate(module = mod)

    n_genes = sum(per_id_counts$n_genes)
    n_ids = length(sub_id)

    list(
      summary = tibble(module = mod, n_genes = n_genes, n_ids = n_ids, n_genes_norm = n_genes / n_ids),
      per_id = per_id_counts,
      genome_counts = genome_counts
    )

    })

  counts_by_module = map_df(module_results, "summary")
  per_id_by_module = map_df(module_results, "per_id")
  genome_counts_all = map_df(module_results, "genome_counts") |>
    distinct(ID, n_genes)

  # Mean/sd of substrate-related genes per genome, grouped by taxonomic family
  family_summary = genome_counts_all |>
    left_join(genome_family, by = "ID") |>
    group_by(family) |>
    summarise(n_genomes = n(),
              mean_genes = mean(n_genes),
              sd_genes = sd(n_genes),
              .groups = "drop") |>
    arrange(desc(mean_genes))

  family_plot = ggplot(family_summary, aes(x = mean_genes, y = reorder(family, mean_genes))) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(xmin = pmax(mean_genes - sd_genes, 0), xmax = mean_genes + sd_genes), height = 0.3) +
    labs(y = "Family", x = "Mean number of genes per genome",
         title = paste(subs, "- genes per genome by family")) +
    theme_minimal(base_size = 9)

  ggsave(fs::path(figures_dir, paste0(subs, "_genes_per_genome_by_family.pdf")),
         plot = family_plot, width = 8, height = max(6, 0.22 * nrow(family_summary)),
         limitsize = FALSE)

  # Density plot: distribution of the per-member gene counts within each module
  density_plot = ggplot(per_id_by_module, aes(x = n_genes)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    facet_wrap(~module, scales = "free") +
    labs(x = "Number of genes (per member)", y = "Probability density") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())

  ggsave(fs::path(figures_dir, paste0(subs, "_gene_density_by_module.pdf")),
         plot = density_plot, width = 12, height = 10)

  return(list(summary = counts_by_module, per_id = per_id_by_module, family_summary = family_summary))

})
names(list_df) = substrates

list_df

saveRDS(object = list_df, file = "results/dbcan_analysis/results_ngenes.RDS")
##################################################
# Same analysis, but based on substrate_prediction.tsv (CGC-level predictions)
# instead of overview.tsv (gene-level annotations)
##################################################

list_df_cgc = purrr::map(substrates, function(subs) {

  subs_terms = if (subs %in% combined_substrates) {
    tolower(stringr::str_extract_all(subs, "[A-Z][a-z]*")[[1]])
  } else {
    tolower(subs)
  }
  subs_pattern = paste(subs_terms, collapse = "|")

  funk_file = fs::dir_ls(fs::path("functionink",subs),recurse = T ,regexp = "Partition.+GT4.txt$")

  funk = read_tsv(file = funk_file, skip = 1,col_names = F)
  colnames(funk) = c("ASV","module")

  subs_table = inner_join(funk,map_asv_genome_f, by = "ASV")

  modules = unique(subs_table[["module"]])

  module_results = map(modules, function(mod) {

    module_table = subs_table |>
      filter(module == mod)

    sub_id = module_table$ID
    unique_id = unique(sub_id)

    subs_df = map_df(unique_id, function(id) {
      # Search for the substrate_prediction.tsv
      pred_file = fs::path(dbcan_path,id,"substrate_prediction.tsv")
      # Column types are pinned explicitly: some genomes have an entirely
      # empty bitscore/score column, which readr would otherwise infer as
      # logical there and double elsewhere, breaking bind_rows across genomes
      suppressMessages(readr::read_tsv(pred_file, col_types = readr::cols(
        .default = readr::col_character(),
        bitscore = readr::col_double(),
        `dbCAN-sub substrate score` = readr::col_double()
      ))) |>
        mutate(ID = id)

    })

    # A CGC counts as a hit if either the dbCAN-PUL or the dbCAN-sub predicted
    # substrate matches the target substrate
    colnames(subs_df)[3] = "PUL_substrate"
    colnames(subs_df)[6] = "sub_substrate"

    # Genomes with no CGC rows at all (empty substrate_prediction.tsv) must
    # still be represented with n_genes = 0, otherwise they're dropped by
    # group_by/summarise and turn into NA (poisoning module-level sums) below
    genome_counts = subs_df |>
      mutate(is_hit = stringr::str_detect(dplyr::coalesce(PUL_substrate,""), subs_pattern) |
                      stringr::str_detect(dplyr::coalesce(sub_substrate,""), subs_pattern)) |>
      group_by(ID) |>
      summarise(n_genes = sum(is_hit), .groups = "drop") |>
      right_join(tibble(ID = unique_id), by = "ID") |>
      mutate(n_genes = tidyr::replace_na(n_genes, 0))

    per_id_counts = tibble(ID = sub_id) |>
      left_join(genome_counts, by = "ID") |>
      mutate(module = mod)

    n_genes = sum(per_id_counts$n_genes)
    n_ids = length(sub_id)

    list(
      summary = tibble(module = mod, n_genes = n_genes, n_ids = n_ids, n_genes_norm = n_genes / n_ids),
      per_id = per_id_counts,
      genome_counts = genome_counts
    )

    })

  counts_by_module = map_df(module_results, "summary")
  per_id_by_module = map_df(module_results, "per_id")
  genome_counts_all = map_df(module_results, "genome_counts") |>
    distinct(ID, n_genes)

  # Mean/sd of substrate-related CGCs per genome, grouped by taxonomic family
  family_summary = genome_counts_all |>
    left_join(genome_family, by = "ID") |>
    group_by(family) |>
    summarise(n_genomes = n(),
              mean_genes = mean(n_genes),
              sd_genes = sd(n_genes),
              .groups = "drop") |>
    arrange(desc(mean_genes))

  family_plot = ggplot(family_summary, aes(x = mean_genes, y = reorder(family, mean_genes))) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(xmin = pmax(mean_genes - sd_genes, 0), xmax = mean_genes + sd_genes), height = 0.3) +
    labs(y = "Family", x = "Mean number of substrate-related CGCs per genome",
         title = paste(subs, "- CGCs per genome by family")) +
    theme_minimal(base_size = 9)

  ggsave(fs::path(figures_dir, paste0(subs, "_CGC_genes_per_genome_by_family.pdf")),
         plot = family_plot, width = 8, height = max(6, 0.22 * nrow(family_summary)),
         limitsize = FALSE)

  # Density plot: distribution of the per-member CGC counts within each module
  density_plot = ggplot(per_id_by_module, aes(x = n_genes)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    facet_wrap(~module, scales = "free") +
    labs(x = "Number of substrate-related CGCs (per member)", y = "Probability density") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())

  ggsave(fs::path(figures_dir, paste0(subs, "_CGC_gene_density_by_module.pdf")),
         plot = density_plot, width = 12, height = 10)

  return(list(summary = counts_by_module, per_id = per_id_by_module, family_summary = family_summary))

})
names(list_df_cgc) = substrates

list_df_cgc
