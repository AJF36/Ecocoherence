##################################################
# cazy_distance_families.R
##################################################
# Compare the functional distance in CAZY subfamilies between the 
# different families
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())

results_dir = "analysis/results/cazy_distance_families.R"
figures_dir = "analysis/figures/cazy_distance_families.R"

for (dir in c(results_dir,figures_dir)) {
if(!(fs::dir_exists(dir))) {
fs::dir_create(dir) }
}

library(tidyverse)
substrates = c("Alginate","Agarose","AgaroseAlginate","Chitin","Carrageenan","AgaroseChitosan")

# dbCAN can assign several substrates to a single gene, ";"-separated
# ("alginate;beta-glucan"). The exact %in% match below is deliberate: those
# multi-substrate hits are not specific to the degradation of our substrate,
# so only genes assigned a single, unambiguous substrate are counted.

# substrate = substrates[1]

# Load the tax table
fileTaxonomy <- file.path("analysis", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
taxonomy = rownames_to_column(taxonomy, var = "ASV")

# Load the asv-genome  mapping
asv_genome_mapping = readr::read_tsv("metabolism/genome_aligment/matched_ESV_id_0.97.tsv",col_names = F)

asv_genome_mapping = asv_genome_mapping |> 
  select(X9,X10)  |> 
  filter(X10 != "*")

  colnames(asv_genome_mapping) <- c("ASV","ID")

clean_id = gsub("[A-Za-z]{2}_([A-Za-z]{3}_[0-9]{9})\\.[0-9]*","\\1",asv_genome_mapping$ID)
asv_genome_mapping[["ID"]] = clean_id

tax_genome_table = inner_join(taxonomy,asv_genome_mapping, by = "ASV")
# Delete the last row, since the id is a sequence
tax_genome_table  = tax_genome_table[-17973,]


families = unique(tax_genome_table$Family)
names(families) = families

id_by_family = map(families, function(family) {
  #
  #  Filter the df 
  f_tax_genome_table = tax_genome_table |> 
    filter(Family == family) |> 
    pull(var = ID)

    return(unique(f_tax_genome_table))
})

id_by_family

# Count the number of genes per id

id_list =  fs::dir_ls(path = "metabolism/dbcan_all/")

enzyme_activity_list = list("PL","GH","AA",c("AA","GH","PL"))

# Proteome size, to normalise the counts. The annotated genomes span 746-8347
# proteins, so a raw "genes per genome" partly tracks genome size (and genome
# size tracks lifestyle), which would confound the family comparison. Counted
# once and cached, it takes ~6 s over the 666 proteomes.
proteome_size_file = fs::path(results_dir, "proteome_size.tsv")

if (fs::file_exists(proteome_size_file)) {
  proteome_size_df = readr::read_tsv(proteome_size_file, show_col_types = FALSE)
} else {
  proteome_files = fs::dir_ls("metabolism/proteomes_all", regexp = "\\.faa$")
  proteome_size_df = tibble(
    ID = fs::path_ext_remove(basename(proteome_files)),
    n_proteins = vapply(proteome_files,
                        function(f) sum(startsWith(readLines(f, warn = FALSE), ">")),
                        numeric(1), USE.NAMES = FALSE))
  readr::write_tsv(proteome_size_df, proteome_size_file)
}

proteome_size = setNames(proteome_size_df$n_proteins, proteome_size_df$ID)


for(enzyme_activity in enzyme_activity_list ) {

  # glue()/paste() vectorise, so collapse the multi-class case to a single
  # string before it reaches any file name or plot label
  activity_label = paste(enzyme_activity, collapse = "+")
  activity_pattern = paste0("(^|\\+)(", paste(enzyme_activity, collapse = "|"), ")[0-9]")

    number_genes_id = map(id_by_family, function(id_vector) {
    
    names(id_vector) = id_vector # So the final vectors has names
    
    # intersect_id = intersect(x =id_vector,y = id_list)
    
    purrr::map_dbl(id_vector, function(id) {

      tryCatch(expr = {
      
      overview_id_path = fs::dir_ls(glue::glue("metabolism/dbcan_all/{id}"),regexp = "overview.tsv")
      
      
      overview_id = suppressMessages( read_tsv(overview_id_path))
      
      colnames(overview_id)[6] = "n_tools"
      
      # Filter the overview file for consensus hits and enzyme activity.
      # DIAMOND holds the CAZy family ("GH16", "PL7", "AA10+CBM2"), never the
      # bare class, so the class has to be matched as a pattern. Anchor it to
      # the start of a "+"-separated token and require the family number:
      # a plain "AA" also matches GenBank accessions such as "AAZ25395.1+GH23".
      # The alternation covers the multi-class case, where a gene carrying more
      # than one class is counted once (union of genes, not sum of counts).
      overview_id_f = overview_id |>
        filter(stringr::str_detect(DIAMOND, pattern = activity_pattern) &
                 n_tools == 3) |>
        filter(Substrate %in% tolower(substrates))

      # Normalise per genome, so each genome contributes a rate rather than a
      # count and the family value below is a mean of rates
      1000 * nrow(overview_id_f) / proteome_size[[id]]

      },
      error = function(e) return(NA_real_)
      )
    
      })
    
    
  }
  ) |> 
  purrr::map(function(vector) {
    vector[!is.na(vector)]
  })  |> 
    discard(function(vector) is_empty(vector) | length(vector) < 3) |>  # Discard empty vectors and families with less than 3 genomes
    map_dbl(function(vector) {
      sum(vector)/length(vector)
    })

  n_genes_df = enframe(number_genes_id,name = "Family",value = "n_genes_per_1k_proteins")

  # Order the family labels by the count so the barplot reads top-down
  n_genes_df_f = n_genes_df |>
    filter(n_genes_per_1k_proteins != 0) |>
    mutate(Family = forcats::fct_reorder(Family, n_genes_per_1k_proteins))

  readr::write_tsv(arrange(n_genes_df_f, desc(n_genes_per_1k_proteins)),
                  fs::path(results_dir, glue::glue("n_genes_per_1k_proteins_by_family_{activity_label}.tsv")))

  n_genes_plot = ggplot(n_genes_df_f, aes(x = n_genes_per_1k_proteins, y = Family)) +
    geom_col(fill = "steelblue") +
    labs(x = glue::glue("{activity_label} genes per 1000 proteins"),
        y = "Family",
        title = glue::glue("{activity_label} genes per 1000 proteins by family"),
        subtitle = glue::glue("consensus hits (3 tools) on marine particle substrates\n",
                              "families with >= 3 annotated genomes")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_blank())

  ggsave(fs::path(figures_dir, glue::glue("n_genes_per_1k_proteins_by_family_{activity_label}.pdf")),
        n_genes_plot, width = 7, height = 11)

  print(n_genes_plot)   # bare expressions do not auto-print inside a for loop

}




## Lets start with prescense/abscense vectors

ids_dirs = fs::dir_ls(path = "metabolism/dbcan_all/",regexp = "[A-Z]{3}") 

id_by_family
test = map(id_by_family, function(ids) {

  #Extract the ids paths
  f_id_dirs = ids_dirs[which(basename(ids_dirs) %in% ids)]

    results_df = map_df(f_id_dirs, function(id_dir) {
      
      overview_path = fs::dir_ls(path = id_dir, regexp = "overview.tsv")
      overview_df = suppressMessages(readr::read_tsv(overview_path))

      colnames(overview_df)[6] = "n_tools"
      colnames(overview_df)[7] = "recommended_results"
      # print("hEY")
      overview_df_f = overview_df |>
        filter((n_tools == 3) & (Substrate %in% tolower(substrates))) |>
        filter(stringr::str_detect(DIAMOND,pattern = "GH")) # To compare degradation
      
        # pull(var = recommended_results)
      

      # print("Ho")
      return(overview_df_f)
    }) 
  

}) 

test_f = keep(test,.p = ~ nrow(.x) > 0)

# Now make the full matrix 
map(.x = test_f,"DIAMOND")
# test_f[[1]][["recommended_results"]]

families = names(test_f)
cazy_sub = unique(unlist(map(.x = test_f,"DIAMOND"),use.names = F))  # Change here the ther matching column of dbcan

f_sub_matrix = matrix(data = 0,nrow = length(cazy_sub), ncol = length(families))

colnames(f_sub_matrix) = families
rownames(f_sub_matrix) = cazy_sub
# Now populate the matrix

for (family in families) {
  
  family_df = test_f[[family]] 
  
  family_enz = family_df[["DIAMOND"]]
  
  # print(family_enz)
  
  for (cazy in cazy_sub) {
    
    # print(cazy %in% family_enz)
    
    if(cazy %in% family_enz) {

      f_sub_matrix[cazy,family] <- 1

    } 
  
  }

}
any(f_sub_matrix == 1)

jaccard_dist = vegan::vegdist(t(f_sub_matrix),method = "jaccard")
class(jaccard_dist)
jaccard_dist
pheatmap::pheatmap(as.matrix(jaccard_dist),show_rownames = T)
plot(hclust(jaccard_dist))
