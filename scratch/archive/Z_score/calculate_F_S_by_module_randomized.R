rm(list = ls())
library(readr)
library(phyloseq)
library(dplyr)
library(grid)
library(gridExtra)
library(this.path)
library(purrr)
library(ggplot2)


work_dir <- this.dir()
# substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")
substrate <- "Alginate"
# taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))
fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
taxonomy <- taxonomy %>% select(-c("Kingdom", "Genus"))

# Función para permutar todas las columnas de la tabla de taxonomía
permute_taxonomy <- function(taxonomy) {
  taxonomy_permuted <- taxonomy %>%
    mutate(across(everything(), ~ sample(.)))
  return(taxonomy_permuted)
}


n_total_families_per_substrate <- list()
results_list <- list()


for (substrate in substrates) {
  ## This bloc is the same as the not randomized
  # This goes inside the substrate loop
  ############################## 
  # functionink table 
  subs_path <- file.path(work_dir,"..","functionink",substrate)
    subs_path <- normalizePath(subs_path)
    setwd(subs_path)
    # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
    ### First we need the complete functionink file
    dir_path <- "functionink_tmp/"
    file <- list.files(
      path = dir_path,
      pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
      full.names = TRUE
    )

    functionink <- read.table(file, sep = "\t", header = TRUE)
    functionink$ESV <- rownames(functionink)
    colnames(functionink) <- c("mod", "ESV")

  # functionink <- as.data.frame(read_tsv(file,skip = 9,col_names = F))
  # colnames(functionink) <- c("ESV","mod")

  #### observed H and S

  ### Add the family info to the functionink df
  functionink$family <- taxonomy$Family[match(functionink$ESV,rownames(taxonomy))]
  families <- unique(functionink$family)
  heatmap_groups <- unique(functionink$mod)

   
  # Crear la matriz guild_matrix
    guild_matrix <- matrix(0, nrow = length(families), ncol = length(heatmap_groups))
    rownames(guild_matrix) <- families
    colnames(guild_matrix) <- heatmap_groups
    
  # Llenar la matriz con el conteo de OTUs por familia y guild
    for (fam in families) {
      for (group in heatmap_groups) {
        guild_matrix[fam, as.character(group)] <- sum(functionink$family == fam & functionink$mod == group)
      }
    }
    
  ### As now we are interested in the modules
    test <- t(guild_matrix)
    # Crear DataFrame con el conteo y clasificaciones de estrategias
    guild_df <- as.data.frame(guild_matrix)
    guild_df$n_members <- rowSums(guild_df)

  ### Store the total number of families
  total_f <- length(unique(functionink$family))
  n_total_families_per_substrate[[substrate]] <- total_f ## To check the effect of the total number of families
  list_n_families_for_each_module <- list()
  ### Now store the number of families in each of the modules
  x <- 1
    for (module in unique(functionink$mod)) {
      
      functionink_mod <- filter(functionink, mod == module)
      families_module <- length(unique(functionink_mod$family))
      list_n_families_for_each_module[[x]] <- families_module
      names(list_n_families_for_each_module)[x] <- module
      x <- x + 1
    }
  total_f

  p_modules <- purrr::map(list_n_families_for_each_module, \(x) (-x/total_f * log(x/total_f)))


  ##########################

  ### Randomized H and S
  ## I see a lot of variability in the results, les try filtering the taxonomy by the families that i have in the modules
  ## of functionink and then permuting
  taxonomy_filtered <- filter(taxonomy, Family %in% unique(functionink$family))
  ### Add the family info to the functionink df
  randomized_results_list <- list()
  mean_randomized_results_list <- list()
  sd_randomized_results_list <- list()
  n_iterations <- 50
  
    for (i in 1:n_iterations){
      print(i)
      functionink_random <- functionink
      permuted_tax <- permute_taxonomy(taxonomy_filtered)
      functionink_random$family <- permuted_tax$Family[match(functionink_random$ESV,rownames(permuted_tax))]

    ### Store the total number of families
      total_f_random <- length(unique(functionink_random$family))

      list_n_families_for_each_module_random <- list()
    ### Now store the number of families in each of the modules
      x <- 1
        for (module in unique(functionink_random$mod)) {
        
        functionink_mod_random <- filter(functionink_random, mod == module)
        families_module_random <- length(unique(functionink_mod_random$family))
        list_n_families_for_each_module_random[[x]] <- families_module_random
        names(list_n_families_for_each_module_random)[x] <- module
        x <- x + 1
        }

        p_modules_random <- purrr::map(list_n_families_for_each_module_random, \(x) (-x/total_f_random * log(x/total_f_random)))
        randomized_results_list[[i]] <- p_modules_random
    }

  
  randomized_results_list
  sum_randomized_results_list <- reduce(map(randomized_results_list,unlist),`+`) ### 
  mean_randomized_results_list <- map(sum_randomized_results_list, \(x) x/n_iterations)

  sd_randomized_results_list <- map_dbl(seq_along(randomized_results_list[[1]]), function(i) {
    sd(map_dbl(randomized_results_list, ~ .x[[i]]))
  })


  ### Calculate the z_score
  z_score_list <- map_dbl(seq_along(p_modules), function(y) {
  (p_modules[[y]] - mean_randomized_results_list[[y]])/sd_randomized_results_list[[y]]
  })

  z_score_list

  names(z_score_list) <- unique(functionink$mod)
  # names(z_score_list)
  # names(z_score_list) <- paste0(substrate,"_", names(z_score_list))
  results_list[[substrate]] <- z_score_list
    
}
results_list[[2]]

test_3 <- imap_dfr(results_list, function(mod_list, substrate) {
  tibble(
    substrate = substrate,
    module = names(mod_list),
    z_score = unlist(mod_list)
  )
})

setwd("/home/ajf/Desktop/CNB/ecocoherence/figures")


ggplot(test_3, aes(x = substrate, y = z_score, color = z_score)) +
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +
  geom_text(aes(label = module), position = position_jitter(width = 0.2, height = 0), vjust = -0.5) +
  scale_color_gradient2(low = "orange", mid = "black", high = "red", midpoint = 0) +
  labs(x = "Módulo", y = "Z-score", color = "Z-score") +
  theme_classic()

