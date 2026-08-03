#############################################################
#
#
#
rm(list = ls())

library(tidyverse)
library(this.path)
library(gt)
library(gtExtras)


script_dir <- this.dir()
figures_dir <- normalizePath(file.path(script_dir, "..", "figures", "9_summarize_strategy_vs_coherence_links"), mustWork = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
family_strategy_dir <- normalizePath(file.path(script_dir, "..", "results", "6g_assign_ecological_strategy_to_families"), mustWork = FALSE)
module_strategy_dir <- normalizePath(file.path(script_dir, "..", "results", "6i_classify_modules_by_ecological_strategy_combined"), mustWork = FALSE)
family_zscore_dir <- normalizePath(file.path(script_dir, "..", "results", "8a_compute_zscore_families"), mustWork = FALSE)
module_zscore_dir <- normalizePath(file.path(script_dir, "..", "results", "8b_compute_zscore_modules"), mustWork = FALSE)
links_table_dir <- normalizePath(file.path(script_dir, "..", "results", "6b_build_family_module_table_with_percentages"), mustWork = FALSE)
setwd(script_dir)



### Function to select columns without stoping the script
   safe_select <- function(df, families_strategy, modules_strategy) {
  tryCatch({
    df[rownames(df) %in% families_strategy, ] %>%
      select(any_of(modules_strategy[!is.na(modules_strategy)]))
  }, error = function(e) {
    message("Warning: selection failed for one table — ", e$message)
    return(data.frame())  # return empty df so code keeps running
  })
}

  
### Load the table with the information of the ecological strategy of the families
strategy_family_table <- read_tsv(file = file.path(family_strategy_dir, "ecological_strategies_families.tsv") ,col_names = T)

### Load the table with the information of the ecological strategy of the modules
strategy_module_table <- read_tsv(file = file.path(module_strategy_dir, "ecological_strategies_modules.tsv") ,col_names = T)
strategy_module_table$module <- gsub("mod_","",strategy_module_table$module)


### Load the table with the information of the coherence of the families
coherence_family_table_path <- normalizePath(file.path(family_zscore_dir,"z_score_families.tsv"))
coherence_family_table <- read_tsv(coherence_family_table_path)

#Filter the statistically significant families
coherence_family_table <- coherence_family_table%>%
  mutate("coherency" = ifelse(abs(Z) >= 2,"coherent","incoherent"))



### Load the table with the information of the coherence of the modules
coherence_module_table_path <- normalizePath(file.path(module_zscore_dir,"z_score_modules.tsv"))
coherence_module_table <- read_tsv(coherence_module_table_path)
coherence_module_table$module <- as.character(coherence_module_table$module)
# coherence_module_table$module <- paste(coherence_module_table$substrate,"_",coherence_module_table$module, sep = "")
coherence_module_table <- coherence_module_table%>%
  mutate("coherency" = ifelse(abs(Z) >= 2,"coherent","incoherent"))
  
full_df_module <- full_join(coherence_module_table, strategy_module_table, by = c("substrate","module"))

full_df_family <- full_join(coherence_family_table, strategy_family_table, by = c("substrate","family"))
# # Filter the statistically significant modules
# coherence_module_table <- coherence_module_table%>%
#   filter(Z >= 2 | Z <= -2)

gt_table_family <- gt(as.data.frame(table(full_df_family$coherency,full_df_family$strategy,useNA = "no",dnn = c("Coherency","E.Strategy")))) %>%
  tab_header(title = md("**By family**"))

gt_table_module <- gt(as.data.frame(table(full_df_module$coherency,full_df_module$strategy,useNA = "no", dnn =c("Coherency","E.Strategy")))) %>%
    tab_header(title = md("**By module**"))
setwd(figures_dir)
gtExtras::gtsave_extra(gt_table_family,"pairwise_comparation_str_cohe_families.pdf")
gtExtras::gtsave_extra(gt_table_module,"pairwise_comparation_str_cohe_module.pdf")
setwd(script_dir)

substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin","Carrageenan")

strategies <- c("selection","generalist","facilitation","attachment","transition")
  # strategies <- unique(substrate_family_table$strategy)
  # strategies <- strategies[!is.na(strategies)]
# sub <- "Alginate"

results_df <- data.frame("s" = character(),"mc" =character(),"fc" = character(),"ms" = character(),"fs" = character(), "n_links" = numeric())

for (sub in substrates){
  print(sub)
  ### Df to store the results of the substrate
  test <- expand.grid("s" = sub,"mc" =c("coherent","incoherent"),"fc" = c("coherent","incoherent"),"ms" =strategies,"fs" = strategies, "n_links" = 0)
  test$mc <- as.character(test$mc)
  test$fc <- as.character(test$fc)
  test$ms <- as.character(test$ms)
  test$fs <- as.character(test$fs)
  ### Load the table with the information of which families goes to which modules
  table_path  <- file.path(links_table_dir,paste("group_analysis_",sub,"_TEST.tsv", sep = ""))
  links_table <- read_tsv(normalizePath(table_path))
  links_table_f <- select(links_table, !matches("members|Family"))
  rownames(links_table_f) <- links_table$Family

  ### Filter the full tables
  substrate_module_table <- filter(full_df_module, substrate == sub)

  substrate_family_table <- filter(full_df_family, substrate == sub)


  ### Lets only count if a family has a link with a substrate
  links_table_f <- as.data.frame(apply(links_table_f,2, \(x)ifelse(x >= 1, 1,0)))
  total_links <- sum(links_table_f == 1)
  


    coherent_modules <- substrate_module_table$module[substrate_module_table$coherency == "coherent"]
    coherent_modules <- coherent_modules[!is.na(coherent_modules)]
    coherent_modules <- coherent_modules[coherent_modules %in% colnames(links_table_f)]

    coherent_families <-substrate_family_table$family[substrate_family_table$coherency == "coherent"]
    coherent_families <- coherent_families[!is.na(coherent_families)]
  
  
    links_table_f_cm <- links_table_f %>%
      select(all_of(coherent_modules)) 

    links_table_f_im <- links_table_f %>%
      select(!all_of(coherent_modules)) 

    links_table_f_cm_cf <- links_table_f %>%
      select(all_of(coherent_modules)) %>%
      filter(rownames(links_table_f) %in% coherent_families)

    links_table_f_cm_if <- links_table_f %>%
      select(all_of(coherent_modules)) %>%
      filter(!rownames(links_table_f) %in% coherent_families)

    links_table_f_im_cf <- links_table_f %>%
      select(!all_of(coherent_modules)) %>%
      filter(rownames(links_table_f) %in% coherent_families)
      
    links_table_f_im_if <- links_table_f %>%
      select(!all_of(coherent_modules)) %>%
      filter(!rownames(links_table_f) %in% coherent_families)


    # sum(links_table_f_cm == 1)/total_links #Nfcoherent modules/Nfamilies
    #     sum(links_table_f_im == 1)/total_links #Nf incoherent modules/Nfamilies

    # ##For coherent modules
    sum(links_table_f_cm_cf == 1)/sum(links_table_f_cm == 1) ### NcoherentF / NtotalF
    sum(links_table_f_cm_if == 1)/sum(links_table_f_cm == 1) ### NincoherentF / NtotalF
    
    ##For incoherent modules
    sum(links_table_f_im_cf == 1)/sum(links_table_f_im == 1) ### NcoherentF / NtotalF
    sum(links_table_f_im_if == 1)/sum(links_table_f_im == 1) ### NincoherentF / NtotalF
    

    # stra <-"facilitation"

    for (mod_stra in strategies){
      print(paste("The module strategy is: ",mod_stra, sep = ""))

      modules_strategy <- substrate_module_table$module[substrate_module_table$strategy == mod_stra]
      
      for(family_stra in strategies){
        print(paste("The family strategy is: ",family_stra, sep = ""))

        families_strategy <- substrate_family_table$family[substrate_family_table$strategy == family_stra]
           
        links_table_f_cm_cf_s <- safe_select(links_table_f_cm_cf, families_strategy, modules_strategy)
        links_table_f_cm_if_s <- safe_select(links_table_f_cm_if, families_strategy, modules_strategy)
        links_table_f_im_cf_s <- safe_select(links_table_f_im_cf, families_strategy, modules_strategy)
        links_table_f_im_if_s <- safe_select(links_table_f_im_if, families_strategy, modules_strategy)

      
        print(sum(links_table_f_cm_cf_s == 1))### NcoherentF / NtotalF
        print(sum(links_table_f_cm_if_s == 1)) ### NincoherentF / NtotalF
    
      ##For incoherent modules
        print(sum(links_table_f_im_cf_s == 1)) ### NcoherentF / NtotalF
        print(sum(links_table_f_im_if_s == 1)) ### NincoherentF / NtotalF
        
        test$n_links[test$mc == "coherent" & test$fc == "coherent" & test$ms == mod_stra & test$fs == family_stra] <- sum(links_table_f_cm_cf_s == 1)
        test$n_links[test$mc == "coherent" & test$fc == "incoherent" & test$ms == mod_stra & test$fs == family_stra] <- sum(links_table_f_cm_if_s == 1)
        test$n_links[test$mc == "incoherent" & test$fc == "coherent" & test$ms == mod_stra & test$fs == family_stra]  <- sum(links_table_f_im_cf_s == 1)
        test$n_links[test$mc == "incoherent" & test$fc == "incoherent" & test$ms == mod_stra & test$fs == family_stra]  <- sum(links_table_f_im_if_s == 1)




      }
      
      
      
  
      # links_table_f_cm_cf_s <- links_table_f_cm_cf[rownames(links_table_f_cm_cf) %in% families_strategy,] %>%
      #   select(any_of(modules_strategy[!is.na(modules_strategy)]))
    
      # links_table_f_cm_if_s <- links_table_f_cm_if[rownames(links_table_f_cm_cf) %in% families_strategy,] %>%
      #   select(modules_strategy[!is.na(modules_strategy)])
      
      # links_table_f_im_cf_s <- links_table_f_im_cf[rownames(links_table_f_cm_cf) %in% families_strategy,] %>%
      #   select(modules_strategy[!is.na(modules_strategy)])
      # links_table_f_im_if_s <- links_table_f_im_if[rownames(links_table_f_cm_cf) %in% families_strategy,] %>%
      #   select(modules_strategy[!is.na(modules_strategy)])





    }
  results_df <- rbind(results_df,test)
}


# str(test)
#                         ###Estos dos van hardcoded           #Cambiar esto for mod_strategy   #cambiar esto por fam_strategy
# test$n_links[test$mc == "coherent" & test$fc == "coherent" & test$ms == mod_stra & test$fs == family_stra] 
# test$n_links[test$mc == "coherent" & test$fc == "incoherent" & test$ms == mod_stra & test$fs == family_stra] 
# test$n_links[test$mc == "incoherent" & test$fc == "coherent" & test$ms == mod_stra & test$fs == family_stra] 
# test$n_links[test$mc == "incoherent" & test$fc == "incoherent" & test$ms == mod_stra & test$fs == family_stra] 

