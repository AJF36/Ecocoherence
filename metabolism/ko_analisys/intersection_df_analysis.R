#################################################
#Script for analyzing the intersection-union df #
#################################################

rm(list = ls())

library(this.path)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

work_dir <- this.dir()
substrates <- c("Alginate","Agarose","AgaroseAlginate","AgaroseCarrageenan","AgaroseChitosan","Carrageenan","Chitin")
# substrates <- "Alginate"
###Load the the matched ESV df
setwd(this.dir())
matched_asv_df <- read_tsv("matched_ESV_id_0.97.tsv",col_names = F)
matched_asv_df <- select(matched_asv_df,c(X9,X10)) %>% 
  filter(X10 != "*")

colnames(matched_asv_df) <- c("ASV","ID")
matched_asv_df$ID <- gsub("[A-Z]{2}_([A-Z]{3}_[0-9]+).[0-9]","\\1",matched_asv_df$ID)


###Load the intersection_union_df

intersection_df <-read_tsv("intersection_union_df_combined.tsv")
intersection_df <- intersection_df %>% 
  separate(col = IDs, into = c("ID1","ID2"),sep = ";") %>%
  separate(col = Familys ,into = c("F1","F2"),sep = "_") 

intersection_df <- intersection_df %>% 
  mutate("union_intersection" = n_union - n_intersection) %>%
  mutate("normalized_difference" =(n_union - n_intersection)/n_union) %>% 
  mutate("overlap" = n_intersection/n_union)
intersection_df$same_family <- ifelse(intersection_df$F1 == intersection_df$F2,1,0)


intersection_df_grouped <- intersection_df %>%
  group_by(same_family) %>%
  summarize("median_overlap" = mean(overlap),"median_normalized_difference" = mean(normalized_difference))


intersection_df_grouped
base_overlap_same_f <- as.numeric(intersection_df_grouped[2,2])
base_overlap_diff_f <- as.numeric(intersection_df_grouped[1,2])
base_difference_same_f <- as.numeric(intersection_df_grouped[2,3])
base_difference_diff_f <- as.numeric(intersection_df_grouped[1,3])



###Distribution plot

intersection_df_same_f <- filter(intersection_df,same_family == 1)
intersection_df_diff_f <- filter(intersection_df,same_family == 0)





# intersection_df_grouped <- intersection_df %>%
#   group_by(same_family) %>%
#   summarize()
colnames(intersection_df)


# "Porphyromonadaceae" %in% intersection_df$F1
# "Porphyromonadaceae" %in% intersection_df$F2
# Get the IDs of all the ASV from funtionink and filter the df
# modules_IDs <- character()

# for(substrate in substrates){
#   print(substrate)
#   #Load the functionink_df 
#   setwd(paste("/home/ajf/Desktop/CNB/ecocoherence/functionink",substrate,"functionink_tmp",sep="/"))
#   file <- list.files(
#     path = ".",
#     pattern = paste("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_",substrate,"_.tsv_guildGT4.txt$",sep=""),
#     full.names = F
#   )
#   # functionink_df <- read_tsv(file,comment = "#",col_names = F)
#   functionink_df <- read.table(file, sep=" " , header = T)
#   functionink_df <- separate(functionink_df,col = guild, into = c("module","ASV"),sep = "\t")
#   colnames(functionink_df) <- c("ASV","module")
#   functionink_df$ID <- matched_asv_df$ID[match(functionink_df$ASV,matched_asv_df$ASV)]
#   ##remove ASV without ID
#   functionink_df <- drop_na(functionink_df)
#   modules_IDs <- c(modules_IDs,functionink_df$ID)
#   length(functionink_df$ID) *7
# 
# }
# functionink_df$ID
# modules_IDs <- unique(modules_IDs)


familys_of_interest <- c("Helicobacteraceae","Nannocystaceae","Campylobacteraceae","Porphyromonadaceae",
                         "Colwelliaceae","Alteromonadaceae","Rhodobacteraceae"
                         ,"Flavobacteriaceae","Oceanospirillaceae","Syntrophaceae","Anaerolineaceae","Desulfobacteraceae","Psychromonadaceae",
                         "Hyphomonadaceae","Phyllobacteriaceae","Clostridiales",
                         "Puniceicoccaceae","Alteromonadales","Pseudomonadales",
                         "Rhodospirillaceae","Cocleimonas","Planctomycetaceae","Verrucomicrobiaceae","Vibrionaceae")

coherency_families_of_interest <- c("early_c","late_c","early_c","late_c","early_c","incoherent","incoherent","incoherent","incoherent","late_c",
                                    "late_c","incoherent","early_c","incoherent","incoherent","incoherent","incoherent","incoherent","incoherent","late_c","early_c",
                                    "late_c","late_c","early_c")
df_for_color <- data.frame("family" = familys_of_interest,"coherency" = coherency_families_of_interest)
intersection_df_f <- filter(intersection_df,F1 %in% familys_of_interest & F2 %in% familys_of_interest)




# intersection_df_f <- filter(intersection_df,F1 == "Porphyromonadaceae")
# any(intersection_df_f$ID1 %in% modules_IDs) & any(intersection_df_f$ID2 %in% modules_IDs)
# intersection_df_f$ID2 %in% modules_IDs
# intersection_df_f$ID1 %in% modules_IDs
# "Porphyromonadaceae" %in% intersection_df_f$F2
intersection_df_f<- intersection_df_f %>%
  mutate("union_intersection" = n_union - n_intersection) %>%
  mutate("normalized_difference" =(n_union - n_intersection)/n_union) %>% 
  mutate("overlap" = n_intersection/n_union)

intersection_df_f$same_family <- ifelse(intersection_df_f$F1 == intersection_df_f$F2, 1, 0)

intersection_df_f_gruped <- intersection_df_f %>%
  group_by(same_family) %>%
  summarize("median_overlap" = median(overlap),"median_normalized_difference" = median(normalized_difference))

###
same_family_df <- filter(intersection_df_f,F1 == F2)
###Add a coherency column


same_family_df$coherency <- ifelse(same_family_df$F1 %in% df_for_color$family,df_for_color$coherency[match(same_family_df$F1,df_for_color$family)],"Other")

family_order <- same_family_df %>%
  group_by(F1) %>%
  # Calculate the most common coherency category or an average score
  arrange(coherency) %>%  # Sort families by coherency score
  pull(F1)

family_order <- unique(family_order)
family_order

  
  
same_family_df$F1 <- factor(same_family_df$F1, 
                                        levels = family_order)




grouped_family_df <- group_by(same_family_df, F1) %>%
  summarize("median_overlap" = median(overlap),"median_difference" = median(normalized_difference))
grouped_family_df$coherence <- ifelse(grouped_family_df$F1 %in% df_for_color$family,df_for_color$coherency[match(grouped_family_df$F1,df_for_color$family)],"Other")







ggplot(grouped_family_df,aes(F1,median_overlap,colour =  coherence)) +
  geom_point(aes(size = 2)) +
  geom_hline(yintercept = base_overlap_same_f, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 10))


ggplot(same_family_df,aes(F1,overlap,colour = coherency)) +
  geom_boxplot() +
  geom_hline(yintercept = base_overlap_same_f, linetype = "dashed", color = "black")
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 10))

ggplot(same_family_df,aes(F1,normalized_difference,colour = coherency)) +
  geom_boxplot() +
  geom_hline(yintercept = base_difference_same_f, linetype = "dashed", color = "black")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 10))

ggplot(grouped_family_df,aes(F1,median_difference,colour =  coherence)) +
  geom_point(aes(size = 2)) +
  geom_hline(yintercept = base_difference_same_f, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 10))



same_family_df$pair <- paste(same_family_df$F1,same_family_df$F2,sep = "_") 


different_f_df <- intersection_df_f %>%
  filter(F1 != F2)
different_f_df$pair <- paste(different_f_df$F1,different_f_df$F2,sep="_")
different_f_df$C1 <- df_for_color$coherency[match(different_f_df$F1,df_for_color$family)]
different_f_df$C2 <- df_for_color$coherency[match(different_f_df$F2,df_for_color$family)]
# different_f_df$Coherency <- paste(different_f_df$C1,different_f_df$C2,sep="_")
different_f_df$Coherency <- apply(different_f_df[, c("C1", "C2")], 1, function(x) {
  sorted_values <- sort(x)
  paste(sorted_values[1], sorted_values[2], sep="_")
})


###Density curves plot 

same_family_coherent_early_df <- filter(same_family_df,coherency == "early_c")
same_family_coherent_late_df <- filter(same_family_df,coherency == "late_c")
same_family_incoherent_df <- filter(same_family_df,coherency == "incoherent")
unique(same_family_df$coherency)


# First, create a combined dataset for density curves
combined_data <- bind_rows(
  data.frame(overlap = intersection_df_diff_f$overlap, group = "Different Family"),
  data.frame(overlap = intersection_df_same_f$overlap, group = "Same Family")
)

# Create named vectors for vertical lines to use in the legend
vline_means <- c(
  "Coherent Early" = mean(same_family_coherent_early_df$overlap),
  "Incoherent" = mean(same_family_incoherent_df$overlap),
  "Coherent Late" = mean(same_family_coherent_late_df$overlap),
  "all_familys_interest" = mean(intersection_df_f$overlap)
)

vline_colors <- c(
  "Coherent Early" = "blue",
  "Incoherent" = "green",
  "Coherent Late" = "red",
  "all_familys_interest" = "purple"
)

# Create the plot
ggplot() +
  # Add density curves with proper grouping
  geom_density(data = combined_data, 
               aes(x = overlap, fill = group),
               alpha = 0.5) +
  # Add vertical lines
  geom_vline(data = data.frame(
    mean = vline_means,
    group = names(vline_means)
  ),
  aes(xintercept = mean, color = group),
  linetype = "dotted", size = 1.2) +
  # Set colors manually
  scale_fill_manual(name = "Groups",
                    values = c("Different Family" = "yellow", 
                               "Same Family" = "brown")) +
  scale_color_manual(name = "Mean Values",
                     values = vline_colors) +
  # Add labels and theme
  theme_bw() +
  xlab("Overlap") +
  # Adjust legend position and formatting
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )






# ggplot() +
#   geom_density(aes(x = intersection_df_diff_f$overlap), fill = "yellow", alpha = 0.5) +
#   geom_density(aes(x = intersection_df_same_f$overlap), fill = "brown", alpha = 0.5) +
#   # geom_density(aes(x = same_family_coherent_df$overlap), fill = "green",alpha = 0.5) +
#   # geom_density(aes(x = same_family_incoherent_df$overlap, fill = "purple",alpha = 0.5))
#   geom_vline(xintercept = mean(same_family_coherent_early_df$overlap),color = "blue", linetype = "dotted") +
#   geom_vline(xintercept = mean(same_family_incoherent_df$overlap), color = "green" , linetype = "dotted") +
#   geom_vline(xintercept = mean(same_family_coherent_late_df$overlap), color = "red" , linetype = "dotted")
#   theme_bw() +
#   xlab("overlap") +
#   theme(legend.location = "top")



different_f_df_grouped <- different_f_df %>%
  group_by(Coherency) %>%
  summarize("median_overlap" = median(overlap),"median_normalized_difference" = median(normalized_difference))

ggplot(different_f_df_grouped,aes(Coherency,median_overlap,colour = Coherency)) +
  geom_point(aes(size = 2))+
  geom_hline(yintercept = base_overlap_diff_f, linetype = "dashed", color = "black") +
  theme_bw()

ggplot(different_f_df_grouped,aes(Coherency,median_normalized_difference,colour = Coherency)) +
  geom_hline(yintercept = base_difference_diff_f, linetype = "dashed", color = "black") +
  geom_point(aes(size = 2))+
  theme_bw()



###Modules analysis

#################

modules_df <- data.frame("substrate" = NA, "overlap_within_same_f" = NA, "overlap_within_diff_f" = NA,
                         "overlap_between_same_f" = NA,"overlap_between_d_f" = NA)
for(substrate in substrates){
  #Load the functionink_df
  setwd(paste("/home/ajf/Desktop/CNB/ecocoherence/functionink",substrate,"functionink_tmp",sep="/"))
  file <- list.files(
    path = ".",
    pattern = paste("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_",substrate,"_.tsv_guildGT4.txt$",sep=""),
    full.names = F
  )
  # functionink_df <- read_tsv(file,comment = "#",col_names = F)
  functionink_df <- read.table(file, sep=" " , header = T)
  functionink_df <- separate(functionink_df,col = guild, into = c("module","ASV"),sep = "\t")
  colnames(functionink_df) <- c("ASV","module")
  ##Add id information to the functionink table
  functionink_df$ID <- matched_asv_df$ID[match(functionink_df$ASV,matched_asv_df$ASV)]
  functionink_df <- drop_na(functionink_df)
  head(functionink_df)
  ###Filter the full df
  substrate_df <- filter(intersection_df, ID1 %in% functionink_df$ID & ID2 %in% functionink_df$ID) 
  # substrate_df <- filter(substrate_df,F1 %in% familys_of_interest & F2 %in% familys_of_interest)
  ggplot(substrate_df,aes(overlap)) +
    geom_density()
  # substrate_df <- substrate_df %>% mutate("union_intersection" = n_union - n_intersection) %>%
  #   mutate("normalized_difference" =(n_union - n_intersection)/n_union) %>%
  #   mutate("overlap" = n_intersection/n_union)
  substrate_df$mod_1 <- functionink_df$module[match(substrate_df$ID1,functionink_df$ID)]
  substrate_df$mod_2 <- functionink_df$module[match(substrate_df$ID2,functionink_df$ID)]
  substrate_df$same_mod <- ifelse(substrate_df$mod_1 == substrate_df$mod_2,1,0)
  substrate_df$pair <- paste(substrate_df$F1,substrate_df$F2,sep = "_")
  paste(min(substrate_df$ID1,substrate_df$ID2),max(substrate_df$ID1,substrate_df$ID2))

  grouped_substrate_df_same_mod_diff_f <- substrate_df %>%
    filter(mod_1 == mod_2 & F1 != F2 ) %>%
    group_by(pair) %>% 
    summarize("median_overlap" = median(overlap), "median_diff" = median(normalized_difference))

  grouped_substrate_df_same_mod_same_f <- substrate_df %>%
    filter(mod_1 == mod_2 & F1 == F2 ) %>%
    group_by(pair) %>% 
    summarize("median_overlap" = median(overlap), "median_diff" = median(normalized_difference))
  # 
  # ggplot(grouped_substrate_df_same_mod,aes(overlap)) +
  #   geom_density()
  
  grouped_substrate_df_diff_mod_diff_f <- substrate_df %>%
    filter(mod_1 != mod_2 & F1 != F2) %>%
    group_by(pair) %>% 
    summarize("median_overlap" = median(overlap), "median_diff" = median(normalized_difference))
  
  grouped_substrate_df_diff_mod_same_f <- substrate_df %>%
    filter(mod_1 != mod_2 & F1 == F2) %>%
    group_by(pair) %>% 
    summarize("median_overlap" = median(overlap), "median_diff" = median(normalized_difference))
  
  same_module_same_f_median_overlap <- median(grouped_substrate_df_same_mod_same_f$median_overlap)
  same_module_diff_f_median_overlap <-median(grouped_substrate_df_same_mod_diff_f$median_overlap)
  diff_module_same_f_median_overlap <-median(grouped_substrate_df_diff_mod_same_f$median_overlap)
  diff_module_diff_f_median_overlap <-median(grouped_substrate_df_diff_mod_diff_f$median_overlap)


  same_module_same_f_median_diff <- median(grouped_substrate_df_same_mod_same_f$median_diff)
  same_module_diff_f_median_diff <-median(grouped_substrate_df_same_mod_diff_f$median_diff)
  diff_module_same_f_median_diff <-median(grouped_substrate_df_diff_mod_same_f$median_diff)
  diff_module_diff_f_median_diff <-median(grouped_substrate_df_diff_mod_diff_f$median_diff)
  

  
  dummy_df <- data.frame("substrate" = substrate, "overlap_within_same_f" = same_module_same_f_median_overlap, "overlap_within_diff_f" = same_module_diff_f_median_overlap,
                         "overlap_between_same_f" = diff_module_same_f_median_overlap,"overlap_between_d_f" = diff_module_diff_f_median_overlap)
                        
  modules_df <- rbind(dummy_df,modules_df)
  # ###Create df with the scenearios 
  # same_module_same_f_df <- filter(substrate_df,mod_1 == mod_2 & F1 == F2)
  # same_module_diff_f_df <- filter(substrate_df,mod_1 == mod_2 & F1 != F2)
  # diff_module_diff_f_df <- filter(substrate_df,mod_1 != mod_2 & F1 != F2)
  # # diff_module_diff_f_df <- filter(substrate_df,mod_1 == mod_2 & F1 != F2)
  # # diff_module_diff_f_df <- filter(substrate_df,mod_1 == mod_2 & F1 != F2)
  # 
  # 
  # ####overlap bewtween same module same family
  # grouped_df_mods_overlap <- group_by(same_module_same_f_df,mod_1) %>%
  #   summarize("median_overlap" = median(overlap))
  # grouped_df_mods_intersection <- group_by(same_module_same_f_df,mod_1) %>%
  #   summarize("median_intersection" = median(normalized_difference))
  # grouped_df_fams_overlap <- group_by(same_module_same_f_df,pair) %>%
  #   summarize("median_overlap" = median(overlap))
  # grouped_df_fams_intersection <- group_by(same_module_same_f_df,pair) %>%
  #   summarize("median_intersection" = median(normalized_difference))
  # ### overlap between same module different family
  # grouped_df_mods2 <-  group_by(same_module_diff_f_df,mod_1) %>%
  #   summarize("median_intersection" = median(normalized_difference))



  # ###plot overlap familys
  # ggplot(same_module_same_f_df,aes(F1,overlap)) +
  #   geom_boxplot()
  # ###
  # 
  # grouped_df_family <- group_by(substrate_df,pair) %>%
  #   summarize("median_normalized_difference" = median(normalized_difference))
  # colnames(substrate_df)


  ###
}

modules_df <- drop_na(modules_df)
modules_df_long <- pivot_longer(modules_df,cols = c(overlap_within_same_f,overlap_within_diff_f,overlap_between_same_f,overlap_between_d_f),
                                values_to = "values",names_to = "variable")
ggplot(modules_df,aes(x = substrate)) +
  geom_point(aes(y=overlap_within_same_f), shape = 1,size = 2) +
  geom_point(aes(y=overlap_within_diff_f),shape = 2, size = 2) +
  geom_point(aes(y=overlap_between_same_f),shape = 3,size = 2) +
  geom_point(aes(y=overlap_between_d_f), shape = 5,size=2)

ggplot(modules_df_long,aes(substrate,values,colour = variable,shape = variable))+
  geom_point(size= 6) +
  theme_bw() +
  ylab("overlap")
# head()

###No separation by substrate
