##### Script to calculate the Z_score by module
rm(list = ls())
library(tidyverse)
library(this.path)

work_dir <- this.dir()
setwd(work_dir)

observed_df <- read_tsv("entropy_by_modules.tsv",col_names = T)
random_df <- read_tsv("randomized_entropy_by_modules.tsv",col_names = T)

## Create the df with all the information for the calculation of the Z_score
full_df <- full_join(observed_df,random_df)

full_df <- full_df %>%
  mutate(Z = (value - mean_X) / sd_X)


### Add the info of the ecological strategy

# load the df
ecological_strategy_df_path <- normalizePath(file.path("../modules_similarity"))
setwd(ecological_strategy_df_path)
list_files <-list.files(path = ".",pattern = "modules_classified")
ecological_df_list <- list()
for ( file in list_files) {
  ecological_strategy_df_substrate <- read_tsv(file)
  substrate <- sub(".*internal(.*)\\.tsv$", "\\1", file)
  ecological_strategy_df_substrate$substrate <- substrate
  ecological_df_list[[substrate]] <- ecological_strategy_df_substrate
}



data_frame_ecological_strategy <- do.call(rbind,ecological_df_list)
colnames(data_frame_ecological_strategy)[[1]] <- "module"
data_frame_ecological_strategy$module <- gsub("mod_","",data_frame_ecological_strategy$module)





ggplot(full_df, aes(x = substrate, y = Z , color = Z)) +
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +
  geom_text(aes(label = module), position = position_jitter(width = 0.2, height = 0), vjust = -0.5) +
  scale_color_gradient2(low = "orange", mid = "black", high = "red", midpoint = 0) +
  labs(x = "Módulo", y = "Z-score", color = "Z-score") +
  theme_classic()

### Lets try to see the difference of the modules strategy assignation

full_df$module <- as.character(full_df$module)
selected_modules_df <- inner_join(full_df,data_frame_ecological_strategy, by = c("module", "substrate"))

fig_path <- normalizePath(file.path("../figures"))

setwd(fig_path)
pdf("coherence_modules_dotplot.pdf",12,12)
ggplot(selected_modules_df, aes(x = substrate, y = Z , color = Strategy)) +
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +
  geom_text(aes(label = module), position = position_jitter(width = 0.2, height = 0), vjust = -0.5) +
  # scale_color_gradient2(low = "orange", mid = "black", high = "red", midpoint = 0) +
  theme_classic() + 
  labs(x = "Substrate", y = "Z-score value", color = "Ecological Strategy") +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size =25),
axis.text.x = element_text(size = 15,angle = 90), axis.text.y = element_text(size = 15),
legend.text = element_text(size = 15), legend.title = element_text(size = 20))
dev.off()


full_df_f <- full_df %>%
  select(substrate,module,Z)

setwd(work_dir)
write_tsv(full_df_f,"z_score_modules.tsv")