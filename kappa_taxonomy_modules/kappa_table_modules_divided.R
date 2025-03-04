# Librerías necesarias
rm(list=ls())
library(readr)
library(dplyr)
library(tidyr)
library(irr)
library(ggplot2)
library(envalysis)
library(viridis)
library(RColorBrewer)
library(this.path)
work_dir <- this.dir()
# Lista de sustratos
substrates <- c("Alginate", "Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan")

# Cargar tabla de taxonomía
# fileTaxonomy <- "/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
# taxa.in <- read.table(fileTaxonomy, sep = ";")
# colnames(taxa.in) <- c("taxa_id", "none", "Kingdom", "sig_Kingdom", "Phylum", "sig_Phylum", "Class", "sig_Class",
#                        "Order", "sig_Order", "Family", "sig_Family", "Genus", "sig_Genus")
# taxonomy <- subset(taxa.in, select = c("taxa_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
# taxonomy$otu <- taxonomy$taxa_id
# taxonomy <- subset(taxonomy, select = c('otu', 'Family', 'Genus', 'Order', 'Class'))
fileTaxonomy <- file.path(work_dir, "..", "data", "marine_particles_source_data","sequence_table.ESV.fasta_RDPclassified.txt")
fileTaxonomy <- normalizePath(fileTaxonomy, mustWork = FALSE)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
taxonomy$otu <- taxonomy$taxa_id
taxonomy <- subset(taxonomy, select = c('otu', 'Family', 'Genus', 'Order', 'Class'))

# rownames(taxonomy)=taxonomy$taxa_id
# taxonomy <- subset(taxonomy, select = c('otu', 'Family', 'Genus', 'Order', 'Class'))

# taxonomy=subset(taxonomy,select=-c(taxa_id))
# tax.pseq = tax_table(as.matrix(taxonomy))


# Cargar las tablas de módulos de cada sustrato en una lista
# otus_modules <- list()
# for (substrate in substrates) {
#   setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
#   dir_path <- "functionink_tmp/"
#   file_name <- list.files(
#     path = dir_path,
#     pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
#     full.names = TRUE
#   )
#   modules_table <- read_tsv(file_name)
#   modules_table <- modules_table %>% separate(col = guild, into = c("otu", "module"), sep = "\t")
#   otus_modules[[substrate]] <- modules_table
# }
# taxonomy$Family[taxonomy$otu == otus_modules$otu]
# otus_modules ###Esta tabla debe incluir una fila que sea module_by_family en la cual sea la union del modulo y la familia de ese otu en este formato mod_2_family
# Paso 1: Obtener todos los OTUs presentes en los módulos de cualquier sustrato
############
otus_modules <- list()
for (substrate in substrates) {
  
  
  subs_path <- file.path(work_dir,"..","functionink",substrate)
  subs_path <- normalizePath(subs_path)
  setwd(subs_path)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  dir_path <- "functionink_tmp/"
  file_name <- list.files(
    path = dir_path,
    pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
    full.names = TRUE
  )
  # # Cambiar directorio
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  # dir_path <- "functionink_tmp/"
  # 
  # # Buscar archivo que coincide con el patrón
  # file_name <- list.files(
  #   path = dir_path,
  #   pattern = paste0("Partition-NL_Average_StopStep-\\d+_interactions_filtered_p0.01_threshold_", substrate, "_.tsv_guildGT4.txt"),
  #   full.names = TRUE
  # )
  # 
  # Leer tabla de módulos
  modules_table <- read_tsv(file_name)
  modules_table <- modules_table %>% separate(col = guild, into = c("otu", "module"), sep = "\t")
  
  # Unir la tabla de módulos con la tabla de taxonomía para agregar la columna 'module_by_family'
  modules_table <- modules_table %>%
    left_join(taxonomy, by = "otu") %>%
    mutate(module_by_family = paste(module, "mod_", Family, "_family", sep = ""))  # Crear la nueva columna combinada
  
  # Añadir la tabla de módulos al listado de módulos por sustrato
  otus_modules[[substrate]] <- modules_table
}



otus_modules



############
all_otus <- unique(unlist(lapply(otus_modules, function(x) x$otu)))

# Paso 2: Generar todas las combinaciones posibles de pares de OTUs
otu_pairs <- expand.grid(otu1 = all_otus, otu2 = all_otus) %>%
  filter(otu1 != otu2) %>%
  distinct()

# Paso 3: Agregar columnas de familia, género, orden y clase
otu_pairs <- otu_pairs %>%
  left_join(taxonomy, by = c("otu1" = "otu")) %>%
  rename(Family1 = Family, Genus1 = Genus, Order1 = Order, Class1 = Class) %>%
  left_join(taxonomy, by = c("otu2" = "otu")) %>%
  rename(Family2 = Family, Genus2 = Genus, Order2 = Order, Class2 = Class)

# Paso 4: Asignar valores para familia, género, orden y clase (1 si coinciden, 0 si no)
otu_pairs <- otu_pairs %>%
  mutate(
    Family = ifelse(Family1 == Family2, 1, 0),
    Genus = ifelse(Genus1 == Genus2, 1, 0),
    Order = ifelse(Order1 == Order2, 1, 0),
    Class = ifelse(Class1 == Class2, 1, 0)
  ) %>%
  select(-Family1, -Family2, -Genus1, -Genus2, -Order1, -Order2, -Class1, -Class2)

# Paso 5: Crear una columna para cada sustrato y asignar los valores de coincidencia en módulos
for (substrate in substrates) {
  module_data <- otus_modules[[substrate]] %>% select(otu, module_by_family)
  module_data
  otu_pairs
  # Unir el módulo de cada OTU en el par
  pairs_with_modules <- otu_pairs %>%
    left_join(module_data, by = c("otu1" = "otu")) %>%
    rename(module1 = module_by_family) %>%
    left_join(module_data, by = c("otu2" = "otu")) %>%
    rename(module2 = module_by_family)
  
  # Asignar los valores según coincidencia de módulos
  otu_pairs[[substrate]] <- ifelse(
    is.na(pairs_with_modules$module1) | is.na(pairs_with_modules$module2), NA,
    ifelse(pairs_with_modules$module1 == pairs_with_modules$module2, 1, 0)
  )
}

# Data frame para almacenar los resultados de kappa
kappa_results_df <- data.frame(Substrate = character(), Taxon = character(), Kappa = numeric(), stringsAsFactors = FALSE)

# Calcular el kappa para cada sustrato y cada nivel taxonómico y almacenar en el data frame
for (substrate in substrates) {
  otu_pairs_filtered <- otu_pairs %>%
    filter(!is.na(.data[[substrate]])) %>% 
    select(Family, Genus, Order, Class, .data[[substrate]]) %>%
    rename(Module = .data[[substrate]])
  # Asegurarnos de que las columnas sean factores
  otu_pairs_filtered$Module <- factor(otu_pairs_filtered$Module, levels = c(0, 1))
  otu_pairs_filtered$Family <- factor(otu_pairs_filtered$Family, levels = c(0, 1))
  otu_pairs_filtered$Genus <- factor(otu_pairs_filtered$Genus, levels = c(0, 1))
  otu_pairs_filtered$Order <- factor(otu_pairs_filtered$Order, levels = c(0, 1))
  otu_pairs_filtered$Class <- factor(otu_pairs_filtered$Class, levels = c(0, 1))
  
  # Calcular y agregar los resultados de kappa al data frame
  kappa_results_df <- rbind(kappa_results_df, data.frame(
    Substrate = substrate,
    Taxon = "Family",
    Kappa = kappa2(otu_pairs_filtered %>% select(Family, Module))$value
  ))
  
  kappa_results_df <- rbind(kappa_results_df, data.frame(
    Substrate = substrate,
    Taxon = "Genus",
    Kappa = kappa2(otu_pairs_filtered %>% select(Genus, Module))$value
  ))
  
  kappa_results_df <- rbind(kappa_results_df, data.frame(
    Substrate = substrate,
    Taxon = "Order",
    Kappa = kappa2(otu_pairs_filtered %>% select(Order, Module))$value
  ))
  
  kappa_results_df <- rbind(kappa_results_df, data.frame(
    Substrate = substrate,
    Taxon = "Class",
    Kappa = kappa2(otu_pairs_filtered %>% select(Class, Module))$value
  ))
}

# Ver la tabla final de resultados de kappa
View(kappa_results_df)
# Ordenar el factor 'Taxon' de acuerdo con la jerarquía taxonómica en orden descendente
kappa_results_df$Taxon <- factor(kappa_results_df$Taxon, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))

# Graficar los valores de kappa con taxones en el eje x y líneas para cada sustrato
ggplot(kappa_results_df, aes(x = Taxon, y = Kappa, color = Substrate, group = Substrate)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
       x = "Taxonomic levels",
       y = "Kappa statistics",
       color = "Substrate") +
  scale_color_brewer(palette = "Set1") +  ###Color blind palette
  theme_bw(base_size = 16)
  theme(axis.text.x = element_text(angle = 45, hjust = 1 , face = "bold")) 

figures_path <- file.path(work_dir,"..","figures")
figures_path <- normalizePath(figures_path)
setwd(figures_path)

# Guardar el gráfico como un archivo de alta calidad para un artículo
ggsave(
  filename = "kappa_statistics_modules_by_family.png",  # Nombre del archivo
  plot = last_plot(),                 # El último gráfico generado por ggplot
  width = 170,                        # Ancho en mm (ajusta según una o dos columnas)
  height = 100,                       # Altura en mm
  units = "mm",                       # Unidades
  dpi = 300                           # Resolución
)
###Lets try manually for alginate and family
otu_pairs_alginate <- otu_pairs %>% select(-c("Agarose", "AgaroseAlginate", "AgaroseCarrageenan", "AgaroseChitosan", "Chitin", "Carrageenan"))
otu_pairs_alginate <- drop_na(otu_pairs_alginate)
anyNA(otu_pairs_alginate)
# ### Define all the parameters
#   Nalginate <- as.numeric(sum(otu_pairs_alginate$Alginate == 1))
#   Nfamily <- as.numeric(sum(otu_pairs_alginate$Family == 1))
#   Number_observations <- as.numeric(length(otu_pairs_alginate$Alginate))
#   Npairs = as.numeric(Number_observations*(Number_observations-1)/2)
#   Nalginate_family <- sum()
#   
#   related_expected <- (Nalginate*Nfamily + (Npairs - Nalginate)*(Npairs - Nfamily))/Npairs
#   related_expected
#   
#   kappa <- (related_observed -related_expected)/(Npairs - related_expected)
#   kappa
#   
#   Nalginate
#   Nfamily
#   NAF
#   NAF <- as.numeric(sum(otu_pairs_alginate$Family == 1 & otu_pairs_alginate$Alginate == 1) + sum(otu_pairs_alginate$Family == 1 & otu_pairs_alginate$Alginate == 1))
#   related_observed <- NAF + (Npairs - Nalginate - Nfamily + NAF)
#   related_observed


# Cálculo de los valores
Nalginate <- as.numeric(sum(otu_pairs_alginate$Alginate == 1))  # Pares relacionados en Alginate
Nfamily <- as.numeric(sum(otu_pairs_alginate$Family == 1))      # Pares relacionados en Family
Npairs <- as.numeric(nrow(otu_pairs_alginate))                  # Número total de observaciones (pares) ya que las filas ya son los pares

# Pares relacionados en ambas clasificaciones (Alginate y Family)
Nalginate_family <- as.numeric(sum(otu_pairs_alginate$Family == 1 & otu_pairs_alginate$Alginate == 1))
Nalginate
# Cálculo del acuerdo esperado (related_expected)
related_expected <- (Nalginate * Nfamily + (Npairs - Nalginate) * (Npairs - Nfamily)) / Npairs
related_expected
# Cálculo del acuerdo observado (related_observed)
related_observed <- Nalginate_family + (Npairs - Nalginate - Nfamily + Nalginate_family)
related_observed
# Cálculo del Kappa
kappa <- (related_observed - related_expected) / (Npairs - related_expected)

# Resultados
Nalginate
Nfamily
Nalginate_family
related_expected
related_observed
kappa ###Mismo kappa que con la funcion

