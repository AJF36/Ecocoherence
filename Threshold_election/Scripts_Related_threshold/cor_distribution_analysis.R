library(readr)
library(ggplot2)
library(phyloseq)
fileSample="/home/ajf/Desktop/CNB/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
sample_metadata = import_qiime_sample_data(fileSample)

head(sample_metadata)


substrates <- unique(sample_metadata$Substrate[2:100])
substrates
plot.list <- list()

for (i in 1:length(substrates)) {
  # Leemos el archivo de datos para cada substrato
  substrate_network <- read_tsv(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", 
                                      substrates[i], 
                                      "/interactions_filtered_0.01._", 
                                      substrates[i], 
                                      "_.tsv", sep = ""))
  
  # Creamos el gráfico de densidad para la correlación
  p <- ggplot(substrate_network, aes(abs(Cor))) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = paste("Correlation distribution", substrates[i], sep = " "))
  
  # Guardamos el gráfico en la lista
  plot.list[[i]] <- p
}


setwd("/home/ajf/Desktop/CNB/ecocoherence_sparcc")
# Crear un archivo PDF para guardar los gráficos
pdf("graficos_distribucion_correlaciones.pdf", width = 8, height = 6)

# Guardar cada gráfico en el PDF
for (i in seq_along(plot.list)) {
  print(plot.list[[i]])
}

# Cerrar el archivo PDF
dev.off()
