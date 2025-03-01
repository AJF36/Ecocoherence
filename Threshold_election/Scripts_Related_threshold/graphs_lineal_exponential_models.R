library(readr)
library(ggplot2)
library(gridExtra)
library(dplyr)
# Lista de sustratos
substrates <- c("Chitosan","Agarose","Alginate","AgaroseAlginate","AgaroseCarrageenan","AgaroseChitosan","Chitin","Carrageenan")

# Loop para cada sustrato
for (substrate in substrates) {
  setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))
  
  # Cargar el archivo CSV para el sustrato actual
  data_path <- paste("Statistics_", substrate,".csv", sep = "")
  data <- read_csv(data_path)
  
  # Cargar el archivo de datos original para generar el ajuste de los modelos
  original_data_path <- paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, "/interactions_filtered_0.01._", substrate, "_.tsv", sep = "")
  original_data <- read_tsv(original_data_path)
  original_data <- original_data %>% 
    arrange(desc(abs(Cor))) %>%
    mutate(x = 1:n())
  
  # Gráfico de Correlación (R-squared) entre modelos
  plot_correlation <- ggplot(data) +
    geom_line(aes(x = lm_counts, y = lm_correlation, color = "Linear Model")) +
    geom_line(aes(x = lme_counts, y = lme_correlation, color = "Exponential Model")) +
    labs(color = "Model Type") +
    theme_minimal() +
    xlab("Counts") +
    ylab("R-squared") +
    ggtitle(paste("Correlation -", substrate))
  
  # Gráfico de MSE entre modelos
  plot_MSE <- ggplot(data) +
    geom_line(aes(x = lm_counts, y = lm_MSE, color = "Linear Model")) +
    geom_line(aes(x = lme_counts, y = lme_MSE, color = "Exponential Model")) +
    labs(color = "Model Type") +
    theme_minimal() +
    xlab("Counts") +
    ylab("MSE") +
    ggtitle(paste("MSE -", substrate))
  
  # Gráfico de p-valor entre modelos
  plot_pvalue <- ggplot(data) +
    geom_line(aes(x = lm_counts, y = lm_pvalue, color = "Linear Model")) +
    geom_line(aes(x = lme_counts, y = lme_pvalue, color = "Exponential Model")) +
    labs(color = "Model Type") +
    theme_minimal() +
    xlab("Counts") +
    ylab("p-value") +
    ggtitle(paste("P-value -", substrate))
  
  # Gráfico de Diferencia de MSE entre modelos
  plot_difference <- ggplot(data) +
    geom_line(aes(x = lm_counts, y = diference), color = "purple") +
    theme_minimal() +
    xlab("Counts") +
    ylab("Difference (LM MSE - LME MSE)") +
    ggtitle(paste("MSE Difference -", substrate))
  
  # Gráfico de Ajuste del Modelo Lineal
  linear_model <- lm(x ~ abs(Cor), data = original_data)
  plot_linear_fit <- ggplot(original_data, aes(x = abs(Cor), y = x)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggtitle(paste("Linear Model Fit -", substrate)) +
    xlab("Absolute Correlation") +
    ylab("Counts") +
    theme_minimal()
  
  # Gráfico de Ajuste del Modelo Exponencial
  exponential_model <- lm(log(x) ~ abs(Cor), data = original_data)
  plot_exponential_fit <- ggplot(original_data, aes(x = abs(Cor), y = x)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = exp(predict(exponential_model, newdata = original_data))), color = "red") +
    ggtitle(paste("Exponential Model Fit -", substrate)) +
    xlab("Absolute Correlation") +
    ylab("Counts") +
    theme_minimal()
  
  # Crear una cuadrícula de los 6 gráficos
  grid_plot <- grid.arrange(
    plot_correlation, plot_MSE, plot_pvalue, plot_difference, plot_linear_fit, plot_exponential_fit, 
    ncol = 2
  )
  
  # Guardar la cuadrícula como un archivo PNG
  ggsave(filename = paste("Model_Comparison_", substrate, ".png", sep = ""), plot = grid_plot, width = 12, height = 15)
}

