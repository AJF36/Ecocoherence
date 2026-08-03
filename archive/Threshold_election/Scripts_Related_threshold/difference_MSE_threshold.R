library(readr)
library(ggplot2)
library(dplyr)

# Función para obtener el p-valor general del modelo
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p) <- NULL
  return(p)
}

# Datos de ejemplo (puedes reemplazarlos con tus datos)
substrates <-  c("Agarose","Alginate","AgaroseAlginate","AgaroseCarrageenan","AgaroseChitosan","Chitin","Carrageenan")
cutoffs <- c()
total_counts <- c()
cor_threshold <- c()

for (i in 1:length(substrates)){
  
  print(substrates[i])
  
  # Inicializamos los vectores para almacenar resultados
  lm_counts <- c()
  lm_MSE <- c()
  lm_correlation <- c()
  lm_pvalue <- c()
  
  # Vectores para el modelo exponencial
  lme_counts <- c()
  lme_MSE <- c()
  lme_correlation <- c()
  lme_pvalue <- c()
  
  # Directorio y lectura de datos
  setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrates[i], sep= ""))
  data <- read_tsv(paste("interactions_filtered_0.01._", substrates[i], "_.tsv", sep = ""))
  
  # Ordenamos los datos y creamos una variable dummy 'x' para los counts
  data <- arrange(data, desc(abs(Cor)))
  data$x <- 1:length(data$Cor)
  data <- subset(data, select = c("Cor", "x"))
  
  # Índice para el bucle de modelos
  z = 1
  
  # Loop que ajusta los modelos y calcula las métricas
  for (i in seq(20, (length(data$x)*0.4), by = 20)) {
    set <- head(data, i)
    
    # Modelos lineal y exponencial
    lm <- lm(x ~ abs(Cor), set)
    lme <- lm(log(x) ~ abs(Cor), set)
    
    # Almacenamos los valores de MSE
    lm_counts[z] <- i
    lm_MSE[z] <- mean(residuals(lm)^2)
    lme_counts[z] <- i
    
    # MSE en la escala original del modelo exponencial
    pred_exponencial <- exp(predict(lme, newdata = set))
    mse_exponencial_original <- mean((set$x - pred_exponencial)^2)
    lme_MSE[z] <- mse_exponencial_original
    
    z <- z + 1
  }
  
  # Crear data frame con resultados
  results_df <- data.frame(lm_counts, lm_MSE, lme_MSE)
  results_df$difference <- results_df$lm_MSE - results_df$lme_MSE
  
  # Ajustar un modelo lineal para la diferencia de MSE en función de lm_counts
  diff_model <- lm(difference ~ lm_counts, data = results_df)
  
  # Evaluar significancia del coeficiente de la pendiente
  summary_diff <- summary(diff_model)
  p_value_slope <- summary_diff$coefficients[2, 4]  # p-valor del coeficiente de la pendiente
  
  # Imprimir resultados del modelo
  print(summary_diff)
  cat("P-valor de la pendiente:", p_value_slope, "\n")
  
  # Guardar el p-valor y el valor de corte en los vectores
  if (p_value_slope > 0.05) {
    cutoff_value <- results_df$lm_counts[which.max(results_df$difference)]
    cutoffs[i] <- cutoff_value
    print(paste("Punto de corte sugerido en counts:", cutoff_value))
  } else {
    print("No se detectó una tendencia significativa en la diferencia de MSE.")
  }
  
  # Graficar la diferencia de MSE y el ajuste del modelo lineal
  ggplot(results_df, aes(x = lm_counts, y = difference)) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(title = paste("Diferencia de MSE entre modelos -", substrates[i]),
         x = "Counts", y = "Diferencia de MSE (LM - LME)") +
    theme_minimal()
}

