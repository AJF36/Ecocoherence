rm(list = ls())
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(this.path)
work_dir <- this.dir()


Cor <- numeric()
mse_lin_forward <- numeric()
mse_exp_forward <- numeric()
mse_lin_backwards <- numeric()
mse_exp_backwards <- numeric()
substrates <- c("Agarose")
g = 1
minimun_values <- numeric()

for (substrate in substrates){
  print(substrate)
  network_path <- file.path(work_dir, "..", "network_creation_sparcc")
  network_path<- normalizePath(network_path, mustWork = FALSE)
  # setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrate, sep = ""))
  setwd(paste(network_path,substrate,sep = "/"))
  data <- read_tsv(paste("interactions_filtered_0.01._", substrate, "_.tsv", sep = ""))
  
  # Asegúrate de tener el dataframe `data` cargado con las columnas `Cor` y `counts`
  # Filtra el rango de `abs(Cor)` de 0 a 0.45
  data <- data %>% 
    filter(abs(Cor) <= 0.45) %>% 
    filter(abs(Cor) >= 0.4) %>%
    arrange(abs(Cor))
  data$counts <- 1:length(data$Cor) 
  z = 1
  for (j in seq(20,length(data$Cor),5)){
    
    set <- head(data, j)
    
    # Linear model and exponential model
    lm <- lm(counts ~ abs(Cor), set)
    lme <- lm(log(counts) ~ abs(Cor), set)
    
    # Store the statistics for the linear model
    Cor[z] <- data$Cor[j]
    mse_lin_forward[z] <- mean(residuals(lm)^2)
    # lm_pvalue[z] <- overall_p(lm)
    # lm_correlation[z] <- summary(lm)$r.squared
    
    # cor_value[z] <- data$Cor[z]
    # Store the statistics for the exponential model
    pred_exponencial <- exp(predict(lme, newdata = set))
    mse_exponencial_original <- mean((set$counts - pred_exponencial)^2)
    
    # lme_counts[z] <- j
    mse_exp_forward[z] <- mse_exponencial_original
    # lme_correlation[z] <- summary(lme)$r.squared
    # lme_pvalue[z] <- overall_p(lme)
    
    z <- z + 1
  
    
  }
  m = 1
  for (i in seq(length(data$Cor),20,-5)){
    
    set <-tail(data, i)
    
    # Linear model and exponential model
    lm <- lm(counts ~ abs(Cor), set)
    lme <- lm(log(counts) ~ abs(Cor), set)
    
    # Store the statistics for the linear model
    mse_lin_backwards[m] <- mean(residuals(lm)^2)
    # lm_pvalue[z] <- overall_p(lm)
    # lm_correlation[z] <- summary(lm)$r.squared
    
    # cor_value[z] <- data$Cor[z]
    # Store the statistics for the exponential model
    pred_exponencial <- exp(predict(lme, newdata = set))
    mse_exponencial_original <- mean((set$counts - pred_exponencial)^2)
    
    # lme_counts[z] <- j
    mse_exp_backwards[m] <- mse_exponencial_original
    # lme_correlation[z] <- summary(lme)$r.squared
    # lme_pvalue[z] <- overall_p(lme)
    
    m <- m + 1
  }
  
  
  #The three combination vectors  
  lineal_lineal <- mse_lin_forward + mse_lin_backwards
  exponential_exponential <- mse_exp_forward + mse_exp_backwards
  lineal_exponential <- mse_lin_forward + mse_exp_backwards
  results_df <- data.frame(Cor,lineal_lineal,lineal_exponential,exponential_exponential)
  
  # Convertimos el dataframe a formato largo para facilitar la asignación de color
  results_long <- results_df %>%
    pivot_longer(cols = c(lineal_lineal, lineal_exponential, exponential_exponential), 
                 names_to = "Modelo", values_to = "MSE")
  
  # Creamos el gráfico con diferentes colores para cada combinación de modelo
  # plot <- ggplot(results_long, aes(x = abs(Cor), y = MSE, color = Modelo)) +
  #   geom_line() +
  #   labs(title = paste("MSE Total vs abs(Cor) para Diferentes Combinaciones de Modelos en " , substrate, sep = ""),
  #        x = "abs(Cor)", y = "MSE Total") +
  #   theme_bw()
  # 
  plot <- ggplot(results_df, aes(x = abs(Cor))) +
    geom_smooth(aes(y = lineal_lineal, color = "lineal_lineal"), se = FALSE) +
    geom_smooth(aes(y = lineal_exponential, color = "lineal_exponential"), se = FALSE) +
    geom_smooth(aes(y = exponential_exponential, color = "exponential_exponential"), se = FALSE) +
    labs(
      title = paste("MSE Total vs abs(Cor) para Diferentes Combinaciones de Modelos en ",substrate, sep = ""),
      x = "abs(Cor)",
      y = "MSE Total",
      color = "Modelo"
    ) +
    theme_bw()
 
  g = g + 1
  min_ll <- results_df$Cor[which.min(results_df$lineal_lineal)]
  min_le <- results_df$Cor[which.min(results_df$lineal_exponential)]
  min_ee <- results_df$Cor[which.min(results_df$exponential_exponential)]
  mins <- c(min_ll,min_le,min_ee)
  minimun_values[g] <- min(mins)
  # Gráfico de MSE total vs abs(Cor) para cada combinación de modelos
  file_path_figures <- file.path(work_dir,"..","figures")
  file_path_figures <- normalizePath(file_path_figures)
  setwd(file_path_figures)
  system(paste("mkdir ",substrate))
  ggsave("double_linear_model.png",plot,height = 5 ,width =10 )
}
minimun_values
minimun_values<- minimun_values[-1]
names(minimun_values) <- substrates
setwd(work_dir)
df <- data_frame(substrates,minimun_values)
df
setwd(work_dir)
write_tsv(df,"minimun_correlations_threshold.tsv")