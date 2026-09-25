library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)


##function to get the overall p-value

overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


substrates <- c("Chitosan")
length(substrates)
cutoffs <- c()
total_counts <- c()
cor_threshold <- c()

for (i in 1:length(substrates)){
  
  print(substrates[i])
  lm_counts <- c()
  lm_MSE <- c()
  lm_correlation <- c()
  lm_pvalue <- c()
  
  # Same for the exponential model
  lme_counts <- c()
  lme_MSE <- c()
  lme_correlation <- c()
  lme_pvalue <- c()
  
  #Set the directory and read the data
  setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrates[i],sep= ""))
  data <- read_tsv(paste("interactions_filtered_0.01._",substrates[i],"_.tsv",sep = ""))
  
  #Order the data by the absolute value of Cor and making a dummy variable(x) which is the number of counts
  data <- arrange(data, desc(abs(Cor)))
  data$x <- 1:length(data$Cor)
  data <- subset(data, select = c("Cor","x"))
  
  
  
  
  
  
  
  # data_test <- data[1000:5000,]
  # data_test
  # 
  # linear_model <- lm(abs(Cor) ~ x,data_test)
  # summary(linear_model)$r.squared
  # mean(residuals(linear_model)^2)
  # 
  # exponential_model <- lm(log(abs(Cor)) ~ x,data_test)
  # summary(exponential_model)$r.squared
  # mean(residuals(exponential_model)^2)
  # 
  # 
  # ggplot(data,aes(((-1*abs(Cor)) ),x)) +
  #   geom_line() +
  #   labs(title = "Agarose") +
  #   theme_bw()
  
  #This variable is an index for the vectors in the next loop
  
  z = 1
  
  #This is the loop which make all the models and stores the statistics in the vectors
  
  for (i in seq(20 ,length(data$x), by= 5)){
    
  
    set <- head(data,i)
    
    #Linear model and the exponential model
    lm <- lm( x ~ abs(Cor) ,set)
    lme <- lm(log(x) ~ abs(Cor), set)
    
    #Store the 3 statistics and the counts for the linear model
    lm_counts[z] <- i
    lm_MSE[z] <- mean(residuals(lm)^2)
    lm_pvalue[z] <- overall_p(lm)
    lm_correlation[z] <- summary(lm)$r.squared
    #Same with the exponential
    lme_counts[z] <- i
    
    # Calcular el MSE en la escala original del modelo exponencial
    pred_exponencial <- exp(predict(lme, newdata = set))  # Transformar las predicciones a la escala original
    
    # MSE en la escala original
    mse_exponencial_original <- mean((set$x - pred_exponencial)^2)
    
    lme_MSE[z] <- mse_exponencial_original
    lme_correlation[z] <-summary(lme)$r.squared
    lme_pvalue[z]<- overall_p(lme)
  
    z <- z +1
  }
  
  #I make the dataframes of the models
  statisticslm_df <- data.frame(lm_counts,lm_MSE,lm_correlation,lm_pvalue)
  statisticslme_df <- data.frame(lme_counts,lme_MSE,lme_correlation,lme_pvalue)
  
  
  #joined_dataframe
  
  join_df <- cbind(statisticslm_df,statisticslme_df)
  join_df$diference <-  join_df$lm_MSE - join_df$lme_MSE
  
  
  # plot(join_df$lm_counts,join_df$diference)
  # 
  # ggplot(data,aes(abs(Cor), log(x))) +
  #   geom_line()
  
  # CORRELATION_PLOT <- ggplot(join_df) +
  #   geom_line(aes(x = lme_counts, y = lme_correlation, color = "LME")) + 
  #   geom_line(aes(x = lm_counts, y = lm_correlation, color = "LM")) +
  #   labs(color = "Model Type") +
  #   theme_minimal() +
  #   xlab("Counts") +
  #   ylab("Correlation") 
  # 
  # 
  # MSE_PLOT <- ggplot(join_df) +
  #   geom_line(aes(x = lme_counts, y = lme_MSE, color = "LME")) + 
  #   geom_line(aes(x = lm_counts, y = lm_MSE, color = "LM")) +
  #   labs(color = "Model Type") +
  #   theme_minimal() +
  #   xlab("Counts") +
  #   ylab("MSE") +
  #   labs(title = "Agarose") +
  #   ylim(0,2.5*10^7)
  #  
  
  
  
  
  threshold <- max(join_df$diference) * 0.01
  cutoff_index <- which(join_df$diference >= threshold)[1]
  cutoff_value <- join_df$lm_MSE[cutoff_index]  # Valor de MSE en el modelo lineal en el punto de corte
  x_cutoff <- lm_counts[cutoff_index]
  
  cutoffs[i] <- x_cutoff
  
  cor_threshold[i] <- abs(data$Cor[x_cutoff])
  
  total_counts[i] <- max(join_df$lm_counts)
  print(total_counts[i])
}

#Make 3 vectors to store the 4 vector to store the number of counts which are taken and the statistics for the lineal model

cor_threshold

total_counts

# difference_plot <- ggplot(join_df,aes(lm_counts,diference)) +
#   geom_line() +
#   abline(h = x_cutoff,col = "purple",lty = 2)
# 
# difference_plot
# sum(abs(data$Cor) > 0.6)
# grid.arrange(MSE_PLOT,CORRELATION_PLOT,ncol = 2)
# which.min(join_df$lm_correlation)
# sum(abs(data$Cor) > 0.6)

