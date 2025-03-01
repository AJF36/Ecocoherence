library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)

## Function to get the overall p-value
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

# Initialize the substrates
substrates <- c("Chitosan","Agarose","Alginate","AgaroseAlginate","AgaroseCarrageenan","AgaroseChitosan","Chitin","Carrageenan")
length(substrates)
length(substrates)
cutoffs <- c()
total_counts <- c()
cor_threshold <- c()
plots <- c()
# Main loop for each substrate
for (i in 1:length(substrates)) {
  
  print(substrates[i])
  
  # Initialize vectors for statistics
  lm_counts <- c()
  lm_MSE <- c()
  lm_correlation <- c()
  lm_pvalue <- c()
  
  lme_counts <- c()
  lme_MSE <- c()
  lme_correlation <- c()
  lme_pvalue <- c()
  
  # Set directory and read data
  setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/", substrates[i], sep = ""))
  data <- read_tsv(paste("interactions_filtered_0.01._", substrates[i], "_.tsv", sep = ""))
  
  # Order data by absolute value of Cor
  data <- arrange(data, desc(abs(Cor)))
  data$x <- 1:length(data$Cor)
  data <- subset(data, select = c("Cor", "x"))
  
  z = 1
  
  # Inner loop for different counts
  for (j in seq(20, length(data$x)* 0.2, by = 5)) {
    
    set <- head(data, j)
    
    # Linear model and exponential model
    lm <- lm(x ~ abs(Cor), set)
    lme <- lm(log(x) ~ abs(Cor), set)
    
    # Store the statistics for the linear model
    lm_counts[z] <- j
    lm_MSE[z] <- mean(residuals(lm)^2)
    lm_pvalue[z] <- overall_p(lm)
    lm_correlation[z] <- summary(lm)$r.squared
    
    # Store the statistics for the exponential model
    pred_exponencial <- exp(predict(lme, newdata = set))  # Transform predictions to original scale
    mse_exponencial_original <- mean((set$x - pred_exponencial)^2)  # MSE in original scale
    
    lme_counts[z] <- j
    lme_MSE[z] <- mse_exponencial_original
    lme_correlation[z] <- summary(lme)$r.squared
    lme_pvalue[z] <- overall_p(lme)
    
    z <- z + 1
  }
  
  # Create data frames for linear and exponential models
  statisticslm_df <- data.frame(lm_counts, lm_MSE, lm_correlation, lm_pvalue)
  statisticslme_df <- data.frame(lme_counts, lme_MSE, lme_correlation, lme_pvalue)
  
  # Join data frames and calculate the difference
  join_df <- cbind(statisticslm_df, statisticslme_df)
  join_df$diference <- join_df$lm_MSE - join_df$lme_MSE
  
  # Calculate the threshold and cutoff points
  threshold <- max(join_df$diference) * 0.01
  cutoff_index <- which(join_df$diference >= threshold)[1]
  cutoff_value <- join_df$lm_MSE[cutoff_index]  # Value of MSE in the linear model at cutoff point
  x_cutoff <- lm_counts[cutoff_index]
  
  # Store cutoff values
  cutoffs[i] <- x_cutoff
  cor_threshold[i] <- abs(data$Cor[x_cutoff])
  total_counts[i] <- max(join_df$lm_counts)
  print(total_counts[i])
}

# Output the vectors to check results
cor_threshold
cutoffs
total_counts
