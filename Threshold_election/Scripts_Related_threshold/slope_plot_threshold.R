library(readr)
library(ggplot2)
library(gridExtra)
##Vector with the substrates
substrates <- c("Chitosan","Agarose","Alginate","AgaroseAlginate","AgaroseCarrageenan","AgaroseChitosan","Chitin","Carrageenan")

##Function that calculates the p-value
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

##Loop over the substrates
for (substrate in substrates){
  setwd(paste("/home/ajf/Desktop/CNB/ecocoherence_sparcc/",substrate,sep = ""))
  data <- read_csv(paste("Statistics_",substrate,".csv",sep = ""))
  
  
  lm_counts <- c()
  lm_MSE <- c()
  lm_pvalue <- c()
  lm_correlation <- c()
  lm_slope <- c()
  
  
  lm_plot <- ggplot(data,aes(lm_counts,diference)) +
    geom_line() +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    theme_bw() +
    ggtitle(substrate)
  
  z = 1
  
  for (j in seq(5,length(data$lm_counts), by = 1)) {
    
    set <- head(data, j)
    # Linear model 
    lm <- lm(diference ~ lm_counts, set)
    summary(lm)
    # Store the statistics for the linear model
    lm_counts[z] <- j
    lm_MSE[z] <- mean(residuals(lm)^2)
    lm_pvalue[z] <- overall_p(lm)
    lm_correlation[z] <- summary(lm)$r.squared
    lm_slope[z] <- summary(lm)$coefficients["lm_counts", "Estimate"]
    
    z <- z + 1
    
    
  }
  
  statistics_data_frame <- data.frame(lm_counts,lm_correlation,lm_MSE,lm_pvalue)
  write_tsv(statistics_data_frame,paste("linear_model_difference",substrate,".tsv"))
  
  correlation_plot <- ggplot(statistics_data_frame,aes(lm_counts * 20,lm_correlation)) +
    geom_line() +
    theme_bw()
  
  pvalue_plot <- ggplot(statistics_data_frame,aes(lm_counts * 20,lm_pvalue)) +
    geom_line() +
    theme_bw()
  
  MSE_plot <- ggplot(statistics_data_frame,aes(lm_counts * 20,lm_MSE)) +
    geom_line() +
    theme_bw()
  
  slope_plot <- ggplot(statistics_data_frame,aes(lm_counts * 20,lm_slope)) +
    geom_line() +
    theme_bw()
  #ggsave(filename = paste("Model_Comparison_", substrate, ".png", sep = ""), plot = grid_plot, width = 12, height = 15)
  
  
  plot_grid <- grid.arrange(lm_plot,correlation_plot,pvalue_plot,MSE_plot,slope_plot, ncol = 2)
  ggsave(filename = paste("difference_linear_model_",substrate,".png",sep = ""),plot_grid, width = 25, height = 15)
}






