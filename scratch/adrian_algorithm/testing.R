library(readr)
library(tidyr)
library(dplyr)
library(yaml)
pak::pak("yaml")
rm(list = ls())
getwd()
setwd("adrian_algorithm/COmmunites")
file <- read.delim("/home/ajf/Desktop/CNB/ecocoherence/adrian_algorithm/COmmunites/cluster_communities_fixed/cluster_1")

file_yaml <- read_yaml(file = "/home/ajf/Desktop/CNB/ecocoherence/adrian_algorithm/COmmunites/cluster_communities_fixed/cluster_1.yaml")
class(file_yaml)
file_yaml 

file_yaml[[1]][[1]][[1]][1]

names(file_yaml)
##### Test Claude
library(yaml)
file_yaml[[1]]$min

### Second try 
convert_yaml_format_2 <- function(input_file, output_file) {
  # Read the original YAML file
  data <- yaml.load_file(input_file)
  
  # Initialize the new structure
  converted_data <- list()
  
  # Process each metabolite (top-level key)
  for (metabolite in names(data)) {
    converted_data[[metabolite]] <- list(min = list())
    
    # Process each community in the min list (there might be multiple)
    for (i in seq_along(data[[metabolite]]$min)) {
      original_min <- data[[metabolite]]$min[[i]]
      
      # Separate Growth variables (community) from other variables (solution)
      community <- list()
      solution <- list()
      
      for (var_name in names(original_min)) {
        if (grepl("^Growth_", var_name)) {
          # Extract model name and set to 1 for community
          # Remove "Growth_" prefix to get model name
          model_name <- sub("^Growth_", "y_", var_name)
          community[[model_name]] <- 1
          
          # ALSO add the Growth variable to solution with original name and value
          solution[[var_name]] <- original_min[[var_name]]
        } else {
          # All other variables go to solution
          solution[[var_name]] <- original_min[[var_name]]
        }
      }
      
      # Create the new structure for this community
      converted_data[[metabolite]]$min[[i]] <- list(
        community = community,
        solution = solution
      )
    }
  }
  
  # Write the converted data to output file
  write_yaml(converted_data, output_file)
  
  cat("Conversion completed successfully!\n")
  cat("Input file:", input_file, "\n")
  cat("Output file:", output_file, "\n")
}


convert_yaml_format_2("/home/ajf/Desktop/CNB/ecocoherence/adrian_algorithm/COmmunites/cluster_communities_fixed/cluster_1.yaml","test_conversion_2.yaml")


###Now lets test if it works
test_cluster <- read_yaml("COmmunites/test_conversion_2.yaml")
adrian_data <- read_yaml("COmmunites/example_data.yaml")
library(ramen)


test_ramen <- ramen::importMisosoup(test_cluster)
test_ramen_f <- test_ramen[[1]] |> dplyr::filter(cons_id == "gal_min_1")


test_adrian_data <- ramen::importMisosoup(adrian_data)
adrian_data_f <- test_adrian_data[[1]] |>
  dplyr::filter(cons_id == "3pg_min_1")


data_for_test_cm <- test_ramen_f
cm_test <- ramen::ConsortiumMetabolism(data_for_test_cm, name = "Test")
cm_adrian <- ramen::ConsortiumMetabolism(adrian_data_f,name = "adrian_data",split_by = "cons_id")
plot(cm_test, type = "EffectiveProduction")
plot(cm_adrian, type = "EffectiveProduction")

getEdges(cm_test)
plot(test_cm)
class(test_ramen)
