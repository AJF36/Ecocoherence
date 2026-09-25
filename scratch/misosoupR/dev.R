
for (i in seq(length(sols))) {
  c_source <- names(sols)[i]
  for (j in seq(length(sols[[i]]))) {
    strain <- names(sols[[i]])[j]
    for (k in seq(length(sols[[i]][[j]]))) {
      x <- sol_to_exchange_network(solution,met_medium)
      name <- paste(c_source, strain, k, sep = "_")
      write.csv(x, file = paste0("test_out/", name, ".csv"), row.names = FALSE)
    }
  }
}
 