rename_net = function(net.tmp){
  # This function simply takes a directed network of metabolites x strains
  # and it simply renames metabolites to have a sufix (either "in" or "out")
  # depending on whether the weight in the network is positive (excretion)
  # or negative (consumption). The intended use is that when we run across 
  # aset of solutions we can more easily extract statistics differentiating
  # both possibilities
  # author: apascualgarcia.github.io
  # date: April 25th, 2022 (Berlin)
  #
  idx_scr=which(net.tmp$weight > 0) # identify secretions
  idx_con=which(net.tmp$weight < 0) # identify consumption
  net.tmp$nodeA[idx_scr]=paste(net.tmp$nodeA[idx_scr],"_out",sep="") # rename secretion
  net.tmp$nodeA[idx_con]=paste(net.tmp$nodeA[idx_con],"_in",sep="") # rename consumption
  return(net.tmp)
}


