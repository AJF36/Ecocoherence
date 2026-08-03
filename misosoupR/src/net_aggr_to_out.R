net_aggr_to_out = function(net.broad){
  # This function works on the data frames containing all 
  # solutions to provide a single aggregated output for 
  # each of them (i.e. means of the fluxes and the proportion
  # of solutions in which each metabolite was observed).
  # author: apascualgarcia.github.io
  # date: April 25th, 2022 (Berlin)
  #
  #browser()
  net.broad[is.na(net.broad)]=0 # convert NA in numeric
  net.names.tmp=net.broad$Row.names # store to use it later
  net.broad=subset(net.broad,select = -c(Row.names)) # and remove
  net.count=net.broad
  net.count[net.count != 0] = 1 # convert observations in counts
  net.broad.mean=as.data.frame(rowMeans(net.broad)) # compute  means
  net.count.mean=as.data.frame(rowMeans(net.count))
  colnames(net.broad.mean)="weight"
  colnames(net.count.mean)="weight"
  net.names=t(sapply(net.names.tmp,FUN=function(x){stri_split(x,fixed="_XX_")[[1]]}))
  net.broad.mean$nodeA=net.names[,1]
  net.broad.mean$nodeB=net.names[,2]
  net.count.mean$nodeA=net.names[,1]
  net.count.mean$nodeB=net.names[,2]
  net.broad.mean=net.broad.mean[,c(2,3,1)] # reorder columns for igraph
  net.count.mean=net.count.mean[,c(2,3,1)] 
  return(list("net.broad.mean"=net.broad.mean,
              "net.count.mean"=net.count.mean))
}

