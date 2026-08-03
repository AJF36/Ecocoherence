output_networks = function(nets,pathOut=".",
                           plot.graph=FALSE,plot.heatmap=FALSE,print.net=TRUE,
                           vec.media=c(),vec.strains=c(),
                           flux=TRUE,filter=0,
                           par.qgraph = list(),par.heatmap = list()){
  # This function prints to file a list (nets) containing a set of networks 
  # in the format provided by sol_to_exchange_network. It optionally plots them
  # as a network (plot.graph) or as a heatmap (plot.heatmap). It can also
  # select subsets of networks specified by their media (vec.media) and strains (vec.strains).
  # flux = Controls how nets and heatmaps are represented depending on whether the
  #        network represent fluxes (TRUE) or proportions (FALSE).
  # filter = [0,1] filters those metabolites that were observed less than a percentage
  #        (determined by the value of filter) of the total flux (or proportion) of the
  #        metabolite with the most abundant sum of fluxes (or proportions), i.e. filter
  #        a metabolite i if sum flux_i <= filter * max_j(sum flux_j)
  # For more details see: https://igraph.org/r/html/latest/plot.common.html
  library(qgraph)
  library(reshape2)
  #library(tkplot) # not used, potentially useful for the future
  #library(igraph) # not used
  
  # --- Set a directory for the output
  current.dir=getwd()
  setwd(pathOut)
  
  # --- Extract the list of solutions of interest
  if(length(vec.media)== 0){ # no list provided by the user
    media=names(nets) # take the list of media
  }
  if(length(vec.strains) == 0){ # no list provided by the user
    strains=c()
    for(medium in media){ # look for the strains in each media
      strains=c(strains,names(nets[[medium]]))
    }
    strains=unique(strains) # remove repeated
  }
  #.... Depending on the dimension in which we aggregate we change the order of
  #     the lists in the for loops below
  lead.list=media
  slave.list=strains
  
  # --- Start computation
  for(lead in lead.list){ # The lead is the media unless it is aggregated for each strain 
    for(slave in slave.list){ # the slave is the strains list
      if(length(nets[[lead]][[slave]]) == 0){
        next
      }
      sols_list=names(nets[[lead]][[slave]])
      for(sol in sols_list){ # for each solution
        net.tmp=nets[[lead]][[slave]][[sol]]
        name.net=paste(lead,slave,sol,sep="_")
        if(plot.graph == "TRUE"){
          net.tmp2=net.tmp
          # we merge metabolism nodes by removing _in and _out. 
          net.tmp2$nodeA=gsub(x=net.tmp2$nodeA,"_in","")
          net.tmp2$nodeA=gsub(x=net.tmp2$nodeA,"_out","")
          # then we convert in a directed graph changing the order of the columns
          id.pos=which(net.tmp2$weight > 0)
          net.tmp2[id.pos,]=net.tmp2[id.pos,c(2,1,3)] # convert in directed
          # set parameters
          default.args=list(maximum = 10)
          final.args=functionArgsList(default.args,par.qgraph)
          # finally plot
          file.igraph=paste("qgraph_",name.net,".pdf",sep="")
          pdf(file=file.igraph,width=12,height=12)
          # qgr = qgraph(net.tmp)
          qgr = do.call("qgraph",c(list(net.tmp2),final.args))
          dev.off()
          # igr = graph_from_data_frame(net_root) # planned future operations
        } # end plot graph
        if(print.net == TRUE){
          if(flux == TRUE){
            type="flux"
          }else{
            type="count"
          }
          file.df=paste("network-",type,"_",name.net,".csv",sep="")
          write.table(net.tmp,file=file.df,sep="\t", 
                      quote = FALSE,row.names = FALSE)
        } # end print network
        if(plot.heatmap == TRUE){
          net.tmp$weight=abs(net.tmp$weight)
          # transform long list into a matrix
          net.matrix=dcast(net.tmp, nodeA ~ nodeB,value.var= "weight")
          rownames(net.matrix)=net.matrix$nodeA
          net.matrix=subset(net.matrix, select = -c(nodeA))
          net.matrix[is.na(net.matrix)]=0
          net.matrix=as.matrix(net.matrix)
          name.net.heat=name.net
          if(filter > 0){
           row.filter=rowSums(net.matrix)
           id.filter=which(row.filter > filter*max(row.filter))
           net.matrix = net.matrix[id.filter,]
           name.net.heat=paste0(name.net,"_filter",filter)
          }
          if(dim(net.matrix)[2] <= 1 ){
            mes=paste("I could't generate a heatmap for a solution
                      with a single strain. Solution: ",lead,slave,sol,sep= "/")
            warning(mes)
            next}
          # ... create a list of common default arguments for heatmap2
          #browser()
          default.args=list(density.info="none",
                            dendrogram="both",
                          trace="none",
                          keysize = 1,
                          cexRow = 0.65,cexCol = 1.25, cex.lab=3,
                          key.par = list(cex.main=2,cex.axis=1.2,
                                         cex.lab=1.5), #usr=c(0,1,0,1)),
                          margins=c(12,12),
                          ylab="Metabolites",xlab="Strains")                  
          if(flux == TRUE){
            file.heatmap=paste("heatmap-flux_",name.net.heat,".pdf",sep="")
            additional.args=list(scale = "column",
                                 col="bluered", #col=cm.colors(255), #,
                                 key.title = "flux")
          }else{
            file.heatmap=paste("heatmap-count_",name.net.heat,".pdf",sep="")
            additional.args=list(key.xlab = "proportion",key.title="")
          } 
          default.args=append(default.args,additional.args)
          #user.args=list(xlab="Other",margins=c(15,15),key.xlab="test")
          final.args=functionArgsList(default.args,par.heatmap)
          # ... finally plot          
          pdf(file=file.heatmap,width=12,height=12)
          do.call("heatmap.2.mod",c(list(net.matrix),final.args))
          dev.off()
        } # end plot  heatmap
      } # end for solutions
    } # end for slave
  } # end for lead
  setwd(current.dir) # come back to the original folder
}

