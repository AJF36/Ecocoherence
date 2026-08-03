
aggregate_networks = function(sols,met_medium,
                              vec.media=c(),vec.strains=c(),
                              vec.sols=c(),vec.suppliers=c(),
                              mode_media="aggr",mode_strains="aggr",mode_sols="aggr"){
  # This function aggregates solutions computing either the mean across a set of
  # solutions or the proportion of observations happening across a set of solutions.
  # An observation is here understood as the secretion or consumption of a metabolite 
  # by a given species across the set of solutions. Note that secretions and consumptions
  # are considered different observations.
  # Aggregation can occur for each specific combination of C source and species
  # or it can additionally be aggregated either fixing a C source and aggregating
  # all solutions for all species or the other way around, fixing a given species
  # and aggregating C sources. It can also aggregate all C sources and species into
  # a single solution (the default). The function also accepts a vector of species or C sources
  # indicating subsets of interest.
  # sols: a list of lists containing the solutions
  # met_medium: a vector? with the metabolites contained in the medium
  # vec.media: (optional) a char vector with the media to analyse from the list of solutions
  # vec.strains: (optional) a char vector with the strains to analyse from the list  of solutions
  # vec.suppliers: (optional) a char vector with strains that should be present as suppliers to consider a solution
  # vec.sols: (optional) a char vector with the solutions ids for each combination of medium and
  #   strains. These ids must be integers increasingly sorted. Note that if a solution does
  #   not exist for a given combination of medium and strain it will silently continue.
  # lead.media: logical, if true results are presented order by media, if false by strains.
  # mode_media: character, "aggr" to aggregate all solutions by media, "split" otherwise
  # mode_strains: character, "aggr" to aggregate all solutions by media, "split" otherwise
  # mode_sols: character, "aggr" to aggregate all solutions. It is the default if any
  #    mode_media or mode_strains is equal to aggr (i.e. it will only split if there is
  #    no aggregation at all)
  #browser()
  
  lead.media = 1 # this is the default, media on the top of hierarchy
  # --- Determine the selected mode to aggregate data
  if((mode_media == "aggr")&(mode_strains=="split")&(mode_sols == "aggr")){
    mode="SAA"
    lead.media = 0 # strains will lead the main for loop
  }else if((mode_media == "split")&(mode_strains=="aggr")&(mode_sols == "aggr")){
    mode="SAA"
  }else if((mode_media == "split")&(mode_strains=="split")&(mode_sols == "aggr")){
    mode="SSA"
  }else if((mode_media == "split")&(mode_strains=="split")&(mode_sols == "split")){
    mode="SSS"
  }else if((mode_media == "aggr")&(mode_strains=="aggr")&(mode_sols == "aggr")){
    mode="AAA" # default mode
  }else{
    mes=paste("The combination of arguments mode_media =",mode_media,", mode
              strains =",mode_strains,", and mode_sols =",mode_sols,", leads
              to an invalid mode. Please see documentation for details.")
    stop(mes)
  }
  
  # --- Determine if all media and strains or a subset should be extracted
  if(length(vec.media)== 0){ # no list provided by the user
    media=names(sols) # take the list of media
  }else{
    media=vec.media
  }
  if(length(vec.strains) == 0){ # no list provided by the user
    strains=c()
    for(medium in media){ # look for the strains in each media
      strains=c(strains,names(sols[[medium]]))
    }
    strains=unique(strains) # remove repeated
  }else{
    strains=vec.strains
  }
  if(length(vec.sols) == 0){ # we want all solutions
    sols.ctrl=FALSE
  }else{
    sols.ctrl=TRUE
  }
  
  # --- Here determine who is lead and who is slave depending on the modes
  if(lead.media == 1){
    lead.list = media
    slave.list=strains
  }else{
    lead.list=strains
    slave.list=media
  }
  
  # --- Main loop, initialize
  net.list=list()
  net.broad=data.frame()
  net.broad.list=list() # initialize output lists
  net.count.list=list()
  key=0
  #browser()
  for(lead in lead.list){ # The lead is the media unless it is aggregated for each strain 
    for(slave in slave.list){ # the slave is the strains list
      if(lead.media == 1){
        if((length(sols[[lead]][[slave]]) == 0)|
           (is.null(sols[[lead]][[slave]]))){
          mes=paste("There was no solution for ",lead,"and",slave)
          warning(mes)
          next
        }
        sol_tmp=sols[[lead]][[slave]]
      }else{
        if((length(sols[[slave]][[lead]]) == 0) |
           (is.null(sols[[slave]][[lead]]))){
          mes=paste("There was no solution for ",lead,"and",slave)
          warning(mes)
          next
        }
        sol_tmp=sols[[slave]][[lead]]
      }
      j=1
      Nsols=length(sol_tmp)
      for(i in 1:Nsols){ # for each solution
        if(sols.ctrl == TRUE){ # if we only want specific solutions
          if(j > length(vec.sols)){ # check that there are still elements in the vector
            break # otherwise finish the for loop 
          }else if(i != vec.sols[j]){ # check if i is the next in the list
            next # skip otherwise
          }else{ # if it is, we continue and the index j is increased 
            j=j+1
          }
        }
        net.tmp.raw=NULL
        net.tmp.raw=sol_to_exchange_network(sol_tmp[[i]],met_medium) # extract the network
        if(is.null(net.tmp.raw)){
          mes=paste("Solution",i,"for",lead,"and",slave,"has no 
                     exchanged reactions, I skip it.")
          warning(mes)
          next()
        }
        if(length(vec.suppliers) > 0){ # if you want to focus on specific suppliers
          # extract species features, we want to discard solutions if certain suppliers are absent
          out.sp=sol_to_sp_feat(sol_tmp[[i]])
          matched=match(vec.suppliers,out.sp$suppliers) # check if all suppliers are present
          if(any(is.na(matched))){ # if there is a supplier absent
            next() # discard solution
          }
        }
        key=key+1
        net.tmp=rename_net(net.tmp.raw) # rename input and output metabolites
        merged.name=paste(net.tmp$nodeA,net.tmp$nodeB,sep="_XX_")
        df.tmp=as.data.frame(net.tmp$weight) # take the values
        sol.name.tmp=paste("Sol",i,sep="_")
        colnames(df.tmp)=paste(lead,slave,sol.name.tmp,sep=".") # name the solution
        #browser()
        df.tmp$Row.names=merged.name # work with merged names for the merge  function below
        #rownames(df.tmp)=merged.name

        if(mode == "SSS" ){ # if we want every single network
          net.broad.list[[lead]][[slave]][[sol.name.tmp]]=net.tmp # just return it
        }else{ # otherwise we aggregate solutions
          if(dim(net.broad)[1] == 0){
            net.broad=df.tmp
          }else{ # we create a df with a solution in each col, filled with NA if the elements are new
            net.broad=merge(net.broad,df.tmp,all=TRUE)#,by.x=0,by.y=0)
            rownames(net.broad)
          }
        } # end if mode SSS
      }  # end for solutions
      if(mode == "SSA"){ # aggregate only solutions
        if(dim(net.broad)[1] == 0){
          mes=paste("No solution for",lead,"and",slave,"with 
                     exchanged reactions, I skip the pair.")
          warning(mes)
          next()
        }
        lists.out=net_aggr_to_out(net.broad) # extract statistics
        net.broad.list[[lead]][[slave]][["sols_aggr"]]=lists.out[["net.broad.mean"]]
        net.count.list[[lead]][[slave]][["sols_aggr"]]=lists.out[["net.count.mean"]]
        net.broad=data.frame() # initialize
      }
      
    }# end for slave
    if(mode == "SAA"){
      if(dim(net.broad)[1] == 0){
        mes=paste("No solution for",lead,"with 
                     exchanged reactions, I skip it.")
        warning(mes)
        next()
      }
      lists.out=net_aggr_to_out(net.broad) 
      if(lead.media == 1){
        net.broad.list[[lead]][["strains_aggr"]][["sols_aggr"]]=lists.out[["net.broad.mean"]]
        net.count.list[[lead]][["strains_aggr"]][["sols_aggr"]]=lists.out[["net.count.mean"]]
      }else{
        net.broad.list[["media_aggr"]][[lead]][["sols_aggr"]]=lists.out[["net.broad.mean"]]
        net.count.list[["media_aggr"]][[lead]][["sols_aggr"]]=lists.out[["net.count.mean"]]
      }
      net.broad=data.frame() # initialize
    }
  } # end for lead
  if(mode == "AAA"){
    lists.out=net_aggr_to_out(net.broad) # here slave and lead are irrelevant, will be overwritten as "root"
    net.broad.list[["media_aggr"]][["strains_aggr"]][["sols_aggr"]]=lists.out[[1]]
    net.count.list[["media_aggr"]][["strains_aggr"]][["sols_aggr"]]=lists.out[[2]]
  }
  if(mode == "SSS"){
    net.count.list=list() # empty, included just to provide a uniform output
  }
  if(key==0){
    mes="No solution was found for the conditions selected"
    stop(mes)
  }
  return(list("nets.flux"=net.broad.list,"nets.count"=net.count.list))
}