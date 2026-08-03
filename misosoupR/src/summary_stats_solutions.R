
summary_stats_solutions = function(sols,met_medium,
                                   vec.media=c(),
                                   vec.strains=c(),
                                   vec.sols=c(),
                                   vec.suppliers=c(),
                             lead.media=TRUE){
  # This function extracts from each solution found in the list, a number
  # of summary statistics (number of species, metabolites, fluxes, etc.) and
  # it returns a data frame containing the different statistics (columns) for
  # each solution (rows)
  # sols: a list of lists containing the solutions
  # met_medium: a vector? with the metabolites contained in the medium
  # vec.media: (optional) a char vector with the media to analyse from the list of solutions
  # vec.strains: (optional) a char vector with the strains to analyse from the list  of solutions
  # vec.suppliers: (optional) a char vector with strains that should be present as suppliers to consider a solution
  # vec.sols: (optional) a char vector with the solutions ids for each combination of medium and
  #   strains. These ids must be integers increasingly sorted. Note that if a solution does
  #   not exist for a given combination of medium and strain it will silently continue.
  # lead.media: logical, if true results are presented order by media, if false by strains.
  ##########
  # apascualgarcia.github.io
  # May 2022, ETH-Zürich
  ###########
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
  mediaToFocal=matrix(0,nrow=length(media),ncol=length(strains)) # this is the only matrix we can initialized in advance
  rownames(mediaToFocal)=media
  colnames(mediaToFocal)=strains
  if(length(vec.sols) == 0){ # we want all solutions
    sols.ctrl=FALSE
  }else{
    sols.ctrl=TRUE
  }
  # --- Here determine who is lead and who is slave depending on the modes
  if(lead.media == TRUE){
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
  mat.list=list()
  #browser()
  # lead=lead.list[1]
  # slave=slave.list[1]
  # i=1
  key=0
  key.suppl=0
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
        focal=slave
      }else{
        if((length(sols[[slave]][[lead]]) == 0) |
           (is.null(sols[[slave]][[lead]]))){
          mes=paste("There was no solution for ",lead,"and",slave)
          warning(mes)
          next
        }
        sol_tmp=sols[[slave]][[lead]]
        focal=lead
      }
      Nsols=length(sol_tmp)
      j=1
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
        # compute statistics, first species features
        out.sp=sol_to_sp_feat(sol_tmp[[i]],focal)
        if(length(vec.suppliers) > 0){ # if you want to focus on specific suppliers
          # extract species features, we want to discard solutions if certain suppliers are absent
          matched=match(vec.suppliers,out.sp$suppliers) # check if all suppliers are present
          if(any(is.na(matched))){ # if there is a supplier absent
            next() # discard solution
          }
        }
        key=key+1
        # then call to function to investigate reactions
        out.stat=sol_to_statistics(net.tmp.raw)
        if(lead.media == TRUE){
          medium=lead
          out=c("medium"=lead,"strain"=slave,"sol"=i,
                "N_eff_sp"=out.sp$N_eff_sp,out.sp$comm_growth,
                out.stat)
        }else{
          medium=slave
          out=c("strain"=lead,"medium"=slave,"sol"=i,
                "N_eff_sp"=out.sp$N_eff_sp,out.sp$comm_growth,
                out.stat)
        }
        out.df.tmp=as.data.frame(t(out))
        if(key == 1){
          out.df.all=(out.df.tmp)
        }else{
          out.df.all=rbind(out.df.all,(out.df.tmp))
        }
        # Update all matrices
        if(out.sp$N_sp == 1){
          mediaToFocal[medium,focal]= -1
        }else{
          key.suppl=key.suppl+1
          mediaToFocal[medium,focal]= mediaToFocal[medium,slave]+1
          mat.list=update_matrices(key.suppl,out.sp,medium,mat.list)
        }
      }  # end for solutions
    }# end for slave
  } # end for lead
  if(key==0){
    mes="No solution was found for the conditions selected"
    stop(mes)
  }
  mat.list$mediaToFocal=mediaToFocal
  diag(mat.list$supplToSuppl)=0
  return(list("stats.df"=out.df.all,"mat.list"=mat.list))
}