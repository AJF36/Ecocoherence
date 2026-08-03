extract_solutions = function(sols,
                             vec.media=c(),vec.strains=c(),
                             vec.sols=c(),vec.suppliers=c()){
  # This function extracts a subset of solutions identified by the media,
  # strains, or numeric id, given by the different vectors. Although the three
  # vectors are indicated as optional, it is expected that at least one
  # vector is not  empty, otherwise no subsetting will be performed and
  # it will return the same list that is provided.
  # sols: a list of lists containing the solutions
  # vec.media: (optional) a char vector with the media to analyse from the list of solutions
  # vec.strains: (optional) a char vector with the strains to analyse from the list  of solutions
  # vec.suppliers: (optional) a char vector with strains that should be present as suppliers to consider a solution
  # vec.sols: (optional) a char vector with the solutions ids for each combination of medium and
  #   strains. These ids must be integers increasingly sorted. Note that if a solution does
  #   not exist for a given combination of medium and strain it will silently continue.
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
  if(length(vec.sols) == 0){ # we want all solutions
    sols.ctrl=FALSE
  }else{
    sols.ctrl=TRUE
  }
  # --- Here determine who is lead and who is slave depending on the modes
  # ..... This structure is unnecessary in this script but I keep it for potential
  #      future developments. It is used in other scripts to present results in a
  #      desired order (either by media or by strain)  
  lead.media=TRUE
  if(lead.media == TRUE){
    lead.list = media
    slave.list=strains
  }else{
    lead.list=strains
    slave.list=media
  }
  
  # --- Main loop, initialize
  sols.sub=list()
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
        focal=slave
        sol_tmp=sols[[lead]][[slave]]
      }else{
        if((length(sols[[slave]][[lead]]) == 0) |
           (is.null(sols[[slave]][[lead]]))){
          mes=paste("There was no solution for ",lead,"and",slave)
          warning(mes)
          next
        }
        focal=lead
        sol_tmp=sols[[slave]][[lead]]
      }
      Nsols=length(sol_tmp)
      j=1 # indexes the vector of desired solutions
      k=0 # indexes final selected solutions
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
        # Finally, check if the desired suppliers are present in the solution
        out.sp=sol_to_sp_feat(sol_tmp[[i]],focal)
        if(length(vec.suppliers) > 0){ # if you want to focus on specific suppliers
          # extract species features, we want to discard solutions if certain suppliers are absent
          matched=match(vec.suppliers,out.sp$suppliers) # check if all suppliers are present
          if(any(is.na(matched))){ # if there is a supplier absent
            next() # discard solution
          }
        }
        k=k+1
        sols.sub[[slave]][[lead]][[k]]=sol_tmp[[i]]
      }
    }
  }
  return(sols.sub)
}