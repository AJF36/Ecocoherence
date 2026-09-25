match_sets_pairwise = function(solsA, solsB, suppress = TRUE){
  # This function takes two sets of solutions following the same structure,
  # i.e. a list with keys sols[[medium]][[strain]][[solution]] and it looks for
  # those solutions in common, meaning that they have the same medium, strain
  # and suppliers, but possibly not the same reactions. 
  # input = 
  # solsA, solsB: lists with solutions
  # suppress: Logical. Suppress warnings when solutions are not found, not
  #           recommended for large datasets
  # output = It returns four lists:
  # * solutions in setA present in both setA and setB
  # * solutions in setB present in both setA and setB
  # * solutions in setA not present in setB
  # * solutions in setB not present in setA
  #################
  # author = apascualgarcia.github.io
  # date = September 7th, 2022 (Lausanne)
  #################
  #browser()
  # --- Initialize lists containing solutions present in only one set
  solsA_not_inB=solsA # we will remove those found in common
  solsB_not_inA=solsB 
  
  # --- Initialize those in common
  solsAB_inA=list()
  solsAB_inB=list()
  
  # --- Extract media and strains in the first set of solutions, will be the reference
  #     for for loops below
  mediaA=names(solsA)
  strainsA=c()
  for(medium in mediaA){ # look for the strains in each media
    strainsA=c(strainsA,names(solsA[[medium]]))
  }
  strainsA=unique(strainsA)
  
  NsetA=0; NsetB=0; NsetA_zero=0; NsetB_zero=0;
  Ncommon_inA=0; Ncommon_inB=0; NsetA_notB=0; NsetB_notA=0
  for(medium in mediaA){ # For each medium in set A
    for(strain in strainsA){ # For each focal strain
      if((length(solsA[[medium]][[strain]]) == 0)| # double check that  the solution is non-empty
         (is.null(solsA[[medium]][[strain]]))){
        if(suppress == FALSE){
          mes=paste("There was no solution for ",medium,"and",strain,"in setA")
          warning(mes)
        }
        next
      }
      sol_tmpA=solsA[[medium]][[strain]]
      if((length(solsB[[medium]][[strain]]) == 0)|
         (is.null(solsB[[medium]][[strain]]))){
        if(suppress == FALSE){
          mes=paste("There was no solution for ",medium,"and",strain,"in setB")
          warning(mes)
        }
        next
      }
      sol_tmpB=solsB[[medium]][[strain]]
      NsolsA=length(sol_tmpA)
      NsolsB=length(sol_tmpB)
      NsetA=NsetA+NsolsA
      NsetB=NsetB+NsolsB
      u=0
      vecSols0A=vector(mode = "numeric", length = NsolsA) # set = 1 if A solutions have no growth
      vecSols0B=vector(mode = "numeric", length = NsolsB) # these vectors are used for counting
      removeA=vector(mode = "numeric", length = NsolsA) # set = 1 if A solutions have no growth or if there is a match
      removeB=vector(mode = "numeric", length = NsolsB) # while these ones are used to remove solutions
      for(i in 1:NsolsA){ # for each solution in A
        key_remove=0
        comm_growth=unlist(sol_tmpA[[i]]["community_growth"])
        if(is.null(comm_growth)){
          vecSols0A[i] = 1
          removeA[i] = 1
          key_remove=1
        }else{
          out.spA=sol_to_sp_feat(sol_tmpA[[i]],strain) # extract suppliers
        }
        for(j in 1:NsolsB){ # for each solution in B
          comm_growth=unlist(sol_tmpB[[j]]["community_growth"])
          if(is.null(comm_growth)){
            vecSols0B[j] = 1
            removeB[j] = 1
            next
          }else if(key_remove == 1){ # if it was zero in A there is nothing else to do
            next
          }
          out.spB=sol_to_sp_feat(sol_tmpB[[j]],strain) # extract suppliers
          matched=match(out.spA$suppliers,out.spB$suppliers) # check if all suppliers in A are present in B
          if(any(is.na(matched))){ # if there is a supplier absent
            next() # discard solution
          }
          matched=match(out.spB$suppliers,out.spA$suppliers) # check that there are no more suppliers in B than in A
          if(any(is.na(matched))){ # if there is a supplier absent
            next() # discard solution
          }
          # reaching this point means that we have a match
          u=u+1 # increase counter
          solsAB_inA[[medium]][[strain]][[u]]=sol_tmpA[[i]]
          solsAB_inB[[medium]][[strain]][[u]]=sol_tmpB[[j]]
          removeA[i] = 1
          removeB[j] = 1
          break # once found go to next solution in A
        } # end solutions in B
      } # end solutions in A
      #browser()
      #removeAtest=removeA[1:3]
      #solsA_not_inB[[medium]][[strain]][removeAtest]=NULL
      id.zeroA = which(vecSols0A == 1)
      id.zeroB = which(vecSols0B == 1)
      id.removeA = which(removeA == 1)
      id.removeB = which(removeB == 1)
      solsA_not_inB[[medium]][[strain]][id.removeA]=NULL # to get those only in A, we remove those in common
      solsB_not_inA[[medium]][[strain]][id.removeB]=NULL # to get those only in B, we remove those in common
      NsetA_zero = NsetA_zero + length(id.zeroA)
      NsetB_zero = NsetB_zero + length(id.zeroB)
      Ncommon_inA = Ncommon_inA + (length(id.removeA) - length(id.zeroA))
      Ncommon_inB = Ncommon_inB + (length(id.removeB) - length(id.zeroB))
    } # end strain
    # remove any empty list within medium (all solutions matched)
    solsA_not_inB[[medium]]=Filter(Negate(function(x){length(x)==0}),
                                   solsA_not_inB[[medium]])
    solsB_not_inA[[medium]]=Filter(Negate(function(x){length(x)==0}),
                                   solsB_not_inA[[medium]])
  } # end medium
  # ... remove any remainder empty list
  solsA_not_inB=Filter(Negate(function(x){length(x)==0}),
                                 solsA_not_inB)
  solsB_not_inA=Filter(Negate(function(x){length(x)==0}),
                                 solsB_not_inA)
  if(length(solsAB_inA)==0){
    mes=paste("There was no solution in common between both sets")
    warning(mes)
  }
  # ... Compute final numbers:
  NsetA_notB = NsetA - Ncommon_inA - NsetA_zero # those only in A, the total minus common and zero
  NsetB_notA = NsetB - Ncommon_inB - NsetB_zero
  #Ncommon_inA=Ncommon_inA-NsetA_zero # those common in A, remove those zero
  #Ncommon_inB=Ncommon_inB-NsetB_zero
  N_elements=c(NsetA,NsetB,
               Ncommon_inA,Ncommon_inB,
               NsetA_zero,NsetB_zero,
               NsetA_notB,NsetB_notA)
  names(N_elements)=c("NsetA","NsetB",
                      "Ncommon_inA","Ncommon_inB",
                      "NsetA_zero","NsetB_zero",
                      "NsetA_notB","NsetB_notA")
  return(list("solsAB_inA" = solsAB_inA,"solsAB_inB" = solsAB_inB,
              "solsA_not_inB" = solsA_not_inB, "solsB_not_inA" = solsB_not_inA,
              "N_elements"=N_elements))
}