match_sets_multiple = function(sets.list){
  # This function takes a list of sets of solutions following the same structure,
  # i.e. a list with keys sols[[medium]][[strain]][[solution]] and it looks for
  # those solutions in common, meaning that they have the same medium, strain
  # and suppliers, but possibly not the same reactions (i.e. the intersection
  # of all common solutions across all datasets). 
  # It returns 2 lists each containing Nsets (with Nsets the number of
  # sets in the list). Each list contain:
  # * List 1 = The solutions in set_i present in the intersection of all datasets 
  # * List 2 = The solutions in set_i not belonging to the intersection
  #################
  # author = apascualgarcia.github.io
  # date = September 7th, 2022 (Lausanne)
  #################
  #browser()
  sets.names = names(sets.list)
  Nsets=length(sets.list)
  
  
  # The first set will be the reference, we will get its solutions after
  # the first round of intersections, and also those of the last one.
  set_lead=sets.names[1]
  setIntLead=sets.list[[set_lead]] 
  set_close=sets.names[length(sets.names)] # the last one will
  
  # --- 
  sets_int=list()
  sets_not_int=list()
  j=2
  for(set_nameB in sets.names[j:Nsets]){
    setB=sets.list[[set_nameB]]
    # match both sets
    setsAB.list=match_sets_pairwise(setIntLead,setB)
    setIntLead=setsAB.list$solsAB_inA
    sets_int[[set_nameB]]=setsAB.list$solsAB_inB # store for 
    if(set_nameB == set_close){ # We can already store those not in the intersection
      sets_not_int[[set_close]]=setsAB.list$solsB_not_inA # for the last set
    }
  }
  sets_int[[set_lead]]=setIntLead
  
  # --- Now, we still need to remove solutions for those sets that are
  #     not the first nor the last one, and we will also retrieve for
  #     each set those solutions that do not belong to the intersection
  #     except for the last one, for which we have all the information already
  
  for(set_nameB in sets.names[1:(Nsets-1)]){ 
    setB=sets.list[[set_nameB]] # we need to compare with the original
    setsAB.list=match_sets_pairwise(setIntLead,setB)
    sets_int[[set_nameB]]=setsAB.list$solsAB_inB # store for 
    sets_not_int[[set_nameB]]=setsAB.list$solsB_not_inA
  }
  return(list("sets_int" = sets_int, "sets_not_int" = sets_not_int))
}