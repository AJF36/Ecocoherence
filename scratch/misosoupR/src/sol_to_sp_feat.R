sol_to_sp_feat = function(sol_tmp,focal){
  # This function takes a solution and the id of the focal
  # species and it extracts features
  # related to species in the community, such as the
  # suppliers ids, their number, individual and community
  # biomass (growth), and their effective number taking into account
  # their relative contribution to community biomass
  id_sp_growth=grep(names(sol_tmp),pattern="Growth") # look for sp growth id, first id is focal sp.
  sp_growth=unlist(sol_tmp[id_sp_growth]) # extract sp growth
  comm_growth=unlist(sol_tmp["community_growth"]) # extract whole community growth
  names(comm_growth)="comm_growth"
  p_sp_growth=sp_growth/comm_growth
  N_sp=length(id_sp_growth)
  N_eff_sp=exp(-sum(p_sp_growth*log(p_sp_growth)))
  #focal=names(sol_tmp)[id_sp_growth[1]] # this was wrong the focal not always is the first  sp
  #focal=sub(focal,pattern = "Growth_",replacement = "") # I include it now as an argument
  if(length(id_sp_growth) > 1){
    suppliers=names(sol_tmp)[id_sp_growth] # get entries containing "growth"
    suppliers=sub(suppliers,pattern = "Growth_",replacement = "") # remove "growth"
    suppliers=suppliers[!suppliers %in% focal] # remove the focal strain
  }else{
    suppliers="none"
  }
  return(list("sp_growth"=sp_growth,"comm_growth"=comm_growth,
              "N_sp"=N_sp,"N_eff_sp"=N_eff_sp,
              "focal"=focal,"suppliers"=suppliers))
}


