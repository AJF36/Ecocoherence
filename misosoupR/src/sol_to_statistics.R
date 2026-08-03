sol_to_statistics=function(net.tmp){
  # This function takes a solution in network format and extracts
  # some statistics
  N_met=length(unique(net.tmp$nodeA)) # number of metabolites
  N_sp=length(unique(net.tmp$nodeB)) # of species
  N_react=dim(net.tmp)[1] # of reactions
  id_pos=which(net.tmp$weight > 0)
  id_neg=which(net.tmp$weight < 0)
  tot_flux=sum(abs(net.tmp$weight))
  tot_flux_pos=sum(abs(net.tmp$weight[id_pos]))
  tot_flux_neg=sum(abs(net.tmp$weight[id_neg]))
  N_react_pos=dim(net.tmp[id_pos,])[1] # of positive reactions
  N_react_neg=dim(net.tmp[id_neg,])[1] # of negative reactions
  p=abs(net.tmp$weight)/tot_flux
  p_pos=abs(net.tmp$weight[id_pos])/tot_flux_pos
  p_neg=abs(net.tmp$weight[id_neg])/tot_flux_neg
  S= -sum(p*log(p))
  S_pos = -sum(p_pos*log(p_pos))
  S_neg = -sum(p_neg*log(p_neg))
  N_eff_react=exp(S) # effective  number  of reactions
  N_eff_react_pos=exp(S_pos) # (among those the number of positive ones)
  N_eff_react_neg=exp(S_neg) # (of negative)
  out=c(N_sp,N_met,N_react,
        N_react_pos,N_react_neg,
        N_eff_react,N_eff_react_pos,N_eff_react_neg,
        tot_flux,tot_flux_pos,tot_flux_neg)
  names(out)=c("N_sp","N_met","N_react",
               "N_react_pos","N_react_neg",
               "N_eff_react","N_eff_react_pos","N_eff_react_neg",
               "tot_flux","tot_flux_pos","tot_flux_neg")
  return(out)
}
