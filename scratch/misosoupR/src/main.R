
rm(list=ls())

library(reshape)
library(yaml)
library(readr)
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirData=paste(this.dir,"/minimal_communities/",sep="") # directory for the solutions and minimal media

# load src

setwd(dirSrc)
source("load_solutions.R")
source("sol_to_exchange_network.R")
source("rename_net.R")
source("net_aggr_to_out.R")
source("aggregate_networks.R")
source("output_networks.R")
source("heatmap.2.mod.R")
source("functionArgsList.R")
source("summary_stats_solutions.R")
source("sol_to_statistics.R")
source("sol_to_sp_feat.R")
source("update_matrices.R")
source("match_sets_pairwise.R")
source("match_sets_multiple.R")


# PRELIMINARIES ---------
# --- Define minimal medium file
#fileMedium="medium_MBM_no_co2_hco3.yaml"
fileMedium="media.yaml"

# .... Load minimal medium
setwd(dirData)
met_medium = read_yaml(file = fileMedium)

# --- Load MisoSoup output
# .... select some examples to work with
# file.list=c("example1.yaml","example2.yaml","example3.yaml")-
file.list=c("misosoup_Agarose_mod_1.yaml")
#file.list=c("misosoup_Agarose.yaml","misosoup_Agarose_mod_48.yaml")


# ... load solutions
sols=load_solutions(file.list)

# Example of structures and basic tools -------------
# ... inspect solutions: The expected structure (and the one generated
#     in downstream computations) is a list with the following  hierarchy 
#     `list[[medium]][[strain]][[solution]]`.
names(sols) # C sources
names(sols[[names(sols)[1]]]) # strains in C source 1
names(sols[[names(sols)[2]]]) # strains in C source 2
# length(sols[["gal"]][["A1R12"]]) # how many solutions in acetate for strain A1R12?

sol_tmp=sols[[1]][[1]][[1]] # extract first solution of first strain and first C source
# str(sol_tmp)
net.tmp=sol_to_exchange_network(sol_tmp,met_medium) # convert into a network
net.tmp.ren=rename_net(net.tmp) # rename network to distinguish excreted and consumed metabolites by name

# "STAT" ANALYSIS ------------------
# Extract summary statistics of each solution 
out_stats=summary_stats_solutions(sols,met_medium)
out.df=out_stats$stats.df
out.mat=out_stats$mat.list

write_tsv(out.df,paste("stats_df_",file.list,".tsv",sep=""))
################ Until here the stats df

#### TO DO --> INCLUDE PLOTS OF THIS MODE

# "AGGR" ANALYSIS -------------------------
# Example 1 ----------
# --- Aggregate only solutions, keeping media and strains explicit
nets.list = aggregate_networks(sols,met_medium,
                               mode_strains="split",
                               mode_media = "split",
                               mode_sols = "aggr")

nets.flux=nets.list[["nets.flux"]]
nets.count=nets.list[["nets.count"]]

# --- Plot the solutions 
# dirOut="media-split_strains-split_sols-aggr"
dirOut=gsub(".yaml","",file.list)
dirOut
dir.create(dirOut)

# ..... First, the network of fluxes, we want both networks and heatmap
nets=nets.flux
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = TRUE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="none", # we remove the zscore for columns, which is the default
                                   col="rainbow"),
                flux=TRUE,
                par.qgraph = list(maximum=10),
                pathOut = dirOut)

# .... for counts, the graph carries less information, so we do not generate it
nets=nets.count
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = FALSE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="none",
                                   col="rainbow"),
                flux=FALSE,
                par.qgraph = list(maximum=7),
                pathOut = dirOut)




# Example 2 ----------
# --- Aggregate all solutions
nets.list = aggregate_networks(sols,met_medium,
                               mode_strains="aggr",
                               mode_media = "aggr",
                               mode_sols = "aggr")

nets.flux=nets.list[["nets.flux"]]
nets.count=nets.list[["nets.count"]]

# --- Plot the solutions 
setwd(dirData)
# dirOut="media-aggr_strains-aggr_sols-aggr"
dirOut <- file.list
dir.create(dirOut)

# ..... First, the network of fluxes, we want both networks and heatmap
nets=nets.flux
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = TRUE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="column", # "none",#
                                   #sepwidth=c(0,0),
                                   cexRow=0.5, # we remove the zscore for columns, which is the default
                                   col="bluered"), # "rainbow"),
                flux=TRUE,
                par.qgraph = list(maximum=10),
                pathOut = dirOut)

# .... for counts, the graph carries less information, so we skip it
nets=nets.count
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = FALSE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="none",
                                   col="rainbow"),
                flux=FALSE,
                par.qgraph = list(maximum=7),
                pathOut = dirOut)

# Example 3 ----------
# --- Strains explicit, aggregate all solutions across all media
nets.list = aggregate_networks(sols,met_medium,
                               mode_strains="split",
                               mode_media = "aggr",
                               mode_sols = "aggr")

nets.flux=nets.list[["nets.flux"]]
nets.count=nets.list[["nets.count"]]

# --- Plot the solutions 
setwd(dirData)
dirOut="media-aggr_strains-split_sols-aggr"
dir.create(dirOut)

# ..... First, the network of fluxes, we want both networks and heatmaps
nets=nets.flux
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = TRUE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="column", # "none",#
                                   #sepwidth=c(0,0),
                                   cexRow=0.5, # we remove the zscore for columns, which is the default
                                   col="bluered"), # "rainbow"),
                flux=TRUE,
                par.qgraph = list(maximum=10),
                pathOut = dirOut)

# .... for counts, the graph carries less information, so we skip it
nets=nets.count
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = FALSE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="none",
                                   col="rainbow"),
                flux=FALSE,
                par.qgraph = list(maximum=7),
                pathOut = dirOut)
# Example 4 ----------
# --- Everything explicit, but only for strain B3M02 (there is only one medium, f6p)
nets.list = aggregate_networks(sols,met_medium,
                               vec.strains = "B3M02", 
                               mode_strains="split",
                               mode_media = "split",
                               mode_sols = "split")

nets.flux=nets.list[["nets.flux"]]
nets.count=nets.list[["nets.count"]]

# --- Plot the solutions 
setwd(dirData)
dirOut="media-f6p_strains-B3M02_sols-split"
dir.create(dirOut)

# ..... First, the network of fluxes, we want both networks and heatmap
nets=nets.flux
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = TRUE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="column", # "none",#
                                   #sepwidth=c(0,0),
                                   cexRow=0.5, # we remove the zscore for columns, which is the default
                                   col="bluered"), # "rainbow"),
                flux=TRUE,
                par.qgraph = list(maximum=10),
                pathOut = dirOut)

# .... for counts, the graph carries less information, so we skip it
nets=nets.count
output_networks(nets,
                plot.heatmap = TRUE,plot.graph = FALSE,
                par.heatmap = list(xlab = "Isolates",
                                   scale="none",
                                   col="rainbow"),
                flux=FALSE,
                par.qgraph = list(maximum=7),
                pathOut = dirOut)

# COMPARISON SETS OF SOLUTIONS -----------

setwd(dirData)

# .... select some examples to work with
#      These are simply the concatenation of the examples used
#      above, e.g. example1-2.yaml contains both the solutions of
#      example1.yaml and example2.yaml
file.list=c("example1-2.yaml","example1-3.yaml","example2-3.yaml","example1-2-3.yaml")

# .... load solutions
sets.list=list()
i=0
for(file.in in file.list){
  sets.list[[file.in]]=read_yaml(file.in)
}

# Example 1 -----
# Match only two sets
# .... Take a couple of sets
setA=sets.list[[1]] # "example1-2.yaml"
setB=sets.list[[2]] # "example1-3.yaml"

# .... Match their solutions
setsAB.list=match_sets_pairwise(setA,setB)

# .... We now have different lists describing  
solsAB_inA=setsAB.list$solsAB_inA # what they have in common (example1)
solsAB_inB=setsAB.list$solsAB_inB
solsA_not_inB=setsAB.list$solsA_not_inB # and those that differ (example 2)
solsB_not_inA=setsAB.list$solsB_not_inA #  (example3)

# Example 2 -----
# Sanity check
# .... Match three set that have no common solutions
sets.list.tmp=sets.list[c(1,2,3)] # 1-2, 1-3 and 2-3

match.list=match_sets_multiple(sets.list.tmp) # warnings should appear

# ... recover the sets of solutions in common
# match.list$sets_int # returns three empty lists, because these three sets have
                    # no solution in common, by  construction.
# ... recover the sets of solutions that do not belong to the intersection
# match.list$sets_not_int # returns three lists, each of them identical to
                    # the input set, because the intersection is empty

# Example 3 -----
# .... Match three sets having the solutions in example1 in common
sets.list.tmp=sets.list[c(1,2,4)] # 1-2, 1-3 and 1-2-3

match.list=match_sets_multiple(sets.list.tmp) # warnings should appear

# ... recover the sets of solutions in common
# match.list$sets_int # returns three identical lists, intersection of 1-2 and 1-3
# ... recover the sets of solutions that do not belong to the intersection
# match.not_int.list=match.list$sets_not_int # for 1-2 its should return the only solution
# in example2, for 1-3 those in example3, for 1-2-3 the union of example2 and example3
