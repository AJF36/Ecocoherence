##################################################
# time_based_clustering.R
##################################################
# In this script I perform ...
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())

results_dir = "results/time_based_clustering"
figures_dir = "figures/time_based_clustering"

results_dir = c(results_dir,figures_dir)

for (dir in results_dir) {
  
  if(!(fs::dir_exists(dir))) {
  
    fs::dir_create(dir) }

}

library(mia)
library(tidyverse)

# Load the main mia_object 
full_mia = readRDS("results/mia_analysis/mia_main_final.RDS")
full_mia$Time = as.factor(full_mia$Time)
# Split it on substrates
list_mia_by_substrate = mia::splitOn(full_mia,"Substrate")
which(names(list_mia_by_substrate) == "Chitosan")
list_mia_by_substrate_f  = list_mia_by_substrate[-4]
# Create the mia_modules objects

funk_dir = "functionink/"

module_mia_list = imap(list_mia_by_substrate_f, function(mia_obj,substrate) {

  substrate_funk_dir = glue::glue("functionink/{substrate}/functionink_tmp")
    
  funk_partition_path = fs::dir_ls(substrate_funk_dir,regexp = "Partition-.+_.tsv$")
    
  funk_partition =  suppressMessages( readr::read_tsv(funk_partition_path,skip = 9,col_names = F))
    
  colnames(funk_partition) = c("ASV","module") 

  # Join this column to the tax table, change it in the mia object and aggregate

  tax = as.data.frame(rowData(mia_obj))
  tax_funk = inner_join(tax,funk_partition, by = "ASV")
  tax_funk[["module"]] = as.factor(tax_funk[["module"]])
  print(nrow(tax_funk))
  # Filter the ASV that have module
  
  mia_obj_f = mia_obj[rownames(mia_obj) %in% tax_funk[["ASV"]],]
  
  rowData(mia_obj_f) = tax_funk
  print(colnames(rowData(mia_obj_f)))
  
  module_mia = mia::agglomerateByVariable(mia_obj_f,by = "rows",group = "module")
  
  # Load substrate 
    

})


# Collapse the replicates: for each time point keep a single "sample" whose
# abundances are the mean of the relative abundances of its replicates.
# (each replicate column sums to 1, so the averaged column sums to 1 as well)

mean_by_time = function(mia_obj) {

  mia_obj = mia::transformAssay(mia_obj,method = "relabundance")

  rel = as.matrix(assay(mia_obj,"relabundance"))
  time = droplevels(as.factor(colData(mia_obj)[["Time"]]))

  # sum the replicates of each time point and divide by how many they are
  sums = rowsum(t(rel), group = time)
  n_rep = as.vector(table(time)[rownames(sums)])
  mean_mat = t(sums / n_rep)

  # Keep the metadata that is constant within a time point (Substrate, SubType,
  # Time, Stage ...), drop the ones that change between replicates (sampleid, Replica)
  cd = as.data.frame(colData(mia_obj))
  first_idx = match(colnames(mean_mat), as.character(time))
  constant_col = vapply(cd, function(x) all(tapply(as.character(x), time, function(v) length(unique(v)) == 1)), logical(1))
  cd_mean = cd[first_idx, constant_col, drop = FALSE]
  cd_mean[["n_replicates"]] = n_rep
  rownames(cd_mean) = colnames(mean_mat)

  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(relabundance = mean_mat),
    rowData = rowData(mia_obj),
    colData = S4Vectors::DataFrame(cd_mean))
}

#PCoA by ASV
iwalk(list_mia_by_substrate, function(mia_obj,substrate) {

  #1.transform to relative abundance and average the replicates of each time point
  mia_obj = mean_by_time(mia_obj)
  #2.Dissimilarity matrix
  # JSD_matrix = mia::getDissimilarity(mia_obj,assay.type = "relabundance", method = "jsd")
  #3.calculate pcoa
  mia_obj = mia::addMDS(mia_obj,FUN = getDissimilarity,assay.type = "relabundance",method = "jsd")

  # Get the variancve explained
  e = attr(reducedDim(mia_obj,"MDS"),"eig")
  rel_eig <- e / sum(e[e > 0])

  p = scater::plotReducedDim(mia_obj,"MDS",colour_by = "Time") + 
    ggtitle(label = glue::glue("Pcoa {substrate} colored by time (replicates averaged)")) +
    labs(
    x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
    y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")) +
    theme_classic(base_size = 16) +
    scale_color_manual(values = pals::glasbey(n = 32))
  
  ggsave(filename = glue::glue("figures/time_based_clustering/PCoA_by_asv_{substrate}.pdf"),plot = p)
})

# Same logic but with the modules, does this improve the separation in the PCoA=

iwalk(module_mia_list, function(mia_obj,substrate) {

  #1.transform to relative abundance and average the replicates of each time point
  mia_obj = mean_by_time(mia_obj)
  #2.Dissimilarity matrix
  # JSD_matrix = mia::getDissimilarity(mia_obj,assay.type = "relabundance", method = "jsd")
  #3.calculate pcoa
  mia_obj = mia::addMDS(mia_obj,FUN = getDissimilarity,assay.type = "relabundance",method = "jsd")

  # Get the variancve explained
  e = attr(reducedDim(mia_obj,"MDS"),"eig")
  rel_eig <- e / sum(e[e > 0])

  p = scater::plotReducedDim(mia_obj,"MDS",colour_by = "Time") + 
    ggtitle(label = glue::glue("Pcoa {substrate} by module colored by time (replicates averaged)")) +
    labs(
    x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
    y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")) +
    theme_classic(base_size = 16) +
    scale_color_manual(values = pals::glasbey(n = 32))
  
    ggsave(filename = glue::glue("figures/time_based_clustering/PCoA_by_module_{substrate}.pdf"),plot = p)

})

##################################################
# Rate of change of the community composition
##################################################
# The PCoAs look like a gradient rather than a set of clusters, so instead of
# forcing groups we ask where along the gradient the community changes fastest:
#   - JSD between consecutive time points -> peaks = transitions
#   - the same divided by the elapsed hours, because the sampling is not evenly
#     spaced (12h at the beginning, up to 48h at the end) and otherwise the last
#     intervals look faster just because they are wider
#   - JSD to t0 -> when does the community stop turning over

rate_of_change = function(mia_obj) {

  mia_obj = mean_by_time(mia_obj)

  d = as.matrix(mia::getDissimilarity(mia_obj,assay.type = "relabundance",method = "jsd"))

  time = as.numeric(as.character(colData(mia_obj)[["Time"]]))
  ord = order(time)
  d = d[ord,ord]
  time = time[ord]
  n = length(time)

  # the dissimilarity of each time point with the next one
  time_from = time[-n]
  time_to = time[-1]
  time_mid = (time_from + time_to) / 2   # a rate belongs to the middle of the interval
  delta_t = time_to - time_from
  jsd_consecutive = d[cbind(seq_len(n - 1), 2:n)]

  bind_rows(
    tibble(time = time_mid, time_from = time_from, time_to = time_to,
           metric = "JSD between consecutive time points", value = jsd_consecutive),
    tibble(time = time_mid, time_from = time_from, time_to = time_to,
           metric = "JSD per hour", value = jsd_consecutive / delta_t),
    tibble(time = time,
           time_from = NA_real_, time_to = NA_real_,
           metric = "JSD to t0", value = d[1,]))
}

metric_levels = c("JSD between consecutive time points","JSD per hour","JSD to t0")

roc = bind_rows(
  imap(list_mia_by_substrate, ~ mutate(rate_of_change(.x), substrate = .y, level = "ASV")) |> bind_rows(),
  imap(module_mia_list,       ~ mutate(rate_of_change(.x), substrate = .y, level = "Module")) |> bind_rows()) |>
  mutate(metric = factor(metric, levels = metric_levels))

readr::write_tsv(roc,"results/time_based_clustering/rate_of_change.tsv")

# The time point where the composition changes fastest, i.e. the main transition
peaks = roc |>
  filter(metric == "JSD per hour") |>
  group_by(substrate,level) |>
  slice_max(value, n = 1) |>
  ungroup()

readr::write_tsv(peaks,"results/time_based_clustering/rate_of_change_peaks.tsv")

# One figure per substrate, ASV and module level together
iwalk(split(roc,roc[["substrate"]]), function(df,substrate) {

  p = ggplot(df,aes(x = time,y = value,colour = level)) +
    geom_vline(data = filter(peaks, substrate == .env$substrate),
               aes(xintercept = time,colour = level),linetype = "dashed",alpha = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ metric,ncol = 1,scales = "free_y") +
    labs(title = glue::glue("Rate of change of the community in {substrate}"),
         subtitle = "dashed line: fastest change (replicates averaged)",
         x = "Time (h)",y = "Jensen-Shannon divergence",colour = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top") +
    scale_colour_manual(values = c(ASV = "#1B6CA8", Module = "#D1495B"))

  ggsave(filename = glue::glue("figures/time_based_clustering/rate_of_change_{substrate}.pdf"),
         plot = p,width = 7,height = 8)
})

# All the substrates together, to see if the transitions happen at the same time
p_all = roc |>
  filter(metric == "JSD per hour") |>
  ggplot(aes(x = time,y = value,colour = substrate)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ level,ncol = 1,scales = "free_y") +
  labs(title = "Rate of change of the community composition",
       x = "Time (h)",y = "JSD per hour",colour = NULL) +
  theme_classic(base_size = 14) +
  scale_colour_manual(values = pals::glasbey(n = 32))

ggsave(filename = "figures/time_based_clustering/rate_of_change_all_substrates.pdf",
       plot = p_all,width = 8,height = 7)

View(as.data.frame(colData(full_mia)))


##################################################
# Assign a stage to each family and to each module
##################################################
# Based in the results of this analysis and previous papers I have decided to divide in
# Three different stages attachment = t0, early = 12-48 and late = 60-204
#
# Same idea as in 6g_assign_ecological_strategy_to_families.R: build an ideal
# vector for each stage (flat inside the stage, zero outside), add them as extra
# rows to the temporal profiles and compute the JSD between all the rows.
#
# The stage of a family/module is the ideal vector it is closest to. The label
# based on the clustering (the ideal vector that falls in its cluster) is kept as
# a second column to compare, but it is not used as the assignment: it depends on
# the linkage and on where the cut stops, and an ideal vector that stays isolated
# leaves its stage empty even when there are features with that shape.
#
# G (flat over all the time points) is there so that a feature that does not
# prefer any stage has a prototype to go to, instead of being forced into the
# closest of A/E/L.

stage_levels = c("A","E","L")            # the three stages in time
ideal_levels = c("A","E","L","G")        # ... plus the generalist vector

stage_of_time = function(time) {
  cut(time,breaks = c(-Inf,0,48,Inf),labels = stage_levels)   # A = t0, E = 12-48, L = 60-204
}

# Temporal profile: each row has to sum 1 over time to be comparable with the
# ideal vectors (this removes the abundance, only the shape in time is left)
temporal_profile = function(mia_obj) {

  mia_obj = mean_by_time(mia_obj)

  mat = as.matrix(assay(mia_obj,"relabundance"))
  time = as.numeric(as.character(colData(mia_obj)[["Time"]]))
  ord = order(time)
  mat = mat[,ord,drop = FALSE]
  time = time[ord]

  # A feature that is absent in this substrate has nothing to say about the stages
  mat = mat[rowSums(mat) > 0,,drop = FALSE]

  list(prof = mat / rowSums(mat),time = time)
}

# Ideal vectors: flat inside the stage and zero outside, plus a flat one (G)
ideal_vectors = function(time) {
  stage = stage_of_time(time)
  ideal = t(vapply(stage_levels,function(s) as.numeric(stage == s) / sum(stage == s),numeric(length(time))))
  rbind(ideal,G = rep(1 / length(time),length(time)))
}

# JSD of every profile against every ideal vector, with the same definition and
# the same pseudocount that BiotypeR::dist.JSD uses, but without building the
# whole n x n matrix (needed at the ASV level, where n is several thousands)
jsd_to_ideal = function(prof,ideal) {

  pseudocount = 10^(round(log10(min(prof[prof > 0])),0) - 1)
  p = prof; p[p == 0] = pseudocount
  q = ideal; q[q == 0] = pseudocount

  out = matrix(0,nrow(p),nrow(q),dimnames = list(rownames(p),rownames(q)))
  for (j in seq_len(nrow(q))) {
    qj = matrix(q[j,],nrow = nrow(p),ncol = ncol(p),byrow = TRUE)
    m = (p + qj) / 2
    out[,j] = sqrt(0.5 * rowSums(p * log(p / m)) + 0.5 * rowSums(qj * log(qj / m)))
  }
  out
}

assign_stage = function(mia_obj) {

  tp = temporal_profile(mia_obj)
  prof = tp[["prof"]]
  time = tp[["time"]]

  ideal = ideal_vectors(time)
  colnames(ideal) = colnames(prof)

  combined = rbind(ideal,prof)

  # BiotypeR expects the observations in the columns
  d = BiotypeR::dist.JSD(t(combined))
  hc = hclust(d,method = "ward.D2")

  # Cut deeper and deeper until each ideal vector ends in a different cluster
  k = 1
  cl = cutree(hc,k = k)
  while (length(unique(cl[ideal_levels])) != length(ideal_levels) && k < length(cl)) {
    k = k + 1
    cl = cutree(hc,k = k)
  }

  # The stage of a cluster is the ideal vector it contains, if any
  cluster_of_stage = setNames(rep(NA_character_,k),seq_len(k))
  cluster_of_stage[as.character(cl[ideal_levels])] = ideal_levels

  dm = as.matrix(d)
  features = setdiff(rownames(combined),ideal_levels)
  d_ideal = dm[features,ideal_levels,drop = FALSE]

  # The assignment: the closest ideal vector. The margin (how much closer the
  # first one is than the second) says how clear the assignment is
  nearest = ideal_levels[max.col(-d_ideal,ties.method = "first")]
  two_best = apply(d_ideal,1,function(v) sort(v)[1:2])

  res = tibble(
    feature = features,
    stage = nearest,
    margin = two_best[2,] - two_best[1,],
    # the label that the clustering would give, only to compare
    cluster = unname(cl[features]),
    cluster_stage = coalesce(unname(cluster_of_stage[as.character(cl[features])]),"no_stage"),
    jsd_A = d_ideal[,"A"], jsd_E = d_ideal[,"E"], jsd_L = d_ideal[,"L"], jsd_G = d_ideal[,"G"],
    n_clusters = k)

  prof_long = as_tibble(prof,rownames = "feature") |>
    pivot_longer(-feature,names_to = "sample",values_to = "profile") |>
    left_join(tibble(sample = colnames(prof),time = time),by = "sample")

  list(assignment = res, hclust = hc, k = k, profiles = prof_long)
}

# Family level (all the ASVs, aggregated by Family) and module level
family_mia_list = map(list_mia_by_substrate,~ suppressWarnings(mia::agglomerateByRank(.x,rank = "Family")))

stage_input = list(Family = family_mia_list, Module = module_mia_list)

stage_res = imap(stage_input,function(mia_list,level) {
  imap(mia_list,function(mia_obj,substrate) {
    out = assign_stage(mia_obj)
    out[["substrate"]] = substrate
    out[["level"]] = level
    out
  })
})

# --- Results table
stage_table = imap(stage_res,function(by_substrate,level) {
  imap(by_substrate,~ mutate(.x[["assignment"]],substrate = .y,level = level)) |> bind_rows()
}) |> bind_rows() |>
  mutate(stage = factor(stage,levels = ideal_levels),
         cluster_stage = factor(cluster_stage,levels = c(ideal_levels,"no_stage")))

readr::write_tsv(stage_table,"results/time_based_clustering/stage_assignment.tsv")

# How many families/modules fall in each stage
stage_summary = count(stage_table,level,substrate,stage)
readr::write_tsv(stage_summary,"results/time_based_clustering/stage_assignment_summary.tsv")

# ... and how much the clustering label would have differed
stage_vs_cluster = count(stage_table,level,substrate,stage,cluster_stage)
readr::write_tsv(stage_vs_cluster,"results/time_based_clustering/stage_assignment_vs_cluster.tsv")

# --- Dendrograms, one per substrate and level
fs::dir_create("figures/time_based_clustering/stage_dendrograms")

stage_colours = c(A = "#1B6CA8", E = "#EDAE49", L = "#D1495B", G = "#4C956C", no_stage = "grey60")

walk(stage_res,function(by_substrate) {
  walk(by_substrate,function(out) {

    dend = as.dendrogram(out[["hclust"]])
    dend = dendextend::set(dend,"branches_k_color",k = out[["k"]])

    lab = labels(dend)
    is_ideal = lab %in% ideal_levels

    # Colour of each leaf: the ideal vector it has been assigned to
    assigned = setNames(as.character(out[["assignment"]][["stage"]]),out[["assignment"]][["feature"]])
    leaf_stage = ifelse(is_ideal,lab,assigned[lab])
    leaf_colour = unname(stage_colours[leaf_stage])

    # The labels are rotated, so what has to fit in the space of each leaf is the
    # height of the line (~12pt at cex 1). Give each leaf a fixed slot and take
    # the cex from it, so that the labels never overlap
    pdf_width = max(15,length(lab) * 0.07)
    slot_pt = pdf_width * 72 / length(lab)
    feature_cex = min(0.5,max(0.1,slot_pt / 13))

    dend = dendextend::set(dend,"labels_cex",ifelse(is_ideal,1.1,feature_cex))
    dend = dendextend::set(dend,"labels_colors",leaf_colour)
    # Make the ideal vectors impossible to miss among the features of their colour
    dend = dendextend::set(dend,"labels",ifelse(is_ideal,paste0("<<",lab,">>"),lab))

    pdf(glue::glue("figures/time_based_clustering/stage_dendrograms/dendrogram_{out$level}_{out$substrate}.pdf"),
        width = pdf_width,height = if (out$level == "Family") 22 else 12)
    par(mar = c(if (out$level == "Family") 14 else 8,4,4,2))
    plot(dend,main = glue::glue("{out$level} clustering by temporal profile - {out$substrate} (k = {out$k})"))
    dendextend::rect.dendrogram(dend,k = out[["k"]],border = 8)

    # A colour strip with the assignment, readable also when there are many leaves.
    # It goes under the labels: they are rotated, so the space they take is their
    # width, which has to be converted from inches to the units of the plot
    usr = par("usr")
    y_per_inch = (usr[4] - usr[3]) / par("pin")[2]
    labels_height = max(strwidth(lab,units = "inches",cex = feature_cex)) * y_per_inch
    bar_top = -labels_height - 0.02 * (usr[4] - usr[3])
    bar_height = 0.03 * (usr[4] - usr[3])

    rect(xleft = seq_along(lab) - 0.5,xright = seq_along(lab) + 0.5,
         ybottom = bar_top - bar_height,ytop = bar_top,
         col = leaf_colour,border = NA,xpd = NA)
    text(x = 0,y = bar_top - bar_height / 2,labels = "assigned  ",adj = 1,cex = 0.9,xpd = NA)
    legend("topright",legend = ideal_levels,text.col = stage_colours[ideal_levels],
           title = "Assigned to",bty = "n",cex = 1.2)
    dev.off()
  })
})

# --- Mean temporal profile of the features assigned to each stage,
#     to check that the assignment makes sense
profile_table = imap(stage_res,function(by_substrate,level) {
  imap(by_substrate,function(out,substrate) {
    left_join(out[["profiles"]],select(out[["assignment"]],feature,stage),by = "feature") |>
      mutate(substrate = substrate,level = level)
  }) |> bind_rows()
}) |> bind_rows()

iwalk(split(profile_table,profile_table[["level"]]),function(df,level) {

  p = df |>
    group_by(level,substrate,stage,time) |>
    summarise(profile = mean(profile),.groups = "drop") |>
    ggplot(aes(x = time,y = profile,colour = stage)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ substrate,scales = "free_y") +
    labs(title = glue::glue("Mean temporal profile of the {level} assigned to each stage"),
         subtitle = "assigned to the closest ideal vector (G = generalist, flat in time)",
         x = "Time (h)",y = "Relative abundance profile (rows sum 1)",colour = "Stage") +
    theme_classic(base_size = 13) +
    scale_colour_manual(values = stage_colours)

  ggsave(filename = glue::glue("figures/time_based_clustering/stage_profiles_{level}.pdf"),
         plot = p,width = 11,height = 7)
})

# How many features per stage, all the substrates together
p_counts = stage_table |>
  ggplot(aes(x = substrate,fill = stage)) +
  geom_bar(position = "fill") +
  facet_wrap(~ level,ncol = 1) +
  labs(title = "Fraction of families / modules assigned to each stage",
       subtitle = "assigned to the closest ideal vector (G = generalist, flat in time)",
       x = NULL,y = "Fraction",fill = "Stage") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values = stage_colours)

ggsave(filename = "figures/time_based_clustering/stage_assignment_counts.pdf",plot = p_counts,width = 8,height = 8)

##################################################
# Coherence of the families with respect to the stages
##################################################
# Same idea as the coherence of the families with respect to the functionink
# modules (6e + 7a + 8a), changing the partition: instead of the modules, the
# four stages. The members are the ASVs of the family present in the substrate,
# each one assigned to the ideal vector its temporal profile is closest to.
#
#   p_i = ASVs of the family assigned to stage i / n_members
#   S   = -sum(p_i log p_i)        Shannon entropy, as in 7a
#   X   = exp(S)                   effective number of stages, between 1 and 4
#   Z   = (X - mean_X_null) / sd_X_null
#
# The null permutes the family label of the ASVs, exactly as 6e does: it keeps
# the size of every family and the global frequency of the stages, so a family
# ends up with n_members ASVs taken at random from the substrate. Z < 0 means
# that the family is split in fewer stages than expected, i.e. it is coherent.

n_perm_stage = 500
min_members = 4      # same threshold as 6e/7a

effective_n_stages = function(counts) {
  p = counts / rowSums(counts)
  S = -rowSums(ifelse(p > 0,p * log(p),0))
  exp(S)
}

# --- Stage of every ASV, substrate by substrate
asv_stage = imap(list_mia_by_substrate,function(mia_obj,substrate) {

  tp = temporal_profile(mia_obj)
  ideal = ideal_vectors(tp[["time"]])
  colnames(ideal) = colnames(tp[["prof"]])

  d = jsd_to_ideal(tp[["prof"]],ideal)

  tibble(substrate = substrate,
         ASV = rownames(d),
         stage = factor(ideal_levels[max.col(-d,ties.method = "first")],levels = ideal_levels))
}) |> bind_rows()

# The family of each ASV, with the same cleaning of the names that 6a/6e do
asv_family = as.data.frame(rowData(full_mia)) |>
  transmute(ASV,Family = gsub("_incertae_sedis|_Incertae Sedis XI","",Family)) |>
  as_tibble()

asv_stage = inner_join(asv_stage,asv_family,by = "ASV")

readr::write_tsv(asv_stage,"results/time_based_clustering/asv_stage_assignment.tsv")

# --- Observed X and the permutation null, per substrate
stage_coherence = imap(split(asv_stage,asv_stage[["substrate"]]),function(df,substrate) {

  counts = table(df[["Family"]],df[["stage"]])
  counts = counts[rowSums(counts) >= min_members,,drop = FALSE]
  families = rownames(counts)

  X = effective_n_stages(counts)

  null_X = replicate(n_perm_stage,{
    counts_perm = table(sample(df[["Family"]]),df[["stage"]])
    effective_n_stages(counts_perm[families,,drop = FALSE])
  })

  mean_X = rowMeans(null_X)
  sd_X = apply(null_X,1,sd)

  tibble(substrate = substrate,
         family = families,
         n_members = as.vector(rowSums(counts)),
         n_A = counts[,"A"],n_E = counts[,"E"],n_L = counts[,"L"],n_G = counts[,"G"],
         S = log(X),
         X = as.vector(X),
         mean_X = mean_X,
         sd_X = sd_X,
         Z = ifelse(sd_X > 0,(X - mean_X) / sd_X,NA_real_))
}) |> bind_rows()

readr::write_tsv(stage_coherence,"results/time_based_clustering/stage_coherence_families.tsv")

# How many families are coherent (Z < -2.5) or anti-coherent (Z > 2.5)
coherence_summary = stage_coherence |>
  group_by(substrate) |>
  summarise(n_families = n(),
            coherent = sum(Z < -2.5,na.rm = TRUE),
            anti_coherent = sum(Z > 2.5,na.rm = TRUE),
            median_Z = median(Z,na.rm = TRUE),.groups = "drop")

readr::write_tsv(coherence_summary,"results/time_based_clustering/stage_coherence_summary.tsv")
print(coherence_summary)

# --- Heatmap of the Z scores, as in 8a
plot_coherence_heatmap = function(df,file,height) {

  fam_order = df |> group_by(family) |> summarise(m = mean(Z,na.rm = TRUE)) |> arrange(m) |> pull(family)
  df = mutate(df,family = factor(family,levels = fam_order))

  p_z = df |>
    mutate(star = ifelse(!is.na(Z) & abs(Z) > 2.5,"*","")) |>
    ggplot(aes(x = substrate,y = family,fill = Z)) +
    geom_tile() +
    geom_text(aes(label = star),colour = "black",size = 2.5,vjust = 0.75) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_gradient2(low = "#1B6CA8",mid = "white",high = "#D1495B",midpoint = 0,na.value = "grey90") +
    labs(title = "Coherence of the families with respect to the stages",
         subtitle = "Z of the effective number of stages\nZ < 0 = fewer stages than expected (coherent). * |Z| > 2.5",
         x = NULL,y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),
          axis.text.y = element_text(size = 5))

  # Annotation: the vector of probabilities that gives the entropy, that is, the
  # fraction of the members of the family that is in each stage
  p_ann = df |>
    select(family,substrate,n_A,n_E,n_L,n_G) |>
    pivot_longer(starts_with("n_"),names_to = "stage",values_to = "n") |>
    mutate(stage = factor(sub("^n_","",stage),levels = ideal_levels)) |>
    group_by(family,substrate) |>
    mutate(p = n / sum(n)) |>
    ungroup() |>
    ggplot(aes(x = p,y = family,fill = stage)) +
    geom_col(width = 1) +
    facet_wrap(~ substrate,nrow = 1) +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(breaks = c(0,1),expand = c(0,0)) +
    scale_fill_manual(values = stage_colours) +
    # The header has to have the same number of lines as the one of the heatmap,
    # otherwise the two panels do not end up aligned
    labs(title = " ",
         subtitle = "Fraction of the members of the family in each stage\n(the p of the entropy)",
         x = NULL,y = NULL,fill = "Stage") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6),
          panel.grid = element_blank(),
          panel.spacing = unit(1.5,"pt"),
          strip.text = element_text(size = 6,angle = 90))

  p = patchwork::wrap_plots(p_z,p_ann,widths = c(1,1.6))

  ggsave(filename = file,plot = p,width = 14,height = height,limitsize = FALSE)
}

plot_coherence_heatmap(stage_coherence,
                       "figures/time_based_clustering/stage_coherence_heatmap.pdf",
                       height = max(7,0.09 * length(unique(stage_coherence[["family"]]))))

# ... and only the families that are significant somewhere, to be able to read it
significant_families = stage_coherence |> filter(abs(Z) > 2.5) |> pull(family) |> unique()

plot_coherence_heatmap(filter(stage_coherence,family %in% significant_families),
                       "figures/time_based_clustering/stage_coherence_heatmap_significant.pdf",
                       height = max(7,0.16 * length(significant_families)))

# --- Observed vs expected, to see the size of the effect
p_xz = stage_coherence |>
  ggplot(aes(x = mean_X,y = X,colour = Z)) +
  geom_abline(slope = 1,intercept = 0,linetype = "dashed",colour = "grey50") +
  geom_point(alpha = 0.8,size = 1.6) +
  facet_wrap(~ substrate) +
  scale_colour_gradient2(low = "#1B6CA8",mid = "grey80",high = "#D1495B",midpoint = 0) +
  labs(title = "Effective number of stages of each family, observed vs randomized",
       x = "Expected (mean of the permutations)",y = "Observed") +
  theme_classic(base_size = 12)

ggsave(filename = "figures/time_based_clustering/stage_coherence_observed_vs_null.pdf",
       plot = p_xz,width = 10,height = 7)

