pak::pak("admarhi/ramen")

library(ramen)
 
test1 <- ms_tb |>
  dplyr::filter(cons_id == "3pg_min_1")

cm_test <- ConsortiumMetabolism(
  data = test1,
  name = "My first object",
  split_by = "cons_id"
)
?ConsortiumMetabolism
plot(cm_test, type = "EffectiveProduction")
x <- getEdges(cm_test)
getSpecies(cm_test)
ms_list <- ms_tb[1:10000, ] |>
  dplyr::group_split(.data$cons_id) |>
  purrr::map(
    .f = \(x) ConsortiumMetabolism(
      data = x,
      name = unique(x$cons_id)
    )
  )
cms_test <- ConsortiumMetabolismSet(
  ms_list,
  name = "test set",
  desc = "this is a super cool test cms with 15 consortia"
)
plot(cms_test)
cms_clust_6 <- cms_test |>
  getCluster(node_id = 6)
getSpecies(cms_clust_6, type = "generalists")
getSpecies(cms_clust_6, type = "specialists")
getFunctionalGroups(cms_clust_6, k = 8)
getEdges(cms_test, type = "pan-cons")
getEdges(cms_test, type = "niche")
getEdges(cms_test, type = "core")
getEdges(cms_test, type = "aux")




















