pak::pak("admarhi/ramen")

ms_raw <- yaml::read_yaml("../ramen-data/misosoup/misosoup_202412.yaml")

ms_imported <- ramen::importMisosoup(ms_raw)


ms_tb <- ms_imported$consortia

cm_test <- ms_tb |> 
  dplyr::filter(cons_id == "3pg_min_1") |> 
  ConsortiumMetabolism(name = "test 1")

plot(cm_test)
