load_solutions = function(file.list){
  # This function takes a character vector in which each entry
  # is a path to a yaml file and it creates a single list with
  # all solutions. The expected structure (and the one generated)  
  # is a list: list[[medium]][[strain]][[solution]]. Note that 
  # the same key (e.g. medium) cannot be repeated in different files.
  # author: apascualgarcia.github.io
  # date: April 25th, 2022 (Berlin)
  #
  require(yaml)
  i=0
  for(fileIn in file.list){
    i=i+1
    sols.tmp=read_yaml(file=fileIn)
    if(i == 1){
      sols = sols.tmp
    }else{
      C.list=names(sols.tmp)
      for(C.tmp in C.list){
        N.C.tmp=length(sols[[C.tmp]]) # check if the C source already is in the list,
        if(N.C.tmp == 0){ #  if not
          sols[[C.tmp]]=sols.tmp[[C.tmp]] # create a new entry
        }else{ # if it is
          N.C.tmp=N.C.tmp+1 # add a new strain
          strain.list=names(sols.tmp[[C.tmp]])
          for(strain.tmp in strain.list){
            sols[[C.tmp]][[strain.tmp]]=sols.tmp[[C.tmp]][[strain.tmp]]
          }
        }
      }
    }
  }
  return(sols)
}



