setup_directory = function(vec.strains,vec.media,vec.sols,vec.suppliers){
  # This function creates a new directory with a name that combines
  # the focal strains, media and/or suppliers desired in a specific analysis.
  # It finally set the directory as the working directory and returns
  # the name, that can be used downstream as a label. Note that if
  # these lists are too long the operating system may have problems.
  spec_strains=length(vec.strains); spec_media=length(vec.media); 
  spec_suppliers=length(vec.suppliers); spec_sols=length(vec.sols)
  dirSpec="spec"
  if(spec_strains > 0){
    dirTmp=paste(vec.strains,collapse="-")
    dirTmp=paste("focal",dirTmp,sep="-")
    dirSpec=paste(dirSpec,dirTmp,sep="_")
  }
  if(spec_media > 0){
    dirTmp=paste(vec.media,collapse="-")
    dirTmp=paste("media",dirTmp,sep="-")
    dirSpec=paste(dirSpec,dirTmp,sep="_")
  }
  if(spec_suppliers > 0){
    dirTmp=paste(vec.suppliers,collapse="-")
    dirTmp=paste("suppliers",dirTmp,sep="-")
    dirSpec=paste(dirSpec,dirTmp,sep="_")
  }
  if(spec_sols > 0){
    dirTmp=paste(vec.sols,collapse="-")
    dirTmp=paste("sols",dirTmp,sep="-")
    dirSpec=paste(dirSpec,dirTmp,sep="_")
  }
  dir.create(dirSpec)
  setwd(dirSpec)
  return(dirSpec)
}