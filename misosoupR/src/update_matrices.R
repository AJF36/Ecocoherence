update_matrices = function(key.suppl,out.sp,medium,mat.list){
  # this function initializes and updates the 
  # medium Vs supplier, focal vs supplier
  # and supplier vs supplier matrices.
  if(key.suppl == 1){# initialize all matrices
    Nsup=length(out.sp$suppliers)
    mediaToSuppl=matrix(1,nrow=1,ncol=Nsup)
    rownames(mediaToSuppl)=medium
    colnames(mediaToSuppl)=out.sp$suppliers
    focalToSuppl=matrix(1,nrow=1,ncol=Nsup)
    rownames(focalToSuppl)=out.sp$focal
    colnames(focalToSuppl)=out.sp$suppliers
    supplToSuppl=matrix(1,nrow=Nsup,ncol=Nsup)
    rownames(supplToSuppl)=out.sp$suppliers
    colnames(supplToSuppl)=out.sp$suppliers
    #diag(supplToSuppl)=0
  }else{
    mediaToSuppl=mat.list$mediaToSuppl
    focalToSuppl=mat.list$focalToSuppl
    supplToSuppl=mat.list$supplToSuppl
    # check if the medium exists, otherwise create a new row
    id.med=match(medium,rownames(mediaToSuppl))
    if(is.na(id.med)){
      newrow=matrix(0,nrow=1,ncol=dim(mediaToSuppl)[2])
      rownames(newrow)=medium
      mediaToSuppl=rbind(mediaToSuppl,newrow)
    }
    # same for focal strains
    id.foc=match(out.sp$focal,rownames(focalToSuppl))
    if(is.na(id.foc)){
      newrow=matrix(0,nrow=1,ncol=dim(focalToSuppl)[2])
      rownames(newrow)=out.sp$focal
      focalToSuppl=rbind(focalToSuppl,newrow)
    }
    # finally suppliers, we work on the three matrices here
    matched=match(out.sp$suppliers,rownames(supplToSuppl)) # identify existing 
    exist.suppl=out.sp$suppliers[!is.na(matched)] 
    new.suppl=out.sp$suppliers[is.na(matched)] # and new suppliers
    all.suppl=c(exist.suppl,new.suppl)
    if(length(new.suppl) > 0){
      # ... extend matrices, cols for media and focal, both rows and cols for suppl
      # ..... media
      newcols=matrix(0,nrow=dim(mediaToSuppl)[1],ncol=length(new.suppl)) 
      colnames(newcols)=new.suppl
      mediaToSuppl=cbind(mediaToSuppl,newcols)
      # ..... focal
      newcols=matrix(0,nrow=dim(focalToSuppl)[1],ncol=length(new.suppl))
      colnames(newcols)=new.suppl
      focalToSuppl=cbind(focalToSuppl,newcols)
      # ..... suppliers
      newcols=matrix(0,nrow=dim(supplToSuppl)[1],ncol=length(new.suppl))
      colnames(newcols)=new.suppl
      supplToSuppl=cbind(supplToSuppl,newcols)
      newrows=matrix(0,ncol=dim(supplToSuppl)[2],nrow=length(new.suppl))
      rownames(newrows)=new.suppl
      supplToSuppl=rbind(supplToSuppl,newrows)
    }
    # We can finally safely update all objects
    mediaToSuppl[medium,all.suppl]=mediaToSuppl[medium,all.suppl]+1
    focalToSuppl[out.sp$focal,all.suppl]=focalToSuppl[out.sp$focal,all.suppl]+1
    supplToSuppl[all.suppl,all.suppl]=supplToSuppl[all.suppl,all.suppl]+1
  }
  return(mat.list=list("mediaToSuppl"=mediaToSuppl,
                       "focalToSuppl"=focalToSuppl,
                       "supplToSuppl"=supplToSuppl))
}