stat_comparison_sets_solutions = function(df.list,addLabel=NULL,
                                          matched.sets=FALSE){
  # This function considers a list containing at least two data frames,
  # each of them obtained with the function "summary_stats_solutions", and it
  # creates a box plot for each of the variables contained in these
  # dataframes to compare their distribution. It also creates a table with
  # the quantiles of each variable and dataset.
  # input: df.list = list with data frames.
  #        addLabel = an additional label for the output directory (defaults to NULL)
  #        matched = logical, indicates if solutions across data.framed are matched, so
  #           the first solution in df1 match first solution in df2, df3, etc. It
  #           is required to generate scatter plots
  # output: A new directory with the box plots and the table of quantiles.
  #######
  # author: apascualgarcia.github.io
  # date: September 7th, 2022 (Lausanne)
  #######
  if(is.null(addLabel)){
    dirOut="stat_comparison_sets_sols"
  }else{
    dirOut=paste("stat_comparison_sets_sols",addLabel,sep="_")
  }
  dir.create(dirOut)
  setwd(dirOut)
  var.vec=colnames(df.list[[1]]) # all dfs have the same colnames
  q=c(0,0.05,0.25,0.5,0.75,0.95,1) # quantiles collected
  q_names=paste0("q_",q)
  key1=0
  for(var in var.vec){
    key2=0
    #if(var == "N_eff_sp"){browser()} # debug
    for(j in 1:length(df.list)){
      vec=df.list[[j]][,var] # extract the vector of values
      if((class(vec) != "numeric")&
         (class(vec) != "integer")){next}
      vec=vec[!is.na(vec)]
      q.vec=quantile(vec,probs = q,na.rm = TRUE) # calculate quantiles
      names(q.vec)=q_names # include colnames
      q.df.tmp=data.frame(t(q.vec)) # convert to df
      row.name.tmp=paste(var,labelIn[j],sep="_") # create row name
      row.names(q.df.tmp)=row.name.tmp
      vals.df.tmp=data.frame(vec,labelIn[j])
      colnames(vals.df.tmp)=c(var,"optimization")
      if(key1 == 0){ # include quantiles in a dataframe
        q.df=q.df.tmp
        key1=1
      }else{
        q.df=rbind(q.df,q.df.tmp)
      }
      if(key2 == 0){ # initialize df in long format
        vals.df0=vals.df.tmp
        vals.df=vals.df0
        colnames(vals.df0)=c(paste0(var,"_no_opt"),"no_opt")
        key2=1
        key3=0 # prepare for wide format (matched == TRUE)
      }else{
        #browser()
        vals.df=rbind(vals.df,vals.df.tmp)
        colnames(vals.df.tmp)=c(paste0(var,"_opt"),"optimization")
        if(matched.sets == TRUE){
          vals.df.wide.tmp=cbind(vals.df0,vals.df.tmp)
          if(key3 == 0){ # initialize wide format
            vals.df.wide=vals.df.wide.tmp
            key3=1
          }else{
            vals.df.wide=rbind(vals.df.wide,vals.df.wide.tmp)
          }
        }
      }
      vals.df[,2]=as.factor(vals.df[,2])
      lev2ord=levels(vals.df[,2]) # reorder levels
      matched=match(labelIn,lev2ord) # according to the vector or directories
      vals.df[,2]=factor(vals.df[,2],levels=levels(vals.df[,2])[matched])
    } # end df.list
    if(key2 == 1){
      # --- Plot bar plot
      file.plot=paste("BarPlot_",var,".pdf")
      pdf(file=file.plot,width = 10,height = 9)
      gg=ggplot(vals.df,aes_string("optimization",y=var))+
        geom_boxplot(aes(colour=optimization))+
        theme_bw()+
        theme(axis.title = element_text(size=28),
              axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
              axis.text = element_text(size=18),
              legend.position = "none")
      
      print(gg)
      dev.off()
      if(matched.sets == TRUE){
        # --- Plot scatter plot
        file.plot=paste("ScatterPlot_",var,".pdf")
        pdf(file=file.plot,width = 12,height = 9)
        xvar=colnames(vals.df.wide)[1]
        yvar=colnames(vals.df.wide)[3]
        # .... determine if log scale should be used
        logscale=FALSE
        max_val=max(c(vals.df.wide[, 1],vals.df.wide[, 3]))
        min_val=min(c(vals.df.wide[, 1],vals.df.wide[, 3]))
        if((max_val - min_val) > 100){logscale=TRUE}
        gg=ggplot(vals.df.wide,aes_string(x=xvar,y=yvar))+
          geom_point(aes(colour=optimization))+
          geom_abline(slope=1, intercept= 0,show.legend = TRUE)+
          theme_bw()+
          theme(axis.title = element_text(size=28),
                axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                axis.text = element_text(size=18),
                legend.text = element_text(size=15),
                legend.title= element_text(size=20))
        if(logscale == TRUE){
         gg=gg+scale_y_continuous(trans = "log10")+
           scale_x_continuous(trans = "log10") 
        }
        print(gg)
        dev.off()
      }
    }
  } # end var
  
  file.table=paste("Table_quantiles_optim-comparison.tsv")
  write.table(q.df,file=file.table,sep="\t",quote=FALSE)
}