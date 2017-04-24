
"kinlist"<-function(...){
  obj <- list(...)

  ### check and modify names of objects ###
  for(i in 1:length(obj)){
    if(!is.list(obj[[i]])){
      Name       <- names(obj)[i]
      if(is.null(Name)){stop(paste("A name must be specified for parameter ",i,".\n",sep=""))}
    }
    if(!is.null(attributes(obj[[i]])$condProb) & length(attributes(obj[[i]])$condProb)==1 & !is.null(names(obj)[i])){
      names(attributes(obj[[i]])$condProb)<-names(obj)[i]
    }
  }

  ### get condProb attribute ###
  cProb <- NULL
  for(i in 1:length(obj)){
    cProb <- append(cProb, attributes(obj[[i]])$condProb)
  }
  
  ### Flatten the list #####
  Res<-list()
  for(i in 1:length(obj)){
    if(!is.list(obj[[i]])){
      Res[[names(obj)[i]]]<-obj[[i]]
    }else{
      Res<-append(Res, obj[[i]])
    }
  }

  if(!is.null(cProb)){attributes(Res)$condProb <- cProb}
  for(k in names(Res)){Res[[k]] <- as(Res[[k]],"matrix")}
  class(Res)<-"kinMatrices"
  
  Res
}