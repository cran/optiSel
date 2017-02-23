
"completeness"<-function(Pedig, keep=NULL, maxd=50, by="Indiv"){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  
  if(is.logical(keep)){keep <- Pedig[keep,1]}
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
  }
  
  Pedig <- prePed(Pedig, keep=keep, addNum=TRUE)
  compl <- rcpp_completeness(as.character(Pedig$Indiv), as.integer(Pedig$numSire), as.integer(Pedig$numDam), as.integer(maxd))
  if(!is.null(keep)){compl<-compl[compl$Indiv %in% keep,]}
  
  if(by=="Indiv"){
    if(PedigAsDataTable){setDT(compl)}
    return(compl)
  }
  
  compl <- merge(compl, Pedig[, c("Indiv", "Sex")], by="Indiv")
  Factors <- list(compl[[by]], compl$Generation)
  names(Factors) <- c(by, "Generation")
  x <- aggregate(compl[,"Completeness",drop=FALSE], Factors, sum)
  x$Completeness <- x$Completeness/mapvalues(x$Sex, from=x[x$Generation==0, by], to=x[x$Generation==0, "Completeness"])
  x[[by]]   <- as.character(x[[by]])
  if(PedigAsDataTable){setDT(x)}
  x
}