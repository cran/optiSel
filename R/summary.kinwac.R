
"summary.kinwac"<-function(object, tlim=NULL, histNe=NULL, base=NULL, df=4, ...){
  object<-object[object[,"used"]==1,-2]
  cohort  <- object[,1]
  Param <- aggregate(object[,-1],list(cohort),mean, na.rm=TRUE)
  colnames(Param)[1]<-"cohort"
  if(is.null(tlim)){tlim<-range(cohort)}
  Param<-Param[Param$cohort>=tlim[1] & Param$cohort<=tlim[2],]
  if("fB" %in% colnames(Param) & "fN" %in% colnames(Param)){
    Param$fD=1-(1-Param[,"fB"])/Param[,"fN"] 
    if(is.numeric(cohort)){
      Param$Ne<-getNe(Param$cohort,1-Param$fD, df=df,I=Param$I)
      I <- mean(Param$I, na.rm=TRUE)
      if(is.null(base)){base<-round(tlim[1]-25*I,0)}
      if(is.null(histNe)){histNe<-round(3*mean(Param$Ne),0)}
      condGD  <- 1-Param$fD
      Param$condGD <- (1-1/(2*histNe))^((tlim[1]-base)/I) * condGD/condGD[1]
      Param$NGE    <- 1/(2*(1-Param$condGD))
    }
  }
  KinNames<-setdiff(colnames(Param),c("condGD", "NGE", "Ne", "I" , "cohort"))
  DivNames<-paste("Div.", KinNames, sep="")
  Param[, DivNames] <- 1-Param[, KinNames]
  
  print(format(round(Param,3), digits=3, scientific=FALSE))
    
  invisible(list(Param=Param,histNe=histNe, base=base, tlim=tlim, df=df))
}