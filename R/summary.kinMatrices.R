
"summary.kinMatrices"<-function(object, Pedig, tlim=NULL, histNe=NULL, base=NULL, df=4, ...){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)} 
  
  Kinwac <- kinwac(K=object, Pedig=Pedig)
  
  condProbMat <- c()
  for( i in names(attributes(Kinwac)$condProb)){
    condProbMat<-c(condProbMat, attributes(Kinwac)$condProb[[i]]["f1"], attributes(Kinwac)$condProb[[i]]["f2"])
  }
  condProb<- attributes(Kinwac)$condProb 
  
  Kinwac <- Kinwac[,-1]
  cohort <- Kinwac[,1]
  Param  <- aggregate(Kinwac[,-1],list(cohort),mean, na.rm=TRUE)
  colnames(Param)[1]<-"cohort"
  if(is.null(tlim)){tlim<-range(cohort)}
  Param   <- Param[Param$cohort>=tlim[1] & Param$cohort<=tlim[2],]
  if(!is.null(condProb)){
    i  <- names(condProb)[1]
    fZ <- condProb[[i]]["f1"]
    fN <- condProb[[i]]["f2"]
    if(fZ %in% colnames(Param) & fN %in% colnames(Param)){
      Param[[i]]=1-(1-Param[,fZ])/Param[,fN] 
      if(is.numeric(cohort)){
        Param$Ne<-getNe(Param$cohort,1-Param[[i]], df=df,I=Param$I)
        I <- mean(Param$I, na.rm=TRUE)
        if(is.null(base)){base<-round(tlim[1]-25*I,0)}
        if(is.null(histNe)){histNe<-round(3*mean(Param$Ne),0)}
        condGD  <- 1-Param[[i]]
        Param$condGD <- (1-1/(2*histNe))^((tlim[1]-base)/I) * condGD/condGD[1]
        Param$NGE    <- 1/(2*(1-Param$condGD))
      }
    }
  }
  Param <- Param[, setdiff(colnames(Param), condProbMat)]
  print(format(round(Param,3), digits=3, scientific=FALSE))
  if(PedigAsDataTable){setDT(Param)}
  attributes(Param)$histNe <- histNe
  attributes(Param)$base   <- base
  attributes(Param)$tlim   <- tlim
  attributes(Param)$df     <- df
  
  invisible(Param)
}