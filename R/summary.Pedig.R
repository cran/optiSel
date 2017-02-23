
"summary.Pedig"<-function(object, keep=NULL, maxd=50, d=4, ...){
  if(is.logical(keep)){keep <- object[keep,1]}
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
  }
  Pedig <- prePed(object, keep=keep, addNum=TRUE)
  Pedig$Inbreeding <- pedInbreeding(Pedig)$Inbr
  compl <- completeness(Pedig, maxd=maxd, by="Indiv")
  setDT(compl)
  setDT(Pedig)
  Completeness <- NULL
  equiGen      <- NULL
  fullGen      <- NULL
  Indiv        <- NULL
  maxGen       <- NULL
  meanCom      <- NULL
  numSire      <- NULL
  numDam       <- NULL
  Inbreeding   <- NULL
  Res   <- compl[, list(maxGen=length(Completeness)-1, equiGen=sum(Completeness)-1), by="Indiv"]
  sRes  <- compl[compl$Completeness==1, list(fullGen=length(Completeness)-1), by="Indiv"]
  Res   <- sRes[Res, list(Indiv, equiGen, fullGen, maxGen), on="Indiv"]
  sRes  <- compl[compl$Generation %in% (0:(d-1)), list(meanCom=sum(Completeness)/(1*d)), by="Indiv"]
  Res   <- sRes[Res, list(Indiv, equiGen, fullGen, maxGen, meanCom), on="Indiv"]
  Res$meanCom[is.na(Res$meanCom)]<-0
  Pedig <- Res[Pedig,list(Indiv, numSire, numDam, equiGen, fullGen, maxGen, meanCom, Inbreeding), on="Indiv"]
  setDF(Pedig)
  Pedig$patIndex <- 0
  Pedig[Pedig$numSire!=0,"patIndex"]<-Pedig[Pedig$numSire[Pedig$numSire!=0],"meanCom"] 
  Pedig$matIndex <- 0
  Pedig[Pedig$numDam!=0,"matIndex"]<-Pedig[Pedig$numDam[Pedig$numDam!=0],"meanCom"] 
  Pedig$PCI <- 2*Pedig$matIndex*Pedig$patIndex/(Pedig$matIndex+Pedig$patIndex)
  Pedig[is.na(Pedig$PCI),"PCI"]<-0
  Pedig <- Pedig[,c("Indiv","equiGen", "fullGen", "maxGen", "PCI","Inbreeding")]
  if(!is.null(keep)){
    Pedig<-Pedig[Pedig$Indiv %in% keep,]
  }
  Pedig
}