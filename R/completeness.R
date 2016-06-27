
"completeness"<-function(Pedig, keep, genNo=-1, type="MEAN"){
  keep <-  setdiff(keep, c(NA))
  Pedig <- prePed(Pedig, keep=keep, addNum=TRUE)
  Ped   <- Pedig[,!(colnames(Pedig)%in% c("numIndiv", "numSire", "numDam"))]
  numKeep  <- Pedig[keep, "numIndiv"]
  Pedig <- Pedig[,c("numIndiv", "numSire", "numDam", "Sex")]
  colnames(Pedig)<-c("ind","father","mother","sex")
  Pedig<-GENLIB::gen.genealogy(Pedig)
  Res<-GENLIB::gen.completeness(Pedig, pro=numKeep, genNo=genNo, type=type)
  if(type=="IND"){
    names(keep)<-numKeep
    Indiv<-keep[str_sub(colnames(Res),5,-1)]
    Res <- t(Res)
    colnames(Res)<-paste("Compl.",0:(ncol(Res)-1), sep="")
    Res<-data.frame(Res, Ped[Indiv,])
    rownames(Res)<-1:nrow(Res)
  }
  Res
}