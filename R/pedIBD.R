
"pedIBD"<-function(Pedig,  keep.only=NULL, keep=keep.only){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  
  if(is.logical(keep)){keep<-Pedig[keep,1]}
  if(!is.null(keep)){keep<-as.character(keep); keep <- setdiff(keep, c(NA, "", " ", "0"))}
  if(!is.null(keep.only)){keep.only<-as.character(keep.only); keep.only <- setdiff(keep.only, c(NA, "", " ", "0"))}
  if(!is.null(keep)){Pedig <- nadiv::prunePed(prePed(Pedig), phenotyped=keep)}
  Indiv<-1; Sire<-2; Dam<-3;
  for(i in c(Indiv, Sire, Dam)){Pedig[,i]<-as.character(Pedig[,i])}
  fA <- 0.5*nadiv::makeA(Pedig[,c(Indiv,Sire,Dam)])
  dimnames(fA)<-list(Pedig[,Indiv], Pedig[,Indiv])
  if(!is.null(keep.only)){fA<-fA[keep.only, keep.only]}
  fA
}






