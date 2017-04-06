
"kinwac"<-function(K, Pedig){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  rownames(Pedig) <- Pedig$Indiv
  cohort <- Pedig[rownames(K[[1]]),"Born"]
 
  Indiv<-1; Sire<-2; Dam<-3;
  for(i in c(Indiv, Sire, Dam)){Pedig[,i]<-as.character(Pedig[,i])}
   
  x <- cbind(Pedig[Pedig[Pedig[,Sire],Dam],  "Born"],
             Pedig[Pedig[Pedig[,Sire],Sire], "Born"],
             Pedig[Pedig[Pedig[,Dam],Dam],   "Born"],
             Pedig[Pedig[Pedig[,Dam],Sire],  "Born"])
  Pedig$I <- (Pedig$Born - apply(x,1,mean,na.rm=TRUE))/2
  
  if(!is.list(K)){K<-list(K=K)}
  cohort[is.na(cohort)]<- -123456789
  Years <- unique(cohort)
  i<-NULL
  j<-NULL
  for(k in Years){
    if(sum(cohort==k)>0){
    i<-c(i, rep(which(cohort==k), sum(cohort==k)))
    j<-c(j, rep(which(cohort==k), each=sum(cohort==k)))
    }
  }
  A<-sparseMatrix(i=i[i!=j],j=j[i!=j],dims=c(length(cohort),length(cohort)))
  A <- as(A, "dgCMatrix")
  cohort[cohort==-123456789]<- NA
  
  kinWithPop <- matrix(NA, nrow=nrow(K[[1]]), ncol=length(K)+2)
  colnames(kinWithPop)<-c("cohort", "I", names(K))
  rownames(kinWithPop)<-rownames(K[[1]])
  N<-apply(A,2,sum)
  
  kinWithPop[,1]<-cohort
  kinWithPop[,2]<-Pedig[rownames(K[[1]]),"I"]
  for(k in 1:length(K)){
    kinWithPop[,k+2]<-apply(K[[k]]*A,2,sum)/N
  }
  kinWithPop <- data.frame(Indiv=rownames(K[[1]]), kinWithPop)
  if(PedigAsDataTable){setDT(kinWithPop)}
  class(kinWithPop)<-c("kinwac", class(kinWithPop))
  attributes(kinWithPop)$condProb <- attributes(K)$condProb 
  kinWithPop
}