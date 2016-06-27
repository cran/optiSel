
"kinwac"<-function(K, Pedig, use=NULL){
  rownames(Pedig) <- Pedig$Indiv
  cohort <- Pedig[rownames(K[[1]]),"Born"]
  if(is.null(use))use<-rep(TRUE,length(cohort))

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
    if(sum((cohort==k) & use)>0){
    i<-c(i,rep(which((cohort==k) & use),sum((cohort==k))))
    j<-c(j,rep(which((cohort==k)),each=sum((cohort==k) & use)))
    }
  }
  A<-sparseMatrix(i=i[i!=j],j=j[i!=j],dims=c(length(cohort),length(cohort)))
  A <- as(A, "dgCMatrix")
  kinWithPop <- matrix(NA, nrow=nrow(K[[1]]), ncol=length(K)+3)
  colnames(kinWithPop)<-c("cohort", "used", "I", names(K))
  rownames(kinWithPop)<-rownames(K[[1]])
  N<-apply(A,2,sum)
  cohort[cohort==-123456789]<- NA
  kinWithPop[,1]<-cohort
  kinWithPop[,2]<-1*use
  kinWithPop[,3]<-Pedig[rownames(K[[1]]),"I"]
  for(k in 1:length(K)){
    kinWithPop[,k+3]<-apply(K[[k]]*A,2,sum)/N
  }
  class(kinWithPop)<-"kinwac"
  kinWithPop
}