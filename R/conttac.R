
"conttac"<-function(cont, cohort, use=rep(TRUE,length(cohort)), mincont=0.05){
  cont<-cont[use==1,]
  cohort<-cohort[use==1]
  ContByYear<-stats::aggregate(cont,by=list(cohort),mean)
  rownames(ContByYear)<-ContByYear$Group.1
  ContByYear<-ContByYear[,-1]
  Migrant<-rowSums(ContByYear)
  ContByYear<-ContByYear[,colMeans(ContByYear)>=mincont]
  if(!("other" %in% colnames(ContByYear)))ContByYear$other<-0
  ContByYear$other<-ContByYear$other+Migrant-rowSums(ContByYear)
  if("unknown"%in% colnames(ContByYear)){
    ContByYear<-ContByYear[,c(setdiff(colnames(ContByYear), c("unknown")),"unknown")]
    }
  
  t(ContByYear)
}

