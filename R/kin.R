
# -----------------------------------------------------------
# --  The following objective functions can be used:       --
# --  fA  1-P(X_N != Y_N)                                  --
# --  fB  1-P(X_N != Y_N AND (X_N %in% F AND Y_N %in% F))  --
# --  fC  1-P(X_N != Y_N AND (X_N %in% F OR  Y_N %in% F))  --
# --  fD  1-P(X_N != Y_N | X_N %in% F AND Y_N %in% F)      --
# ----------------------------------------------------------- 

"kin"<-function(Pedig, thisBreed=NA,  method=c('A','B','C','D'), keep.only=NULL, keep=keep.only, nGen=6){
  if(!is.null(keep)){keep <- setdiff(keep, c(NA))}
  if(!is.null(keep.only)){keep.only <- setdiff(keep.only, c(NA))}
  if(!is.null(keep)){Pedig <- nadiv::prunePed(prePed(Pedig), phenotyped=keep)}
  condProb <- list()
  if("D" %in% method){method <- union(method, c("B", "N"))}
  if("B" %in% method & "N" %in% method){condProb$fD=c(f1="fB",f2="fN")}
  
  Indiv<-1; Sire<-2; Dam<-3; Sex<-4; Breed<-5;
  if(is.na(thisBreed)){thisBreed<-"unkown"}
  if(length(colnames(Pedig))==4){Pedig$Breed<-"unkown"}
  
  for(i in c(Indiv, Sire, Dam, Breed)){Pedig[,i]<-as.character(Pedig[,i])}

  Pedig2<-Pedig

  Rassen<-setdiff(names(table(Pedig[,Breed])),c(thisBreed))
  Selfing <- data.frame(
    paste(rep(c('Founder','Migrant'),20), rep(20:1,each=2),sep=''),
    c(NA, NA, paste(rep(c('Founder','Migrant'),19),rep(20:2,each=2),sep='')),
    c(NA, NA, paste(rep(c('Founder','Migrant'),19),rep(20:2,each=2),sep='')),
    0, "Dummy", stringsAsFactors=FALSE)
  colnames(Selfing)<-colnames(Pedig)[1:5]
  
  Pedig<-rbind.data.frame(Selfing, Pedig[,1:5])

  #Pedig[Pedig[,Breed] %in% Rassen,Sire]<-NA
  #Pedig[Pedig[,Breed] %in% Rassen, Dam]<-NA
  fOO <- 0.5*nadiv::makeA(Pedig[-(1:40),c(Indiv,Sire,Dam)])
  dimnames(fOO)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  Pedig[Pedig[,Breed] %in% Rassen, Sire] <- 'Migrant1'
  Pedig[Pedig[,Breed] %in% Rassen,  Dam] <- 'Migrant1'
  suppressWarnings(fOI <- 0.5*nadiv::makeA(Pedig[,c(Indiv,Sire,Dam)])[- (1:40),- (1:40)])
  dimnames(fOI)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  Pedig[is.na(Pedig[,Sire])| Pedig[,Sire]=="0",Sire] <- 'Founder1'
  Pedig[is.na(Pedig[,Dam]) | Pedig[,Dam]=="0",  Dam] <- 'Founder1'
  Pedig[1:2, c(Sire, Dam)] <- NA
  suppressWarnings(fII <- 0.5*nadiv::makeA(Pedig[,c(Indiv,Sire,Dam)])[- (1:40),- (1:40)])
  dimnames(fII)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  
  Res<-list()
  if("A" %in% method){Res$fA<- fOO}
  if("B" %in% method|"D" %in% method){Res$fB<- fOI + matrix(1,nrow=nrow(fII),ncol=ncol(fII)) - fII}
  if("C" %in% method){Res$fC<- fOI} 
  if("N" %in% method){
    Pedig[is.na(Pedig[,2]),2]<-"0"
    Pedig[is.na(Pedig[,3]),3]<-"0"
    Cont <- genecont(Pedig[,1], Pedig[,2], Pedig[,3], NAncestors=40)[- (1:40),"Migrant1"]
    Res$fN<- 1 - 0.5*(matrix(Cont, nrow=nrow(fII),ncol=ncol(fII),byrow=TRUE) +  matrix(Cont,nrow=nrow(fII),ncol=ncol(fII),byrow=FALSE)) - 0.5*(1-fII)
    dimnames(Res$fN)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  }
  
  if("B" %in% method & "N" %in% method & !is.null(keep)){
    nativeNe <- round(nativeNe(Pedig=Pedig2, Kin=Res, keep=keep, thisBreed=thisBreed, nGen=nGen),1)
    cat("Native Ne = ", nativeNe," (estimated from ",nGen," previous generations)\n", sep="")
    attr(Res,"nativeNe") <- nativeNe
  }
  if(!is.null(keep.only)){for(i in names(Res))Res[[i]]<-Res[[i]][keep.only, keep.only]}
    
  class(Res)<-"kinMatrices"
  attr(Res,"condProb")<- condProb
  Res
}



nativeNe <- function(Pedig, Kin, keep, thisBreed, nGen=3){
  U  <- upper.tri(Kin$fN[keep, keep])
  fN <- mean(Kin$fN[keep, keep][U])
  fB <- mean(Kin$fB[keep, keep][U])
  fD <- 1-(1-fB)/fN
  ID <- keep
#ID <- setdiff(unlist(Pedig[ID,c("Sire","Dam")]),c(NA,"0"))
 for(i in 1:nGen){ID <- setdiff(unlist(Pedig[ID,c("Sire","Dam")]),c(NA,"0"))}
#  for(i in 2:nGen){ID <- unlist(Pedig[ID,c("Sire","Dam")]); ID<-ID[!is.na(ID) & Pedig[ID,5]==thisBreed]}
 ID <- intersect(Pedig[Pedig[,5]==thisBreed,1],ID)
  U <-  upper.tri(Kin$fN[ID, ID])
  fN <- mean(Kin$fN[ID, ID][U])
  fB <- mean(Kin$fB[ID, ID][U])
  fD0<- 1-(1-fB)/fN
  1/(2*(1-((1-fD)/(1-fD0))^(1/nGen)))
}


