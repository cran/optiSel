
# -----------------------------------------------------------
# --  The following objective functions can be used:       --
# --  fA  1-P(X_N != Y_N)                                  --
# --  fB  1-P(X_N != Y_N AND (X_N %in% F AND Y_N %in% F))  --
# --  fC  1-P(X_N != Y_N AND (X_N %in% F OR  Y_N %in% F))  --
# --  fD  1-P(X_N != Y_N | X_N %in% F AND Y_N %in% F)      --
# ----------------------------------------------------------- 

"pedIBDatN"<-function(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only, nGen=6){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  
  if(is.logical(keep)){keep<-Pedig[keep,1]}
  if(!is.null(keep)){ keep<-as.character(keep); keep <- setdiff(keep, c(NA, "", " ", "0"))}
  if(!is.null(keep.only)){keep.only<-as.character(keep.only); keep.only <- setdiff(keep.only, c(NA, "", " ", "0"))}
  if(!is.null(keep)){Pedig <- nadiv::prunePed(prePed(Pedig), phenotyped=keep)}
  condProb <- list()
  condProb$pedIBDofN <- c(f1="pedZ",f2="pedN")
  
  Indiv<-1; Sire<-2; Dam<-3; Sex<-4; Breed<-5;
  if(is.na(thisBreed)){stop("The name of this breed is not specified.\n")}
  if(length(colnames(Pedig))==4){stop("Column breed is missing.\n")}
  
  for(i in c(Indiv, Sire, Dam, Breed)){Pedig[,i]<-as.character(Pedig[,i])}

  Pedig2  <- Pedig
  Rassen  <- setdiff(names(table(Pedig[,Breed])),c(thisBreed))
  Selfing <- data.frame(
    paste(rep(c('Founder','Migrant'),20), rep(20:1,each=2),sep=''),
    c(NA, NA, paste(rep(c('Founder','Migrant'),19),rep(20:2,each=2),sep='')),
    c(NA, NA, paste(rep(c('Founder','Migrant'),19),rep(20:2,each=2),sep='')),
    0, "Dummy", stringsAsFactors=FALSE)
  colnames(Selfing)<-colnames(Pedig)[1:5]
  
  Pedig <- rbind.data.frame(Selfing, Pedig[,1:5])
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

  Res$pedZ <- as(fOI + matrix(1,nrow=nrow(fII),ncol=ncol(fII)) - fII, "matrix")

  Pedig[is.na(Pedig[,2]),2]<-"0"
  Pedig[is.na(Pedig[,3]),3]<-"0"
  Cont   <- genecont(Pedig[,1], Pedig[,2], Pedig[,3], NAncestors=40)[- (1:40),"Migrant1"]
  Res$pedN <- as(1 - 0.5*(matrix(Cont, nrow=nrow(fII),ncol=ncol(fII),byrow=TRUE) +  matrix(Cont,nrow=nrow(fII),ncol=ncol(fII),byrow=FALSE)) - 0.5*(1-fII),"matrix")
  dimnames(Res$pedN)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  
  if(!is.null(keep)){
    nativeNe <- round(nativeNe(Pedig=Pedig2, Kin=Res, keep=keep, thisBreed=thisBreed, nGen=nGen),1)
    cat("Native Ne = ", nativeNe," (estimated from ",nGen," previous generations)\n", sep="")
    attr(Res,"nativeNe") <- nativeNe
  }
  if(!is.null(keep.only)){for(i in names(Res))Res[[i]]<-Res[[i]][keep.only, keep.only]}
  
  attr(Res,"meanpedIBDatN") <- 1-(1-mean(Res$pedZ))/mean(Res$pedN)
  if(!is.null(keep)){
  cat("Mean kinship at native alleles: ", round(attr(Res,"meanpedIBDatN"), 4), "\n")
  }
  class(Res)<-"kinMatrices"
  attr(Res,"condProb")<- condProb
  Res
}





nativeNe <- function(Pedig, Kin, keep, thisBreed, nGen=3){
  U  <- upper.tri(Kin$pedN[keep, keep])
  fN <- mean(Kin$pedN[keep, keep][U])
  fB <- mean(Kin$pedZ[keep, keep][U])
  fD <- 1-(1-fB)/fN
  ID <- keep
  #ID <- setdiff(unlist(Pedig[ID,c("Sire","Dam")]),c(NA,"0"))
  for(i in 1:nGen){ID <- setdiff(unlist(Pedig[ID,c("Sire","Dam")]),c(NA,"0"))}
  #  for(i in 2:nGen){ID <- unlist(Pedig[ID,c("Sire","Dam")]); ID<-ID[!is.na(ID) & Pedig[ID,5]==thisBreed]}
  ID <- intersect(Pedig[Pedig[,5]==thisBreed,1],ID)
  U <-  upper.tri(Kin$pedN[ID, ID])
  fN <- mean(Kin$pedN[ID, ID][U])
  fB <- mean(Kin$pedZ[ID, ID][U])
  fD0<- 1-(1-fB)/fN
  1/(2*(1-((1-fD)/(1-fD0))^(1/nGen)))
}


