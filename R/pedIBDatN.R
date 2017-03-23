
# -----------------------------------------------------------
# --  The following objective functions can be used:       --
# --  fA  1-P(X_N != Y_N)                                  --
# --  fB  1-P(X_N != Y_N AND (X_N %in% F AND Y_N %in% F))  --
# --  fC  1-P(X_N != Y_N AND (X_N %in% F OR  Y_N %in% F))  --
# --  fD  1-P(X_N != Y_N | X_N %in% F AND Y_N %in% F)      --
# ----------------------------------------------------------- 

"pedIBDatN"<-function(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only, nGen=NA){
  getNe <- !is.na(nGen)
  if("data.table" %in% class(Pedig)){
    Pedig <- as.data.frame(Pedig)
    setDF(Pedig)
    }
  ids <- as.character(Pedig[[1]])
  
  if(is.logical(keep)){keep <- as.character(Pedig[keep,1])}
  if(!is.null(keep)){keep <- as.character(keep)}
  if(is.na(thisBreed)){stop("The name of this breed is not specified.\n")}
  if(!("Breed" %in% colnames(Pedig))){stop("Column breed is missing.\n")}
  
  Pedig <- prePed(Pedig, keep=keep, lastNative=1234567, thisBreed=thisBreed)
  hasWrongBreed <- (!is.na(Pedig$Sire) & !is.na(Pedig$Dam) & !(Pedig$Breed %in% thisBreed) & ((Pedig[Pedig$Sire, "Breed"] %in% thisBreed) | (Pedig[Pedig$Dam, "Breed"] %in% thisBreed)) )
  Pedig[hasWrongBreed, "Breed"] <- thisBreed
  Pedig[Pedig$Breed != thisBreed, "Sire"]<-NA
  Pedig[Pedig$Breed != thisBreed, "Dam"]<-NA
  Pedig <- prePed(Pedig, keep=keep)
  
  if(!is.null(keep)){
    keep <- Pedig$Indiv[Pedig$Indiv %in% keep]
    }
  
  if(is.null(keep.only)){
    keep.only <- ids
  }else{
    keep.only <- as.character(keep.only)
    keep.only <- ids[ids %in% keep.only]
  }
  if(getNe){
    Selection <- Pedig$Indiv
  }else{
    Selection <- keep.only
  }
  condProb <- list()
  condProb$pedIBDatN <- c(f1="pedZ",f2="pedN")
  
  Rassen     <- setdiff(names(table(Pedig$Breed)), c(thisBreed))
  MigFounder <- Pedig$Indiv[(is.na(Pedig$Sire)|is.na(Pedig$Dam)) &   Pedig$Breed %in% Rassen]
  NatFounder <- Pedig$Indiv[(is.na(Pedig$Sire)|is.na(Pedig$Dam)) & !(Pedig$Breed %in% Rassen)]
  nMig <- length(MigFounder)
  nNat <- length(NatFounder)
  AMig <- matrix(2, nMig, nMig, dimnames=list(MigFounder, MigFounder))
  ANat <- matrix(2, nNat, nNat, dimnames=list(NatFounder, NatFounder))
  cat(paste0("Number of Migrant Founders: ", nrow(AMig), "\n"))
  cat(paste0("Number of Native  Founders: ", nrow(ANat), "\n"))
  cat(paste0("Individuals in Pedigree   : ", nrow(Pedig), "\n"))
  GB <- ((nrow(AMig)+nrow(ANat))^2 + length(Selection)^2 + nrow(Pedig)^2)*(7.45058066987776e-09)*1.1
  if(GB>1){cat(paste0("Ensure that you have more than ", round(GB, 1), " GB memory available.\n"))}
  if(GB>1){cat("Computing fOI ...")}
  fOI  <- 0.5*makeA(Pedig[,1:3], AFounder=AMig)[Selection, Selection]
  gc()
  if(GB>1){cat("finished\nComputing fII ...")}
  fII  <- 0.5*makeA(Pedig[,1:3], AFounder=adiag(AMig, ANat))[Selection, Selection]
  if(GB>1){cat("finished\nCombining results ...")}
  rm(AMig)
  rm(ANat)
  gc()
  
  Res<-list()
  Res$pedZ <- fOI + 1 - fII
  rm(fOI)
  
  Cont <- pedBreedComp(Pedig, thisBreed=thisBreed)
  Cont <- 1 - Cont[Selection, "native"]
  Res$pedN <- 1 - 0.5*(matrix(Cont, nrow=nrow(fII), ncol=ncol(fII), byrow=TRUE) +  matrix(Cont, nrow=nrow(fII), ncol=ncol(fII), byrow=FALSE)) - 0.5*(1-fII)
  dimnames(Res$pedN) <- list(Selection, Selection)
  rm(fII)
  
  Res$pedIBDandN <- Res$pedZ + Res$pedN - 1
  dimnames(Res$pedIBDandN) <- dimnames(Res$pedN)
  if(GB>1){cat("finished\n")}
  if(!is.null(keep) & getNe){
    nativeNe <- round(nativeNe(Pedig=Pedig, Kin=Res, keep=keep, thisBreed=thisBreed, nGen=nGen),1)
    cat("Native Ne = ", nativeNe, " (estimated from ", nGen, " previous generations)\n", sep="")
    attr(Res,"nativeNe") <- nativeNe
  }
  if(getNe){
    for(i in names(Res)){
      Res[[i]]<-Res[[i]][keep.only, keep.only]
    }
  }
  
  attr(Res,"meanpedIBDatN") <- 1-(1-mean(Res$pedZ))/mean(Res$pedN)
  if(!is.null(keep)){
    cat("Mean kinship at native alleles: ", round(attr(Res,"meanpedIBDatN"), 4), "\n")
  }
  class(Res) <- "kinMatrices"
  attr(Res,"condProb") <- condProb
  
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


