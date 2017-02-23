


"opticont4mb"<-function(method, K, phen, bc, thisBreed=names(bc)[1], con=list(), solver="cccp", quiet=FALSE, make.definite=solver=="csdp", ...){
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  
  phen$Breed <- as.character(phen$Breed)
  phen$Indiv <- as.character(phen$Indiv)
  rownames(phen)<-phen$Indiv
    
  BreedsInAnalysis <- names(table(phen$Breed))
  bc[setdiff(BreedsInAnalysis, names(bc))] <- 0
  bc[bc<0.001]<-0.001
  bc <- bc/sum(bc)
  if(!all(names(bc)%in%BreedsInAnalysis)){
    cat("Error: For some breeds included in the contribution vector bc no phenotypes or genotypes are provided.\n")
    return(NULL)
  }
  ####################################################################################
  ###  Create vector 'condProbMat' containing the names of all matrices needed     ###
  ###  for computing conditional kinships                                          ###
  ####################################################################################
  condProbMat <- c()
  for( i in names(attributes(K)$condProb)){
    condProbMat<-c(condProbMat, attributes(K)$condProb[[i]]["f1"], attributes(K)$condProb[[i]]["f2"])
  }
  condProb<-attributes(K)$condProb
  
  ###############################################################################
  ###  Enlarge matrices needed for computing conditional kinships             ###
  ###  by a diagonal matrix                                                   ###
  ###############################################################################
  indiv <- phen$Indiv[phen$Breed==thisBreed]
  for(i in condProbMat){
    G <- diag(1,nrow=nrow(phen),ncol=nrow(phen))
    rownames(G) <- phen$Indiv
    colnames(G) <- phen$Indiv
    G[indiv,indiv] <- K[[i]][indiv,indiv]
    K[[i]] <- G
  } 

  ###############################################################################
  ### Double and rename matrices suitable for computing across breed kinships ###  
  ###############################################################################
  Names <- setdiff(names(K), condProbMat)
  for(i in Names){
    BreedsKi <- names(table(phen[rownames(K[[i]]),"Breed"]))
    if(length(BreedsKi)>1){
      newi <- paste(i,".acrossBreeds",sep="")
      K[[newi]] <- K[[i]][phen$Indiv, phen$Indiv]
    }
  }

  
  #####################################################################################
  ###  - For within breed kinships the corresponding matrix is set equal            ###
  ###    to an identity matrix except for individuals not from the breed.           ###
  ###  - Create vector 'wBreedKin' containing names of all kinship matrices for     ###
  ###    which the mean  kinship within breed can be computed in the summary.       ###  
  #####################################################################################

  Names     <- setdiff(names(K), condProbMat)
  wBreedKin <- Names[!str_detect(Names, "acrossBreeds")]
  for(i in wBreedKin){
    G <- diag(1,nrow=nrow(phen),ncol=nrow(phen))
    rownames(G) <- phen$Indiv
    colnames(G) <- phen$Indiv
    G[indiv,indiv] <- K[[i]][indiv,indiv]
    K[[i]] <- G
  }
  
  ####################################################################################
  ###  Create vector 'cwBreedKin' containing names of all kinship matrices for     ###
  ###             which the mean kinship within breed should be constrained.       ###
  ####################################################################################
  cwBreedKin <- wBreedKin[paste("ub.", wBreedKin,sep="") %in% names(con)]
  
  ##########################################################################
  ### Use simulated sexes for breeds whose contributions                 ###
  ###  are not optimized                                                 ### 
  ##########################################################################
  Sex   <- 1+(1:nrow(phen))%%2
  phen$Sex[phen$Breed!=thisBreed] <- Sex[phen$Breed!=thisBreed]
  
  ##########################################################################
  ###  set phenotypes from other breeds equal to 0                       ###
  ##########################################################################
  Traits  <- colnames(phen)[-1]
  x<-NULL;for(i in 2:ncol(phen)){x<-c(x,is.numeric(phen[,i]))}
  Traits  <- setdiff(Traits[x], c("lb", "ub","oc", "Sex"))
  phen[phen$Breed!=thisBreed, Traits] <- 0
  
  
  ###########################################################################
  ### Define ub and lb:                                                   ###
  ### For animals from all other breeds ub=lb is defined such that males  ###
  ### and females from each breed have equal contributions and breed      ###
  ### contributions are as in vector bc.                                  ###
  ###########################################################################
  
  if(!("ub" %in% names(con))){con$ub <- c(M=NA,F=NA)}
  if(!("lb" %in% names(con))){con$lb <- c(M=0, F=0 )}
  
  nIndiv  <- table(phen$Breed, phen$Sex)
  nBreeds <- nrow(nIndiv)
  nIndiv  <- ifelse(phen$Sex==1, nIndiv[phen$Breed, 1], nIndiv[phen$Breed, 2])
  ub      <- bc[phen$Breed]/(2*nIndiv)
  names(ub)<-NULL
  lb      <- ub
  
  if("M" %in% names(con$ub) & "F" %in% names(con$ub)){
    equalMaleCont   <- con$ub["M"] %in% c(-1)
    equalFemaleCont <- con$ub["F"] %in% c(-1)
    if(!equalMaleCont){
      ub[phen$Breed == thisBreed & phen$Sex == 1] <- bc[thisBreed]*con$ub["M"]
      lb[phen$Breed == thisBreed & phen$Sex == 1] <- bc[thisBreed]*con$lb["M"]
      }
    if(!equalFemaleCont){
      ub[phen$Breed == thisBreed & phen$Sex == 2] <- bc[thisBreed]*con$ub["F"]
      lb[phen$Breed == thisBreed & phen$Sex == 2] <- bc[thisBreed]*con$lb["F"]
      }
  }else{
    ub[names(con$ub)] <- bc[thisBreed]*con$ub
    lb[names(con$lb)] <- bc[thisBreed]*con$lb
  } 
  con$ub <- ub
  con$lb <- lb
  
  ###########################################################################
  ### The bounds for the constraints of kinships within breed need to be  ###
  ### adjusted because the kinship matrices are enlarged by identity      ###
  ### matrices:                                                           ###
  ###########################################################################
  con2 <- con
  for(i in cwBreedKin){
    j <- paste("ub.",i,sep="")
    con2[[j]] <- con[[j]]*bc[thisBreed]^2+sum(con$ub[phen$Breed!=thisBreed]^2)
    }

  ###########################################################################
  ### Adjust parameters for conditional kinships:                         ###
  ### The functions for computing conditional kinship needs to be modified###
  ### because the kinship matrices are enlarged by identity               ###
  ### matrices to cover individuals from all breeds.                      ###
  ###########################################################################
    
  if(!is.null(condProb)){
    for(i in names(condProb)){
      attributes(K)$condProb[[i]]["a"] <- sum(con$ub[phen$Breed!=thisBreed]^2) - 2*sum(con$ub[phen$Breed!=thisBreed]) + sum(con$ub[phen$Breed!=thisBreed])^2
      attributes(K)$condProb[[i]]["b"] <- sum(con$ub[phen$Breed!=thisBreed]^2)
      }
  }
  
  ############################################################################
  ### Bounds for linear constraints need to be adjusted.                   ###
  ### The new bound refers to animals from all breeds                      ###
  ############################################################################
  
  const <- getconst(con, Traits)
  for(i in rownames(const)){con2[[i]]<-con[[i]]*bc[thisBreed];names(con2[[i]])<-NULL}
  
  ############################################################################
  ###  perform optimization                                                ###
  ############################################################################
  gc()
  Res <- opticont(method=method, K=K, phen=phen, con=con2, solver=solver, quiet=quiet, make.definite=make.definite, ...)
  
  ######## readjust bounds for constraints  #######
  
  Res$con <- con
  Res$con$ub<-Res$con$ub[Res$parent$Breed==thisBreed]
  Res$con$lb<-Res$con$lb[Res$parent$Breed==thisBreed]
  Res$con$ub[is.na(Res$con$ub)]<-0.5
  Res$con$lb[is.na(Res$con$lb)]<-0.0
  
  for(i in cwBreedKin){
    Res$quadcon[i] <-  con[[paste("ub.",i,sep="")]]
  }
  
  if(nrow(const)>0){
    for(i in rownames(const)){
      Res$lincon[i,"val"]<-con[[i]]
    }
  }
  
  Res$parent    <- Res$parent[Res$parent$Breed==thisBreed,]
  Res$parent$lb <- Res$parent$lb/bc[thisBreed]
  Res$parent$oc <- Res$parent$oc/bc[thisBreed]
  Res$parent$ub <- Res$parent$ub/bc[thisBreed]
  Res$parent$ub[Res$parent$ub>0.5]<-0.5
  
  for(i in wBreedKin){
    indiv <- phen$Indiv[phen$Breed==thisBreed]
    Res$meanKin[i]<- t(Res$parent$oc)%*%(K[[i]][indiv,indiv])%*%(Res$parent$oc)
  }
  
  if(!is.null(condProb)){
    indiv <- phen$Indiv[phen$Breed==thisBreed]
    for(i in names(condProb)){
      f1 <- condProb[[i]]["f1"]
      f2 <- condProb[[i]]["f2"]
      Res$meanKin[i]<- 1-(1-t(Res$parent$oc)%*%(K[[f1]][indiv,indiv])%*%(Res$parent$oc))/(t(Res$parent$oc)%*%(K[[f2]][indiv,indiv])%*%(Res$parent$oc))
     # Res$meanKin<-Res$meanKin[setdiff(names(Res$meanKin),c(f1,f2,paste(c(f1,f2),"withinBreed",sep="")))]
    }
  }
  if(phenAsDataTable){setDT(Res$parent)}  
  Res
}