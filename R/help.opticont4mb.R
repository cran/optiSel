
"help.opticont4mb"<-function(K, phen){
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  
  phen$Breed <- as.character(phen$Breed)
  phen$Indiv <- as.character(phen$Indiv)
  rownames(phen) <- phen$Indiv
  
  condProb <- attr(K,"condProb")
  condProbMat <- c()
  for( i in names(condProb)){
    condProbMat<-c(condProbMat, condProb[[i]]["f1"], condProb[[i]]["f2"])
  }
  Names1 <- setdiff(names(K), condProbMat)
  NamesAcross <- character(0)
  for(i in Names1){
    if(length(table(phen[rownames(K[[i]]),"Breed"]))>1){
      NamesAcross <- c(NamesAcross, paste(i, ".acrossBreeds", sep=""))
    }
  }
  Names2 <- names(attributes(K)$condProb)
  
  Names <- c(NamesAcross, Names1, Names2)
  
  Traits  <- colnames(phen)[c(-1,-2)]
  if(ncol(phen)>2){
    x<-NULL;for(i in 3:ncol(phen)){x<-c(x,is.numeric(phen[,i]))}
    Traits  <- setdiff(Traits[x], c("lb", "ub","oc", "Born","Sex"))
  }
  cat("Available objective functions: \n")
  cat("     ", paste("min.", Names, sep=""),"\n", sep="  ")
  if(length(Traits)>0){cat("     ", paste("min.", Traits, sep="") ,"\n", sep="  ")}
  if(length(Traits)>0){cat("     ", paste("max.", Traits, sep="") ,"\n\n", sep="  ")}
  cat("Available constraints: \n")
  cat("     ", paste("ub.", Names, sep=""),"\n", sep="  ")
  if(length(Traits)>0){cat("     ", paste("lb.", Traits, sep="") ,"\n", sep="  ")}
  if(length(Traits)>0){cat("     ", paste("ub.", Traits, sep="") ,"\n", sep="  ")}
  if(length(Traits)>0){cat("     ", paste("eq.", Traits, sep="") ,"\n", sep="  ")}
  cat("       ub  lb\n\n")
  cat("Attention: Probably not all these objective functions and\n")
  cat("           constraints make sense in animal breeding\n\n")
}
