
"help.opticont"<-function(K, phen){
  condProb <- attr(K,"condProb")
  condProbMat <- c()
  for( i in names(condProb)){
    condProbMat<-c(condProbMat, condProb[[i]]["f1"], condProb[[i]]["f2"])
  }
  Traits  <- colnames(phen)[c(-1,-2)]
  if(ncol(phen)>2){
    x<-NULL;for(i in 3:ncol(phen)){x<-c(x,is.numeric(phen[[i]]))}
    Traits  <- setdiff(Traits[x], c("lb", "ub", "oc", "Born", "Sex"))
  }
  cat("Available objective functions: \n")
  cat("     ", paste("min.", c(setdiff(names(K), condProbMat), names(condProb)), sep=""),"\n", sep="  ")
  if(length(Traits)>0){cat("     ", paste("min.", Traits, sep="") ,"\n", sep="  ")}
  if(length(Traits)>0){cat("     ", paste("max.", Traits, sep="") ,"\n\n", sep="  ")}
  cat("Available constraints: \n")
  cat("     ", paste("ub.", c(setdiff(names(K), condProbMat), names(attr(K,"condProb"))), sep=""),"\n", sep="  ")
  if(length(Traits)>0){cat("     ", paste("lb.", Traits, sep="") ,"\n", sep="  ")}
  if(length(Traits)>0){cat("     ", paste("ub.", Traits, sep="") ,"\n", sep="  ")}
  if(length(Traits)>0){cat("     ", paste("eq.", Traits, sep="") ,"\n", sep="  ")}
  cat("       ub  lb\n\n")
  cat("Attention: Probably not all these objective functions and\n")
  cat("           constraints make sense in animal breeding\n\n")
}
