

"matings" <- function(cand, Kin, N=2*sum(cand$Sex=="female"), alpha=1, ub.nOff=NA, max=FALSE){
  if(!("Indiv" %in% colnames(cand))){stop("Column 'Indiv' is missing in cand.")}
  if(!("Sex"   %in% colnames(cand))){stop("Column 'Sex' is missing in cand.")}
  if(!("oc"    %in% colnames(cand))){stop("Column 'oc' is missing in cand.")}
  
  ub.herd <- alpha * table(cand$herd)
  cand$nOff <- noffspring(cand, N=N)$nOff
  cand <- cand[cand$nOff>0,]
  Sire <- cand$Indiv[cand$Sex=="male"]
  Dam  <- cand$Indiv[cand$Sex=="female"]

  if(is.null(rownames(Kin))||is.null(rownames(Kin))){
    stop("Row names and column names of Matrix 'Kin' must be the individual IDs.")
  }

  if(!all(Sire %in% rownames(Kin))){
    stop("All male candidates must appear in the row names of matrix Kin.")
  }
  
  if(!all(Dam %in% colnames(Kin))){
    stop("All female candidates must appear in the column names of matrix Kin.")
  }
  
  if(alpha<1){
    if("herd" %in% colnames(cand)){
      cand$herd <- as.character(cand$herd)
    }else{
      stop("Column 'herd' is needed if alpha<1.")
    }
  }

  Kin <- Kin[Sire, Dam]
  Zeros <- 0*Kin
  
  rhsM <- cand$nOff[cand$Indiv %in% Sire]
  ConM <- NULL
  for(k in 1:length(Sire)){
    Con <- Zeros
    Con[k,] <- 1
    ConM <- rbind(ConM, c(Con))
  }
  
  rhsF <- cand$nOff[cand$Indiv %in% Dam]
  ConF <- NULL
  for(k in 1:length(Dam)){
    Con <- Zeros
    Con[,k] <- 1
    ConF <- rbind(ConF, c(Con))
  }

  Mat <- rbind(ConF, ConM)
  Rhs <- c(rhsF, rhsM)
  Dir <- rep("==", length(rhsF)+length(rhsM))
  
  if(alpha<1){
    herds   <- setdiff(unique(cand$herd), NA)
    rhsH    <- rep(ub.herd[herds], each=length(Sire))
    ConH  <- NULL
    for(l in herds){
      for(k in 1:length(Sire)){
        Con <- Zeros
        Con[k,] <- cand[Dam,"herd"] %in% l  
        ConH <- rbind(ConH, c(Con))
      }
    }
    Mat <- rbind(Mat, ConH)
    Rhs <- c(Rhs, rhsH)
    Dir <- c(Dir, rep("<=", length(rhsH)))
  }
  
  if(is.na(ub.nOff)){
    bounds <- NULL
  }else{
    bounds <- list(upper=list(ind=1:(length(Sire)*length(Dam)), val=rep(ub.nOff, length(Sire)*length(Dam))))
  }
  
  res <- Rsymphony_solve_LP(obj=c(Kin), mat=Mat, dir=Dir, rhs=Rhs, types="I", bounds=bounds, max=max)
  
  Mating <- matrix(res$solution, nrow=nrow(Zeros), ncol=ncol(Zeros))
  Mating <- as.data.table(Mating)
  colnames(Mating) <- Dam
  Mating$Sire <- Sire
  Matings <- melt(Mating, id.vars="Sire", variable.name="Dam", value.name="nOff")
  Matings$Dam <- as.character(Matings$Dam)
  Matings <- Matings[Matings$nOff>0,]
  setDF(Matings)
  if(alpha<1){
    herds <- cand$herd
    names(herds) <- cand$Indiv
    Matings$herd <- herds[Matings$Dam]
    #Matings$ub.herd <- floor(c(ub.herd[Matings$herd]))
    }
  attributes(Matings)$objval <- res$objval/sum(res$solution)
  Matings
}