

"matings" <- function(cand, Kin, N=2*sum(cand$Sex=="female"), alpha=1, ub.nOff=NA, max=FALSE, ...){
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
  }
  
  nVar <- nrow(Kin)*ncol(Kin)

  G <- -diag(nVar)
  h <- rep(0, nVar)
  
  if(!is.na(ub.nOff)){
    G <- rbind(G, diag(nVar))
    h <- c(h, rep(ub.nOff, nVar))
  }
  
  if(alpha<1){
    G <- rbind(G,  ConH)
    h <- c(h, rhsH)
  }
  
  opt <- list(...)
  if("maxit"        %in% names(opt)){opt$maxit        <- as.integer(opt$maxit)}
  if("verbose"      %in% names(opt)){opt$verbose      <- as.integer(opt$verbose)}
  if("mi_max_iters" %in% names(opt)){opt$mi_max_iters <- as.integer(opt$mi_max_iters)}
  opt <- do.call(ecos.control, opt)
    
  dims <- list(l=length(h), q=NULL, e=0L)
  A <- Matrix(rbind(ConF, ConM), sparse=TRUE)
  G <- Matrix(G, sparse=TRUE)
  sig <- ifelse(max,-1,1)
  res <- ECOS_csolve(c=sig*c(Kin), G=G, h=h, dims=dims, A=A, b=c(rhsF, rhsM), int_vars=as.integer(1:nVar), control=opt)
  res$x <- round(res$x, 0)

  Mating <- matrix(res$x, nrow=nrow(Zeros), ncol=ncol(Zeros))
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
    }
  attributes(Matings)$objval <- sum(c(Kin)*res$x)/sum(res$x)
  attributes(Matings)$info <- res$infostring
  cat(paste(res$infostring,"\n"))
  Matings
}