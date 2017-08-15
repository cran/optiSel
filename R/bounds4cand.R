
"bounds4cand" <- function(phen, uniform, lb, ub, eqConSexes, quiet=FALSE){
  
  id     <- phen$Indiv
  Breeds <- names(eqConSexes)
  
  lbval  <- setNames(rep(0,  length(id)), id)
  ubval  <- setNames(rep(NA, length(id)), id)
  
  if(!is.null(uniform)){
    if(identical(uniform, "female")){uniform <- paste0(Breeds,".female")}
    if(identical(uniform,   "male")){uniform <- paste0(Breeds,".male")}
    
    #### Check parameter 'uniform' ###
    validUnif <- c(Breeds, paste0(Breeds[eqConSexes], ".female"), paste0(Breeds[eqConSexes], ".male"))
    if(!is.vector(uniform)){stop("Constraint 'uniform' must be a character vector.\n")}
    uniform <- as.character(uniform)
    if(any( ((paste0(Breeds, ".female") %in% uniform) && !eqConSexes))){
      stop("'BREEDNAME.female' in constraint 'uniform' is only permitted if the breed does not contain NA in column Sex of cand$phen.\n")
    }
    if(any( ((paste0(Breeds, ".male") %in% uniform) && !eqConSexes))){
      stop("'BREEDNAME.male' in constraint 'uniform' is only permitted if the breed does not contain NA in column Sex of cand$phen.\n")
    }
    if(!all(uniform %in% validUnif)){
      uniform <- uniform[uniform %in% validUnif]
      if(!quiet){cat(paste0("Note: uniform=c('", paste(uniform, collapse="', '"),"') is used.\n"))}
    }
    ### Define upper and lower bounds ###########
    ### for groups with uniform contributions ###
    for(b in Breeds){
      isBreedb <- phen$Breed==b
      
      if(eqConSexes[b] && ((b %in% uniform)||(paste0(b,".female")%in% uniform))){
        Nf <- sum(isBreedb & phen$Sex=="female")
        lbval[isBreedb & phen$Sex=="female"] <- 1/(2*Nf)
        ubval[isBreedb & phen$Sex=="female"] <- 1/(2*Nf)
      }
      if(eqConSexes[b] && ((b %in% uniform)||(paste0(b,".male")%in% uniform))){
        Nm <- sum(isBreedb & phen$Sex=="male")
        lbval[isBreedb & phen$Sex=="male"]   <- 1/(2*Nm)
        ubval[isBreedb & phen$Sex=="male"]   <- 1/(2*Nm)
      }
      if((!eqConSexes[b]) && (b %in% uniform)){
        lbval[isBreedb] <- 1/sum(isBreedb)
        ubval[isBreedb] <- 1/sum(isBreedb)
      }
    }
  }
  
  if(!is.null(ub)){
    if(!is.vector(ub)){stop("Constraint ub must be a named numeric vector.\n")}
    if(is.null(names(ub))){stop("Constraint ub must be a named numeric vector.\n")}
    if(!all(names(ub) %in% id)){
      stop(paste0("The following individuals from constraint 'ub' do not appear in data frame cand$phen:", setdiff(names(ub), id),".\n"))
    }
    if(any(ub<0)){stop("Some upper bounds are negative.\n")}
    
    ubval[names(ub)] <- ub 
  }
  
  if(!is.null(lb)){
    if(!is.vector(lb)){stop("Constraint lb must be a named numeric vector.\n")}
    if(is.null(names(lb))){stop("Constraint lb must be a named numeric vector.\n")}
    if(!all(names(lb) %in% id)){
      stop(paste0("The following individuals from constraint 'lb' do not appear in data frame cand$phen:", setdiff(names(lb),id),".\n"))
    }
    if(any(lb<0)){stop("Some lower bounds are negative.\n")}
    
    lbval[names(lb)] <- lb 
  }
  if(any(!is.na(ubval) & !is.na(lbval) & lbval>ubval)){
    stop("Some lower bounds are larger than the upper bounds.\n")
  }
  
  
  return(list(upper=ubval, lower=lbval))
}