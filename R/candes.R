
"candes"<-function(phen, N=250, quiet=FALSE, bc=NA, ...){
  ### Check data table phen ###################################
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  if(!("data.frame" %in% class(phen))){
    stop("Argument 'phen' is not a data frame.\n")
    }
  if(!("Sex" %in% colnames(phen))){
    phen$Sex <- NA
  }
  if(!("Indiv" %in% colnames(phen))){
    stop("Column 'Indiv' with IDs of individuals is missing in data frame 'phen'.\n")
  }
  if(!("Breed" %in% colnames(phen))){
    phen$Breed <- "missing"
  }
  
  phen$Indiv <- as.character(phen$Indiv)
  phen$Sex   <- as.character(phen$Sex)
  phen$Breed <- as.character(phen$Breed)
  
  if(any(is.na(phen$Breed))){
    stop("Column Breed of data table phen contains 'NA'.\n")
    }
  if(any(is.na(phen$Indiv))){
    stop("Some Individual IDs are NA in data frame phen.\n")
    }
  if(any(duplicated(phen$Indiv))){
    stop("Some Individuals appear twice in data frame phen.\n")
    }
  rownames(phen)<-phen$Indiv
  
  if(!all(phen$Sex %in% c("male", "female", NA))){
    stop("Some sexes are not coded as 'male' and 'female'.\n")
  }
  if(all(phen$Sex %in% c("male"))){
    stop("Females are missing in data frame 'phen'.\n")
    }
  if(all(phen$Sex %in% c("female"))){
    stop("Males are missing in data frame 'phen'.\n")
  }
  
  ### Make assumption about the number N of selection ###
  ######## candidates in the next generation ############
  Breeds      <- unique(phen$Breed)
  singleBreed <- length(Breeds)==1

  N <- as.numeric(N)
  if(!is.vector(N) || length(N)>1 || any(is.na(N)) || N<=0){
    stop("N must be positive or Inf.\n")
  }
  
  if(!quiet){
    #cat(paste0("Assuming N=",N," selection candidates per generation.\n"))
  }
  
 
  ### Convert kinships to class quadFun or ratioFun ###########
  obj <- list(...)
  
  for(i in seq_along(obj)){
    if((class(obj[[i]])!="quadFun") && (class(obj[[i]])!="ratioFun") && (!is.matrix(obj[[i]]))){
      stop("All '...' arguments must have class 'quadFun', 'ratioFun', or 'matrix'.\n")
    }
    
    Name <- names(obj)[i]
    if(is.null(Name)){stop(paste("A name must be specified for all additional parameters.\n",sep=""))}
    
    if(class(obj[[i]])=="quadFun" || class(obj[[i]])=="ratioFun"){
      obj[[i]]$name <- Name
    }
    if(is.matrix(obj[[i]])){
      if(is.null(rownames(obj[[i]]))){
        stop(paste("Rownames must be specified for the matrix at position ",i,".\n",sep=""))
      }
      obj[[i]]   <- optiSolve::quadfun(Q=obj[[i]], d=0, id=rownames(obj[[i]]), name=Name)
    }
  }

  ### Append kinship matrices for single breeds ###
  Seq <- seq_along(obj)
  for(i in Seq){
    if((class(obj[[i]])=="quadFun") && all(phen$Indiv %in% obj[[i]]$id) && (length(Breeds)>1)){
      Name <- obj[[i]]$name
      bname <- str_extract(Name, paste(Breeds,collapse="|"))
      if(!is.na(bname)){stop(paste0("Kinship" , Name, " contains kinships from more than one breed, so the parameter name should not contain a breed name.\n"))}
      for(b in Breeds){
        Kname <- paste(Name, b, sep=".")
        if(!(Kname %in% names(obj))){
          use <- obj[[i]]$id %in% phen$Indiv[phen$Breed==b]
          obj[[Kname]] <- optiSolve::quadfun(Q=obj[[i]]$Q[use, use], a=obj[[i]]$a[use],  d=0, id=obj[[i]]$id[use], name=Kname)
        }
      }
    }
  }
  
  
  ### Get initial number of individuals ###
  initN <- rep(0, length(obj))
  for(i in seq_along(obj)){
    initN[i] <- length(obj[[i]]$id)
  }
  
  ### Check if kinship matrices contain all required individuals ###
  for(i in seq_along(obj)){
    if(singleBreed){
      ### Check if all individuals from 'phen' are included in ###############
      ### all kinship matrices ###############################################
      if(!all(phen$Indiv %in% obj[[i]]$id)){
        stop(paste0("Data frame 'phen' contains individuals not included in argument ", obj[[i]]$name, ".\n"))
      }
      obj[[i]]$breed <- Breeds
    }else{
      ### Check if all kinship matrices contain either all individuals or ###
      ### all individuals from exactly one breed ############################
      if(identical(sort(phen$Indiv), sort(obj[[i]]$id))){
        if(class(obj[[i]])=="ratioFun"){stop("Kinships at native alleles must include individuals from only one breed.\n")}
        obj[[i]]$breed <- "across breeds"
      }else{
        thisBreed <- unique(phen[intersect(phen$Indiv, obj[[i]]$id),"Breed"])
        if(length(thisBreed)>1){
          stop(paste0("Argument ",obj[[i]]$name, " must include either all individuals from 'phen', or all individuals from exactly one breed.\n"))
        }
        if(!all(phen$Indiv[phen$Breed==thisBreed] %in% obj[[i]]$id)){
          stop(paste0("Some individuals from breed ", thisBreed, " are missing in argument", obj[[i]]$name, ".\n"))
        }
        obj[[i]]$breed <- thisBreed
      }
    }
  }
  
  ### Adjust drift terms d, d1, d2 ###########

  for(i in seq_along(obj)){
    if(obj[[i]]$breed=="across breeds"){
      if(class(obj[[i]])=="quadFun"){
        obj[[i]]$d <- 0
      }else{
        obj[[i]]$d1 <- 0
        obj[[i]]$d2 <- 0
      }
    }else{
      if(is.na(N)){stop("Parameter N is NA.\n")}
      if(class(obj[[i]])=="quadFun"){
        ids        <- obj[[i]]$id
        diagQ      <- diag(obj[[i]]$Q)
        if(any(is.na(phen$Sex))){
          obj[[i]]$d <- c(1-mean(diagQ[ids %in% phen$Indiv]))/(2*N)
        }else{
          isfemale   <- (ids %in% phen$Indiv) & (phen[ids,"Sex"] %in% "female")
          ismale     <- (ids %in% phen$Indiv) & (phen[ids,"Sex"] %in% "male")
          obj[[i]]$d <- c(1-(0.5*mean(diagQ[isfemale]) + 0.5*mean(diagQ[ismale])))/(2*N)
        }
      }else{
        obj[[i]]$d1 <- obj[[i]]$d1*initN[i]/N
        obj[[i]]$d2 <- obj[[i]]$d2*initN[i]/N
      }
    }
  }
  
  ### Remove individuals not included in data frame 'phen' from kinships ##############
  ### Enlarge kinship matrices for one breed to include individuals from all breeds ###
  ### and sort individuls in kinship matrices according to 'phen' #####################
  
  for(i in seq_along(obj)){
    obj[[i]] <- adjust(obj[[i]], phen$Indiv)
  }
  
  
  ### Check that trait values are provided for exactly ###
  ### one breed or for all breeds and define an ##########
  ### individual trait for each breed  ###################
  Traits <- gettraits(phen, quiet=quiet)
  for(trait in Traits){
    if(all(!is.na(phen[[trait]]))&&(!singleBreed)){
      phenT <- phen[[trait]]
      phen[[trait]] <- NULL
      phen[[trait]] <- phenT
      for(b in Breeds){
        phen[[paste(trait,b,sep=".")]] <- ifelse(phen$Breed==b, phenT, NA)
      }
    }
  }
  Traits <- gettraits(phen, quiet=TRUE)

  ### Define BreedforVar  #######################################
  
  if(singleBreed){
    BreedforVar <- setNames(rep(Breeds,length(Traits)), Traits)
  }else{
    BreedforVar <- setNames(rep("across breeds", length(Traits)), Traits)
    for(trait in Traits){
      if(any(is.na(phen[[trait]]))){
        BreedforVar[trait] <- unique(phen$Breed[!is.na(phen[[trait]])])
      }
    } 
  }
  
  BreedforVar <- setNames(c(BreedforVar, unlist(lapply(obj,function(x){x$breed}))), c(Traits, names(obj)))
  
  ### Get kinship bcKin for estimating the optimum breed composition ##
  bcKin <- NA
  if(!singleBreed){
    if(is.character(bc)){
      bcKin <- bc
      bc    <- NA
      if((!(bcKin %in% names(obj)))|| (class(obj[[i]])!="quadFun")){
        stop(paste0("bc=", bcKin, " is not a valid name of a kinship.\n"))
      }
      if(!identical(obj[[i]]$id, phen$Indiv)){
        stop(paste0("Kinship ",bcKin," must contain all individuals in data frame phen.\n"))
      }
    }else{
      if(any(is.na(bc))){
        iKin <- NA
        for(i in rev(seq_along(obj))){
          if(class(obj[[i]])=="quadFun"){
            if(identical(sort(phen$Indiv), sort(obj[[i]]$id))){
              iKin <- i
            }
          }
        }
        if(is.na(iKin)){
          stop("No kinship is suitable for estimating the optimum breed composition. Please provide parameter 'bc',\n")
        }else{
          bcKin <- names(obj)[iKin]
          if(!quiet){cat(paste0("Breed contributions minimize ",bcKin," across breeds.\n"))}
        }
      }
    }
  }
  
  ### Get optimum breed composition ##
  if(singleBreed){
    bc <- setNames(c(1), Breeds)
  }else{
    if(any(is.na(bc)) && !is.na(bcKin)){
      lbbc <- setNames(rep(0.10/length(Breeds),length(Breeds)), Breeds)
      bc   <- opticomp(obj[[bcKin]]$Q, phen, lb=lbbc)$bc
      bc   <- checkbc(phen$Breed, bc)
    }else{
      bc <- checkbc(phen$Breed, bc)
    }
  }
  


  
  
  ### Table current values of kinships and traits ###
  
  cur <- data.frame(
    Type     = c(rep("trait", length(Traits)), unlist(lapply(obj,function(x){class(x)})), rep("breed cont",length(bc))),
    Var       = c(Traits, names(obj),           paste0("bc.",names(bc))), 
    Breed     = c(BreedforVar,                  names(bc)), 
    Val       = c(rep(NA, length(BreedforVar)), bc),
    obj.fun.and.constraints = '',
    stringsAsFactors=FALSE)
  
  cur$Type <- mapvalues(cur$Type, c("quadFun", "ratioFun"), c("kinship", "nat. kin."), warn_missing = FALSE)
  
  for(i in 1:nrow(cur)){
    if(cur$Type[i] %in% c("kinship", "nat. kin.")){
      cur[i,"obj.fun.and.constraints"] <- paste0(c("min.", "ub."), cur$Var[i], collapse=", ")
    }
    if(cur$Type[i] %in% c("trait")){
      cur[i,"obj.fun.and.constraints"] <- paste0(c("min.", "max.", "lb.","eq.", "ub."), cur$Var[i], collapse=", ")
    }
  }
  
  
  if(singleBreed){cur <- cur[cur$Var!=paste0("bc.",names(bc)),]}
  cur$Name <- str_replace(cur$Var, paste0(paste0(".",Breeds),collapse="|"),"")
  cur$Name <- str_replace(cur$Name, "bc", "Breed Contribution")
  
  Seq <- which(cur$Type %in% c("trait", "kinship", "nat. kin."))
  for(i in Seq){
    v <- cur$Var[i]
    u <- getu(phen, cur$Breed[i], bc)
    if(!is.null(u)){
      if(cur$Type[i]=="trait"){
        cur[i,"Val"] <- sum(u*phen[[v]], na.rm=TRUE)
      }
      if(cur$Type[i]=="kinship"){
        Q <- obj[[v]]$Q
        #uTu <- sum(u*u)
        if(singleBreed){
          Nf <- sum(phen[rownames(Q),"Sex"]=="female")
          Nm <- sum(phen[rownames(Q),"Sex"]=="male")
          uf <- ifelse(phen[rownames(Q),"Sex"]=="female", 1/Nf, 0)
          um <- ifelse(phen[rownames(Q),"Sex"]=="male",   1/Nm, 0)
          u  <- (1/2)*um + (1/2)*uf
          s1 <- ((N-2)/(4*N))*c((Nf/(Nf-1))*t(uf)%*%Q%*%uf - (1/(Nf-1))*t(uf)%*%diag(Q))
          s2 <- ((N-2)/(4*N))*c((Nm/(Nm-1))*t(um)%*%Q%*%um - (1/(Nm-1))*t(um)%*%diag(Q))
          cur[i,"Val"] <- s1 + s2 + (1/4)*t(um)%*%Q%*%uf + (1/4)*t(uf)%*%Q%*%um + (1/N)*t(u)%*%diag(Q)
        }else{
          cur[i,"Val"] <- c(t(u)%*%Q%*%u)
        }
       }
      if(cur$Type[i]=="nat. kin."){
        
        
        if(singleBreed){
          Nf <- sum(phen[obj[[v]]$id,"Sex"]=="female")
          Nm <- sum(phen[obj[[v]]$id,"Sex"]=="male")
          uf <- ifelse(phen[obj[[v]]$id,"Sex"]=="female", 1/Nf, 0)
          um <- ifelse(phen[obj[[v]]$id,"Sex"]=="male",   1/Nm, 0)
          u  <- (1/2)*um + (1/2)*uf
          
          Q  <- obj[[v]]$Q1
          s1 <- ((N-2)/(4*N))*c((Nf/(Nf-1))*t(uf)%*%Q%*%uf - (1/(Nf-1))*t(uf)%*%diag(Q))
          s2 <- ((N-2)/(4*N))*c((Nm/(Nm-1))*t(um)%*%Q%*%um - (1/(Nm-1))*t(um)%*%diag(Q))
          x1 <- s1 + s2 + (1/4)*t(um)%*%Q%*%uf + (1/4)*t(uf)%*%Q%*%um + (1/N)*t(u)%*%diag(Q)

          Q  <- obj[[v]]$Q2
          s1 <- ((N-2)/(4*N))*c((Nf/(Nf-1))*t(uf)%*%Q%*%uf - (1/(Nf-1))*t(uf)%*%diag(Q))
          s2 <- ((N-2)/(4*N))*c((Nm/(Nm-1))*t(um)%*%Q%*%um - (1/(Nm-1))*t(um)%*%diag(Q))
          x2 <- s1 + s2 + (1/4)*t(um)%*%Q%*%uf + (1/4)*t(uf)%*%Q%*%um + (1/N)*t(u)%*%diag(Q)
        }else{
          Q  <- obj[[v]]$Q1
          x1 <- c(t(u)%*%Q%*%u)
          Q  <- obj[[v]]$Q2
          x2 <- c(t(u)%*%Q%*%u)
        }
        cur[i,"Val"] <- x1/x2
      }
    }
  }
  
  if(!quiet){
    cur2      <- cur[cur$Type %in% c("trait", "kinship", "nat. kin."), c("Type", "Var",  "Breed", "obj.fun.and.constraints","Val")]
    cur2$Val  <- sprintf("%8.4f", cur2$Val)
    cur2$Breed[cur2$Breed %in% Breeds] <- paste0("in ", cur2$Breed[cur2$Breed %in% Breeds])
    cur2$Type <- paste0("for ", cur2$Type)
    cur2$Var  <- paste0(" '", cur2$Var, "' ")
    cur2$Type <- str_pad(as.character(cur2$Type),  max(nchar(cur2$Type)), "right")
    cur2$Var  <- str_pad(as.character(cur2$Var),   max(nchar(cur2$Var)), "right")
    cur2$Breed<- str_pad(as.character(cur2$Breed), max(nchar(cur2$Breed)), "right")
    if(identical(Breeds, "missing")){cur2$Breed<-''}
    cur2$Type <- paste0(cur2$Type, cur2$Var, cur2$Breed, ":")
    cur2      <- cur2[,c("Type", "obj.fun.and.constraints", "Val")]
    cat("\n")
    colnames(cur2) <- c("Mean values of the parameters are:", "1", "  Value")
    print(cur2[,c(1,3),drop=FALSE], right = FALSE, row.names=FALSE)
    cat("\n")
    
    colnames(cur2) <- c("Available objective functions and constraints:", "2", "  Value")
    cur2[[1]] <- paste(cur2[[1]], cur2[[2]])
    print(cur2[,1,drop=FALSE], right = FALSE, row.names=FALSE, col.names=FALSE)
    cat("\n")
    cat(" ub  lb uniform\n\n")
  }
  
  obj$phen    <- phen  
  obj$current <- cur
  obj$mean    <- as.data.frame(as.list(setNames(cur$Val, cur$Var)),stringsAsFactors=FALSE)
  obj$bc      <- bc
  class(obj)  <- "candes"
  
  attr(obj,"singleBreed") <- singleBreed
  
  return(obj)
}