

"opticont"<-function(method, cand, con, bc=NA,  solver="default", quiet=FALSE, make.definite=solver=="csdp", ...){
  
  ### Check parameters #############################
  if(class(cand)!="candes"){
    stop("Parameter 'cand' must be created with function candes.")
  }
  
  #if(!("Breed" %in% colnames(cand$phen))){stop("Column 'Breed' is missing in data frame cand$phen.\n")}
  
  namesK   <- setdiff(names(cand), c("phen", "current", "mean", "bc"))
  Traits   <- gettraits(cand$phen)
  validCon <- c(paste0(c("lb.","ub."),  rep(Traits,each=length(Traits))), paste0("ub.",  namesK), "ub", "lb", "uniform")
  validMet <- c(paste0(c("min.","max."),rep(Traits,each=length(Traits))), paste0("min.", namesK))
  
  if(!(method %in% validMet)){
    stop(paste0("The method must be one of ", paste(validMet, collapse=", "), "\n."))
  }
  if(!all(names(con) %in% validCon)){
    stop(paste0("The following constraints do not have valid names: ", paste(setdiff(names(con), validCon), collapse=", "), "\n."))
  }

  ### Find optimum breed composition ###  
  if(any(is.na(bc))){bc <- cand$bc}
  bc <- checkbc(cand$phen$Breed, bc)
  if(any(is.na(bc))){stop("Argument 'bc' is not defined.\n")}

  ### Convert constraints to data frame #######
  
  condf <- con2df(con, Traits)
  
  #####################################
  ### Define prelim. cop including  ###
  ### - objective function,         ###
  ### - quadratic constraints,      ###
  ### - rational constraints,       ###
  ### - max                         ###
  #####################################
  
  mycop     <- qcrc4cand(cand, condf, namesK, bc)
  mycop$f   <- f4cand(cand, method, Traits, bc)
  mycop$max <- str_sub(method, 1, 3)=="max"
  
  ####################################################
  ### Determine for which breeds the contributions ###
  ###    of males and females should be equal.     ###
  ####################################################
  
  Breeds    <- unique(cand$phen$Breed)
  eqConSexes <- setNames(rep(NA, length(Breeds)), Breeds)
  
  for(b in Breeds){
    eqConSexes[b] <- all(!is.na(cand$phen$Sex[cand$phen$Breed==b]))
  }
  if((!all(eqConSexes)) && (!quiet)){
    cat("Assuming: Sexes do not have equal contributions for \n")
    cat(paste0("          ", paste0(Breeds[!eqConSexes], collapse=", "), ", because some sexes are NA:\n"))
  }

  #####################################
  ### Define upper and lower bounds ###
  #####################################
  
  bound <- bounds4cand(cand$phen, con[["uniform"]], con[["lb"]], con[["ub"]], eqConSexes,quiet=quiet)
  
  mycop$lb <- lbcon(val=bound$lower*bc[cand$phen$Breed], id=cand$phen$Indiv)
  mycop$ub <- ubcon(val=bound$upper*bc[cand$phen$Breed], id=cand$phen$Indiv)

  #####################################
  ### Define linear constraints     ###
  #####################################
  
  mycop$lc <- lc4cand(cand$phen, bc, condf, Traits, eqConSexes)
  
  
  ######################################
  ### Solve the optimization problem ###
  ######################################
  
  mycop  <- do.call(cop, mycop)
  res    <- solvecop(mycop, solver=solver, make.definite=make.definite, quiet=quiet, ...)
  sy     <- validate(mycop, res, quiet=TRUE)
  res$summary <- sy$summary
  res$info    <- sy$info
  bound$upper[is.na(bound$upper)] <- 0.5
  res$parent  <- data.frame(cand$phen, lb=bound$lower, oc=res$x, ub=bound$upper)
  
  ####### Add breed name to summary ######## 
  
  BreedforVar <- setNames(cand$current$Breed, cand$current$Var)
  res$summary$Breed <- BreedforVar[res$summary$Var]
  Breeds <- unique(cand$phen$Breed)
  x <- str_extract(res$summary$Var, paste(Breeds, collapse="|"))
  res$summary$Breed[!is.na(x)] <- x[!is.na(x)]
  res$summary$Breed[is.na(res$summary$Breed)] <- "" 
  
  res$summary$Name <- str_replace(res$summary$Var, "bc.", "breed contribution.")
  res$summary$Name <- str_replace(res$summary$Name,"scd.", "sex contrib. diff..")
  res$summary$Name <- str_replace(res$summary$Name, paste(paste0("\\.",Breeds), collapse="|"),"")
  
  ####### Add original threshold values to summary ######## 
  
  cnames <- res$summary$Var[res$summary$Var %in% condf$var]
  res$summary[res$summary$Var %in% condf$var, "Bound"] <- condf[cnames, "val"]

  
  ####### Compute costraint values for summary ######## 
  
  for(i in 1:nrow(res$summary)){
    var <- res$summary$Var[i]
    brd <- res$summary$Breed[i]
    
    if(!(var %in% colnames(cand$phen))){
      redvar <- str_replace(var, paste(paste0(".",Breeds), collapse="|"),"")
      if(redvar %in% colnames(cand$phen)){
        var <- redvar
      }
    }
    
    if(var %in% colnames(cand$phen)){
      TraitVal <- cand$phen[[var]]
      if(brd!="across breeds"){
        TraitVal[cand$phen$Breed != brd] <- NA
        }
      thisbc <- sum((!is.na(TraitVal))*res$x)
      res$summary[i,"Val"] <- sum(TraitVal*res$x, na.rm=TRUE)/thisbc
    }else{
      if((var %in% namesK) && (class(cand[[var]])=="quadFun")){
        nameB <- cand[[var]]$breed
        if(nameB == "across breeds"){
          thisbc <- 1
        }else{
          thisbc <- bc[nameB]
        }
        res$summary[i,"Val"] <- c(t(res$x)%*%(cand[[var]]$Q)%*%(res$x)/(thisbc^2)+ cand[[var]]$d)
      }
    }
  }
  
  
  if(!quiet){
    Sy     <- res$summary
    if(identical(Breeds, "missing")){Sy$Breed[Sy$Breed=="missing"] <- ""}
    Sy$Var <- paste0(str_pad(Sy$Name, max(nchar(Sy$Name))+1, "right"), Sy$Breed)
    Sy <- list(summary=Sy, info=res$info)
    class(Sy)<- "copValidation"
    print.copValidation(Sy)
  }
  
  res$x <- NULL  
  Selection <- (!(res$summary$Var %in% c("lower bounds", "upper bounds"))) 
  Selection <- Selection & (!(str_detect(res$summary$Var,"scd.")))
  Selection <- Selection & (!(str_detect(res$summary$Var,"bc.")))
  Selection[1] <- FALSE
  res$mean  <- setNames(res$summary$Val[Selection],res$summary$Var[Selection])
  res$mean  <- as.data.frame(as.list(res$mean), stringsAsFactors = FALSE)
  res$bc    <- as.data.frame(as.list(bc), stringsAsFactors = FALSE)
  res$solver <- NULL
  res$status <- NULL
  res$obj.fun <- setNames(res$summary$Val[1], res$summary$Var[1])
  if(length(Breeds)==1){res$mean<-res$mean[,colnames(res$mean)!=paste0("bc.",Breeds)]}
  for(b in Breeds){
    use <- res$parent$Breed==b
    res$parent$oc[use] <- res$parent$oc[use]/bc[b]
  }
   
  
  return(res)
  
}