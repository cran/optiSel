
"summary.opticont"<-function(object, ...){
  x <- object
  phen   <- x$parent
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  
  oc  <- phen[,"oc"]
  con <- x$con
  sex <- as.integer(mapvalues(phen[,"Sex"], from=c("male","female"), to=c(1,2)))
  obj.var <- substr(x$method,5,nchar(x$method))
  min.Kin <- (x$method %in% paste("min.",names(x$meanKin),sep=""))

  VarName<-deparse(substitute(object))[1]
  if(str_detect(VarName,"\\(")){VarName<-""}
  
  Res<-data.frame(VarName = VarName, stringsAsFactors = FALSE)
  rownames(Res)<-Res$VarName
  Res$method= x$method
  
#  for(i in names(x$cKin)){
#    f1 <- x$cKin[[i]][1]
#    f2 <- x$cKin[[i]][2]
#    x$meanKin[i]<- 1-c((1-x$meanKin[f1])/(x$meanKin[f2]))
#    }
  
  Res$obj.fun = NA
  if(obj.var %in% names(x$meanKin)){
    Res$obj.fun<-c(x$meanKin[obj.var])
  }else{
    Res$obj.fun<-c(t(oc)%*%phen[,obj.var])
  }
  Res$valid   = FALSE  

  
  for(i in names(x$meanKin)){
    Res[,i]<-x$meanKin[i]
    Res[, paste("ub.",i,sep="")]<-x$quadcon[i]
  }

#  for(i in names(x$meanKin)){
#    Res[, paste("Div.",i,sep="")]<-1-Res[,i]
#  }

  if(length(con$ub==2)&!is.null(names(con$ub))){
    Res$ubM = con$ub[1]
    Res$ubF = con$ub[2]
    con$ub <- con$ub[sex]
    con$ub[is.na(con$ub)]<-0.5
  }

  Res$ContMales    = sum(oc[sex==1])
  Res$ContFemales  = sum(oc[sex==2])
  Res$minCont      = min(oc)
  Res$maxContMale  = max(oc[sex==1])
  Res$maxContFemale= max(oc[sex==2])
  Res$solver    = x$solver

  valid<-TRUE
  
  cat("Checking constraints:\n")
  isOK  <- min(oc)+0.0001>=0
  isOK <- !is.na(isOK) & isOK
  valid <- valid & isOK
  cat("  min(oc) >= 0           : ", isOK, "\n", sep="")
  isOK  <- sum(oc[sex==1])>0.4999 & sum(oc[sex==1])<0.5001
  isOK <- !is.na(isOK) & isOK
  valid <- valid & isOK
  cat("  total male cont   = 0.5: ", isOK, "\n", sep="")
  isOK  <- sum(oc[sex==2])>0.4999 & sum(oc[sex==2])<0.5001
  isOK <- !is.na(isOK) & isOK
  valid <- valid & isOK
  cat("  total female cont = 0.5: ", isOK, "\n", sep="")
  equalFemaleCont <- "F" %in% names(con$ub) & (con$ub["F"] == -1)
  equalMaleCont   <- "M" %in% names(con$ub) & (con$ub["M"] == -1)
  equalFemaleCont <- equalFemaleCont | (sum(is.na(phen$ub[sex==2]))==0 & all(phen$lb[sex==2]==phen$ub[sex==2]) & sd(phen$ub[sex==2])==0)
  equalMaleCont   <- equalMaleCont   | (sum(is.na(phen$ub[sex==1]))==0 & all(phen$lb[sex==1]==phen$ub[sex==1]) & sd(phen$ub[sex==1])==0)
  if(equalMaleCont){
    isOK <- sd(oc[sex==1])==0
    isOK <- !is.na(isOK) & isOK
    valid <- valid & isOK
    cat("  males have equal cont  : ", isOK, "\n", sep="")
  }
  if(equalFemaleCont){
    isOK <- sd(oc[sex==2])==0
    isOK <- !is.na(isOK) & isOK
    valid <- valid & isOK
    cat("  females have equal cont: ", isOK, "\n", sep="")
  }
  if(!equalMaleCont){
    isOK <- all(oc[sex==1]-0.0001<=con$ub[sex==1])
    isOK <- !is.na(isOK) & isOK
    valid <- valid & isOK
    cat("  all male cont <= ub    : ", isOK, "\n", sep="")
  } 
  if(!equalFemaleCont){
    isOK <- all(oc[sex==2]-0.0001<=con$ub[sex==2])
    isOK <- !is.na(isOK) & isOK
    valid <- valid & isOK
    cat("  all female cont <= ub  : ", isOK, "\n", sep="")
  } 
  for(i in names(x$quadcon)){
    isOK  <- x$meanKin[i]-0.0001 <= x$quadcon[i]
    isOK <- !is.na(isOK) & isOK
    valid <- valid & isOK
    cat("  mean ",i," <= ub.",i,"       : ", isOK, "\n", sep="")    
  }
  
  Traits  <- colnames(phen)[-(1:2)]
  x2<-NULL;for(i in 3:ncol(phen)){x2<-c(x2,is.numeric(phen[,i]))}
  Traits  <- setdiff(Traits[x2], c("lb", "ub", "oc", "Born", "Sex"))
   for(vari in Traits){
    lb.var<- paste("lb.", vari, sep="")
    ub.var<- paste("ub.", vari, sep="")
    eq.var<- paste("eq.", vari, sep="")
    Res[, lb.var]                   = as.numeric(x$lincon[lb.var,"val"])
    Res[, paste("mean",vari,sep="")] = c(t(oc)%*%phen[,vari])
    Res[, ub.var]                   = as.numeric(x$lincon[ub.var,"val"])
    if(is.na(Res[, lb.var]))Res[, lb.var]<-as.numeric(x$lincon[eq.var,"val"])
    if(is.na(Res[, ub.var]))Res[, ub.var]<-as.numeric(x$lincon[eq.var,"val"])
    if(!is.na(Res[, lb.var])){
      isOK <- Res[, lb.var]<= c(t(oc)%*%phen[,vari])+0.0001
      isOK <- !is.na(isOK) & isOK
      valid <- valid & isOK
      cat("  mean ",vari," >= lb.",vari,"       : ", isOK, "\n", sep="")   
    }
    if(!is.na(Res[, ub.var])){
      isOK  <- Res[, ub.var]>= c(t(oc)%*%phen[,vari])-0.0001
      isOK <- !is.na(isOK) & isOK
      valid <- valid & isOK
      cat("  mean ",vari," <= ub.",vari,"       : ", isOK, "\n", sep="")   
    }
  }
  cat(" \n")
  if(is.na(valid)){valid<-FALSE}
  Res$valid <- valid
  
  if(phenAsDataTable){setDT(Res)}
  invisible(Res)
}

