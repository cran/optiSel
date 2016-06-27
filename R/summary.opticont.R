
"summary.opticont"<-function(object, ...){
  x <- object
  phen   <- x$parent
  oc  <- phen[,"oc"]
  con <- x$con
  sex <- phen[,1]
  obj.var <- substr(x$method,5,nchar(x$method))
  min.Kin <- (x$method %in% paste("min.",names(x$meanKin),sep=""))

  VarName<-deparse(substitute(object))[1]
  if(str_detect(VarName,"\\(")){VarName<-""}
  
  Res<-data.frame(VarName = VarName, stringsAsFactors = FALSE)
  rownames(Res)<-Res$VarName
  Res$method= x$method
  
  for(i in names(x$cKin)){
    f1 <- x$cKin[[i]][1]
    f2 <- x$cKin[[i]][2]
    x$meanKin[i]<- 1-c((1-x$meanKin[f1])/(x$meanKin[f2]))
    }
  
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

  for(i in names(x$meanKin)){
    Res[, paste("Div.",i,sep="")]<-1-Res[,i]
  }

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
  valid <- valid & isOK
  cat("  min(oc) >= 0           : ", isOK, "\n", sep="")
  isOK  <- sum(oc[sex==1])>0.4999 & sum(oc[sex==1])<0.5001
  valid <- valid & isOK
  cat("  total male cont   = 0.5: ", isOK, "\n", sep="")
  isOK  <- sum(oc[sex==2])>0.4999 & sum(oc[sex==2])<0.5001
  valid <- valid & isOK
  cat("  total female cont = 0.5: ", isOK, "\n", sep="")
  equalMaleCont   <- (sum(is.na(con$ub[sex==1]))==0 & sum(con$ub[sex==1])<0.5)
  equalFemaleCont <- (sum(is.na(con$ub[sex==2]))==0 & sum(con$ub[sex==2])<0.5)
  if(equalMaleCont){
    isOK <- sd(oc[sex==1])==0
    valid <- valid & isOK
    cat("  males have equal cont  : ", isOK, "\n", sep="")
  }
  if(equalFemaleCont){
    isOK <- sd(oc[sex==2])==0
    valid <- valid & isOK
    cat("  females have equal cont: ", isOK, "\n", sep="")
  }
  if(!equalMaleCont){
    isOK <- all(oc[sex==1]-0.0001<=con$ub[sex==1])
    valid <- valid & isOK
    cat("  all male cont <= ub    : ", isOK, "\n", sep="")
  } 
  if(!equalFemaleCont){
    isOK <- all(oc[sex==2]-0.0001<=con$ub[sex==2])
    valid <- valid & isOK
    cat("  all female cont <= ub  : ", isOK, "\n", sep="")
  } 
  for(i in names(x$quadcon)){
    isOK  <- x$meanKin[i]-0.0001 <= x$quadcon[i]
    valid <- valid & isOK
    cat("  mean ",i," <= ub.",i,"       : ", isOK, "\n", sep="")    
  }
  
  x2<-apply(phen[,-1],2,is.numeric)
  Traits <- setdiff(names(x2)[x2], c("lb", "ub","oc"))
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
      isOK  <- Res[, lb.var]<= c(t(oc)%*%phen[,vari])+0.0001
      valid <- valid & isOK
      cat("  mean ",vari," >= lb.",vari,"       : ", isOK, "\n", sep="")   
    }
    if(!is.na(Res[, ub.var])){
      isOK  <- Res[, ub.var]>= c(t(oc)%*%phen[,vari])-0.0001
      valid <- valid & isOK
      cat("  mean ",vari," <= ub.",vari,"       : ", isOK, "\n", sep="")   
    }
  }
  cat(" \n")
  Res$valid <- valid
  invisible(Res)
}

