
"asDefinite"<-function(K, quiet=FALSE){
  isdsCMatrix<- class(K)=="dsCMatrix"
  K <- as(K, "matrix")
  K <- 0.5*(K + t(K))
  #cat(dim(K),"\n")
  x<-base::eigen(K,only.values=TRUE, EISPACK=TRUE)$values
  if(min(x)<0.0001){
    if(!quiet){cat("Matrix has negative or zero eigenvalues. Looking for approximate solution...\n")}
    x<-base::eigen(K,only.values=FALSE, EISPACK=TRUE)
    v<-x$values
    v[v<0.0001]<-0.0001
    K<-x$vectors%*%diag(v)%*%t(x$vectors)
  }
  if(isdsCMatrix)K<-as(K,"dsCMatrix")
  K
}

"getconst"<-function(con, Traits){
  x<-names(con)[gsub("([ul]b\\.)|(eq\\.)", "", names(con)) %in% Traits]
  const<-data.frame(var=gsub("([ul]b\\.)|(eq\\.)", "", x), con=str_extract(substr(x,1,2),"[ul]b"), val=as.numeric(unlist(con[x])),stringsAsFactors = FALSE)
  rownames(const)<-x
  const$con[is.na(const$con)]<-"eq"
  const<-const[!is.na(const$val),]
  const
}


"opticont"<-function(method, K, phen, con=list(), solver="cccp", quiet=FALSE, make.definite=solver=="csdp", ...){
  gc()
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  
  if(!("ub" %in% names(con))){con$ub<-c(M=NA, F=NA)}
  if(!("lb" %in% names(con))){con$lb<-c(M=0,  F=0)}
  cKin  <- attr(K,"condProb")
  phen2 <- phen
  rownames(phen)<-phen[,1]
  phen  <- phen[,-1]
  for(i in 1:length(K)){K[[i]]<-K[[i]][rownames(phen), rownames(phen)]}
  Traits  <- colnames(phen)[-1]
  x<-NULL;for(i in 2:ncol(phen)){x<-c(x,is.numeric(phen[,i]))}
  Traits  <- setdiff(Traits[x], c("lb", "ub","oc"))

  condProbMat <- c()
  for( i in names(attributes(K)$condProb)){
    condProbMat<-c(condProbMat, attributes(K)$condProb[[i]]["f1"], attributes(K)$condProb[[i]]["f2"])
  }
  
  #cat(x,"\n")
  #cat(Traits,"\n")
  min.Kin    <- (method %in% paste("min.",names(K),    sep=""))
  min.cKin   <- (method %in% paste("min.",names(cKin), sep=""))
  min.Trait  <- (method %in% paste("min.",Traits,      sep=""))
  max.Trait  <- (method %in% paste("max.",Traits,      sep=""))
  const   <- getconst(con, Traits)
  quadcon <- getconst(con,c(names(K),names(cKin)))
  quadcon <- setNames(quadcon$val[quadcon$con=="ub"],quadcon$var[quadcon$con=="ub"])

  if(!(solver %in% c("cccp","cccp2","slsqp","alabama","csdp"))){cat("Solver not available.\n"); return(NULL)}
  if(solver=="csdp" & min.cKin){cat("Solver not suitable\n"); return(NULL)}
  
  opt<-list(...)
  if(solver %in% c("cccp", "cccp2")){
    if(!("abstol"  %in% names(opt))){opt$abstol  = (1e-06)*(10^length(quadcon)) }
    if(!("feastol" %in% names(opt))){opt$feastol = 1e-05}
    if(!("trace"   %in% names(opt))){opt$trace   = TRUE}
    if(!("stepadj" %in% names(opt))){if(min.cKin){opt$stepadj = 0.40}else{opt$stepadj = 0.90}}
  }
  
  lb     <- con$lb
  ub     <- con$ub
  if("Sex" %in% colnames(phen)){
    sex    <- phen[,"Sex"]
  }else{
    sex    <- phen[,1]
  }
  Res     <- list(parent=phen, con=con, method=method, solver=solver, quadcon=quadcon, lincon=const, cKin=cKin)
  relax   <- FALSE#length(quadcon)==0
  
  
  for(i in intersect(names(quadcon), names(cKin))){
    ad <- 0
    bd <- 0
    if("a" %in% names(cKin[[i]])){ad <- as.numeric(cKin[[i]]["a"])}
    if("b" %in% names(cKin[[i]])){bd <- as.numeric(cKin[[i]]["b"])}
    #cat("a=",ad,"\n")
    #cat("b=",bd,"\n")
    K[[i]] <- (1-quadcon[i])*(K[[cKin[[i]][2]]]-bd) + K[[cKin[[i]][1]]]
    quadcon[i] <- 1+ad
  }
  
  if(!(min.Kin | min.cKin | min.Trait | max.Trait)){cat("Method not implemented.\n");return(NULL)}
  obj.var <- substr(method,5,nchar(method))
  fun.sig <- ifelse(substr(method,1,3)=="max",-1,1)
  if(min.Kin){
    if(!quiet)cat("Objective: minimizing mean kinship ", obj.var, " in the offspring.\n",sep="")
  }else{
    if(min.cKin){
      if(!quiet)cat("Objective: minimizing conditional kinship ", obj.var, " in the offspring.\n",sep="")
    }else{
      if(!quiet)cat("Objective: ",substr(method,1,3),"imizing mean ", obj.var, " of the offspring.\n",sep="")
    }
  }
 
  if(!quiet)cat("Constraints:\n")
  for(i in rownames(const)){
    if(const[i,"con"]=="ub" & !quiet)cat("  Mean ",const[i,"var"]," in the offspring is not exceeding ", i,".\n",sep="")
    if(const[i,"con"]=="lb" & !quiet)cat("  Mean ",const[i,"var"]," in the offspring is at least ", i,".\n",sep="")
    if(const[i,"con"]=="eq" & !quiet)cat("  Mean ",const[i,"var"]," in the offspring is equal to ", i,".\n",sep="")
  }
  for(i in setdiff(names(Res$quadcon), names(cKin))){
    if(!quiet)cat("  Mean kinship ", i, " in the offspring is not exceeding ub.", i,".\n",sep="")
  }
  for(i in intersect(names(Res$quadcon), names(cKin))){
    if(!quiet)cat("  Conditional kinship ", i, " in the offspring is not exceeding ub.", i,".\n",sep="")
  }
  
  if(length(lb)==2 & !is.null(names(lb))){
    if(!quiet)cat("  Minimum contribution of males to offspring is defined.\n",sep="")                
    if(!quiet)cat("  Minimum contribution of females to offspring is defined\n",sep="")
    lb<-lb[sex]
  } else{
    if(!quiet)cat("  Minimum contributions of animals to the offspring were defined\n",sep="")
  }
    
  if(length(ub)==2 & !is.null(names(ub))){
    if(is.na(ub[1])){
      if(!quiet)cat("  Number of offspring of males is not limited.\n")              
      if(!quiet)cat("    (Thus, the maximum contribution per male to the offspring is 0.5).\n")              
    }else{
      if(ub[1]==-1){
        if(!quiet)cat("  All males have equal contributions to the offspring.\n")
        if(!quiet)cat("    (Thus, no optimization is done for the males).\n")                      
      }else{
        if(!quiet)cat("  Maximum contribution of males to offspring is provided.\n",sep="")                
      }
    }
    if(is.na(ub[2])){
      if(!quiet)cat("  Number of offspring of femmales is not limited.\n")              
      if(!quiet)cat("    (Thus, the maximum contribution per female to the offspring is 0.5).\n")              
    }else{
      if(ub[2]==-1){
        if(!quiet)cat("  All females have equal contributions to the offspring.\n")
        if(!quiet)cat("    (Thus, no optimization is done for the females).\n")                      
      }else{
        if(!quiet)cat("  Maximum contribution of females to offspring is provided.\n",sep="")
      }
    }
    ub<-ub[sex]
    equalMaleCont   <- sum(is.na(ub[sex==1]))==0 & (all(ub[sex==1]==-1))
    equalFemaleCont <- sum(is.na(ub[sex==2]))==0 & (all(ub[sex==2]==-1))
  }else{
    equalMaleCont   <- sum(is.na(ub[sex==1]))==0 & (all(ub[sex==1]==-1))
    equalFemaleCont <- sum(is.na(ub[sex==2]))==0 & (all(ub[sex==2]==-1))
    if(equalMaleCont){
      if(!quiet)cat("  All males have equal contributions to the offspring.\n")
      if(!quiet)cat("    (Thus, no optimization is done for the males).\n")                      
    }else{
      if(!quiet)cat("  Maximum contributions of males to the offspring were provided.\n",sep="")
    } 
    if(equalFemaleCont){
      if(!quiet)cat("  All females have equal contributions to the offspring.\n")
      if(!quiet)cat("    (Thus, no optimization is done for the females).\n")                      
    }else{
      if(!quiet)cat("  Maximum contributions of females to the offspring were provided.\n",sep="")
    }
  }
  if(!quiet)cat("  The total genetic contribution of   males to offspring is 0.5.\n")
  if(!quiet)cat("  The total genetic contribution of females to offspring is 0.5.\n")
  
  lb[is.na(lb)]<-0.0
  ub[is.na(ub)]<-0.5

  for(v in unique(const$var)){
    m <- mean(phen[,v])
    s <- sd(phen[,v])
    phen[,v] <- (phen[,v]-m)/s
    const[const$var==v,"val"]<-(const[const$var==v,"val"]-m)/s
  }
  
  useChol <- (length(quadcon)>0) & solver == "cccp"
  
  if(equalMaleCont){
    ub[sex==1] <- rep(1/(2*sum(sex==1)),sum(sex==1))
    lb[sex==1] <- ub[sex==1]
    }
  if(equalFemaleCont){
    ub[sex==2] <- rep(1/(2*sum(sex==2)),sum(sex==2))
    lb[sex==2] <- ub[sex==2]
  }
  Res$parent$lb <- lb
  Res$parent$oc <- lb
  Res$parent$ub <- ub
  isC <- lb==ub
  isV <- lb!=ub
  nV  <- sum(isV)
  nC  <- sum(isC)
  if(nV==0){
    Res$parent$oc  <- lb
    kinNames <- setdiff(names(K),union(names(cKin),condProbMat))
    Res$meanKin <- setNames(rep(NA,length(kinNames)), kinNames)
    for(i in kinNames){
      Res$meanKin[i]<-c(as(t(Res$parent$oc)%*%K[[i]]%*%Res$parent$oc,"matrix"))
    }
    if(!is.null(cKin)){
      for(i in names(cKin)){
        f1 <- cKin[[i]]["f1"]
        f2 <- cKin[[i]]["f2"]
        Res$meanKin[i]<- c(as(1-(1-t(Res$parent$oc)%*%(K[[f1]])%*%(Res$parent$oc))/(t(Res$parent$oc)%*%(K[[f2]])%*%(Res$parent$oc)),"matrix"))
      }
    }
    class(Res)<-"opticont"
    Res$parent <- data.frame(Indiv=rownames(Res$parent), Res$parent, stringsAsFactors = FALSE)
    if(phenAsDataTable){setDT(Res$parent)}
    return(Res)
  }
  b  <- matrix(c(0.5-sum(lb[isC & sex==1]), 0.5 - sum(lb[isC & sex==2])), ncol=1)
  cC <- lb[isC]
  X  <- ub
  X[sex==1 & isV] <- b[1,1]/(sum(sex==1 & isV)) 
  X[sex==2 & isV] <- b[2,1]/(sum(sex==2 & isV))
  X <- X[isV]
  if(length(table(sex[isV]))==1){
    A  <- matrix(1, nrow=1, ncol=nV)
    b  <- b[b[,1]>0.0000001,, drop=FALSE]
  }else{
    A  <- t(model.matrix(~as.factor(sex[isV])-1))
  }
  
  P <- NULL
  q <- NULL
  if(min.Kin){
    P <- as(K[[obj.var]][isV, isV], "matrix")
    if(make.definite){P <- asDefinite(P, quiet=quiet)}
    q <- rep(0,nV)
    if(nC>0){q <-c(as(K[[obj.var]][isV, isC],"matrix")%*%cC)}
  }

  f0 <- NULL
  g0 <- NULL
  h0 <- NULL
  if(min.cKin){
    f1 <- cKin[[obj.var]][1]
    f2 <- cKin[[obj.var]][2]
    ad <- 0
    bd <- 0
    if("a" %in% names(cKin[[obj.var]])){ad <- as.numeric(cKin[[obj.var]]["a"])}
    if("b" %in% names(cKin[[obj.var]])){bd <- as.numeric(cKin[[obj.var]]["b"])}
    
  #fU <- cKin[[obj.var]][3]
    AB <- as((K[[f1]][isV, isV]),"matrix")#asDefinite
    AN <- as((K[[f2]][isV, isV]),"matrix")#asDefinite
 #   if(make.definite){
#      AB <- asDefinite(AB)
#      AN <- asDefinite(AN)
#    }
    uB <- rep(0,nV)
    uN <- rep(0,nV)
    CB <- 0 - ad
    CN <- 0 - bd
    if(nC!=0){
      uB <- 2*as(K[[f1]][isV, isC],"matrix")%*%cC
      uN <- 2*as(K[[f2]][isV, isC],"matrix")%*%cC
      CB <- t(cC)%*%as(K[[f1]][isC, isC],"matrix")%*%cC - ad
      CN <- t(cC)%*%as(K[[f2]][isC, isC],"matrix")%*%cC - bd
    }
    f0 <- function(x){c(1+(t(x)%*%AB%*%x+t(uB)%*%x+CB-1)/(t(x)%*%AN%*%x+t(uN)%*%x+CN))}
    g0 <- function(x){
      hx <- c(t(x)%*%AN%*%x+t(uN)%*%x+CN)
      gx <- c(t(x)%*%AB%*%x+t(uB)%*%x+CB)
      ax <- 1/hx
      bx <- (1-gx)/(hx^2)
      c(ax*(2*AB%*%x+uB) + bx*(2*AN%*%x+uN))
    }
    h0 <- function(x){
      hx <- c(t(x)%*%AN%*%x+t(uN)%*%x+CN)
      gx <- c(t(x)%*%AB%*%x+t(uB)%*%x+CB)
      ax <- 1/hx
      bx <- (1-gx)/(hx^2)
      dax <- - (2*AN%*%x+uN)/hx^2
      dbx <- - (2*AB%*%x+uB)/hx^2 - 2*(1-gx)*(2*AN%*%x+uN)/hx^3
      ax*(2*AB)+dax%*%t(2*AB%*%x+uB)+bx*(2*AN)+dbx%*%t(2*AN%*%x+uN)
    }
    con2<-Res$con
    names(con2)[names(con2)=="ub.MC"]<-"lb.MC"
    for(i in intersect(names(con2), paste("ub.", names(K), sep="")))con2[[i]]<-con2[[i]]-0.001
    
    con2$lb <- Res$parent$lb
    con2$ub <- Res$parent$ub
    x  <- 1-diag(as(K[[f2]],"matrix"))
    sm <- sum(sex==1 & con2$lb < con2$ub)
    sf <- sum(sex==2 & con2$lb < con2$ub)
    gm <- min(150, sm%/%2)
    gf <- min(150, sf%/%2)
    xm <- sort(x[sex==1 & con2$lb < con2$ub])[gm]
    xf <- sort(x[sex==2 & con2$lb < con2$ub])[gf]
    weg <- rep(FALSE, length(sex))
    if(sm>5){weg <- weg | (sex==1 & con2$lb < con2$ub & x>xm)}
    if(sf>5){weg <- weg | (sex==2 & con2$lb < con2$ub & x>xf)}
    con2$lb[weg]<-0
    con2$ub[weg]<-0
    #X <- opticont(paste("min.", f1, sep=""), K=K, phen=Res$parent, con=con2, solver=Res$solver, abstol=abstol, feastol = feastol, quiet=TRUE, make.definite=TRUE, ...)$parent$oc[isV]
    X <-  opticont(paste("min.", f1, sep=""), K=K, phen=phen2, con=con2, solver=Res$solver, quiet=TRUE, make.definite=TRUE, ...)$parent$oc[isV]
    b <- A%*%X
    }
  
  
  g <- setNames(vector("list", length(quadcon)), names(quadcon))
  F <- setNames(vector("list", length(quadcon)), names(quadcon)) 
  for(i in names(quadcon)){
    if(useChol){
      cat("Computing Cholesky decomposition for ",i,"...")
      F[[i]] <- K[[i]][isV, isV]
      F[[i]] <- asDefinite(F[[i]], quiet=quiet)
      F[[i]] <- chol(F[[i]])
      F[[i]] <- as(F[[i]],"matrix")
      cat("finished\n")
      g[[i]] <- rep(0,nV)
      if(nC>0){
        g[[i]]     <- as(solve(t(F[[i]]))%*%(K[[i]][isV,isC])%*%cC,"matrix")
        quadcon[i] <- quadcon[i] - c(t(cC)%*%as(K[[i]][isC, isC],"matrix")%*%cC)+t(g[[i]])%*%g[[i]]
      }
    }else{
      F[[i]] <- as(K[[i]][isV, isV],"matrix")
      if(make.definite){F[[i]] <- asDefinite(F[[i]], quiet=quiet)}
      g[[i]] <- rep(0,nV)
      if(nC>0){
        g[[i]]     <- as(K[[i]][isV, isC],"matrix")%*%cC
        quadcon[i] <- quadcon[i] - c(t(cC)%*%as(K[[i]][isC, isC],"matrix")%*%cC)
      }
    }
  }
 
  lb <- lb[isV]
  ub <- ub[isV]
  
  G <- NULL
  h <- NULL
  for(i in rownames(const)){
    if(nC>0){const[i,"val"] <- (const[i,"val"]-cC%*%phen[isC, const[i,"var"]])/sum(b)}
    if(const[i,"con"]=="lb"){
      G <- rbind(G,-(phen[isV, const[i,"var"]]-const[i,"val"]))
      h <- c(h, 0)    
    }
    if(const[i,"con"]=="ub"){
      G <- rbind(G, phen[isV, const[i,"var"]]-const[i,"val"])
      h <- c(h, 0)    
    }
  }
  
  if(relax){
    G <- rbind(G,-A)
    h <- c(h,-b)
    A <- NULL
    b <- NULL
  }
  
  for(i in rownames(const)[const$con=="eq"]){
    A <- rbind(A, phen[isV, const[i,"var"]]-const[i,"val"])
    b <- c(b, 0)    
  }

  if(min.Kin){
    Res$parent$oc[isV] <- opticontx(X=X, P=as(P,"matrix"), q=q,                lb=lb, ub=ub, A=A, b=b, G=G, h=h,  F=F, g=g, d=NULL, quadcon=quadcon, isChol=useChol, solver=solver, quiet=quiet, opt=opt)
  }else{
    if(min.cKin){
      Res$parent$oc[isV] <- opticontx(X=X, f0=f0, g0=g0,  h0=h0,               lb=lb, ub=ub, A=A, b=b, G=G, h=h,  F=F, g=g, d=NULL, quadcon=quadcon, isChol=useChol, solver=solver, quiet=quiet, opt=opt)
    }else{
      Res$parent$oc[isV] <- opticontx(X=X, P=NULL, q=fun.sig*phen[isV,obj.var],lb=lb, ub=ub, A=A, b=b, G=G, h=h,  F=F, g=g, d=NULL, quadcon=quadcon, isChol=useChol, solver=solver, quiet=quiet, opt=opt)
    }
  }
  kinNames<-setdiff(names(K),union(names(cKin),condProbMat))
  Res$meanKin<-setNames(rep(NA,length(kinNames)),kinNames)
  for(i in kinNames){
    Res$meanKin[i]<-c(as(t(Res$parent$oc)%*%K[[i]]%*%Res$parent$oc,"matrix"))
  }
  if(!is.null(cKin)){
    for(i in names(cKin)){
      f1 <- cKin[[i]]["f1"]
      f2 <- cKin[[i]]["f2"]
      Res$meanKin[i]<- c(as(1-(1-t(Res$parent$oc)%*%(K[[f1]])%*%(Res$parent$oc))/(t(Res$parent$oc)%*%(K[[f2]])%*%(Res$parent$oc)),"matrix"))
    }
  }
  Res$parent <- data.frame(Indiv=rownames(Res$parent), Res$parent, stringsAsFactors = FALSE)
  if(phenAsDataTable){setDT(Res$parent)}
  class(Res)<-"opticont"
  Res
}