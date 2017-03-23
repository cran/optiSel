
"opticomp"<-function(f, Breed, obj.fun="NGD", lb=NULL, ub=NULL, ...){
  X <- aggregate(f, list(Breed), mean)
  rownames(X)<-X[,"Group.1"]
  X <- t(as.matrix(X[,-1]))
  X <- aggregate(X, list(Breed), mean)
  rownames(X)<-X[,"Group.1"]
  f <- as.matrix(X[,-1])
  
  l <- rep(0,ncol(f))
  u <- rep(1,ncol(f))
  names(l)<-colnames(f)
  names(u)<-colnames(f)
  if(!is.null(lb)){l[names(lb)]<-lb}
  if(!is.null(ub)){u[names(ub)]<-ub}
  H <- matrix(c(f),ncol=ncol(f),nrow=nrow(f))
  gc()
  # maximize neutral gene diversity
  if(obj.fun=="NGD"){
    Res   <- solve.QP(Dmat=2*H,dvec=rep(0,nrow(H)),Amat=cbind(diag(nrow(H)),-diag(nrow(H)),1),bvec=c(l,-u,1), ...)
    bc    <- setNames(Res$solution, colnames(f))
    value <- (1-Res$value)
  }
  #maximize neutral trait diversity
  if(obj.fun=="NTD"){
    F    <- diag(f)
    eins <- rep(1,nrow(f))
    Dmat <- 2*(eins%*%t(eins)-(F%*%t(eins)-2*f+eins%*%t(F)))
    dvec <- -F
    Res  <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=cbind(diag(nrow(H)),-diag(nrow(H)),1),bvec=c(l,-u,1), ...)
    bc   <- setNames(Res$solution, colnames(f))
    value<- t(bc)%*%(eins-F)+t(bc)%*%(F%*%t(eins)-2*f+eins%*%t(F))%*%bc
  }
   
  #Compute genetic distance between breeds
  A <- matrix(diag(f),ncol=ncol(f),nrow=nrow(f),byrow=TRUE)
  Dist <- sqrt((A+t(A))/2-f)
  
  list(bc=bc, value=value, f=f,  Dist=Dist)
}

