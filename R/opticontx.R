
"opticontx"<-function(X=NULL, P=NULL, q=NULL, f0=NULL, g0=NULL, h0=NULL, lb=NULL, ub=NULL, A=NULL, b=NULL, G=NULL, h=NULL, F=NULL, g=NULL, d=NULL, quadcon=NULL, isChol=FALSE, solver="cccp", opt=NULL, quiet=FALSE){
  
  if(solver %in% c("cccp", "cccp2")){
    if(! "maxiters" %in% names(opt)){opt$maxiters = 100L}
    if(! "abstol"   %in% names(opt)){opt$abstol  = 1e-07}
    if(! "reltol"   %in% names(opt)){opt$reltol  = 1e-06}
    if(! "feastol"  %in% names(opt)){opt$feastol = 1e-07}
    if(! "stepadj"  %in% names(opt)){opt$stepadj = 0.95}
    if(! "beta"     %in% names(opt)){opt$beta    = 0.5}
    if(! "trace"    %in% names(opt)){opt$trace   = TRUE}
    params.opt <- c("maxiters", "abstol", "reltol", "feastol", "stepadj", "beta", "trace")
    opt    <- opt[intersect(names(opt), params.opt)]
    optctrl<- do.call(cccp::ctrl, opt)
  }

  if(solver == "alabama"){
    params.outer  <- c("lam0", "sig0", "eps", "itmax", "ilack.max", "trace", "method", "NMinit", "i.scale", "e.scale", "kkt2.check")
    params.optim  <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    optctrl.outer <- opt[intersect(names(opt), params.outer)]
    optctrl.optim <- opt[intersect(names(opt), params.optim)]
    if(! "lam0"     %in% names(opt)){opt$lam0     = "default"}
    if(! "sig0"     %in% names(opt)){opt$sig0     = "default"}
    if(! "eps"      %in% names(opt)){opt$eps      = "default"}
    if(! "itmax"    %in% names(opt)){opt$itmax    = "default"}
    if(! "ilack.max"%in% names(opt)){opt$ilack.max= "default"}
    if(! "trace"    %in% names(opt)){opt$trace    = "default"}
    if(! "method"   %in% names(opt)){opt$method   = "default"}
    if(! "NMinit"   %in% names(opt)){opt$NMinit   = "default"}
    if(! "i.scale"   %in% names(opt)){opt$i.scale   = "default"}
    if(! "e.scale"   %in% names(opt)){opt$e.scale   = "default"}
    if(! "kkt2.check"   %in% names(opt)){opt$kkt2.check   = "default"}
    if(! "fnscale"  %in% names(opt)){opt$fnscale  = "default"}
    if(! "parscale" %in% names(opt)){opt$parscale = "default"}
    if(! "ndeps"    %in% names(opt)){opt$ndeps    = "default"}
    if(! "maxit"    %in% names(opt)){opt$maxit    = "default"}
    if(! "abstol"   %in% names(opt)){opt$abstol   = "default"}
    if(! "reltol"   %in% names(opt)){opt$reltol   = "default"}
    if(! "alpha"    %in% names(opt)){opt$alpha    = "default"}
    if(! "beta"     %in% names(opt)){opt$beta     = "default"}
    if(! "gamma"    %in% names(opt)){opt$gamma    = "default"}
    if(! "REPORT"   %in% names(opt)){opt$REPORT   = "default"}
    if(! "type"     %in% names(opt)){opt$type     = "default"}
    if(! "lmm"      %in% names(opt)){opt$lmm      = "default"}
    if(! "factr"    %in% names(opt)){opt$factr    = "default"}
    if(! "pgtol"    %in% names(opt)){opt$pgtol    = "default"}
    if(! "temp"     %in% names(opt)){opt$temp     = "default"}
    if(! "tmax"     %in% names(opt)){opt$tmax     = "default"}
    opt <- opt[unique(c(params.outer, params.optim))]
    }

  if(solver == "slsqp"){
    if(! "nl.info"  %in% names(opt)){opt$nl.info  <- FALSE}
    if(! "stopval"  %in% names(opt)){opt$stopval  <- -Inf}
    if(! "xtol_rel" %in% names(opt)){opt$xtol_rel <- 1e-06}
    if(! "maxeval"  %in% names(opt)){opt$maxeval  <- 1000}
    if(! "ftol_rel" %in% names(opt)){opt$ftol_rel <- 0.0}
    if(! "ftol_abs" %in% names(opt)){opt$ftol_abs <- 0.0}
    if(! "check_derivatives" %in% names(opt)){opt$check_derivatives <- FALSE}
    optlist <- opt[intersect(names(opt), c("stopval", "xtol_rel", "maxeval", "ftol_rel", "ftol_abs", "check_derivatives"))]
  }
  
  if(solver=="csdp"){
    if(! "axtol"       %in% names(opt)){opt$axtol = 1e-08}
    if(! "atytol"      %in% names(opt)){opt$atytol = 1e-08}
    if(! "objtol"      %in% names(opt)){opt$objtol = 1e-05}
    if(! "pinftol"     %in% names(opt)){opt$pinftol = 1e+20}
    if(! "dinftol"     %in% names(opt)){opt$dinftol = 1e+20}
    if(! "maxiter"     %in% names(opt)){opt$maxiter = 100}
    if(! "minstepfrac" %in% names(opt)){opt$minstepfrac = 0.8}
    if(! "maxstepfrac" %in% names(opt)){opt$maxstepfrac = 0.97}
    if(! "minstepp"    %in% names(opt)){opt$minstepp = 1e-20}
    if(! "minstepd"    %in% names(opt)){opt$minstepd = 1e-20}
    if(! "usexzgap"    %in% names(opt)){opt$usexzgap = 1}
    if(! "tweakgap"    %in% names(opt)){opt$tweakgap = 0}
    if(! "affine"      %in% names(opt)){opt$affine = 0}
    if(! "printlevel"  %in% names(opt)){opt$printlevel = 1} 
    if(! "perturbobj"  %in% names(opt)){opt$perturbobj = 1}
    if(! "fastmode"    %in% names(opt)){opt$fastmode = 0}
    params.opt <- c("axtol", "atytol", "objtol", "pinftol", "dinftol", "maxiter", "minstepfrac", "maxstepfrac", "minstepp", "minstepd", "usexzgap", "tweakgap", "affine", "printlevel", "perturbobj", "fastmode")
    opt<-opt[intersect(names(opt), params.opt )]
    }
  
  
  if(!quiet)cat("\nUsing solver ", solver, " with parameters: ",  paste(names(rapply(opt,c)),rapply(opt,c),sep="=", collapse=", "), "\n\n", sep="")
  
  noBounds <- (solver %in% c("cccp","cccp2","alabama", "csdp"))
  if(noBounds){
    N <- length(lb)
    G <- rbind(-diag(N),G)
    h <- c(-lb,     h)
    if(sum(ub<0.5)>0){
      G <- rbind(diag(N)[ub<0.5,], G)
      h <- c(ub[ub<0.5], h)
    }  
  }

  if(solver=="csdp"){
    N <- length(lb)
    nblocks<-length(quadcon)+(!is.null(A))+(!is.null(G))+(!is.null(P))
    Y<-vector("list", N+(!is.null(P)))
    C<-vector("list", nblocks)
    K<-list(type=rep('', nblocks), size=rep(0,nblocks))
    is<-1
    if(!is.null(A)){C[[is]] <- c(c(b),c(-b));K$type[is]<-"l"; K$size[is]<-length(C[[is]]); is<-is+1}
    if(!is.null(G)){C[[is]] <- c(-h); K$type[is]<-"l"; K$size[is]<-length(C[[is]]); is<-is+1}
    for(j in names(quadcon)){
      C[[is]] <- as(Matrix::bdiag((-1)*ginv(F[[j]]), (-1)*quadcon[j]),"matrix")
      K$type[is]<-"s"; 
      K$size[is]<-N+1;
      is <- is+1
    }
    if(!is.null(P)){
      C[[is]] <- as(Matrix::bdiag((-1)*ginv(P),0),"matrix")
      K$type[is]<-"s"; 
      K$size[is]<-N+1;
      is <- is+1 
    }
    for(i in 1:N){
      Y[[i]]<-vector("list", nblocks)
      is<-1
      if(!is.null(A)){Y[[c(i,is)]] <- c(A[,i],-A[,i]); is<-is+1}
      if(!is.null(G)){Y[[c(i,is)]] <- -G[,i]; is<-is+1}
      for(j in names(quadcon)){
        Y[[c(i,is)]]<-simple_triplet_sym_matrix(i=c(i,N+1),j=c(N+1,N+1),v=c(1,-2*g[[j]][i]))
        is<-is+1
      }
      if(!is.null(P)){
        Y[[c(i,is)]]<-simple_triplet_sym_matrix(i=c(i),j=c(N+1),v=c(1))
        is<-is+1
      }
    }
    if(!is.null(P)){
      Y[[N+1]]<-vector("list", nblocks)
      is<-1
      if(!is.null(A)){Y[[c(N+1,is)]] <- 0*c(A[,1],-A[,1]); is<-is+1}
      if(!is.null(G)){Y[[c(N+1,is)]] <- 0*G[,1]; is<-is+1}
      for(j in names(quadcon)){
        Y[[c(N+1,is)]]<-simple_triplet_sym_matrix(i=c(1,N+1),j=c(N+1,1),v=c(0,0))
        is<-is+1
      }
      Y[[c(N+1,is)]]<-simple_triplet_sym_matrix(i=c(N+1),j=c(N+1),v=c(1))
      is<-is+1
    }
    if(is.null(P)){
      B<-q
    }else{
      B<-c(rep(0,N),1)
    }
    return(Rcsdp::csdp(C,Y,B,K,control= opt)$y[1:N])
  }
  
  
  if(solver %in% c("alabama", "slsqp")){
    if(is.null(f0)){
      if(is.null(P)){
        f0<-function(x){c(q%*%x)}
        g0<-function(x){c(q)}        
      }else{
        f0<-function(x){c(0.5*t(x)%*%P%*%x)+c(t(q)%*%x)}
        g0<-function(x){c(P%*%x+q)}        
      }
    }
    hin    <- NULL
    hin.jac<- NULL
    if(length(quadcon)+length(h)>0){
      hin     <- function(x){}
      hin.jac <- function(x){}
      Body     <- NULL
      Body.jac <- NULL
      if(length(h)>0){
        Body     <- "c(G%*%x-h)"
        Body.jac <- "G"
      }
      for(i in names(quadcon)){
        next.hin <- paste("c(t(x)%*%F[['",i,"']]%*%x+2*t(g[['",i,"']])%*%x-quadcon['",i,"'])", sep="")
        next.jac <- paste("c(2*t(x)%*%F[['",i ,"']]+2*t(g[['",i,"']]))", sep="")
        if(is.null(Body)){    Body     <- next.hin}else{Body     <- paste(Body,     next.hin, sep=",")}
        if(is.null(Body.jac)){Body.jac <- next.jac}else{Body.jac <- paste(Body.jac, next.jac, sep=",")}
      }
      Body         <- paste("c(",     Body,     ")", sep="")
      Body.jac     <- paste("rbind(", Body.jac, ")", sep="")
      Body         <- paste("(-1)*",     Body, sep="")
      Body.jac     <- paste("(-1)*", Body.jac, sep="")
      
      body(hin)    <- parse(text=Body)
      body(hin.jac)<- parse(text=Body.jac)
    }
    heq    <- NULL
    heq.jac<- NULL
    if(length(b)>0){
      heq     <- function(x){c(A%*%x-b)}
      heq.jac <- function(x){A}
    }
    if(solver=="alabama"){
      return(alabama::auglag(par=X, fn=f0, gr=g0, hin=hin, hin.jac=hin.jac, heq=heq, heq.jac=heq.jac, control.outer=optctrl.outer, control.optim=optctrl.optim)$par)
    }
    if(solver=="slsqp"){
      return(nloptr::slsqp(x0=X, fn=f0, gr=g0, lower=lb, upper=ub, hin=hin, hinjac=hin.jac, heq=heq, heqjac=heq.jac, nl.info=opt$nl.info, control=optlist)$par)
    }
  }
    
  if(solver=="cccp"){
    cList<-setNames(vector("list", length(quadcon)+1), c("nnoc1", names(quadcon)))
    cList[["nnoc1"]] <- nnoc(G=G, h=c(h))
    for(i in names(quadcon)){
      cList[[i]]<-socc(F=F[[i]], g=c(g[[i]]), d=0*c(g[[i]]), f=c(sqrt(quadcon[i])))
    }
    if(is.null(f0)){
      return(getx(cccp(P=P, q=c(q), A=A, b=c(b), cList=cList, optctrl=optctrl))) 
    }else{
      return(getx(cccp(f0=f0, g0=g0, h0=h0, x0=X, A=A, b=c(b), cList=cList, optctrl=optctrl))) 
    }
  }
  
  if(solver=="cccp2"){
    nnoc1<-nnoc(G=G, h=c(h))   
    if(length(quadcon)>0){
      nlfList<-setNames(vector("list", length(quadcon)), names(quadcon))
      nlgList<-setNames(vector("list", length(quadcon)), names(quadcon))
      nlhList<-setNames(vector("list", length(quadcon)), names(quadcon))
      for(i in names(quadcon)){
        nlfList[[i]] <- function(x){}
        nlgList[[i]] <- function(x){}
        nlhList[[i]] <- function(x){}
        body(nlfList[[i]]) <- parse(text=paste("c(t(x)%*%F[['",i,"']]%*%x+2*t(g[['",i,"']])%*%x-quadcon['",i,"'])",sep=""))
        body(nlgList[[i]]) <- parse(text=paste("c(2*t(x)%*%F[['",i,"']]+2*t(g[['",i,"']]))",sep=""))
        body(nlhList[[i]]) <- parse(text=paste("2*F[['",i,"']]",sep=""))
      }
      gc()
      myf<-list()
      myg<-list()
      myh<-list()
      for(i in names(quadcon)){
        if(nlfList[[i]](X)>=0 & length(myf)==0){X<-1.0*(0.00001*X+0.99999*getx(cccp(P=F[[i]], q=c(g[[i]]), A=NULL, b=NULL, cList=list(nnoc(G=rbind(G,-A), h=c(h,-b))), optctrl=optctrl))[,1])}
        if(nlfList[[i]](X)>=0 & length(myf)>0){ X<-1.0*(0.00001*X+0.99999*getx(cccp(P=F[[i]], q=c(g[[i]]), A=A, b=c(b), cList=list(nnoc(G=G, h=h)), nlfList=myf, nlgList=myg, nlhList=myh, x0=X, optctrl=optctrl))[,1])}
        myf[[i]]<-nlfList[[i]]
        myg[[i]]<-nlgList[[i]]
        myh[[i]]<-nlhList[[i]]
        gc()
        #cat("Initial values have f.",i,"(X)=",nlfList[[i]](X),"\n",sep="")
        if(nlfList[[i]](X)>0){cat("No solution exists.\n");return(X)}
      }
      gc()
      if(is.null(f0)){
        return(getx(cccp(P=P, q=c(q),         A=A, b=c(b), cList=list(nnoc1), nlfList=nlfList, nlgList=nlgList, nlhList=nlhList,x0=X,optctrl=optctrl))[,1])
      }else{
        return(getx(cccp(f0=f0, g0=g0, h0=h0, A=A, b=c(b), cList=list(nnoc1), nlfList=nlfList, nlgList=nlgList, nlhList=nlhList,x0=X,optctrl=optctrl))[,1])
      }
    }else{
      gc()
      if(is.null(f0)){
        return(getx(cccp(P=P, q=c(q),         A=A, b=c(b), cList=list(nnoc1),       optctrl=optctrl))[,1])
      }else{
        return(getx(cccp(f0=f0, g0=g0, h0=h0, A=A, b=c(b), cList=list(nnoc1), x0=X, optctrl=optctrl))[,1])
      }
    }      
  }
  

}