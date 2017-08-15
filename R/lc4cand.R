
"lc4cand" <- function(phen, bc, condf, Traits, eqConSexes){
  phen[is.na(phen$Sex),"Sex"]<-"undefined"
  Breeds <- names(eqConSexes)
  
  if(length(Breeds)==1){
    A <- matrix(1,nrow=1,ncol=nrow(phen))
    rownames(A) <- Breeds
  }else{
    A <- t(model.matrix(~Breed-1, data=phen))
    rownames(A) <- str_replace(rownames(A), "Breed", "")
  }
  
  rnamesA     <- rownames(A)
  rownames(A) <- paste0("bc.", rnamesA)
  
  eqConSexes <- eqConSexes[rnamesA]
  
  use  <- rep(TRUE, nrow(A))
  dir  <- rep("==", nrow(A))
  val  <- bc[rnamesA]
  name <- rownames(A)
  
  if(!all(phen$Sex=="undefined")){
    if(length(Breeds)==1){
      A2 <- t(model.matrix(~Sex-1, data=phen))
      rownames(A2) <- str_replace(rownames(A2), "Sex",   "")
      rownames(A2) <- paste(rnamesA, rownames(A2), sep=":")
    }else{
      A2 <- t(model.matrix(~Breed:Sex-1, data=phen))
      rownames(A2) <- str_replace(rownames(A2), "Breed", "")
      rownames(A2) <- str_replace(rownames(A2), "Sex",   "")
    }

    A2 <- A2[paste0(rnamesA,":female"),,drop=FALSE] - A2[paste0(rnamesA,":male"),,drop=FALSE]
    rownames(A2) <- paste0("scd.",str_replace(rownames(A2), ":female", ""))
    A2 <- A2[paste0("scd.", rnamesA[eqConSexes]),,drop=FALSE]

    if(nrow(A2)>0){
      A    <- rbind(A, A2)
      use  <- c(use, rep(TRUE, nrow(A2)))
      dir  <- c(dir, rep("==", nrow(A2)))
      val  <- c(val, rep(0,    nrow(A2)))
      name <- c(name, rownames(A2))
    }
  }
  
  alpha <- ifelse(Traits %in% condf$var, condf[Traits, "val"], -Inf)
  A2 <- t(phen[, Traits] - matrix(alpha, nrow=nrow(phen), ncol=length(alpha), byrow=TRUE))
  A2[is.na(A2)] <- 0
  
  if(length(Traits)>0){
    A    <- rbind(A, A2)
    use  <- c(use,  Traits %in% condf$var)
    dir  <- c(dir,  ifelse(Traits %in% condf$var, condf[Traits, "dir"], ">="))
    val  <- c(val,  rep(0,length(Traits)))
    name <- c(name, Traits)
  }
  
  lc <- lincon(A=A,dir=dir, val=val, id=phen$Indiv, use=use, name=name)
  
}