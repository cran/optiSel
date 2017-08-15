
"kinwac"<-function(cand){
  if(class(cand)!="candes"){stop("Argument cand must be created with function candes.\n")}
  if(!("Born" %in% colnames(cand$phen))){stop("Column 'Born' with Year-of-Birth is missing.\n")}
  if(!is.numeric(cand$phen$Born)){       stop("Column 'Born' is not numeric.\n")}  
  if(all(is.na(cand$phen$Born))){        stop("Column 'Born' contains only NA.\n")}  


  ### Define Indicator Matrix A with ########################
  ### A[i,j]=1 if i and j are from the same birth cohort. ###
  
  cohort <- cand$phen$Born
  cohort[is.na(cohort)]<- -123456789
  Years <- unique(cohort)
  i<-NULL
  j<-NULL
  for(k in Years){
    if(sum(cohort==k)>0){
    i<-c(i, rep(which(cohort==k), sum(cohort==k)))
    j<-c(j, rep(which(cohort==k), each=sum(cohort==k)))
    }
  }
  A <- sparseMatrix(i=i[i!=j],j=j[i!=j],dims=c(length(cohort),length(cohort)))
  A <- as(A, "dgCMatrix")
  cohort[cohort==-123456789]<- NA

  ### Compute mean kinship of individuals with all ##########
  ### individuals from the same birth cohort. ###############
  
  N <- apply(A, 2, sum)
  
  for(i in seq_along(cand)){
    if(class(cand[[i]])%in% c("quadFun","ratioFun")){
      cand$phen[[cand[[i]]$name]] <- NA
    }
  }
  
  
  for(i in seq_along(cand)){
    if(class(cand[[i]])=="quadFun"){
      name <- cand[[i]]$name
      cand$phen[[name]] <- apply((cand[[i]]$Q)*A, 2, sum)/N
      cand$phen[is.na(cohort), name] <- NA
    }
    if(class(cand[[i]])=="ratioFun"){
      name  <- cand[[i]]$name
      val1 <- apply((cand[[i]]$Q1)*A, 2, sum)/N
      val2 <- apply((cand[[i]]$Q2)*A, 2, sum)/N
      cand$phen[[name]]  <- val1/val2
      cand$phen[is.na(cohort), name]  <- NA
    }
  }
  
  return(cand$phen)
}