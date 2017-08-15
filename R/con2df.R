
"con2df" <- function(con, Traits){
  cnames <- setdiff(names(con),c("lb", "ub", "uniform"))
  if(length(cnames)>0){
    condf  <- data.frame(dir=str_sub(cnames,1,2), var=str_sub(cnames,4,-1), stringsAsFactors = FALSE)
    if(any(duplicated(condf$var))){
      stop("Some variables appear in more than one constraint.\n")
    }
    rownames(condf) <- condf$var
    condf$dir   <- mapvalues(condf$dir, c("ub", "lb", "eq"), c("<=", ">=", "=="), warn_missing = FALSE)
    condf$val   <- unlist(con[cnames])
    condf$isLin <- condf$var %in% Traits
  }else{
    condf <- data.frame(dir=character(0),var=character(0),val=numeric(0), isLin=logical(0))
  }
  condf
}