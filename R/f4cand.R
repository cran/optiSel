

"f4cand" <- function(cand, method, Traits, bc=NULL){
  optiParam <- str_sub(method, 5, -1)
  if(optiParam %in% Traits){
    a <- cand$phen[[optiParam]]
    a[is.na(a)] <- 0
    f <- linfun(a=a, id=cand$phen$Indiv, name=optiParam)
  }else{
    f <- cand[[optiParam]]
    if(class(f)=="ratioFun" && (f$breed!="across breeds") && !is.null(bc)){
      f$d1 <- f$d1*(bc[f$breed])^2
      f$d2 <- f$d2*(bc[f$breed])^2
    }
  }
  return(f)
}