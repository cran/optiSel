
makeA <- function(Pedig, keep.only=NULL, keep=keep.only, AFounder=NULL){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  if(PedigAsDataTable){
    Pedig <- as.data.frame(Pedig)
    setDF(Pedig)
    }
  
  ids   <- as.character(Pedig[[1]])
  Pedig <- prePed(Pedig[,1:3], keep=keep, addNum=TRUE)
  
  if(is.null(keep.only)){
    keep.only <- ids
  }else{
    keep.only <- as.character(keep.only)
    keep.only <- setdiff(keep.only, c(NA, "", " ", "0"))
    keep.only <- ids[ids %in% keep.only]
  }

  if(is.null(AFounder)){
    AFounder   <- matrix(1, 0, 0)
    numFounder <- integer(0)
  }else{
    idF        <- Pedig$Indiv[Pedig$Indiv %in% rownames(AFounder)]
    AFounder   <- AFounder[idF, idF]
    numFounder <- Pedig[idF, "numIndiv"]
  }
  pedKin <- rcpp_makeA(as.integer(Pedig$numSire), as.integer(Pedig$numDam), AFounder, numFounder-1, as.character(Pedig$Indiv))
  if(identical(Pedig$Indiv, keep.only)){return(pedKin)}
  pedKin[keep.only, keep.only]
}
