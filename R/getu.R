
"getu" <- function(phen, breed, bc){
  if(breed=="across breeds"){breed <- unique(phen$Breed)}
  u <- setNames(rep(0, nrow(phen)), phen$Indiv)
  for(b in breed){
    thisBreed <- phen$Breed==b
    if(any(is.na(phen$Sex[thisBreed]))){
      u[thisBreed] <- 1/sum(thisBreed)
    }else{
      isMale       <- phen$Sex=="male"    & thisBreed
      isFemale     <- phen$Sex=="female"  & thisBreed
      u[isMale]    <- 1/(2*sum(isMale))
      u[isFemale]  <- 1/(2*sum(isFemale))
    }
  }
  if(length(breed)>1){
    u <- bc[phen$Breed]*u
    }
  return(u)
}