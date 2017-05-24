
"segIBDatN" <- function(files, phen, map, thisBreed, refBreeds="others", ubFreq=0.01, minSNP=20, unitP="Mb", minL=1.0, unitL="Mb", a=0.0, keep=NULL, lowMem=TRUE, skip=NA, cskip=NA, cores=1){
  ##################################################
  # Convert data tables to data frames             #
  ##################################################
  mapAsDataTable <- "data.table" %in% class(map)
  map <- as.data.frame(map)
  if(mapAsDataTable){setDF(map)}
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(mapAsDataTable){setDF(phen)}
  
  ##################################################
  #       Data preparation                         #
  ##################################################
  if(!("Indiv" %in% colnames(phen))){stop("Column 'Indiv' is missing in phen.")}
  if(!("Breed" %in% colnames(phen))){stop("Column 'Breed' is missing in phen.")}
  if(!("Name"  %in% colnames(map))){ stop("Column 'Name' is missing in map.")}
  if(!("Chr"   %in% colnames(map))){ stop("Column 'Chr'  is missing in map.")}
  phen$Indiv <- as.character(phen$Indiv)
  phen$Breed <- as.character(phen$Breed)
  if(is.logical(keep)){keep <- phen[keep,"Indiv"]}
  if(is.null(keep)){keep <- phen$Indiv}
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
  }
  phen <- phen[!duplicated(phen$Indiv),]
  phen[is.na(phen$Breed),"Breed"] <- "notSpecified"
  map$Name      <- as.character(map$Name)
  map$Chr       <- as.character(map$Chr)
  rownames(map) <- map$Name
  if(is.character(files)){
    files <- list(hap.thisBreed=files)
  }
  if(is.na(skip)){  skip <- getskip(files$hap.thisBreed[1])}
  if(is.na(cskip)){cskip <- getcskip(files$hap.thisBreed[1], skip)}
  
  ##################################################
  #       Main part                                #
  ##################################################
  cat("Identifying native alleles...\n")
  if(lowMem){
    wdir <- file.path(tempdir(), "tempHapGFFdBvcw")
    dir.create(wdir, showWarnings=FALSE)
    Res <- haplofreq(files=files, phen=phen, map=map, thisBreed=thisBreed, refBreeds=refBreeds, minSNP=minSNP, minL=minL, unitL=unitL, ubFreq=ubFreq, keep=keep, skip=skip, cskip=cskip, what="match", w.dir=wdir, cores=cores)
    Native <- Res$match
  }else{
    Native <- haplofreq(files=files, phen=phen, map=map, thisBreed=thisBreed, refBreeds=refBreeds, minSNP=minSNP, minL=minL, unitL=unitL, keep=keep, skip=skip, cskip=cskip, what="freq", cores=cores)$freq < ubFreq
  }
  Res  <- list()
  cat("Computing probabilities for segments to be shared and native...\n")
  Res$segIBDandN <- segIBDandN(files=files$hap.thisBreed, Native=Native, map=map, minSNP=minSNP, unitP=unitP, minL=minL, unitL=unitL, a=a, keep=keep, skip=skip, cskip=cskip, cores=cores)
  cat("Computing probabilities for segments to be native...\n")
  Res$segN       <- segN(Native=Native, map=map, unitP=unitP, keep=keep, cores=cores)
  Res$segZ       <- 1+Res$segIBDandN-Res$segN
  attributes(Res)$condProb <- list(segIBDatN=c(f1="segZ", f2="segN"))
  class(Res)     <- "kinMatrices"
  attributes(Res)$meanIBDatN <- mean(Res$segIBDandN)/mean(Res$segN)
  cat("Mean kinship at native segments:", round(attributes(Res)$meanIBDatN,4), "\n")
  if(lowMem){
    file.remove(Native)
    unlink(wdir, recursive=TRUE)
  }
  Res
}
