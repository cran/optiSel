

"segN"<-function(Native, map, unitP="kb", keep=NULL){
  ##################################################
  # Convert data tables to data frames             #
  ##################################################
  mapAsDataTable <- "data.table" %in% class(map)
  map <- as.data.frame(map)
  if(mapAsDataTable){setDF(map)}
  
  ##################################################
  #       Data preparation                         #
  ##################################################
  if(!("Name"  %in% colnames(map))){ stop("Column 'Name' is missing in map.")}
  if(!("Chr"   %in% colnames(map))){ stop("Column 'Chr'  is missing in map.")}  
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
  }  
  checkMap(map, unitP, unitP)
  map$Name <- as.character(map$Name)
  map$Chr       <- as.character(map$Chr)
  rownames(map) <- map$Name
  if(!(unitP %in% c("SNP", colnames(map)))){
    stop("Parameter unitP must be either 'SNP' or the name of a column in the marker map.\n")
  }
  
  if(is.vector(Native)){
    names(Native)<-str_sub(str_extract(basename(Native), "Chr[0-9,a-z]*"),4,-1)
    if(!all(names(Native) %in% map$Chr)){
      stop("Some chromosomes are not in the map.")
    }
    if(any(duplicated(names(Native)))){
      stop("For some chromosomes different files were provided.")
    }
    IndivFileN <- scan(Native[1], nlines=1, what="character",quiet=TRUE)[-1]
  }else{
    IndivFileN <- colnames(Native)
  }
  if(is.null(keep)){
    indexN <- which(IndivFileN %in% IndivFileN)
  }else{
    indexN <- which(IndivFileN %in% keep)
  }
  NFileN <- length(IndivFileN)
  Indiv  <- IndivFileN[indexN]
  NC     <- length(indexN)

  if(is.vector(Native)){
    Res  <- matrix(0,nrow=NC, ncol=NC)
    GenL <- 0
    for(chr in names(Native)){
      cat(paste("Reading chromosome ", chr, "...  "))
      submap <- map[map$Chr==chr, ]
      M      <- nrow(submap)
      submap$SNP <- 1:M
      x   <- submap[, unitP]
      kb  <- (c(0,x)+c(x,x[length(x)]+x[1]))/2
      Nkb <- kb[2:length(kb)]-kb[1:(length(kb)-1)]
      Res <- Res + rcpp_segN(as.character(Native[chr]), as.integer(NFileN), as.integer(NC), as.integer(indexN-1), as.integer(M), as.double(Nkb))
      GenL <- GenL + sum(Nkb)
      }
  }else{
    map <- map[rownames(Native),]
    map$SNP <- NA
    map$nUnits <- NA
    for(chr in unique(map$Chr)){
      map[map$Chr==chr, "SNP"] <- (1:sum(map$Chr==chr))
      x   <- map[map$Chr==chr, unitP]
      kb  <- (c(0,x)+c(x,x[length(x)]+x[1]))/2
      Nkb <- kb[2:length(kb)]-kb[1:(length(kb)-1)]
      map[map$Chr==chr, "nUnits"] <- Nkb
    }
    Units  <- matrix(map$nUnits,ncol=NC, nrow=nrow(Native), byrow=FALSE)
    Native <- Native[, indexN]
    Res    <- t(Native)%*%(Native*Units)
    GenL   <- sum(map$nUnits)
  }
  
  N   <- ncol(Res)/2
  Res <- ((Res[2*(1:N)-1,2*(1:N)-1] + Res[2*(1:N)-1,2*(1:N)] + Res[2*(1:N),2*(1:N)-1] + Res[2*(1:N),2*(1:N)])/GenL)/4
  rownames(Res) <- Indiv[2*(1:N)]
  colnames(Res) <- Indiv[2*(1:N)]
  Res
}