## ------------------------------------------------------------------------
library("optiSel")
data(Cattle)
head(Cattle)

## ------------------------------------------------------------------------
data(map)
dir     <- system.file("extdata", package="optiSel")
GTfiles <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
head(map)

## ---- results="hide"-----------------------------------------------------
Kin  <- kinlist(sKin=segIBD(GTfiles, map))

## ------------------------------------------------------------------------
phen <- Cattle[Cattle$Breed=="Angler",]

## ------------------------------------------------------------------------
help.opticont(Kin, phen)

## ------------------------------------------------------------------------
con <- list()

## ------------------------------------------------------------------------
con$ub <- c(M=NA, F=-1)

## ------------------------------------------------------------------------
Ne <- 100
meanKin     <- mean(Kin$sKin[phen$Indiv, phen$Indiv])
con$ub.sKin <- meanKin + (1-meanKin)*(1/(2*Ne))

## ------------------------------------------------------------------------
maxBV <- opticont(method="max.BV", K=Kin, phen=phen, con=con, trace=FALSE)

## ------------------------------------------------------------------------
maxBV.s <- summary(maxBV)
maxBV.s$obj.fun

## ------------------------------------------------------------------------
Candidate <- maxBV$parent[,  c("Indiv", "Sex", "oc")]
head(Candidate[rev(order(Candidate$oc)),])

## ------------------------------------------------------------------------
Candidate$nOff <- noffspring(Candidate, N=250)$nOff
head(Candidate[rev(order(Candidate$oc)),])

## ---- results="hide"-----------------------------------------------------
Kin  <- kinlist(
            sKin    = segIBD(GTfiles, map, minSNP=20, minL=1000), 
            sKinatN = segIBDatN(GTfiles, Cattle, map, thisBreed="Angler", ubFreq=0.01, minL=1000)
            )

## ---- results="hide"-----------------------------------------------------
wdir  <- file.path(tempdir(), "HaplotypeEval")
wfile <- haplofreq(GTfiles, Cattle, map, thisBreed="Angler", minSNP=20, minL=1000, w.dir=wdir)
Comp  <- segBreedComp(wfile$match, map)
Cattle[rownames(Comp), "MC"] <- 1 - Comp$native

## ------------------------------------------------------------------------
head(Cattle[,-1])

## ------------------------------------------------------------------------
help.opticont(Kin, phen=Cattle)

## ------------------------------------------------------------------------
help.opticont4mb(Kin, phen=Cattle)

## ------------------------------------------------------------------------
phen        <- Cattle[Cattle$Breed=="Angler",]
con         <- list(ub=c(M=NA, F=-1))
meanKin     <- mean(Kin$sKin[phen$Indiv, phen$Indiv])
con$ub.sKin <- meanKin + (1-meanKin)*(1/(2*Ne))

## ---- results="hide"-----------------------------------------------------
maxBV   <- opticont(method="max.BV", K=Kin, phen=phen, con=con, trace=FALSE)
maxBV.s <- summary(maxBV)

## ------------------------------------------------------------------------
maxBV.s[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN")]

## ---- results="hide"-----------------------------------------------------
meanKinatN     <- mean(Kin$segIBDandN)/mean(Kin$segN)
con$ub.sKinatN <- meanKinatN +(1-meanKinatN)*(1/(2*Ne))
con$ub.MC      <- mean(phen$MC)
maxBV2         <- opticont(method="max.BV", K=Kin, phen=phen, con=con, solver="slsqp")
maxBV2.s       <- summary(maxBV2)

## ------------------------------------------------------------------------
Results <- rbind(maxBV.s, maxBV2.s)
Results[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN")]

## ------------------------------------------------------------------------
phen <- Cattle[Cattle$Breed=="Angler",]
con  <- list(ub=c(M=NA, F=-1))

## ---- results="hide"-----------------------------------------------------
minKin   <- opticont(method="min.sKin", K=Kin, phen=phen, con=con, trace=FALSE)
minKin.s <- summary(minKin)

## ------------------------------------------------------------------------
minKin.s[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN")]

## ---- results="hide"-----------------------------------------------------
con$ub.MC      <- 0.50
meanKinatN     <- mean(Kin$segIBDandN)/mean(Kin$segN)
con$ub.sKinatN <- meanKinatN +(1-meanKinatN)*(1/(2*Ne))
minKin2        <- opticont(method="min.sKin", K=Kin, phen=phen, con=con, solver="slsqp", trace=FALSE)
minKin2.s      <- summary(minKin2)

## ------------------------------------------------------------------------
Results <- rbind(minKin.s, minKin2.s)
Results[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN")]

## ---- results="hide"-----------------------------------------------------
phen           <- Cattle[Cattle$Breed=="Angler",]
con            <- list(ub=c(M=NA, F=-1))
meanKin        <- mean(Kin$sKin[phen$Indiv, phen$Indiv])
meanKinatN     <- mean(Kin$segIBDandN)/mean(Kin$segN)
con$ub.sKin    <- meanKin    + (1-meanKin)*(1/(2*Ne))
con$ub.sKinatN <- meanKinatN + (1-meanKinatN)*(1/(2*Ne))
minMC   <- opticont(method="min.MC", K=Kin, phen=phen, con=con, solver="slsqp", trace=FALSE)
minMC.s <- summary(minMC)

## ------------------------------------------------------------------------
minMC.s[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN")]

## ---- results="hide"-----------------------------------------------------
con$lb.BV <- mean(phen$BV)
minMC2    <- opticont(method="min.MC", K=Kin, phen=phen, con=con, solver="cccp", trace=FALSE)
minMC2.s  <- summary(minMC2)

## ------------------------------------------------------------------------
Results <- rbind(minMC.s, minMC2.s)
Results[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN")]

## ------------------------------------------------------------------------
CoreSet <- opticomp(Kin$sKin, Cattle$Breed, lb=c(Angler=0.1))
CoreSet$bc

## ------------------------------------------------------------------------
con         <- list(ub=c(M=NA, F=-1))
meanKin     <- mean(Kin$sKin[phen$Indiv, phen$Indiv])
con$ub.sKin <- meanKin  + (1-meanKin)*(1/(2*Ne))

## ---- results="hide"-----------------------------------------------------
minKin4mb <- opticont4mb("min.sKin.acrossBreeds", Kin, phen=Cattle, CoreSet$bc, 
                       thisBreed="Angler", con=con, trace=FALSE)
minKin4mb.s <- summary(minKin4mb)

## ------------------------------------------------------------------------
minKin4mb.s[,c("valid", "meanBV", "meanMC", "sKin", "sKinatN", "sKin.acrossBreeds")]

## ---- results="hide"-----------------------------------------------------
data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors, thisBreed="Hinterwaelder", lastNative=1970)

## ------------------------------------------------------------------------
head(Pedig)

## ------------------------------------------------------------------------
Summary <- summary(Pedig, keep=Pedig$Born %in% (2006:2007) & !is.na(Pedig$BV))
keep    <- Summary[Summary$equiGen>=5.0, "Indiv"]
table(Pedig[keep, "Sex"])

## ------------------------------------------------------------------------
Kin <- kinlist(
    pKin    = as(pedIBD(Pedig, keep.only=keep), "matrix"),
    pKinatN = pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep)
)

## ------------------------------------------------------------------------
cont     <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")
Pedig$MC <- 1-cont$native
head(cont[keep, 2:6])

## ------------------------------------------------------------------------
phen <- Pedig[Pedig$Indiv %in% keep, c("Indiv", "Sex", "Breed", "Born", "BV", "MC")]
phen$BV <- phen$BV + 4*(phen$MC-mean(phen$MC))
head(phen[,-1])

## ------------------------------------------------------------------------
help.opticont(Kin, phen)

## ------------------------------------------------------------------------
con         <- list(ub=c(M=NA, F=-1))
meanKin     <- mean(Kin$pKin[phen$Indiv, phen$Indiv])
con$ub.pKin <- meanKin + (1-meanKin)*(1/(2*Ne))

## ---- results="hide"-----------------------------------------------------
maxBV   <- opticont(method="max.BV", K=Kin, phen=phen, con=con, trace=FALSE)
maxBV.s <- summary(maxBV)

## ------------------------------------------------------------------------
maxBV.s[,c("valid", "meanBV", "meanMC", "pKin", "pKinatN")]

## ---- results="hide"-----------------------------------------------------
meanKinatN     <- mean(Kin$pedIBDandN)/mean(Kin$pedN)
con$ub.pKinatN <- meanKinatN +(1-meanKinatN)*(1/(2*Ne))
con$ub.MC      <- mean(phen$MC)
maxBV2         <- opticont(method="max.BV", K=Kin, phen=phen, con=con, solver="slsqp")
maxBV2.s       <- summary(maxBV2)

## ------------------------------------------------------------------------
Results <- rbind(maxBV.s, maxBV2.s)
Results[,c("valid", "meanBV", "meanMC", "pKin", "pKinatN")]

## ------------------------------------------------------------------------
con  <- list(ub=c(M=NA, F=-1))

## ---- results="hide"-----------------------------------------------------
minKin   <- opticont(method="min.pKin", K=Kin, phen=phen, con=con, trace=FALSE)
minKin.s <- summary(minKin)

## ------------------------------------------------------------------------
minKin.s[,c("valid", "meanBV", "meanMC", "pKin", "pKinatN")]

## ---- results="hide"-----------------------------------------------------
con  <- list(ub=c(M=NA, F=-1))
con$ub.MC   <- 0.50
meanKin     <- mean(Kin$pKin[phen$Indiv, phen$Indiv])
con$ub.pKin <- meanKin + (1-meanKin)*(1/(2*Ne))

minKin2     <- opticont(method="min.pKinatN", K=Kin, phen=phen, con=con, solver="slsqp")
minKin2.s   <- summary(minKin2)

## ------------------------------------------------------------------------
Results <- rbind(minKin.s, minKin2.s)
Results[,c("valid", "meanBV", "meanMC", "pKin", "pKinatN")]

## ---- results="hide"-----------------------------------------------------
con            <- list(ub=c(M=NA, F=-1))
meanKin        <- mean(Kin$pKin[phen$Indiv, phen$Indiv])
meanKinatN     <- mean(Kin$pedIBDandN)/mean(Kin$pedN)
con$ub.pKin    <- meanKin    + (1-meanKin)*(1/(2*Ne))
con$ub.pKinatN <- meanKinatN + (1-meanKinatN)*(1/(2*Ne))
minMC   <- opticont(method="min.MC", K=Kin, phen=phen, con=con, trace=FALSE)
minMC.s <- summary(minMC)

## ------------------------------------------------------------------------
minMC.s[,c("valid", "meanBV", "meanMC", "pKin", "pKinatN")]

## ---- results="hide"-----------------------------------------------------
con$lb.BV <- mean(phen$BV)
minMC2    <- opticont(method="min.MC", K=Kin, phen=phen, con=con, trace=FALSE)
minMC2.s  <- summary(minMC2)

## ------------------------------------------------------------------------
Results <- rbind(minMC.s, minMC2.s)
Results[,c("valid", "meanBV", "meanMC", "pKin", "pKinatN")]

