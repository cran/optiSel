## ------------------------------------------------------------------------
library("optiSel")
data(Cattle)
head(Cattle)

## ------------------------------------------------------------------------
data(map)
dir     <- system.file("extdata", package="optiSel")
GTfiles <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
head(map)

## ------------------------------------------------------------------------
phen <- Cattle[Cattle$Breed=="Angler",]

## ------------------------------------------------------------------------
sKin <- segIBD(GTfiles, map)

## ------------------------------------------------------------------------
cand  <- candes(phen=phen, sKin=sKin)

## ------------------------------------------------------------------------
con <- list()

## ------------------------------------------------------------------------
con$uniform <- "female"

## ------------------------------------------------------------------------
cand$mean

## ------------------------------------------------------------------------
Ne <- 100

con$ub.sKin <- cand$mean$sKin + (1-cand$mean$sKin)*(1/(2*Ne))

## ------------------------------------------------------------------------
Offspring <- opticont("max.BV", cand, con, trace=FALSE)

## ------------------------------------------------------------------------
Offspring$info

## ------------------------------------------------------------------------
Offspring$obj.fun

## ------------------------------------------------------------------------
Offspring$mean

## ------------------------------------------------------------------------
Candidate <- Offspring$parent[,  c("Indiv", "Sex", "oc")]
head(Candidate[rev(order(Candidate$oc)),])

## ------------------------------------------------------------------------
Candidate$nOff <- noffspring(Candidate, N=250)$nOff
head(Candidate[rev(order(Candidate$oc)),])

## ------------------------------------------------------------------------
Mating <- matings(Candidate, Kin=sKin)
head(Mating)

## ------------------------------------------------------------------------
attributes(Mating)$objval

## ---- results="hide"-----------------------------------------------------
wdir  <- file.path(tempdir(), "HaplotypeEval")
wfile <- haplofreq(GTfiles, Cattle, map, thisBreed="Angler", minL=1.0, w.dir=wdir)
Comp  <- segBreedComp(wfile$match, map)
Cattle[rownames(Comp), "NC"] <- Comp$native

## ------------------------------------------------------------------------
head(Cattle[,-1])

## ---- results="hide"-----------------------------------------------------
phen    <- Cattle[Cattle$Breed=="Angler",]
sKin    <- segIBD(GTfiles, map, minL=1.0)
sKinatN <- segIBDatN(GTfiles, Cattle, map, thisBreed="Angler", minL=1.0)

## ------------------------------------------------------------------------
cand  <- candes(phen=phen, sKin=sKin, sKinatN=sKinatN)

## ------------------------------------------------------------------------
cand$mean

## ------------------------------------------------------------------------
con         <- list(uniform="female")
con$ub.sKin <- cand$mean$sKin + (1-cand$mean$sKin)*(1/(2*Ne))

## ---- results="hide"-----------------------------------------------------
Offspring   <- opticont("max.BV", cand, con, trace=FALSE)

## ------------------------------------------------------------------------
Offspring$info

## ------------------------------------------------------------------------
Offspring$mean

## ---- results="hide"-----------------------------------------------------
con$ub.sKinatN <- cand$mean$sKinatN +(1-cand$mean$sKinatN)*(1/(2*Ne))
con$lb.NC      <- cand$mean$NC
Offspring2     <- opticont("max.BV", cand, con)

## ------------------------------------------------------------------------
rbind(Ref=cand$mean, maxBV=Offspring$mean, maxBV2=Offspring2$mean)

## ------------------------------------------------------------------------
con  <- list(uniform="female")

## ---- results="hide"-----------------------------------------------------
Offspring   <- opticont("min.sKin", cand, con)

## ------------------------------------------------------------------------
Offspring$mean

## ---- results="hide"-----------------------------------------------------
con$lb.NC      <- 0.50
con$ub.sKinatN <- cand$mean$sKinatN +(1-cand$mean$sKinatN)*(1/(2*Ne))
Offspring2     <- opticont("min.sKin", cand, con)

## ------------------------------------------------------------------------
rbind(Ref=cand$mean, minKin=Offspring$mean, minKin2=Offspring2$mean)

## ---- results="hide"-----------------------------------------------------
con            <- list(uniform="female")
con$ub.sKin    <- cand$mean$sKin    + (1-cand$mean$sKin   )*(1/(2*Ne))
con$ub.sKinatN <- cand$mean$sKinatN + (1-cand$mean$sKinatN)*(1/(2*Ne))
Offspring <- opticont("max.NC", cand, con)

## ------------------------------------------------------------------------
Offspring$info

## ------------------------------------------------------------------------
Offspring$mean

## ---- results="hide"-----------------------------------------------------
con$lb.BV   <- cand$mean$BV
Offspring2  <- opticont("max.NC", cand, con)

## ------------------------------------------------------------------------
rbind(Ref=cand$mean, maxNC=Offspring$mean, maxNC2=Offspring2$mean)

## ------------------------------------------------------------------------
cand <- candes(phen=Cattle, sKin=sKin, sKinatN.Angler=sKinatN, bc="sKin")
Unif <- c("Angler.female", "Fleckvieh", "Holstein", "Rotbunt")
mKin <- cand$mean$sKinatN.Angler
con  <- list(uniform = Unif, ub.sKinatN.Angler = mKin + (1-mKin)*(1/(2*Ne)))

## ------------------------------------------------------------------------
Offspring <- opticont("min.sKin", cand, con, trace=FALSE)

## ------------------------------------------------------------------------
Offspring$info

## ------------------------------------------------------------------------
Offspring$mean

## ------------------------------------------------------------------------
head(Offspring$parent[,c("Breed","lb","oc","ub")])

## ---- results="hide"-----------------------------------------------------
data("PedigWithErrors")
Pedig <- prePed(PedigWithErrors, thisBreed="Hinterwaelder", lastNative=1970)

## ------------------------------------------------------------------------
head(Pedig[,-1])

## ------------------------------------------------------------------------
cont     <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")
Pedig$NC <- cont$native
tail(cont[, 2:5])

## ------------------------------------------------------------------------
data("Phen")
Summary <- summary(Pedig, keep=Pedig$Indiv %in% Phen$Indiv)
keep    <- Summary$Indiv[Summary$equiGen>=3.0]
table(Pedig[keep, "Sex"])

## ------------------------------------------------------------------------
Phen <- merge(Pedig, Phen[,c("Indiv", "BV")], by="Indiv")
Phen <- Phen[Phen$Indiv %in% keep, c("Indiv", "Sex","Breed", "Born", "BV", "NC")]
head(Phen)

## ---- results="hide"-----------------------------------------------------
phen    <- Phen[Phen$Indiv %in% keep, ]
pKin    <- pedIBD(Pedig, keep.only=keep)
pKinatN <- pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep)

## ------------------------------------------------------------------------
cand <- candes(phen = phen, pKin = pKin, pKinatN = pKinatN)

## ------------------------------------------------------------------------
cand$mean

## ------------------------------------------------------------------------
con         <- list(uniform="female")
con$ub.pKin <- cand$mean$pKin + (1-cand$mean$pKin)*(1/(2*Ne))

## ---- results="hide"-----------------------------------------------------
Offspring <- opticont("max.BV", cand, con)

## ------------------------------------------------------------------------
Offspring$mean

## ---- results="hide"-----------------------------------------------------
con$ub.pKinatN <- cand$mean$pKinatN +(1-cand$mean$pKinatN)*(1/(2*Ne))
con$lb.NC      <- cand$mean$NC
Offspring2     <- opticont("max.BV", cand, con)

## ------------------------------------------------------------------------
rbind(Ref=cand$mean, maxBV=Offspring$mean, maxBV2=Offspring2$mean)

## ------------------------------------------------------------------------
con  <- list(uniform="female")

## ---- results="hide"-----------------------------------------------------
Offspring <- opticont("min.pKin", cand, con)

## ------------------------------------------------------------------------
Offspring$mean

## ---- results="hide"-----------------------------------------------------
con  <- list(uniform="female")
con$lb.NC   <- 1.05*cand$mean$NC
con$ub.pKin <- cand$mean$pKin + (1-cand$mean$pKin)*(1/(2*Ne))

Offspring2  <- opticont("min.pKinatN", cand, con)

## ------------------------------------------------------------------------
rbind(Ref=cand$mean, minKin=Offspring$mean, minKin2=Offspring2$mean)

## ---- results="hide"-----------------------------------------------------
con            <- list(uniform="female")
con$ub.pKin    <- cand$mean$pKin    + (1-cand$mean$pKin   )*(1/(2*Ne))
con$ub.pKinatN <- cand$mean$pKinatN + (1-cand$mean$pKinatN)*(1/(2*Ne))
Offspring <- opticont("max.NC", cand, con)

## ------------------------------------------------------------------------
Offspring$mean

## ---- results="hide"-----------------------------------------------------
con$lb.BV  <- cand$mean$BV
Offspring2 <- opticont("max.NC", cand, con)

## ------------------------------------------------------------------------
rbind(Ref=cand$mean, maxNC=Offspring$mean, maxNC2=Offspring2$mean)

