## ------------------------------------------------------------------------
Pedigree <- data.frame(
  ID   = c("Iffes",     "Peter",    "Anna-Lena", "Kevin",    "Horst"),
  dad  = c("Kevin",     "Kevin",     NA,          0,         "Horst"),
  mom  = c("Chantalle", "Angelika", "Chantalle", "",           NA),
  Breed= c("Angler",    "Angler",   "Angler",    "Holstein", "Angler"),
  Born = c(2015,        2016,       2011,       2010,        2015)
  )
Pedigree

## ------------------------------------------------------------------------
library("optiSel")
Pedig <- prePed(Pedigree)

## ------------------------------------------------------------------------
Pedig

## ---- warning=FALSE------------------------------------------------------
sPed <- subPed(Pedig, keep = c("Chantalle","Angelika"), prevGen = 2, succGen = 1)
pedplot(sPed, label = c("Indiv", "Born", "Breed"), cex = 0.55)

## ---- results="hide"-----------------------------------------------------
data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)

## ------------------------------------------------------------------------
compl <- completeness(Pedig, keep=Pedig$Born %in% (2006:2007), by="Indiv")
head(compl)

## ---- warning=FALSE------------------------------------------------------
compl <- completeness(Pedig, keep=Pedig$Born %in% (2006:2007), by="Sex")
library("ggplot2")
ggplot(compl, aes(x=Generation, y=Completeness, col=Sex)) + geom_line()

## ------------------------------------------------------------------------
keep <- Pedig$Indiv[Pedig$Born %in% (2006:2007)]
Summary <- summary(Pedig, keep.only=keep)
head(Summary[Summary$equiGen>3.0, -1])

## ------------------------------------------------------------------------
Animal <- pedInbreeding(Pedig)
mean(Animal$Inbr[Animal$Indiv %in% keep])

## ------------------------------------------------------------------------
pedKIN <- pedIBD(Pedig, keep.only=keep)
use    <- Pedig$Sex==1 & Pedig$Indiv %in% keep & summary(Pedig)$equiGen>5 & Pedig$BV>1.7
Males  <- Pedig$Indiv[use]
pedKIN[rownames(pedKIN) %in% Males, "276000813677683", drop=FALSE]

## ---- results="hide"-----------------------------------------------------
Pedig <- prePed(PedigWithErrors, thisBreed="Hinterwaelder", lastNative=1970)

## ------------------------------------------------------------------------
cont  <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")
Pedig$MC <- 1-cont$native
head(cont[rev(keep), 2:6])

## ---- results="hide"-----------------------------------------------------
fD  <- pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep)

## ------------------------------------------------------------------------
pedKINatN <- fD$pedIBDandN/fD$pedN
use <- rownames(fD$pedN)[diag(fD$pedN)>0.2]
mean(pedKINatN[use,use])

## ------------------------------------------------------------------------
sK    <- pedKIN[use, use]
sKatN <- pedKINatN[use, use]
cor(sK[upper.tri(sK)], sKatN[upper.tri(sKatN)], use="complete.obs")

## ------------------------------------------------------------------------
use   <- rownames(fD$pedN)[diag(fD$pedN)>0.01]
sK    <- pedKIN[use, use]
sKatN <- pedKINatN[use, use]
cor(sK[upper.tri(sK)], sKatN[upper.tri(sKatN)], use="complete.obs")

## ------------------------------------------------------------------------
1-mean(pedKIN[keep, keep])

## ------------------------------------------------------------------------
mean(fD$pedIBDandN)/mean(fD$pedN)

## ------------------------------------------------------------------------
1- mean(fD$pedIBDandN)/mean(fD$pedN)

## ------------------------------------------------------------------------
fD  <- pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep, nGen=6)

## ------------------------------------------------------------------------
attributes(fD)$nativeNe

## ------------------------------------------------------------------------
id     <- Summary$Indiv[Summary$equiGen>=4 & Summary$Indiv %in% keep]
g      <- Summary[id, "equiGen"]
N      <- length(g)
n      <- (matrix(g, N, N, byrow=TRUE) + matrix(g, N, N, byrow=FALSE))/2
deltaC <- 1 - (1-pedKIN[id,id])^(1/n)
Ne     <- 1/(2*mean(deltaC))
Ne

## ------------------------------------------------------------------------
data(ExamplePed)
Pedig    <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
Kinships <- kinlist(pedIBD=pedIBD(Pedig), pedIBDatN=pedIBDatN(Pedig, thisBreed="Hinterwaelder"))
Kin      <- kinwac(Kinships, Pedig=Pedig, use=Pedig$Breed=="Hinterwaelder")
sy       <- summary(Kin, tlim=c(1970, 1995), histNe=150, base=1800)

