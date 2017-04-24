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
data("PedigWithErrors")
data("Phen")
PedigWithErrors <- merge(PedigWithErrors, Phen[,c("Indiv", "BV")], by="Indiv", all="TRUE")
Pedig <- prePed(PedigWithErrors)

## ------------------------------------------------------------------------
compl <- completeness(Pedig, keep=Phen$Indiv, by="Indiv")
head(compl)

## ---- warning=FALSE------------------------------------------------------
compl <- completeness(Pedig, keep=Phen$Indiv, by="Sex")
library("ggplot2")
ggplot(compl, aes(x=Generation, y=Completeness, col=Sex)) + geom_line()

## ------------------------------------------------------------------------
keep <- Phen$Indiv
Summary <- summary(Pedig, keep.only=keep)
head(Summary[Summary$equiGen>3.0, -1])

## ------------------------------------------------------------------------
Animal <- pedInbreeding(Pedig)
mean(Animal$Inbr[Animal$Indiv %in% keep])

## ------------------------------------------------------------------------
pedKIN <- pedIBD(Pedig, keep.only=keep)
use    <- Pedig$Sex=="male" & Pedig$Indiv %in% keep & summary(Pedig)$equiGen>5 & Pedig$BV>1.0
Males  <- Pedig$Indiv[use]
pedKIN[rownames(pedKIN) %in% Males, "276000812750188", drop=FALSE]

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

## ---- results="hide"-----------------------------------------------------
data("PedigWithErrors")
set.seed(1)
Pedig <- prePed(PedigWithErrors, thisBreed="Hinterwaelder", lastNative=1970)
use   <- Pedig$Breed=="Hinterwaelder" & Pedig$Born %in% (1970:2000)
IDs   <- split(Pedig$Indiv[use], Pedig$Born[use])
keep  <- unlist(mapply(sample, x=IDs, size=pmin(lapply(IDs,length), 50), SIMPLIFY=FALSE))

Kinships <- kinlist(pedIBD   = pedIBD(Pedig, keep.only=keep), 
                    pedIBDatN= pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep))

## ---- fig.show='hold'----------------------------------------------------
sy <- summary(Kinships, Pedig, tlim=c(1970, 2000), histNe=150, base=1800, df=4)
ggplot(sy, aes(x=cohort, y=Ne)) + geom_line() + ylim(c(0,100))
ggplot(sy, aes(x=cohort, y=NGE)) + geom_line() + ylim(c(0,7))

## ---- fig.width = 5, fig.height = 3--------------------------------------
use  <- Pedig$Breed=="Hinterwaelder" & Pedig$Born %in% (1950:1995)
cont <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")
contByYear <- conttac(cont, cohort=Pedig$Born, use=use, mincont = 0.01)
ggplot(contByYear, aes(x=Year, y=Contribution, fill=Breed)) + geom_area(color="black")

