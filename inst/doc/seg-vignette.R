## ------------------------------------------------------------------------
library("optiSel")
data(Cattle)
phen <- Cattle
head(phen)

## ------------------------------------------------------------------------
table(phen$Breed)

## ------------------------------------------------------------------------
data(map)
head(map)

## ------------------------------------------------------------------------
tapply(map$kb, map$Chr, max)

## ------------------------------------------------------------------------
dir     <- system.file("extdata", package="optiSel")
GTfiles <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))

## ---- results="hide"-----------------------------------------------------
Animal <- segInbreeding(GTfiles, map, minSNP=20, minL=1000)

## ------------------------------------------------------------------------
head(Animal)

## ---- results="hide"-----------------------------------------------------
segKIN <- segIBD(GTfiles, map, minSNP=20, minL=1000)

## ------------------------------------------------------------------------
segKIN[1:3,1:3]

## ------------------------------------------------------------------------
Males  <- phen$Indiv[phen$Sex==1 & phen$Breed=="Angler" & phen$BV>2.0]
segKIN[rownames(segKIN) %in% Males, "276000102372349", drop=FALSE]

## ------------------------------------------------------------------------
D     <- sim2dis(segKIN, a=6.0, baseF=0.03, method=1)
color <- c(Angler="red", Rotbunt="green", Fleckvieh="blue", Holstein="black")
col   <- color[phen[rownames(D), "Breed"]]
Res   <- cmdscale(D)
plot(Res, pch=18, col=col, main="Multidimensional Scaling", cex=0.5, xlab="",ylab="", asp=1)

## ---- fig.width = 5, results="hide"--------------------------------------
Haplo <- haplofreq(GTfiles, phen, map, thisBreed="Angler", refBreeds="Rotbunt",   minSNP=20, minL=1000)
plot(Haplo, ID="276000101676415", hap=2)

## ---- fig.width = 5, results="hide"--------------------------------------
Haplo <- freqlist(
  haplofreq(GTfiles, phen, map, thisBreed="Angler", refBreeds="Rotbunt",   minSNP=20, minL=1000),
  haplofreq(GTfiles, phen, map, thisBreed="Angler", refBreeds="Holstein",  minSNP=20, minL=1000),
  haplofreq(GTfiles, phen, map, thisBreed="Angler", refBreeds="Fleckvieh", minSNP=20, minL=1000)
  )

plot(Haplo, ID=1, hap=2, refBreed="Rotbunt")

## ---- results="hide"-----------------------------------------------------
Haplo <- haplofreq(GTfiles, phen, map, thisBreed="Angler", refBreeds="others", ubFreq=0.01, minL=2500)

## ------------------------------------------------------------------------
Haplo$freq[1:10,1:3]

## ------------------------------------------------------------------------
Haplo$match[1:10,1:3]

## ---- results="hide"-----------------------------------------------------
wdir  <- file.path(tempdir(), "HaplotypeEval")
wfile <- haplofreq(GTfiles, phen, map, thisBreed="Angler", minSNP=20, minL=1000, w.dir=wdir)

## ------------------------------------------------------------------------
Comp  <- segBreedComp(Haplo$match, map)
head(Comp[,-1])

## ------------------------------------------------------------------------
Average <- apply(Comp[,-1],2,mean)
round(Average, 3)

## ---- results="hide"-----------------------------------------------------
fD <- segIBDatN(GTfiles, phen, map, thisBreed="Angler", ubFreq=0.01, minL=1000)

## ------------------------------------------------------------------------
segKINatN <- fD$segIBDandN/fD$segN
segKINatN[c(2,4,5), c(2,4,5)]

## ------------------------------------------------------------------------
keep <- phen$Indiv[phen$Breed=="Angler"]
1 - mean(segKIN[keep, keep])

## ------------------------------------------------------------------------
mean(fD$segIBDandN)/mean(fD$segN)

## ------------------------------------------------------------------------
1 - mean(fD$segIBDandN)/mean(fD$segN)

## ---- results="hide"-----------------------------------------------------
segKIN  <- segIBD(GTfiles, map, minSNP=20, minL=1000)

## ------------------------------------------------------------------------
Breed   <- phen[rownames(segKIN),"Breed"]
CoreSet <- opticomp(segKIN, Breed)
round(CoreSet$f, 3)

## ------------------------------------------------------------------------
round(CoreSet$Dist, 3)

## ------------------------------------------------------------------------
CoreSet <- opticomp(segKIN, Breed)
CoreSet$bc

## ------------------------------------------------------------------------
CoreSet$value

## ------------------------------------------------------------------------
CoreSet <- opticomp(segKIN, Breed, ub=c(Angler=0))
CoreSet$bc

## ------------------------------------------------------------------------
CoreSet$value

