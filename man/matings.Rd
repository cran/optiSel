\name{matings}
\Rdversion{1.1}
\alias{matings}
\title{Mate Allocation}

\description{
Males and females are allocated for mating such that all breeding animals have the desired number of offspring. The mean inbreeding coefficient in the offspring is minimized.
}
\usage{
matings(cand, Kin, N=2*sum(cand$Sex=="female"), alpha=1, ub.nOff=NA, max=FALSE)}

\arguments{
\item{cand}{Data frame with optimum contributions (column \code{oc}), sexes (column \code{Sex}), and IDs (column \code{Indiv}) of the selection candidates. The data frame may also contain column \code{herd} containing the name of the herd to which each female belongs (\code{NA} for males).}
\item{Kin}{Kinship matrix (or an other similarity matrix) for selection candidates.}
\item{N}{Desired number of offspring that should be available as selection candidates in the next generation.}
\item{alpha}{If \code{alpha<1} then the proportion of females mated with the same male is at most \code{alpha} in each herd. A value \code{alpha<1} increases genetic connectedness between herds and enables to estimate more accurate breeding values.}
\item{ub.nOff}{Maximum number of offspring per mating to be used as selection candidates. Without this constraint (i.e. \code{ub.nOff=NA}), mate allocation can result in many full sib matings.}
\item{max}{The default \code{max=FALSE} means that the objective function is minimized.}
}

\details{
Males and females are allocated for mating such that all breeding animals have the desired number of offspring. If \code{Kin} is a kinship matrix, then the mean inbreeding coefficient in the offspring is minimized. The maximum number of offspring per mating can be constrained. For each herd, the proportion \code{alpha} of females mated with the same male can be constrained as well. In practice, the number of offspring per mating may be higher since not all offspring will be suitable as selection candidates in the next generation. 
}

\value{
Data frame with columns \code{Sire}, \code{Dam}, \code{nOff}, and \code{herd}, whereby column \code{nOff} contains the desired number of offspring from each mating, and column \code{herd} contains the herd of the dam.
}


\examples{

data(map)
dir   <- system.file("extdata", package = "optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")
sKin  <- segIBD(files, map, minSNP=20, minL=2.0)
Kin   <- kinlist(sKin = sKin)

data(Cattle)
Phen  <- Cattle[Cattle$Breed=="Angler", ]
head(Phen)

help.opticont(Kin, Phen)

con   <- list(ub.sKin = 0.057)
maxBV <- opticont("max.BV", K=Kin, phen=Phen, con=con, solver="cccp2", trace=FALSE)
(summary(maxBV))

#Minimize inbreeding
Mating <- matings(maxBV$parent, Kin=sKin,  alpha=0.3, ub.nOff=5)
Mating 
attributes(Mating)$objval


\dontshow{
### Check results ###

library(data.table)

Candidate      <- maxBV$parent
Candidate$nOff <- noffspring(Candidate, N=2*sum(Candidate$Sex=="female"))$nOff

setDT(Mating)
DF <- as.data.frame(Mating[,.(nOff=sum(nOff)),by=Sire])
rownames(DF)<-DF$Sire
range(DF$nOff- Candidate[rownames(DF),"nOff"])

DF <- as.data.frame(Mating[,.(nOff=sum(nOff)),by=Dam])
rownames(DF)<-DF$Dam
range(DF$nOff- Candidate[rownames(DF),"nOff"])

mean(sKin)
attributes(Mating)$objval #Expected mean inbreeding in the offspring

nOffperHerd <- Mating[,.(nOff=sum(nOff)),by=.(herd,Sire)]
setkey(nOffperHerd, herd, nOff)
nOffperHerd

}

}


\author{Robin Wellmann}
