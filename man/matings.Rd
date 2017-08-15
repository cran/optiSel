\name{matings}
\Rdversion{1.1}
\alias{matings}
\title{Mate Allocation}

\description{
Males and females are allocated for mating such that all breeding animals have the desired number of offspring. The mean inbreeding coefficient in the offspring is minimized if matrix \code{Kin} contains pairwise kinships of the selection candidates.
}
\usage{
matings(cand, Kin, N=2*sum(cand$Sex=="female"), alpha=1, 
    ub.nOff=NA, max=FALSE, solver="default", ...)}

\arguments{
\item{cand}{Data frame with optimum contributions (column \code{oc}), sexes (column \code{Sex}), and IDs (column \code{Indiv}) of the selection candidates. The data frame may also contain column \code{herd} containing the names of the herds to which the females belong (\code{NA} for males).}
\item{Kin}{Kinship matrix (or an other similarity matrix) for selection candidates.}
\item{N}{Desired number of offspring that should be available as selection candidates in the next generation.}
\item{alpha}{If \code{alpha<1} then the proportion of females mated with the same male is at most \code{alpha} in each herd. A value \code{alpha<1} increases genetic connectedness between herds and enables to estimate more accurate breeding values.}
\item{ub.nOff}{Maximum number of offspring per mating. Without this constraint (i.e. \code{ub.nOff=NA}), some superior animals may always be mated to the same inferior animal, so their offspring would likely not be suitable for breeding.}
\item{max}{The default \code{max=FALSE} means that the objective function is minimized.}
\item{solver}{Either \code{solver="default"}, or \code{solver=Rsymphony_solve_LP}. The latter is possible only if package \code{Rsymbhony} is loaded, which is not available for all platforms. }
\item{...}{Further optimization parameters. By default, they are passed to function \link[ECOSolveR]{ecos.control}}.
}

\details{
Males and females are allocated for mating such that all breeding animals have the desired number of offspring. If \code{Kin} is a kinship matrix, then the mean inbreeding coefficient in the offspring is minimized. In general, the mean similarity of the parents is minimized.

The maximum number of offspring per mating can be constrained. For each herd, the proportion \code{alpha} of females mated with the same male can be constrained as well, but this increases computation time. 

Note that the desired number of offspring per mating is the number of offspring that should be used as selection candidates in the next generation, which is not necessarily the total number of offspring.
}

\value{
Data frame with columns \code{Sire}, \code{Dam}, \code{nOff}, and possibly \code{herd}, whereby column \code{nOff} contains the desired number of offspring from each mating, and column \code{herd} contains the herd of the dam.

The data frame has attributes \code{objval} with the value of the objective function (usually the mean inbreeding coefficient), and attribute \code{info} describing the solution as reported by the solver.

}


\examples{

data("map")
data("Cattle")
dir   <- system.file("extdata", package = "optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")
sKin  <- segIBD(files, map, minSNP=20, minL=2.0)
Phen  <- Cattle[Cattle$Breed=="Angler", ]
cand  <- candes(phen=Phen, sKin = sKin)
con   <- list(ub.sKin = 0.057)
Offspring <- opticont("max.BV", cand, con, trace=FALSE)


#####  Minimize inbreeding   #####
Candidate <- Offspring$parent
Mating <- matings(Candidate, sKin, ub.nOff=5, maxit=50L)
Mating 
attributes(Mating)$objval
#[1] 0.007270397

\dontrun{
## This is faster but not available on all platforms:
library("Rsymphony")
Mating <- matings(Candidate, sKin, ub.nOff=5, alpha=0.13, solver=Rsymphony_solve_LP)
Mating
attributes(Mating)$objval
#[1] 0.009648924

attributes(Mating)$info
#[1] "Optimum solution found"
}

\dontshow{
### Check results ###

library(data.table)

Candidate      <- Offspring$parent
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

#nOffperHerd <- Mating[,.(nOff=sum(nOff)),by=.(herd,Sire)]
#setkey(nOffperHerd, herd, nOff)
#nOffperHerd

}

}


\author{Robin Wellmann}
