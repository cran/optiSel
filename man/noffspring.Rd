\name{noffspring}
\Rdversion{1.1}
\alias{noffspring}
\title{Calculates Optimum Numbers of Offspring}

\description{
Calculates the optimum numbers of offspring from optimum contributions of selection candidates.
}
\usage{
noffspring(cand, N)}

\arguments{
\item{cand}{Data frame with optimum contributions (column \code{oc}), sexes (column \code{Sex}), and IDs (column \code{Indiv}) of the selection candidates.}
\item{N}{Desired number of individuals in the offspring population.}
}

\details{
The function calculates the optimum numbers of offspring of the selection candidates from the optimum contributions \code{cand$oc} and the size \code{N} of the offspring population.
}

\value{

Data frame with column \code{Indiv} containing the individual IDs and  integer column \code{nOff} containing the optimum numbers of offspring. This column 
 is approximately \code{2*N*cand$oc} with \code{sum(noff[cand$Sex==1])=N} and \code{sum(noff[cand$Sex==2])=N}.
}


\examples{
data(PedigWithErrors)
data(Phen)

keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
Kin   <- kinlist(fA=pedIBD(Pedig, keep.only=keep))
con   <- list(ub.fA=0.03, ub=c(M=NA, F=-1))
maxBV <- opticont("max.BV", K=Kin, phen=Phen, con=con, trace = FALSE)
summary(maxBV)

N <- 150
cand <- maxBV$parent
cand$nOff <- noffspring(cand, N)$nOff

sum(cand$nOff[cand$Sex==1])
#[1] 150

sum(cand$nOff[cand$Sex==2])
#[1] 150

round(2*N*cand$oc-cand$nOff, 2)
}

\author{Robin Wellmann}
