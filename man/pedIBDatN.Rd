\name{pedIBDatN}
\Rdversion{1.1}
\alias{pedIBDatN}
\title{Calculates the Pedigree Based Kinship at Native Alleles
}
\description{
Calculates the \bold{ped}igree based probability of alleles to be \bold{IBD} \bold{at N}ative alleles.
}
\usage{
pedIBDatN(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only, nGen=6)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, (4) Sex, and (5) Breed. Missing parents are coded as NA. Both parents must either be missing or present. If this is not the case use \link{prePed}.}
\item{thisBreed}{Name of the breed in column (5) of the pedigree for which the kinships are to be computed.}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
\item{nGen}{Number of generations taken into account for estimating the native effective size.}
}

\details{
Calculates a list containing matrices needed to compute pedigree based kinships at native alleles, defined as the conditional probability that two randomly chosen alleles are IBD, given that both originate from native founders.
A native founder is an individual with unkown parents belonging to \code{thisBreed}.
}

\value{
A list with the following components:
\item{pedZ}{matrix with pedZij = Bij = Probability that two alleles chosen from individuals i and j are IBD or at least one of them is a migrant allele.}
\item{pedN}{matrix with pedNij = Probability that two alleles chosen from individuals i and j are both native alleles.}

The list has class attribute \code{"kinMatrices"} and the additional attribute \code{condProb}.
}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
fD    <- pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep)

#Native Ne = 49.5 (estimated from 6 previous generations)
#Mean kinship of native alleles:  0.0777
}



\author{Robin Wellmann}
