\name{kin}
\Rdversion{1.1}
\alias{kin}
\title{Computation of kinships according to different definitions
}
\description{
Computation of pedigree based kinship and modified kinships taking allele origin into account
}
\usage{
kin(Pedig, thisBreed=NA, method=c("A","B","C","D"), keep.only, keep=keep.only, nGen=6)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, (4) Sex, and (5) Breed. Missing parents are coded as NA. Both parents must either be missing or present. If this is not the case use \link{prePed}.}
\item{thisBreed}{Name of the breed in column (5) of the pedigree for which the kinships are to be computed.}
\item{method}{Kinship matrices which are to be computed. For methods "A", "B", "C", and "N" the matrices fA, fB, fC, and fN are computed, respectively. For method "D", the matrices fB and fN are computed since they are needed to compute the conditional kinship fD.
}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
\item{nGen}{For method "D" an estimate of the native effective size computed from \code{nGen} generations is provided as additional information if vector \code{keep} contains individuals from a birth cohort. }

}

\details{
Computation of pedigree based kinship matrix fA and modified kinship matrices fB, fC, fN taking allele origin into account. They are defined as
\tabular{ll}{
A) \tab fAij = Probability that two alleles chosen from i and j are IBD, \cr
B) \tab fBij = Probability that two alleles chosen from i and j are IBD or at least one of them is a migrant allele, \cr
C) \tab fCij = Probability that two alleles chosen from i and j are IBD or both are migrant alleles, \cr
N) \tab fNij = Probability that two alleles chosen from i and j are both native alleles. \cr
}
Matrices fB and fN may be needed to compute the conditional kinship fD, defined as the conditional probability that two randomly chosen alleles are IBD given that both originate from native founders.
A native founder is an individual with unkown parents belonging to \code{thisBreed}.
}

\value{
A list with the following components:
\item{fA}{matrix with kinships fAij (if "A" is in \code{method})}
\item{fB}{matrix with kinships fBij (if "B" or "D" are in \code{method})}
\item{fC}{matrix with kinships fCij (if "C" is in \code{method})}
\item{fN}{matrix with kinships fNij (if "N" or "D" are in \code{method})}
The list has class attribute \code{"kinMatrices"}. The additional attribute \code{condProb} specifies for conditional kinships the matrices needed for computing them.
}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
Kin   <- kin(Pedig, thisBreed="Hinterwaelder", method=c("A","D"), keep.only=keep)

#Native Ne = 49.5 (estimated from 6 previous generations)

#for(i in names(Kin))Kin[[i]]<-as(Kin[[i]],"matrix")
	
}



\author{Robin Wellmann}
