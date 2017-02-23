\name{pedIBD}
\Rdversion{1.1}
\alias{pedIBD}
\title{Calculates the Pedigree Based Kinship Matrix 
}
\description{                                                 
Calculates the \bold{ped}igree based probability of alleles to be \bold{IBD}. This
pedigree based kinship matrix is also called coancestry matrix and is half the additive relationship matrix.
}
\usage{
pedIBD(Pedig, keep.only=NULL, keep=keep.only)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA. Both parents must either be missing or present. If this is not the case use \link{prePed}.}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
}

\details{
Computation of pedigree based kinship matrix f which is half the additive relationship matrix. For individuals i and j it is defined as
\tabular{l}{
 fij = Probability that two alleles chosen from individuals i and j are IBD. \cr
}
This function is a wrapper function for \code{\link[nadiv:makeA]{makeA}} from package  \code{\link[nadiv:nadiv-package]{nadiv}} 
}

\value{
A sparse matrix.
}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
pedA  <- pedIBD(Pedig, keep.only=keep)
pedA  <- as(pedA,"matrix")
}



\author{Robin Wellmann}
