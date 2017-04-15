\name{Kin}
\docType{data}
\alias{Kin}
\title{A list of class kinMatrices containing kinship matrices}
\description{
 A list of class kinMatrices containing kinship matrix fA and the matrices and the matrices needed for computing fD for the animals included in datat set \link{Phen}. 
}
\usage{Kin}
\format{A list of class kinMatrices.}
\keyword{datasets}
\examples{
data(PedigWithErrors)
data(Phen)
Phen[1:5,]

keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
fA    <- pedIBD(Pedig, keep.only=keep)
fD    <- pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=keep)
Kin2   <- kinlist(fA=fA, fD=fD)

data(Kin)

identical(Kin, Kin2)
#[1] TRUE

}