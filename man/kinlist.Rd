\name{kinlist}
\alias{kinlist}

\title{Combines Kinship Matrices into a List}

\description{
Combines matrices or objects of class kinMatrices into a list}
\usage{
kinlist(...)
}

\arguments{
 \item{...}{One or more matrices or objects of class kinMatrices}
}
\value{
List of class kinMatrices
}
\details{
Combines matrices or objects of class kinMatrices into a list
}
\examples{
data(PedigWithErrors)
data(Phen)
Phen  <- Phen[80:140,]
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
fA    <- pedIBD(Pedig, keep.only=keep)
fD    <- pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=keep)
x     <- pedIBDorM(Pedig, thisBreed="Hinterwaelder", keep.only=keep)
Kin   <- kinlist(fA=fA, fB=x$pedIBDorM, fC=x$pedIBDorMM, fD=fD)
help.opticont(Kin, Phen)
}

\author{Robin Wellmann}

