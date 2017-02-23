\name{pedBreedComp}
\Rdversion{1.1}
\alias{pedBreedComp}
\title{Calculates the Pedigree Based Breed Composition of Individuals
}
\description{
Computes for every individual the genetic contribution from native founders and from other breeds according to the pedigree.
}
\usage{
pedBreedComp(Pedig, thisBreed)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, (4) Sex, and (5) Breed name. Missing parents are coded as \code{NA}. All animals have no parent or both parents missing. If this is not the case use \link{prePed}.}
\item{thisBreed}{Name of this breed as denoted in column (5) of the pedigree.}
}

\details{
For every individual the genetic contribution from native founders and from other breeds is computed. It is the fraction of genes that originate from the respective breed.
}

\value{
Data frame with one row for each individual and the following columns
\item{Indiv}{IDs of the individuals}
\item{native}{Native Contribution: The genetic contribution from native founders.}
\item{...}{Genetic contributions from other breeds, one column for each breed. The columns are ordered, so that the most influential breeds come first.}
}


\examples{
data(ExamplePed)
Pedig    <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
cont     <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")
Pedig$MC <- 1-cont$native
cont[1000:1010,2:5]

contByYear <- conttac(cont, Pedig$Born, use=Pedig$Breed=="Hinterwaelder", mincont=0.04, long=FALSE)
round(contByYear,2)

barplot(contByYear,ylim=c(0,1), col=1:10, ylab="genetic contribution",
        legend=TRUE, args.legend=list(x="bottomleft",cex=0.6))
}



\author{Robin Wellmann}
