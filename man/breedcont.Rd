\name{breedcont}
\Rdversion{1.1}
\alias{breedcont}
\title{breed contributions 
}
\description{
Computes for every individual the genetic contribution from other breeds
}
\usage{
breedcont(Pedig, thisBreed)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, (4) Sex, and (5) Breed name. Missing parents are coded as \code{NA}. All animals have no parent or both parents missing. If this is not the case use \link{prePed}.}
\item{thisBreed}{Name of this breed as denoted in column (5) of the pedigree.}
}

\details{
For every individual the genetic contribution from other breeds is computed. It is the fraction of genes that originate from the respective breed.
}

\value{
Data frame containing the genetic contribution from every breed to every individual.
}


\examples{
data(ExamplePed)
Pedig    <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
cont     <- breedcont(Pedig, thisBreed="Hinterwaelder")
Pedig$MC <- rowSums(cont)
cont[1:10,]

contByYear <- conttac(cont, cohort=Pedig$Born, use=Pedig$Breed=="Hinterwaelder", mincont=0.04)
round(contByYear,2)

barplot(contByYear, ylim=c(0,1), col=1:10, ylab="genetic contribution",
        legend=TRUE, args.legend=list(x="topleft",cex=0.6))
}



\author{Robin Wellmann}
