\name{kinwac}
\Rdversion{1.1}
\alias{kinwac}
\title{Calculates Kinships With Age Cohorts}

\description{
Calculates for every individual it's mean kinship with the age cohort to which it belongs (excluding the kinship with itself). 
}
\usage{
kinwac(K, Pedig)}

\arguments{
\item{K}{Named list containing kinship matrices of the individuals to be used for computing mean kinships. Typically only individuals from the breed of interest are included.}
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, and column \code{Born} containing for every individual the age cohort to which it belongs (typically year of birth).}
}

\details{
The function computes for every individual it's mean kinship with the age cohort to which it belongs. This can be computed for different notions of kinship simultaneously. If an individual has high kinship with it's age cohort then it should not be extensively used for breeding if the aim is to conserve genetic diversity.
}

\value{
A data frame with column names \code{c("cohort",  "I", names(K))} containing for every individual and every matrix provided in \code{K} the average value of the individual with the birth cohort to which it belongs.


The data frame has class attribute "kinwac".
}


\examples{
data(ExamplePed)
Pedig    <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
Kinships <- kinlist(pedIBD=pedIBD(Pedig), pedIBDatN=pedIBDatN(Pedig, thisBreed="Hinterwaelder"))
Kin      <- kinwac(Kinships, Pedig=Pedig)
Kin[1001:1010,]
}



\author{Robin Wellmann}
