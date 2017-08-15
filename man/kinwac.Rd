\name{kinwac}
\Rdversion{1.1}
\alias{kinwac}
\title{Calculates Kinships With Age Cohorts}

\description{
Calculates for every individual it's mean kinship with the age cohort to which it belongs (excluding the kinship with itself). 
}
\usage{
kinwac(cand)}

\arguments{
\item{cand}{R-Object created with function candes containing phenotypes and kinship information on the individuals. Data frame \code{cand$phen} must include additional column \code{Born} with Year-of-Birth.}
}

\details{
The function computes for every individual it's mean kinship with the age cohort to which it belongs. This can be computed for different notions of kinship simultaneously. If an individual has high kinship with it's age cohort then it should not be extensively used for breeding if the aim is to conserve genetic diversity.
}

\value{
The data frame \code{cand$phen} with an additional column for each kinship and each native kinship, containig the average kinship of the individual with individuals from the same age-cohort.
}

\examples{
data(ExamplePed)
Pedig   <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
phen    <- Pedig[Pedig$Breed=="Hinterwaelder",]
pKin    <- pedIBD(Pedig)
pKinatN <- pedIBDatN(Pedig, thisBreed="Hinterwaelder")
pop     <- candes(phen=phen, pKin=pKin, pKinatN=pKinatN, quiet=TRUE, N=Inf)
Kin     <- kinwac(pop)
Kin[1001:1010,]
}



\author{Robin Wellmann}
