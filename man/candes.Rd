\name{candes}
\alias{candes}

\title{Candidate Description}

\description{
An R-Object is created containing all information describing the individuals, which are usually the selection candidates. Average kinships and trait values, and the available objective functions and constraints for optimum contribution selection are reported.
}

\usage{
candes(phen, N=250, quiet=FALSE, bc=NA, ...)
}

\arguments{
 \item{phen}{Data frame with column \code{Indiv} containing animal IDs and \code{Sex} containng sexes, coded as \code{'male'} and \code{'female'}. The other columns may contain traits, e.g. breeding values or native contributions, and column \code{Breed} with breed names for multi-breed evaluations.}
 \item{N}{Total number of selection candidates in the next generation. This value affects the increase in kinship per generation due to genetic drift. The default is \code{N=250}.}
  \item{quiet}{Should the report be suppressed?}
\item{bc}{Only needed if multi-breed data is provided. Named vector with breed contributions, with component names being the names of the breeds in \code{phen}.  It contains the proportion of each breed to a hypothetical multi-breed population for which the diversity across breeds should be managed. Alternatively, \code{bc} can be a character string containig the name of a kinship. In this case, optimum contributions of the breeds are determined automatically so that the mean kinship across breeds is minimized.
}
 \item{...}{One or more objects of class \code{'matrix'}, \code{'quadFun'}, or \code{'ratioFun'}.}
}
\value{
List of class candes with the following components:
\item{phen}{Data frame containing phenotypes and individual IDs.}
\item{mean}{Data frame containing estimates of the 'current' mean values of the parameters in a population consisting of \code{N} individuals for which the individuals in \code{cand$phen} are representative.}
\item{current}{Data frame containing current mean values of the parameters and the permitted names for the objective functions and constraints.}
\item{bc}{Character vector with optimum breed contributions (see above).}
\item{...}{Objects of class \code{'quadFun'}, or \code{'ratioFun'}, one for each additional parameter. These objects define the functions needed to estimate the mean kinships and mean native kinships in the offspring.}
}

\details{
The arguments are checked for plausibility. An R-Object is created containing all information describing the candidates. Average kinships and trait values are estimated and reported. Males and females contribute equally to these estimates. The available objective functions and constraints for optimum contribution selection are reported.
}

\examples{
data(PedigWithErrors)
data(Phen)
Phen  <- Phen
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
pKin    <- pedIBD(Pedig, keep.only=keep)
pKinatN <- pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=keep)
cand    <- candes(Phen, pKin=pKin, pKinatN=pKinatN)

cand$current[,c("Type", "Var", "Breed", "Val", "Name")]

#             Type     Var   Breed           Val    Name
#BV          trait      BV missing -4.903304e-17      BV
#NC          trait      NC missing  4.235493e-01      NC
#pKin      kinship    pKin missing  2.197712e-02    pKin
#pKinatN nat. kin. pKinatN missing  7.531245e-02 pKinatN

### Shorthand:

cand$mean
#             BV        NC       pKin    pKinatN
#1 -4.903304e-17 0.4235493 0.02197712 0.07531245

cand$mean$pKin
#[1] 0.02197712

#Object cand can now be used as an argument of function opticont.
}

\author{Robin Wellmann}

