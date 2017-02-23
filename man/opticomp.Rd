\name{opticomp}
\Rdversion{1.1}
\alias{opticomp}
\title{Calculates the Optimum Breed Composition}

\description{
Calculates optimum contributions of breeds to a hypothetical synthetic population with maximum diversity. Additionally the average kinship within and
between breeds and the genetic distances between breeds are computed.
}
\usage{
opticomp(f, Breed, obj.fun="NGD", lb=NULL, ub=NULL, ...)
}

\arguments{
\item{f}{Kinship matrix (e.g. a segment based kinship matrix, or a IBS based kinship matrix).}
\item{Breed}{Vector containing the breed name of each individual.}
\item{obj.fun}{The objective function to be maximized. For \code{"NGD"} the diversity 1-\bold{c}'\bold{f}\bold{c} is maximized, or equivalently, the mean kinship \bold{c}'\bold{f}\bold{c} is minimized, where \bold{f} is the matrix containing the mean kinships within and between breeds.
For \code{"NTD"} the term \bold{c}'(\bold{1}-\bold{F})+\bold{c}'(\bold{F}\bold{1}' - 2\bold{f} + \bold{1}\bold{F}')\bold{c} is maximized, where \bold{F}=diag(\bold{f}). This puts more weight on between population diversity.
}
\item{lb}{Named vector providing lower bounds for the contributions of the breeds can be provided. The names of the components are the breed names. The default \code{lb=NULL} means that the lower bound is 0 for all breeds.}
\item{ub}{Named vector providing upper bounds for the contributions of the breeds can be provided. The names of the components are the breed names. The default \code{ub=NULL} means that the upper bound is 1 for all breeds.}
\item{...}{Further parameters passed to the solver \link[quadprog]{solve.QP} of R package \code{quadprog}. }
}

\details{
Calculates optimum contributions of breeds to a hypothetical synthetic population with maximum diversity. Additionally the average kinship within and
between breeds and the genetic distances between breeds are computed.
}

\value{
A list with the following components:
\item{bc}{Vector with optimum contributions of breeds to a synthetic population with maximum  diversity}
\item{value}{The value of the objective function, i.e. the maximum diversity that can be achieved.}
\item{f}{Matrix containing the mean kinships within and between breeds.}
\item{Dist}{Genetic distances between breeds.}
}

\author{Robin Wellmann}

\references{
Wellmann, R., Bennewitz, J., Meuwissen, T.H.E. (2014) A unified approach to characterize and conserve adaptive and neutral genetic diversity in subdivided populations. {Genetics Selection Evolution}. 69, e16
}


\examples{
library(optiSel)
data(map) 
data(Cattle)
dir   <- system.file("extdata", package = "optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")

#####################################################################
#   Find the optimum breed composition using segment based kinship  #
#####################################################################
IBD <- segIBD(files, minSNP=15, map=map, minL=1000, unitL="kb", unitP="kb")
mb  <- opticomp(IBD, Breed=Cattle$Breed, obj.fun="NGD")

#### Optimum breed composition: ###
round(mb$bc,3)
#    Angler Fleckvieh  Holstein   Rotbunt
#     0.512     0.404     0.084     0.000 

#### Average kinships within and between breeds: ###
round(mb$f,4)
#          Angler Fleckvieh Holstein Rotbunt
#Angler    0.0709    0.0180   0.0629  0.0664
#Fleckvieh 0.0180    0.0950   0.0147  0.0158
#Holstein  0.0629    0.0147   0.1269  0.1130
#Rotbunt   0.0664    0.0158   0.1130  0.1327

#### Genetic distances between breeds: ###
round(mb$Dist,4)
#          Angler Fleckvieh Holstein Rotbunt
#Angler    0.0000    0.2549   0.1896  0.1882
#Fleckvieh 0.2549    0.0000   0.3102  0.3132
#Holstein  0.1896    0.3102   0.0000  0.1296
#Rotbunt   0.1882    0.3132   0.1296  0.0000

#####################################################################
#   The optimum breed composition depends on the kinship matrix     #
#   and the objective function:                                     #
#####################################################################

bc <- opticomp(IBD, Breed=Cattle$Breed, obj.fun="NTD")$bc
round(bc,3)
#    Angler Fleckvieh  Holstein   Rotbunt
#     0.248     0.445     0.176     0.131 

}

