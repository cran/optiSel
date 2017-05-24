\name{opticont4mb}
\Rdversion{1.1}
\alias{opticont4mb}
\title{Calculates Optimum Contributions of Selection Candidates using Multi-Breed Genotype Data
}
\description{
Calculates optimum genetic contributions for selection candidates from one breed using multi-breed genotype data. Genotype data from multiple breeds may be used in order to increase the genetic distance between the breed of interest (\code{thisBreed}) and other breeds. 
}
\usage{
opticont4mb(method, K, phen, bc, thisBreed=names(bc)[1], con=list(), 
    solver="cccp", quiet=FALSE, make.definite=solver=="csdp", ...)}

\arguments{
\item{method}{Possible values are \code{"min.VAR"}, and \code{"max.VAR"}, where \code{VAR} is the name of a column in data frame \code{phen}, or \code{"min.KIN"}, or \code{"min.KIN.acrossBreeds"}, where \code{KIN} is the name of a kinship as defined by function \link{kinlist}. Use \link{help.opticont4mb} to see the available objective functions. If kinship \code{KIN} is available for all animals from the multi-breed population, then \code{"min.KIN.acrossBreeds"} minimizes the kinship in the multi-breed population by optimizing the contributions of selection candidates from the breed of interest.
}
\item{K}{List created by function \link{kinlist} containing kinships of genotyped individuals.}
\item{phen}{Data frame with  one row for each animal from the multi-breed population that is to be included in the analysis.  The animal IDs is in column 1 (named \code{Indiv}) and the sex is in column 2 (named \code{Sex}).  The sex is coded as \code{'male'} and \code{'female'}. Column \code{Breed} contains the breed name of every genotyped animal. Further columns contain e.g. breeding values or migrant contributions that may be used for defining linear constraints.
}
\item{bc}{Named vector containing the proportion of every genotyped breed in the hypothetical multi-breed offspring population. The names  of the components are the breed names. Note that only the contributions of selection candidates from \code{thisBreed} will be optimized. Animals from other breeds have fixed contributions. }
\item{thisBreed}{The breed to which the selection candidates belong. }
\item{con}{List defining the constraints. The components are described in the Details section. If a component is missing, then the respective constraint is not applied. Use \link{help.opticont4mb} to see the available constraints. 
}

\item{solver}{Name of the algorithm for optimization. Available solvers are  \code{"alabama"}, \code{"cccp"}, \code{"cccp2"}, \code{"csdp"},   and \code{"slsqp"}. The default is \code{"cccp"}. The solvers are described in the Details section.
}

\item{quiet}{If \code{quiet=FALSE} then detailed information is shown.}
\item{make.definite}{If \code{make.definite=TRUE} then all non-definite matrices are approximated by positive definite matrices before optimization.
This is the default setting for the solver \code{csdp}.
}
\item{...}{Tuning parameters of the solver. The available parameters depend on the solver and will be printed when function \code{opticont} is used with default values. An overview is given in the Details section. 
}
}
\details{
Computation of optimum genetic contributions for the selection candidates from one breed using multi-breed marker data. Marker data from multiple breeds may be used in order to increase the genetic distance between the breed of interest (\code{thisBreed}) and the other breeds. 

In this case a hypothetical subdivided population is considered consisting of purebred offspring of genotyped individuals. That is, the offspring population consists of several breeds with specified breed proportions (e.g.
20\% Angler cattle, 40\% Holstein cattle, and 40\% Fleckvieh cattle). Only the contributions of the selection candidates from \code{thisBreed} will be optimized. Animals from other breeds have equal contributions.

The aim is to reduce the average genomic relationship in this multi-breed  population since this causes the genetic distance between \code{thisBreed} and other breeds to increase. This may increase the conservation value of the breed.

If managing diversity across breeds is not intended then function \link{opticont} could be used instead.

\bold{Constraints}

A list possibly containing the following components providing the constraints:

\bold{ub.KIN}: Upper bound for the mean kinship in the offspring, where \code{KIN} must be replaced by the name of a kinship as defined by function \link{kinlist}. Use \link{help.opticont4mb} to see available methods.

\bold{ub.KIN.acrossBreeds}:  Upper bound for the mean kinship in the next generation of the multi-breed population, where \code{KIN} must be replaced by the name of a kinship as defined by function \link{kinlist}. Use \link{help.opticont4mb} to see available methods.

\bold{lb}: Either a named vecor of the form \code{c(M=a, F=b)} containing lower bounds for the contributions of males (\code{a}) and females (\code{b}) from \code{thisBreed}, or a named vector containing  the minimum permissible contribution of each selection candidate. The default is \code{c(M=0, F=0)}.

\bold{ub}: Either a named vecor of the form \code{c(M=a, F=b)} containing upper bounds for the contributions of males (\code{a}) and females (\code{b}) from \code{thisBreed}, or a named vector containing  the maximum permissible contribution of each selection candidate. For \code{M=-1} (\code{F=-1}) it is assumed that all males (females) have equal contributions to the offspring.  If a number is \code{NA} then the number of offspring for that sex/individual is not bounded. The default is \code{c(M=NA, F=NA)}.

\bold{lb.VAR}: Lower bound for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring from \code{thisBreed}. For example \code{lb.BV=a} defines a lower bound for the mean breeding value in the offspring from \code{thisBreed} to be \code{a} if data frame \code{phen} has column \code{BV} with breeding values of the parents. Lower bounds for an arbitrary number of variables can be defined.

\bold{ub.VAR}: Upper bound for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring from \code{thisBreed}.
For example \code{ub.MC=a} defines the upper bound for the genetic contributions from migrant breeds in the offspring of \code{thisBreed}  to be \code{a} if data frame \code{phen} has column \code{MC} with migrant contributions for the parents. Upper bounds for an arbitrary number of variables can be defined.

\bold{eq.VAR}: Equality constraint for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring from \code{thisBreed}.
For example \code{eq.MC=a} forces the genetic contribution from migrant breeds in the offspring from \code{thisBreed} to be \code{a} if data frame \code{phen} has column \code{MC} with migrant contributions for the parents. Equality constraints for an arbitrary number of variables can be defined.



\bold{Solver}


 \code{"alabama"}: The augmented lagrangian minimization algorithm \link[alabama]{auglag} from package \code{alabama} is used.
 That is, the method combines the objective function and a penalty for each constraint  into a single function. This modified objective function is then passed to another optimization algorithm with no constraints. If the constraints are violated by the solution of this sub-problem, then the size of the penalties is increased and the process is repeated. The default methods for the uncontrained optimization in the inner loop is the quasi-Newton method called \code{BFGS}. The available parameters used for the outer loop are described in the details section of the help page of function \link[alabama]{auglag}. The available parameters used for the inner loop are described in the details section of the help page of function \link[stats]{optim}.
 
\code{"cccp", "cccp2"}: Function  \link[cccp]{cccp} from package \code{cccp} for solving cone constrained convex programs is used. For \code{cccp} quadratic constraints are defined as second order cone constraints. This solver is not suitable if computation of the Cholesky decomposition fails.  For \code{cccp2} quadratic constraints are defined by functions. The implemented algorithms are partially ported from CVXOPT. The parameters are those from function \link[cccp]{ctrl}. They are among others the maximum count of iterations as an integer value (\code{maxiters}), the feasible level of convergence to be achieved (\code{feastol}) and whether the solver's progress during the iterations is shown (\code{trace}). If numerical problems are encountered increase the optimization parameter \code{feastol} or reduce parameter \code{stepadj}.

 \code{"csdp"}: The problem is reformulated as a semidefinite programming problem and solved with the CSDP library.
 Non-definite matrices are approximated by positive definite matrices. This solver is not suitable when the objective is to minimize kinship at native alleles. Available parameters are described in the CSDP User's Guide: \code{https://projects.coin-or.org/Csdp/export/49/trunk/doc/csdpuser.pdf} .
 
 \code{"slsqp"}: The sequential (least-squares) quadratic programming (SQP) algorithm \link[nloptr]{slsqp} for  gradient-based optimization from package \code{nloptr} is used. The algorithm optimizes successive second-order (quadratic/least-squares) approximations of the objective function, with first-order (affine) approximations of the constraints. Available parameters are described in \link[nloptr]{nl.opts}.

\bold{Remark}

If the function does not provide a valid result due to numerical problems then try the following modifications:

\tabular{ll}{
\code{*} \tab modify the optimization parameters,\cr
\code{*} \tab use another \code{solver},\cr
\code{*} \tab change the order of the kinship constraints if more than one kinship is constrained,\cr
\code{*} \tab define upper or lower bounds instead of equality constraints.\cr
\code{*} \tab increase the upper bounds for the kinships.\cr
}
Validity of the result can be checked with function \link{summary.opticont}. Use
\link{help.opticont4mb} to see available objective functions and constraints.
}

\value{
A list with class \code{"opticont"}
 which has component \code{parent}.  This is the data frame \code{phen} but includes ony the animals from the breed of interest. It has the additional column \code{oc} containing the optimum genetic contribution of each selection candidate to the next generation, \code{lb} containing the lower bounds of the optimum contributions, and \code{ub} containing the upper bounds.

}


\examples{
data(map) 
data(Cattle)
dir  <- system.file("extdata", package = "optiSel")
files<- file.path(dir, paste("Chr", 1:2, ".phased", sep=""))

### Compute genomic kinship and genomic kinship at native segments
G    <- segIBD(files, map, minSNP=20, minL=3.0)
GN   <- segIBDatN(files, Cattle, map, thisBreed="Angler", refBreeds="others", 
           ubFreq=0.02, minSNP=20, minL=3.0, lowMem=TRUE)
Kin  <- kinlist(G=G, GN=GN)

### Compute migrant contributions of selection candidates 
Haplo<- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="others",
           minSNP=20, minL=3.0, ubFreq=0.02, what="match")
Comp <- segBreedComp(Haplo$match, map)
Cattle$MC <- NA
Cattle[rownames(Comp), "MC"] <- 1-Comp$native
apply(Comp[,-1],2,mean)
#     native           F           H           R 
#0.551844104 0.009739393 0.202216271 0.236200232 


########################################
#  Find optimum breed contributions    #
########################################
lb <- c(Angler=0.10, Holstein=0.20, Fleckvieh=0.20)
bc <- opticomp(G, Breed=Cattle$Breed, obj.fun="NGD", lb=lb)$bc
round(bc,3)
#   Angler Fleckvieh  Holstein   Rotbunt 
#    0.355     0.445     0.200     0.000 

########################################
#  Check available objective functions #
#  and constraints                     #
########################################

help.opticont4mb(Kin, Cattle)


##################################################################
#   Compute the minimum segment based kinship achievable         #
#  across breeds while constraining it within the breed          #
##################################################################

con  <- list(ub.G=0.05, ub=c(M=NA, F=-1))
minG <- opticont4mb("min.G.acrossBreeds", Kin, Cattle, bc, thisBreed="Angler", con=con, trace=FALSE)
minG.s <- summary(minG)
minG.s[,c("G.acrossBreeds")]
#[1] 0.02289039

##################################################################
# Compute the genetic progress achievable while constraining     #
#      segment based kinship  within and across breeds           #
#                     and migrant contributions                  #
##################################################################

con          <- list(ub.G=0.05, ub.G.acrossBreeds=0.026, ub.MC=0.32,  ub=c(M=NA, F=-1))
maxBV.G.MC   <- opticont4mb("max.BV", Kin, Cattle, bc, thisBreed="Angler", con=con, trace=FALSE)
maxBV.G.MC.s <- summary(maxBV.G.MC)
maxBV.G.MC.s$meanBV 
# [1] 0.2851194

##################################################################
#    Compute the minimum achievable kinship at native alleles    #
#    while constraining kinship within and across breeds         #
#    and migrant contributions                                   #
##################################################################

con   <- list(ub.G=0.05, ub.G.acrossBreeds=0.026, ub.MC=0.32, ub=c(M=NA, F=-1))
minGN <- opticont4mb("min.GN", Kin, Cattle, bc, thisBreed="Angler", con=con, solver="slsqp")
minGN.s <- summary(minGN)
minGN.s$GN
#[1] 0.04114953


##################################################################
# Summary statistics from different optimizations                #
# can be combined in a data frame. The most important parameters #
# are printed for comparison:                                    #
##################################################################
Res <- rbind(minG.s, maxBV.G.MC.s, minGN.s)
format(Res[,c("valid","meanBV", "meanMC", "G.acrossBreeds", "G", "GN")],digits=4)

#           valid  meanBV meanMC G.acrossBreeds       G      GN
#minG        TRUE -0.4329 0.3379        0.02289 0.03451 0.03986
#maxBV.G.MC  TRUE  0.2851 0.3200        0.02475 0.05000 0.06830
#minGN       TRUE -0.5321 0.3200        0.02312 0.03600 0.04115

cor(cbind(minG$parent$oc, maxBV.G.MC$parent$oc, minGN$parent$oc))
#          [,1]      [,2]      [,3]
#[1,] 1.0000000 0.2741736 0.8540842
#[2,] 0.2741736 1.0000000 0.3608900
#[3,] 0.8540842 0.3608900 1.0000000
}

\references{
Borchers, B. (1999). CSDP, A C Library for Semidefinite Programming Optimization Methods and Software 11(1):613-623
\code{http://euler.nmt.edu/~brian/csdppaper.pdf}

Kraft, D. (1988). A software package for sequential quadratic programming, Technical Report DFVLR-FB 88-28, Institut fuer Dynamik der Flugsysteme, Oberpfaffenhofen, July 1988.

Lange K, Optimization, 2004, Springer.

Madsen K, Nielsen HB, Tingleff O, Optimization With Constraints, 2004, IMM, Technical University of Denmark.
}


\author{Robin Wellmann}
