\name{opticont}
\Rdversion{1.1}
\alias{opticont}
\title{Calculates Optimum Contributions of Selection Candidates
}
\description{
Calculates optimum genetic contributions of selection candidates to the next generation. 
}
\usage{
opticont(method, K, phen, con=list(), solver="auto", quiet=FALSE, 
    make.definite=solver=="csdp", ...)}

\arguments{
\item{method}{Possible values are \code{"min.VAR"}, and \code{"max.VAR"}, where \code{VAR} is the name of a column in data frame \code{phen}, or \code{"min.KIN"}, where \code{KIN} is the name of a kinship as defined by function \link{kinlist}. Use \link{help.opticont} to see the available objective functions.
}
\item{K}{List created by function \link{kinlist}, containing e.g. kinship matrices of the selection candidates.}
\item{phen}{Data frame with IDs of selection candidates in column 1, and sex in column 2. The other columns may contain breeding values or migrant contributions. The sex is coded as \code{'male'} and \code{'female'}.}
\item{con}{List defining the constraints. The compomnents are described in the Details section. If a component is missing, then the respective constraint is not applied. Use \link{help.opticont} to see the available constraints. 
}

\item{solver}{Name of the algorithm for optimization. Available solvers are  \code{"alabama"}, \code{"cccp"}, \code{"cccp2"}, \code{"csdp"},   and \code{"slsqp"}. The default \code{"auto"} means that the solver is chosen automatically. The solvers are described in the Details section.
}

\item{quiet}{If \code{quiet=FALSE} then detailed information is shown.}
\item{make.definite}{If \code{make.definite=TRUE} then all non-definite matrices are approximated by positive definite matrices before optimization.
This is the default setting for the solver \code{csdp}.
}
\item{...}{Tuning parameters of the solver. The available parameters depend on the solver and will be printed when function \code{opticont} is used with default values. An overview is given in the Details section. 
}
}
\details{
Computation of optimum genetic contributions of selection candidates. 

\bold{Constraints}

Constraints are defined in argument \code{con}, which is a list whose components may have the following names:

\bold{ub.KIN}: Upper bound for the mean kinship  in the offspring, where \code{KIN} must be replaced by the name of a kinship as defined by function \link{kinlist}. Upper bounds for an arbitrary number of different kinships may be provided. 

\bold{lb}: Either a named vecor of the form \code{c(M=a, F=b)} containing lower bounds for the contributions of males (\code{a}) and females (\code{b}), or a N-vector containing  the minimum permissible contribution of each selection candidate. The default is \code{c(M=0, F=0)}.

\bold{ub}: Either a named vecor of the form \code{c(M=a, F=b)} containing upper bounds for the contributions of males (\code{a}) and females (\code{b}), or a N-vector containing  the maximum permissible contribution of each selection candidate. For \code{M=-1} (\code{F=-1}) it is assumed that all males (females) have equal contributions to the offspring.  If a number is \code{NA} then the number of offspring for that sex/individual is not bounded. The default is \code{c(M=NA, F=NA)}.

\bold{lb.VAR}: Lower bound for the expected mean value of variable \code{VAR} from data frame \code{phen} in the offspring. For example \code{lb.BV=a} defines a lower bound for the mean breeding value in the offspring to be \code{a} if data frame \code{phen} has column \code{BV} with breeding values of the parents. Lower bounds for an arbitrary number of variables can be defined.

\bold{ub.VAR}:   Upper bound for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring.
For example \code{ub.MC=a} defines the upper bound for the genetic contributions from migrant breeds in the offspring to be \code{a} if data frame \code{phen} has column \code{MC} with migrant contributions for the parents. Upper bounds for an arbitrary number of variables can be defined.

\bold{eq.VAR}: Equality constraint for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring. Equality constraints for an arbitrary number of variables can be defined.


\bold{Solver}


 \code{"alabama"}: The augmented lagrangian minimization algorithm \link[alabama]{auglag} from package \code{alabama}.
 The method combines the objective function and a penalty for each constraint  into a single function. This modified objective function is then passed to another optimization algorithm with no constraints. If the constraints are violated by the solution of this sub-problem, then the size of the penalties is increased and the process is repeated. The default methods for the uncontrained optimization in the inner loop is the quasi-Newton method called \code{BFGS}.
The available parameters used for the outer loop are described in the details section of the help page of function \link[alabama]{auglag}. The available parameters used for the inner loop are described in the details section of the help page of function \link[stats]{optim}.
 
\code{"cccp", "cccp2"}: Function  \link[cccp]{cccp} from package \code{cccp} for solving cone constrained convex programs. For \code{cccp} quadratic constraints are defined as second order cone constraints, whereas for \code{cccp2} quadratic constraints are defined by functions. The implemented algorithms are partially ported from CVXOPT. The parameters are those from function \link[cccp]{ctrl}. They are among others the maximum count of iterations as an integer value (\code{maxiters}), the feasible level of convergence to be achieved (\code{feastol}) and whether the solver's progress during the iterations is shown (\code{trace}). If numerical problems are encountered increase the optimization parameter \code{feastol} or reduce parameter \code{stepadj}.

 \code{"csdp"}: The problem is reformulated as a semidefinite programming problem and solved with the CSDP library.
 Non-definite matrices are approximated by positive definite matrices. This solver is not suitable when the objective is to minimize kinship at native alleles.
 Available parameters are described in the CSDP User's Guide: \code{https://projects.coin-or.org/Csdp/export/49/trunk/doc/csdpuser.pdf}.
 
 \code{"slsqp"}: The sequential (least-squares) quadratic programming (SQP) algorithm \link[nloptr]{slsqp} for  gradient-based optimization from package \code{nloptr}. The algorithm optimizes successive second-order (quadratic/least-squares) approximations of the objective function, with first-order (affine) approximations of the constraints.  Available parameters are described in \link[nloptr]{nl.opts}


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
\link{help.opticont} to see available objective functions and constraints.
}

\value{
A list with class \code{"opticont"}
 which has component \code{parent}.  This is the data frame \code{phen} with the additional column \code{oc} containing the optimum genetic contribution of each selection candidate to the next generation, \code{lb} containing the lower bounds, and \code{ub} containing the upper bounds.


%Note: In windows R3.2.4, if there are a lot / large objects in memory, eigen or chol, which are called from this function, may return NULL and subsequently no result is obtained. This is a bug in R.
}


\examples{
#################################################
# Example 1: Advanced OCS using pedigree data   #
# Objective: maximize genetic gain              #
# Constraints:                                  #
#   - mean kinship                              #
#   - mean kinship at native alleles            #
#   - genetic contributions from other breeds   #
#################################################
data(PedigWithErrors)
data(Phen)

keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
Pedig$MC <- 1-pedBreedComp(Pedig, thisBreed="Hinterwaelder")$native
Phen  <- merge(Pedig, Phen[,c("Indiv", "BV")], by="Indiv")
Kin   <- kinlist(pKin    = pedIBD(Pedig, keep.only=keep), 
                 pKinatN = pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=keep))

head(Phen)
cor(Phen$MC, Phen$BV)
help.opticont(Kin, Phen)

# Compute offspring parameters for unselected population
noSel      <- opticont(method="min.pKin", K=Kin, phen=Phen, con=list(ub=c(M=-1, F=-1)))
noSel.s    <- summary(noSel)
round(noSel.s[,c("pKin","pKinatN", "meanMC", "meanBV")],4)
meanMC     <- noSel.s$meanMC
meanKin    <- noSel.s$pKin
meanKinatN <- noSel.s$pKinatN
meanBV     <- noSel.s$meanBV

# Define Constraints
Ne <- 100

con            <- list(ub=c(M=NA, F=-1))
con$ub.pKin    <- meanKin    + (1-meanKin   )*(1/(2*Ne))
con$ub.pKinatN <- meanKinatN + (1-meanKinatN)*(1/(2*Ne))
con$ub.MC      <- 0.97*meanMC

# Compute the genetic progress achievable
maxBV   <- opticont("max.BV", K=Kin, phen=Phen, con=con)
maxBV.s <- summary(maxBV)
maxBV.s$meanBV 
# [1] 0.5428925

# Get optimum contributions of sires
Sire <- maxBV$parent[maxBV$parent$Sex=="male",]
ord  <- order(Sire$oc, decreasing=TRUE)
head(Sire[ord,])

#################################################
# Example 2: Advanced OCS using genotype data   #
# Objective: minimize inbreeding                #
# Constraints:                                  #
#   - breeding values                           #
#   - mean kinship at native alleles            #
#   - genetic contributions from other breeds   #
#################################################

data(map)
data(Cattle)
dir   <- system.file("extdata", package = "optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")
Kin   <- kinlist(sKin    = segIBD(files, map, minSNP=20, minL=2.0), 
                 sKinatN = segIBDatN(files, Cattle, map, thisBreed="Angler", 
                                     ubFreq=0.01, minSNP=20, minL=2.0))
Haplo <- haplofreq(files, Cattle, map, thisBreed="Angler", 
                   minSNP=20, minL=2.0, ubFreq=0.01, what="match")
Comp  <- segBreedComp(Haplo$match, map)
Comp$MC <- 1-Comp$native
Phen  <- merge(Cattle, Comp[,c("Indiv", "MC")], by="Indiv", all=FALSE)

help.opticont(Kin, Phen)

cor(Phen$MC, Phen$BV,use="complete.obs")
#[1] 0.5033714

# Compute offspring parameters for unselected population
noSel      <- opticont(method="min.sKin", K=Kin, phen=Phen, con=list(ub=c(M=-1, F=-1)))
noSel.s    <- summary(noSel)
round(noSel.s[,c("sKin","sKinatN", "meanMC", "meanBV")],4)
meanMC     <- noSel.s$meanMC
meanKin    <- noSel.s$sKin
meanKinatN <- noSel.s$sKinatN
meanBV     <- noSel.s$meanBV

# Define Constraints
Ne <- 100

con            <- list(ub=c(M=NA, F=-1))
con$ub.sKinatN <- meanKinatN + (1-meanKinatN)*(1/(2*Ne))
con$ub.MC      <- 0.97*meanMC
con$lb.BV      <- meanBV

# Compute the smallest mean kinship achievable
minKin   <- opticont("min.sKin", K=Kin, phen=Phen, con=con)
minKin.s <- summary(minKin)
minKin.s$sKin
# [1] 0.03881304

# Get optimum contributions of sires
Sire <- minKin$parent[minKin$parent$Sex=="male",]
ord  <- order(Sire$oc, decreasing=TRUE)
head(Sire[ord,])
}

\references{
Borchers, B. (1999). CSDP, A C Library for Semidefinite Programming Optimization Methods and Software 11(1):613-623
\code{http://euler.nmt.edu/~brian/csdppaper.pdf}

Kraft, D. (1988). A software package for sequential quadratic programming, Technical Report DFVLR-FB 88-28, Institut fuer Dynamik der Flugsysteme, Oberpfaffenhofen, July 1988.

Lange K, Optimization, 2004, Springer.

Madsen K, Nielsen HB, Tingleff O, Optimization With Constraints, 2004, IMM, Technical University of Denmark.
}


\author{Robin Wellmann}
