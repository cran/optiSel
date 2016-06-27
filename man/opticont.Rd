\name{opticont}
\Rdversion{1.1}
\alias{opticont}
\title{Computation of optimum genetic contributions for the breeding animals
}
\description{
Computation of optimum genetic contributions for the breeding animals
}
\usage{
opticont(method, K, phen, con=list(), solver="cccp", quiet=FALSE, make.definite, ...)}

\arguments{
\item{method}{Possible values are \code{"min.VAR"}, and \code{"max.VAR"}, where \code{VAR} must be replaced by the name of a column in data frame \code{phen} or by the name of a kinship matrix contained in list \code{K} or by the name of a conditional kinship contained in the attribute \code{condProb} of list \code{K}.

For example, if \code{"fA"} is in \code{names(K)} then method \code{"min.fA"} minimizes the mean kinship \code{fA} in the offspring.
If \code{phen} has e.g. column \code{BV} with breeding values and column \code{MC} with migrant contributions then \code{"max.BV"} maximizes the mean breeding value of the offspring, and  \code{"min.MC"} minimizes the contribution of migrant breeds to the offspring. 
If \code{"fD"} is in \code{names(attr(K,"condProb"))} then \code{"min.fD"} minimizes the conditional kinship \code{fD} in the offspring. Optimization is done under the side constraints given in parameter \code{con}. 
}
\item{K}{A named list containing one or more NxN Kinship matrices for the N breeding animals, e.g. the output of function \link{kin}. Row names and column names of the matrices are the animal IDs.}
\item{phen}{Data frame with N rows and various columns, with animal IDs in column 1 and sex in column 2. The other columns may contain breeding values or migrant contributions. The sex is 1 for males and 2 for females.}
\item{con}{List possibly containing the following components providing the constraints:
\tabular{ll}{
\code{ub.KIN} \tab Upper bound for the mean kinship  in the offspring, where \code{KIN} must be replaced by the name of a kinship matrix included in list \code{K} or by the name of a conditional kinship included in \code{names(attr(K,"condProb"))}. Upper bounds for an arbitrary number of different kinships may be provided. 

For example, if \code{"fD"} is in \code{names(attr(K,"condProb"))}
then \code{ub.fD=a} defines the upper bound for the conditional kinship \code{fD} in the offspring to be \code{a}.
\cr
\code{lb}  \tab Either a named vecor of the form \code{c(M=a, F=b)} containing lower bounds for the contributions of males (\code{a}) and females (\code{b}), or a N-vector containing  the minimum permissible contribution of each breeding animal. The default is \code{c(M=0, F=0)}.\cr
\code{ub}  \tab Either a named vecor of the form \code{c(M=a, F=b)} containing upper bounds for the contributions of males (\code{a}) and females (\code{b}), or a N-vector containing  the maximum permissible contribution of each breeding animal. For \code{M=-1} (\code{F=-1}) it is assumed that all males (females) have equal contributions to the offspring.  If a number is \code{NA} then the number of offspring for that sex/individual is not bounded. The default is \code{c(M=NA, F=NA)}.
\cr
\code{lb.VAR}  \tab Lower bound for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring. For example \code{lb.BV=a} defines a lower bound for the mean breeding value in the offspring to be \code{a} if data frame \code{phen} has column \code{BV} with breeding values of the parents. Lower bounds for an arbitrary number of variables can be defined.
\cr
\code{ub.VAR}  \tab Upper bound for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring.
For example \code{ub.MC=a} defines the upper bound for the genetic contributions from migrant breeds in the offspring to be \code{a} if data frame \code{phen} has column \code{MC} with migrant contributions for the parents. Upper bounds for an arbitrary number of variables can be defined.
\cr
\code{eq.VAR}  \tab Equality constraint for the mean value of variable \code{VAR} from data frame \code{phen} in the offspring.
For example \code{eq.MC=a} forces the genetic contribution from migrant breeds in the offspring to be \code{a} if data frame \code{phen} has column \code{MC} with migrant contributions for the parents. Equality constraints for an arbitrary number of variables can be defined.
\cr
}
If a component is missing, then the respective constraint is not applied.
}

\item{solver}{The name of the algorithm chosen for optimization. 
All solvers should provide almost the same results but sometimes one of them fails due to numerical problems, or provides an approximate solution because a non-definite matrix had to approximated by a positive definite matrix. Available solvers are 

\tabular{ll}{
 \code{alabama} \tab The augmented lagrangian minimization algorithm \link[alabama]{auglag} from package \code{alabama} is used.
 That is, the method combines the objective function and a penalty for each constraint  into a single function. This modified objective function is then passed to another optimization algorithm with no constraints. If the constraints are violated by the solution of this sub-problem, then the size of the penalties is increased and the process is repeated. The default methods for the uncontrained optimization in the inner loop is the quasi-Newton method called \code{BFGS}.
 \cr
\code{cccp, cccp2} \tab Function  \link[cccp]{cccp} from package \code{cccp} for solving cone constrained convex programs is used. For \code{cccp} quadratic constraints are defined as second order cone constraints, whereas for \code{cccp2}
quadratic constraints are defined by functions.
The implemented algorithms are partially ported from CVXOPT.\cr
 \code{csdp} \tab The problem is reformulated as a semidefinite programming problem and solved with the CSDP library.
 Non-definite matrices are approximated by positive definite matrices. This solver is not suitable when the objective is to minimize a conditional kinship.
 \cr
 \code{slsqp} \tab The sequential (least-squares) quadratic programming (SQP) algorithm \link[nloptr]{slsqp} for  gradient-based optimization from package \code{nloptr} is used. The algorithm optimizes successive second-order (quadratic/least-squares) approximations of the objective function, with first-order (affine) approximations of the constraints.\cr
}
}

\item{quiet}{If \code{quiet=FALSE} then detailed information is shown.}
\item{make.definite}{If \code{make.definite=TRUE} then all non-definite matrices are approximated by positive definite matrices before optimization.
This is the default setting for the solver \code{csdp}.
}
\item{...}{Tuning parameters of the solver. The available parameters depend on the solver and will be printed when function \code{opticont} is used with default values:
\tabular{ll}{
 \code{alabama} \tab The available parameters used for the outer loop are described in the details section of the help page of function \link[alabama]{auglag}. The available parameters used for the inner loop are described in the details section of the help page of function \link[stats]{optim}. \cr
\code{cccp, cccp2} \tab The parameters are those from function \link[cccp]{ctrl}. They are among others the maximum count of iterations as an integer value (\code{maxiters}), the feasible level of convergence to be achieved (\code{feastol}) and whether the solver's progress during the iterations is shown (\code{trace}). If numerical problems are encountered increase the optimization parameter \code{feastol} or reduce parameter \code{stepadj}.\cr
 \code{csdp} \tab Available parameters are described in the CSDP User's Guide: \code{https://projects.coin-or.org/Csdp/export/49/trunk/doc/csdpuser.pdf} . \cr
 \code{slsqp} \tab Available parameters are described in \link[nloptr]{nl.opts}\cr
}
}
}
\details{
Computation of optimum genetic contributions for the breeding animals. If the function does not provide a valid result due to numerical problems then try the following modifications:

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
 which has component \code{parent}.  This is the data frame \code{phen} with the additional column \code{oc} containing the optimum genetic contribution of each breeding animal to the next generation, \code{lb} containing the lower bounds of the optimum contributions, and \code{ub} containing the upper bounds.


Note: In windows R3.2.4, if there are a lot / large objects in memory, eigen or chol, which are called from this function, may return NULL and subsequently no result is obtained. This is a bug in R.
}


\examples{
data(Kin)  #see ?kin for creating this object
data(Phen)
Phen[1:5,]

########################################
#  Check available objective functions #
#  and constraints                     #
########################################

help.opticont(Kin, Phen)

#########################
#   Check data          #
#########################

plot(Phen$MC, Phen$BV)
points(Phen$MC[Phen$Sex==1], Phen$BV[Phen$Sex==1], col="red", pch=18)

#################################################################
#       Compute the mean kinship fA that would be achieved      #
#          in the offspring without selection                   #
#################################################################
con     <- list(ub=c(M=-1, F=-1))
noSel   <- opticont(K=Kin, phen=Phen, con=con)
noSel.s <- summary(noSel)
noSel.s$fA
# 0.02284205
#===============================================================#
# => Allow the kinship in the next generation                   #
#    to be slightly larger, e.g. 0.03                           #
#===============================================================#


##################################################################
#                   Compute the genetic progress                 #
#            achievable while constraining only fA               #
##################################################################
con     <- list(ub.fA=0.03, ub=c(M=NA, F=-1))
maxBV   <- opticont("max.BV", K=Kin, phen=Phen, con=con, trace = FALSE)
maxBV.s <- summary(maxBV)
maxBV.s$meanBV - noSel.s$meanBV
# [1] 1.159723


#################################################################
#  Compute the minimum migrant contributions achievable while   #
#        constraining mean kinship fA in the offspring.         #
#################################################################
con     <- list(ub.fA=0.03, ub=c(M=NA, F=-1))
minMC   <- opticont("min.MC", K=Kin, phen=Phen, con=con, trace = FALSE)
minMC.s <- summary(minMC)
c(minMC.s$meanMC, noSel.s$meanMC)
#[1] 0.4958524 0.5764507
#===============================================================#
# => choose an upper bound for the migrant contribution MC in   #
#    the offspring between 0.492 and 0.574, e.g. 0.55           #
#===============================================================#


##################################################################
#                   Compute the genetic progress                 #
#            achievable while constraining fA and MC             #
##################################################################
con     <- list(ub.fA=0.03, ub.MC=0.55, ub=c(M=NA, F=-1))
maxBV.MC<- opticont("max.BV", K=Kin, phen=Phen, con=con, trace = FALSE)
maxBV.MC.s <- summary(maxBV.MC)
maxBV.MC.s$meanBV - noSel.s$meanBV
# [1] 0.5809614


##################################################################
# compute the minimum conditional kinship fD achievable while    #
# constraining fA and MC:                                        #
##################################################################
con     <- list(ub.fA=0.03, ub.MC=0.55, ub=c(M=NA, F=-1))
minfD   <- opticont("min.fD", K=Kin, phen=Phen, con=con, trace = FALSE)
minfD.s <- summary(minfD)
c(minfD.s$fD, noSel.s$fD)
#[1]  0.06660931 0.07698101

#================================================================#
# => Choose constraint for fD e.g.  0.073                        #
#================================================================#


##################################################################
#                   Compute the genetic progress                 #
#            achievable while constraining MC, fA and fD         #
##################################################################
con        <- list(ub.fA=0.03, ub.fD=0.073, ub.MC=0.55, ub=c(M=NA, F=-1))
maxBV.fD   <- opticont("max.BV", K=Kin, phen=Phen, con=con, trace = FALSE)
maxBV.fD.s <- summary(maxBV.fD)
maxBV.fD.s$meanBV - noSel.s$meanBV
#[1]  0.3593038
#================================================================#
# => Compared to no selection, the BV can be increased           #
#    by  0.3593038                                               #
#================================================================#



##################################################################
# Summary statistics from different optimizations                #
# can be combined in a data frame. The most important parameters #
# are printed for comparison:                                    #
##################################################################
Res<-rbind(noSel.s,  maxBV.fD.s, maxBV.MC.s, maxBV.s)
format(Res[,c("valid","meanBV", "meanMC", "fA", "fD")],digits=3)
#         valid meanBV meanMC     fA     fD
#noSel     TRUE 0.0985  0.576 0.0228 0.0770
#maxBV.fD  TRUE 0.4578  0.550 0.0241 0.0730
#maxBV.MC  TRUE 0.6794  0.550 0.0300 0.0908
#maxBV     TRUE 1.2582  0.760 0.0300 0.0908

#================================================================#
# => Genetic progress could be reduced considerably when migrant #
# contributions and diversity of native alleles are taken into   #
# account in optimum contribution selection.                     #
#================================================================#


##################################################################
# Compare the optimum contributions with and without controlling #
# diversity of native alleles by using different constraints:    #
##################################################################
plot(maxBV.fD$parent$oc, maxBV.MC$parent$oc, pch=18)
plot(maxBV$parent$oc,    maxBV.fD$parent$oc, pch=18)
x<-list(maxBV.fD, maxBV.MC, maxBV)
cor(sapply(x,function(x){x$parent$oc}))
#           [,1]      [,2]       [,3]
#[1,] 1.00000000 0.7895464 0.07438411
#[2,] 0.78954643 1.0000000 0.20823346
#[3,] 0.07438411 0.2082335 1.00000000

}

\references{
Borchers, B. (1999). CSDP, A C Library for Semidefinite Programming Optimization Methods and Software 11(1):613-623
\code{http://euler.nmt.edu/~brian/csdppaper.pdf}

Kraft, D. (1988). A software package for sequential quadratic programming, Technical Report DFVLR-FB 88-28, Institut fuer Dynamik der Flugsysteme, Oberpfaffenhofen, July 1988.

Lange K, Optimization, 2004, Springer.

Madsen K, Nielsen HB, Tingleff O, Optimization With Constraints, 2004, IMM, Technical University of Denmark.
}


\author{Robin Wellmann}
