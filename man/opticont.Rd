\name{opticont}
\Rdversion{1.1}
\alias{opticont}
\title{Optimum Contributions of Selection Candidates}

\description{
The optimum contributions of selection candidates to the offspring are calculated. The optimization procedure can take into account conflicting breeding goals, which are to achieve genetic gain, to reduce the rate of inbreeding, and to recover the original genetic background of a breed. 

Optimization can be done for several breeds simultaneously, which is adviseable if the aim is to increase diversity or genetic distance between breeds.
}
\usage{
opticont(method, cand, con, bc=NA,  solver="default", quiet=FALSE, 
         make.definite=solver=="csdp", ...)
}

\arguments{
\item{method}{Character string \code{"min.VAR"}, or \code{"max.VAR"}, whereby \code{VAR} is the name of the variable to be minimized or maximized. Available methods are reported by function \link{candes}.}

\item{cand}{An R-Object containing all information describing the individuals (phenotypes and kinships), which are the selection candidates, and possibly random samples of individuals from other breeds. It can be created with function \link{candes}. If column \code{Sex} of data frame \code{cand$phen} contains \code{NA} for one breed, then the constraint stating that contributions of both sexes must be equal is not applied.}

\item{con}{List defining threshold values for constraints. The components are described in the Details section. If one is missing, then the respective constraint is not applied. Permitted constraint names are reported by function \link{candes}.}

\item{bc}{Named numeric vector with breed contributions, which is only needed if \code{cand$phen} contains individuals from different breeds. It contains the proportion of each breed in a hypothetical multi-breed population for which the diversity across breeds should be managed. The names of the components are the breed names. }

\item{solver}{Name of the solver used for optimization. Available solvers are  \code{"alabama"}, \code{"cccp"}, \code{"cccp2"}, \code{"csdp"}, and \code{"slsqp"}. By default, the solver is chosen automatically. The solvers are the same as for function \link[optiSolve]{solvecop} from package \code{optiSolve}.}

\item{quiet}{If \code{quiet=FALSE} then detailed information is shown.}
\item{make.definite}{Logical variable indicating whether non-positive-semidefinite matrices should be approximated by positive-definite matrices. This is always done for solvers that are known not to convergue otherwise.}

\item{...}{Tuning parameters of the solver. The available parameters depend on the solver and will be printed when function \code{opticont} is used with default values. Definitions of the tuning parameters can be found for \code{alabama}  in \link[alabama]{auglag} and \link[stats]{optim}, for \code{cccp} and \code{cccp2} in \link[cccp]{ctrl}, for \code{csdp} in \url{https://projects.coin-or.org/Csdp/export/49/trunk/doc/csdpuser.pdf}, and for \code{slsqp}  in \link[nloptr]{nl.opts}.}

}
\details{

The optimum contributions of selection candidates to the offspring are calculated. The proportion of offspring that should have a particular selection candidate as parent is twice its optimum contribution.

\bold{Constraints}

Argument \code{con} is a list defining the constraints. Permitted names for the components are displayed by  function \link{candes}. Their meaning is as follows:

\bold{uniform}: Character vector specifying the breeds or sexes for which the contributions are not to be optimized. Within each of these groups it is assumed that all individuals have equal (uniform) contributions. Character string \code{"BREED.female"} means that all females from breed \code{BREED} have equal contributions and thus equal numbers of offspring.

\bold{lb}: Named numeric vector containing lower bounds for the contributions of the individuals. The component names are their IDs. By default the lower bound is 0 for all individuals.

\bold{ub}: Named numeric vector containing upper bounds for the contributions of the individuals. Their component names are the IDs. By default no upper bound is specified.

\bold{ub.VAR}: Upper bound for the expected mean value of kinship or trait \bold{VAR} in the offspring. Upper bounds for an arbitrary number of different kinships and traits may be provided. If data frame \code{cand$phen} contains individuals from several breeds, the bound refers to the mean value of the kinship or trait in the multi-breed offspring population.

\bold{ub.VAR.BREED}: Upper bound for the expected mean value of kinship or trait \bold{VAR} in the offspring from breed \bold{BREED}. Upper bounds for an arbitrary number of different kinships and traits may be provided. 

Note that \bold{VAR} must be replaced by the name of the variable and \bold{BREED} by the name of the breed. For traits, lower bounds can be defined as \bold{lb.VAR} or \bold{lb.VAR.BREED}. Equality constraints can be defined as  \bold{eq.VAR} or \bold{eq.VAR.BREED}.


\bold{Application to multi-breed data}

Optimization can be done for several breeds simultaneously, which is adviseable if the aim is to increase genetic diversity in a multi-breed population, or to increase the genetic distance between the breed of interest and other breeds. However, for computing the kinship of individuals from different breeds, marker data is needed. 

The multi-breed population referred above is a hypothetical subdivided population consisting of purebred animals from the breeds included in column \code{Breed} of \code{cand$phen}. The proportion of individuals from a given breed in this population is its breed contribution specified in argument \code{bc}. It is not the proportion of individuals of this breed in data frame \code{cand$phen}.

The aim is to minimize or to constrain the average genomic kinship in this multi-breed population. This causes the genetic distance between the breed of interest and other breeds to increase, and thus may increase the conservation value of the breed.


\bold{Remark}

If the function does not provide a valid result due to numerical problems then try to use another solver, use other optimization parameters, define upper or lower bounds instead of equality constraints, or relax the constraints to ensure that the optimization problem is solvable. 
}

\value{
A list with the following components

\item{parent}{Data frame \code{cand$phen} with some appended columns. Column \code{oc} contains the optimum contributions of the selection candidates, column \code{lb} the lower bounds, and \code{ub} the upper bounds for the contributions.}
\item{info}{Data frame with component \code{valid} indicating if all constraints are fulfilled, component \code{solver} containing the name of the solver used for optimization, and component \code{status} describing the solution as reported by the solver.}
\item{mean}{Data frame containing the expected mean value of each kinship and trait in the offspring population.}
\item{bc}{Data frame with breed contributions in the hypothetical multi-breed population used for computing the average kinship across breeds.}
\item{obj.fun}{Named numeric value with value and name of the objective function.}
\item{summary}{Data frame containing one row for each constraint with the value of the constraint in column \code{Val}, and the bound for the constraint in column \code{Bound}. Column \code{OK} states if the constraint is fulfilled, and column \code{Breed} contains the name of the breed to which the constraint applies.  The value of the objective function is shown in the first row.  Additional rows contain the values of traits and kinships in the offspring which are not constrained.}
}


\examples{
######################################################
# Example 1: Advanced OCS using pedigree data        #
######################################################

### Prepare pedigree data
data(PedigWithErrors)
data(Phen)

Pedig    <- prePed(PedigWithErrors, keep=Phen$Indiv, 
                   thisBreed="Hinterwaelder", lastNative=1970)
Pedig$NC <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")$native
Phen     <- merge(Pedig, Phen[,c("Indiv", "BV")], by="Indiv")
pKin     <- pedIBD(Pedig, keep.only=Phen$Indiv)
pKinatN  <- pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=Phen$Indiv)
cand     <- candes(phen=Phen, pKin=pKin, pKinatN=pKinatN)
head(cand$phen)

######################################################
# Objective: Maximize genetic gain in Angler cattle  #
# Constraints:                                       #
#   - mean kinship                                   #
#   - mean kinship at native alleles                 #
#   - genetic contributions from other breeds        #
######################################################

# Define constraints for the next generation
Ne             <- 100
con            <- list(uniform="female")
con$ub.pKin    <- cand$mean$pKin    + (1-cand$mean$pKin   )*(1/(2*Ne))
con$ub.pKinatN <- cand$mean$pKinatN + (1-cand$mean$pKinatN)*(1/(2*Ne))
con$lb.NC      <- 1.03*cand$mean$NC

# Compute the genetic progress achievable
Offspring  <- opticont("max.BV", cand, con, trace=FALSE)

# Check if the optimization problem is solved 
Offspring$info           
#  valid solver  status
#1  TRUE  cccp2 optimal

# Average values of traits and kinships 
Offspring$mean   
#         NC        BV       pKin    pKinatN
#1 0.4362558 0.5140699 0.02573526 0.07993624

Offspring$mean$BV

# Value of the objective function 
Offspring$obj.fun
#       BV 
#0.5140699 

# Data frame with optimum contributions
head(Offspring$parent)   

# See how much the values deviate from the constraints
Offspring$summary        

######################################################
# Further Examples: Advanced OCS using genotype data #
######################################################
\dontrun{
### Prepare genotype data

library("optiSel")
data(map) 
data(Cattle)
dir  <- system.file("extdata", package = "optiSel")
files<- file.path(dir, paste("Chr", 1:2, ".phased", sep=""))

### Compute genomic kinship and genomic kinship at native segments
G  <- segIBD(files, map, minL=2.0)
GN <- segIBDatN(files, Cattle, map, thisBreed="Angler", refBreeds="others", minL=2.0)

### Compute migrant contributions of selection candidates 
Haplo<- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="others",
                  minL=2.0, what="match")
Comp <- segBreedComp(Haplo$match, map)
Cattle$NC <- NA
Cattle[rownames(Comp), "NC"] <- Comp$native
apply(Comp[,-1],2,mean)
#    native          F          H          R 
#0.43798759 0.02997242 0.24603493 0.28600506 

######################################################
# Objective: Minimize inbreeding in Angler cattle    #
# Constraints:                                       #
#   - breeding values                                #
#   - mean kinship at native alleles                 #
#   - genetic contributions from other breeds        #
######################################################

cand <- candes(phen=Cattle[Cattle$Breed=="Angler",], sKin=G, sKinatN=GN)

# Define Constraints for the next generation
Ne <- 100
con            <- list(uniform="female")
con$ub.sKinatN <- cand$mean$sKinatN + (1-cand$mean$sKinatN)*(1/(2*Ne))
con$lb.NC      <- 1.03*cand$mean$NC
con$lb.BV      <- cand$mean$BV

# Compute optimum contributions; the objective is to minimize mean kinship 
Offspring   <- opticont("min.sKin", cand, con=con)

# Check if the optimization problem is solved 
Offspring$info           

# Average values of traits and kinships 
Offspring$mean           

# Value of the objective function 
Offspring$obj.fun
#      sKin 
#0.03977612 


######################################################
# Objective: Maximize breeding value in Angler       #
# Constraints:                                       #
#   - kinship at native alleles                      #
#   - kinship across breeds                          #
######################################################

cand <- candes(phen=Cattle, sKin=G, sKinatN.Angler=GN)

Unif      <- c("Angler.female", "Fleckvieh", "Holstein", "Rotbunt")
con       <- list(uniform=Unif, ub.sKin=0.027, ub.sKinatN.Angler=0.06)
Offspring <- opticont("max.BV", cand, con, trace=FALSE)

Offspring$mean

Offspring$obj.fun
#       BV 
#0.4986473 

######################################################
# Objective: Minimize kinship across breeds          #
# Constraints:                                       #
#   - native contributions                           #
#   - breeding values                                #
# by optimizing contributions of Angler males        #
######################################################

Unif      <- c("Angler.female", "Fleckvieh", "Holstein", "Rotbunt")
con       <- list(uniform=Unif, ub.sKinatN.Angler=0.06, lb.NC=0.7, lb.BV=0.5)
Offspring <- opticont("min.sKin", cand, con, trace=FALSE)

Offspring$mean

Offspring$obj.fun
#      sKin 
#0.02953495 

######################################################
# Objective:  Maximize breeding values in all breeds #
# Constraints:                                       #
#   - average kinships within each breed             #
#   - average kinships across breeds                 #
#   - average native kinship of Angler               #
#   - average native contribution in Angler          #
# by optimizing contributions of males from all breeds#
######################################################

set.seed(1)
Phen <- Cattle
Phen[Phen$Breed!="Angler","BV"] <- rnorm(sum(Phen$Breed!="Angler")) #simulate BV
cand <- candes(phen=Phen, sKin=G, sKinatN.Angler=GN)

Unif <- paste0(c("Angler", "Fleckvieh", "Holstein", "Rotbunt"), ".female")
con  <- list(uniform          = Unif, 
             ub.sKin          = cand$mean$sKin - 0.002,
             ub.sKin.Angler   = cand$mean$sKin.Angler   + 0.005,
             ub.sKin.Holstein = cand$mean$sKin.Holstein + 0.005,
             ub.sKin.Rotbunt  = cand$mean$sKin.Rotbunt  + 0.005,
             ub.sKin.Fleckvieh= cand$mean$sKin.Fleckvieh+ 0.005,
             ub.sKinatN.Angler= cand$mean$sKinatN.Angler+ 0.005, 
             lb.NC            = cand$mean$NC + 0.05)
            
Offspring <- opticont("max.BV", cand, con, trace=FALSE)

Offspring$mean

Offspring$obj.fun
#       BV 
#0.6440849  
}

}

\references{
Wellmann, R. (2017). Optimum Contribution Selection and Mate
Allocation for Breeding: The R Package optiSel. submitted
}


\author{Robin Wellmann}
