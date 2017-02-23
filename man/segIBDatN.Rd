\name{segIBDatN}
\alias{segIBDatN}

\title{Calculates Segment Based Kinship at Native Alleles.}

\description{
\bold{Seg}ment based probability of alleles to be \bold{IBD at N}ative haplotype segments ("kinship at native segments").}
\usage{
segIBDatN(files, phen, map, thisBreed, refBreeds="others", ubFreq=0.01, minSNP=20,
  unitP="kb", minL=1000, unitL="kb", a=0.0, keep=NULL, lowMem=TRUE, skip=NA, cskip=NA)
}

\arguments{

 \item{files}{This can be a character vector with names of the phased marker files, one file for each chromosome. 
 Alternatively \code{files} can be a list with the following two components:

a) \code{hap.thisBreed}: A character vector with names of the phased marker files for the individuals from \code{thisBreed}, one file for each chromosome.
 
b) \code{hap.refBreeds}: A character vector with names of the phased marker files for the individuals from the reference breeds (\code{refBreeds}), one file for each chromosome. If this component is missing, then it is assumed that the haplotypes of these animals are also included in \code{hap.thisBreed}.

File names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}. The required format of the marker files is described under \code{Details}.
}

\item{phen}{Data frame containing the ID (column \code{"Indiv"}) and the breed name (column \code{"Breed"}) of each individual.}
\item{map}{Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in kilo base pairs \code{'kb'}, and the position in centimorgan \code{'cM'}. (The position in base pairs could result in an integer overflow). The order of the markers must bethe same as in the files.}
\item{thisBreed}{Breed name: Results will be computed for individuals from \code{thisBreed}.}
\item{refBreeds}{Vector containing names of genotyped breeds. A segment is considered native if its frequency is smaller than \code{ubFreq}  in all \code{refBreeds}. The default \code{"others"} means that all genotyped breeds except \code{thisBreed} are considered.}
\item{ubFreq}{A segment is considered native if its frequency is smaller than \code{ubFreq} in all reference breeds.}
\item{minSNP}{Minimum number of marker SNPs included in a segment.}
\item{unitP}{The unit for measuring the proportion of the genome included in native segments.
  Possible units are the number of marker SNPs included in shared segments (\code{'SNP'}), the number of kilo base pairs (\code{'kb'}), and the total length of the shared segments in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
\item{minL}{Minimum length of a segment in \code{unitL} (e.g. in cM).}
\item{unitL}{The unit for measuring the length of a segment. Possible units are the number of marker SNPs included in the segment (\code{'SNP'}), the number of kilo base pairs (\code{'kb'}), and the genetic  distances between the first and the last marker in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
\item{a}{The function providing the weighting factor for each segment is w(x)=x*x/(a+x*x). The parameter of the function is the length of the segment in \code{unitL}. The default value \code{a=0.0} implies no weighting, whereas \code{a>0.0} implies that old inbreeding has less influence on the result than new inbreeding.}
 \item{keep}{
 Subset of the IDs of the individuals from data frame \code{phen} (including individuals from other breeds) or a logical vector indicating the animals in data frame \code{phen} that should be used. By default all individuals included in \code{phen} will be used.}
\item{lowMem}{If \code{lowMem=TRUE} then temporary files will be created and deleted.}
     \item{skip}{Take line \code{skip+1} of the genotype files as the row with column names. By default, the number is determined automatically.}
 \item{cskip}{Take column \code{cskip+1} of the genotype files as the first column with genotypes. By default, the number is determined automatically.}
}
\value{
A list containing matrices needed for computing the segment based pobability of alleles to be IBD at native segments. The list has components
\item{segN}{This matrix contains for each pair of individuals the probability   that two SNPs taken at random position from randomly chosen haplotypes both belong to native segments.}
\item{segIBDandN}{This matrix contains for each pair of individuals the probability   that two SNPs taken at random position from randomly chosen haplotypes belong to a shared native segment.}
\item{segZ}{1+segIBDandN-segN.}
The list has attribute \code{meanIBDatN} providing the probability of randomly chosen alleles to be IBD at native haplotype segments. Note that \code{1-meanIBDatN} is the genetic diversity at native segments within the genotyped individuals from \code{thisBreed}.
}
\details{
Computation of the segment based probability of alleles to be IBD at native haplotype segments. 

\bold{Genotype file format}: Each file containing phased genotypes has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. The first \code{cskip} columns are ignored. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces. The name of each file must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}.
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package = "optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")
Res   <- segIBDatN(files, Cattle, map, thisBreed="Angler", ubFreq=0.01, minL=1000, lowMem=FALSE)
               
## Mean kinship at native segments:
attributes(Res)$meanIBDatN
#[1] 0.06742614

## Results for individuals:
kin <-Res$segIBDandN/Res$segN
use <-upper.tri(kin) & Res$segN>0.2
boxplot(kin[use], ylim=c(0,1))

## Using temporary files to reduce working memory:
\dontrun{
Res <- segIBDatN(files, Cattle, map, thisBreed="Angler", ubFreq=0.01, minL=1000)
               
## Mean kinship at native segments:
attributes(Res)$meanIBDatN
#[1] 0.06742614
}
}
\author{Robin Wellmann}

