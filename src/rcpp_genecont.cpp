
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_genecont(const arma::ivec& numSire, const arma::ivec& numDam, const arma::ivec&  numAnc, const Rcpp::CharacterVector rNames, const Rcpp::CharacterVector cNames){
  int i, j, nSire, nDam, nP;
  double pCont;
  int N    = numSire.n_elem;
  int NAnc = numAnc.n_elem;
  
  Rcpp::NumericMatrix rGeneCont(N, NAnc);
  arma::mat GeneCont(rGeneCont.begin(), rGeneCont.nrow(), rGeneCont.ncol(), false);
  GeneCont.zeros();
  
  for(i=0; i<NAnc;i++){
    GeneCont.at(numAnc.at(i), i) = 1.0;
  }
  
 
  for(i=0; i<N;i++){
    nSire = numSire.at(i);
    nDam  = numDam.at(i);
    /*Rprintf("i=%d\n",i);*/
    
    if(nSire | nDam){
      nP = ((nSire<nDam)?nDam:nSire);
      for(j=0; (numAnc.at(j)<=nP)&(j<NAnc); j++){
        pCont = ((nSire)?(GeneCont.at(nSire-1, j)):0.0) + ((nDam)?(GeneCont.at(nDam-1,  j)):0.0);
        if(pCont>0){GeneCont.at(i, j) = 0.5*pCont;}
      }
    }
  }
  
  rGeneCont.attr("dimnames") = Rcpp::List::create(rNames, cNames);
  
  return rGeneCont;
}
