
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector rcpp_nativecont(std::string pathNative, int NFileN, int NC, const arma::ivec& ArmaIndexN, int M, const arma::vec& ArmaNkb) {
  int m, i;
  char str2[100];
  FILE *fN;
  Rcpp::NumericVector ArmaNatCont(NC);

  size_t bufsize = 2*NFileN;  
  char* Line = (char*)malloc(bufsize*sizeof(char));
  if(Line == NULL){error_return("Memory allocation failed.");};
  
  
  int* indexN     = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  double* Nkb     = (double*)calloc(ArmaNkb.n_elem, sizeof(double)); /*    M+1 - vector */
  double* NatCont = (double*)calloc(NC, sizeof(double));             /*     NC - vector */
  if(indexN  == NULL){error_return("Memory allocation failed.");};
  if(Nkb     == NULL){error_return("Memory allocation failed.");};
  if(NatCont == NULL){error_return("Memory allocation failed.");};
  
  for(m=0;m< M;m++){Nkb[m]    = ArmaNkb.at(m);}
  for(i=0;i<NC;i++){indexN[i] = ArmaIndexN.at(i);}
  
  /* ******* Main part ******** */
  fN = fopen(pathNative.c_str(),"r");
  if(fN == NULL){error_return("File opening failed.");};	 
  while(fgetc(fN)!='\n'){}
  
  m=0;
  while(fscanf(fN, "%s ", str2)>0){
    fgets(Line, 2*NFileN, fN);
    for(i=0; i<NC;i++){
      if(Line[2*indexN[i]]=='1'){
        NatCont[i] += Nkb[m];
      }
    }
    m=m+1;
  }
  fclose(fN);
  
  Rprintf("M=%d\n",m);
  
  for(i=0; i<NC;i++){
      ArmaNatCont.at(i) = NatCont[i];
  }

  free(NatCont);
  free(Nkb);
  free(indexN);
  free(Line);
  
  return ArmaNatCont;
}
