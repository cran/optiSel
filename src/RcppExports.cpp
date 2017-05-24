// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_completeness
Rcpp::DataFrame rcpp_completeness(Rcpp::StringVector Indiv, const arma::ivec& ArmanumSire, const arma::ivec& ArmanumDam, int maxd);
RcppExport SEXP optiSel_rcpp_completeness(SEXP IndivSEXP, SEXP ArmanumSireSEXP, SEXP ArmanumDamSEXP, SEXP maxdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type Indiv(IndivSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmanumSire(ArmanumSireSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmanumDam(ArmanumDamSEXP);
    Rcpp::traits::input_parameter< int >::type maxd(maxdSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_completeness(Indiv, ArmanumSire, ArmanumDam, maxd));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_genecont
Rcpp::NumericMatrix rcpp_genecont(const arma::ivec& numSire, const arma::ivec& numDam, const arma::ivec& numAnc, const arma::ivec& numKeep, const arma::ivec& ainKeep, const Rcpp::CharacterVector rNames, const Rcpp::CharacterVector cNames, const arma::ivec& anOff);
RcppExport SEXP optiSel_rcpp_genecont(SEXP numSireSEXP, SEXP numDamSEXP, SEXP numAncSEXP, SEXP numKeepSEXP, SEXP ainKeepSEXP, SEXP rNamesSEXP, SEXP cNamesSEXP, SEXP anOffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type numSire(numSireSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numDam(numDamSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numAnc(numAncSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numKeep(numKeepSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ainKeep(ainKeepSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type rNames(rNamesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type cNames(cNamesSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type anOff(anOffSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_genecont(numSire, numDam, numAnc, numKeep, ainKeep, rNames, cNames, anOff));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_haplofreq
Rcpp::List rcpp_haplofreq(std::string pathThisBreed, std::string pathRefBreeds, std::string pathFreq, std::string pathOrig, std::vector< std::string > MarkerName, std::string stdBreedSymbol, const arma::ivec& ArmaIndexC, const arma::imat& ArmaIndexR, int NFileC, int NFileR, int NC, const arma::ivec& ArmaNR, int minSNP, double minL, double ubFreq, const arma::vec& ArmaPos, std::string stdsymB, int skip, int cskip, int getFreq, int getOrig);
RcppExport SEXP optiSel_rcpp_haplofreq(SEXP pathThisBreedSEXP, SEXP pathRefBreedsSEXP, SEXP pathFreqSEXP, SEXP pathOrigSEXP, SEXP MarkerNameSEXP, SEXP stdBreedSymbolSEXP, SEXP ArmaIndexCSEXP, SEXP ArmaIndexRSEXP, SEXP NFileCSEXP, SEXP NFileRSEXP, SEXP NCSEXP, SEXP ArmaNRSEXP, SEXP minSNPSEXP, SEXP minLSEXP, SEXP ubFreqSEXP, SEXP ArmaPosSEXP, SEXP stdsymBSEXP, SEXP skipSEXP, SEXP cskipSEXP, SEXP getFreqSEXP, SEXP getOrigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pathThisBreed(pathThisBreedSEXP);
    Rcpp::traits::input_parameter< std::string >::type pathRefBreeds(pathRefBreedsSEXP);
    Rcpp::traits::input_parameter< std::string >::type pathFreq(pathFreqSEXP);
    Rcpp::traits::input_parameter< std::string >::type pathOrig(pathOrigSEXP);
    Rcpp::traits::input_parameter< std::vector< std::string > >::type MarkerName(MarkerNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type stdBreedSymbol(stdBreedSymbolSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexC(ArmaIndexCSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type ArmaIndexR(ArmaIndexRSEXP);
    Rcpp::traits::input_parameter< int >::type NFileC(NFileCSEXP);
    Rcpp::traits::input_parameter< int >::type NFileR(NFileRSEXP);
    Rcpp::traits::input_parameter< int >::type NC(NCSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaNR(ArmaNRSEXP);
    Rcpp::traits::input_parameter< int >::type minSNP(minSNPSEXP);
    Rcpp::traits::input_parameter< double >::type minL(minLSEXP);
    Rcpp::traits::input_parameter< double >::type ubFreq(ubFreqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmaPos(ArmaPosSEXP);
    Rcpp::traits::input_parameter< std::string >::type stdsymB(stdsymBSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type cskip(cskipSEXP);
    Rcpp::traits::input_parameter< int >::type getFreq(getFreqSEXP);
    Rcpp::traits::input_parameter< int >::type getOrig(getOrigSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_haplofreq(pathThisBreed, pathRefBreeds, pathFreq, pathOrig, MarkerName, stdBreedSymbol, ArmaIndexC, ArmaIndexR, NFileC, NFileR, NC, ArmaNR, minSNP, minL, ubFreq, ArmaPos, stdsymB, skip, cskip, getFreq, getOrig));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_makeA
Rcpp::NumericMatrix rcpp_makeA(const arma::ivec& numSire, const arma::ivec& numDam, const arma::mat& AFounder, const arma::ivec& numFounder, const Rcpp::CharacterVector IndivName);
RcppExport SEXP optiSel_rcpp_makeA(SEXP numSireSEXP, SEXP numDamSEXP, SEXP AFounderSEXP, SEXP numFounderSEXP, SEXP IndivNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type numSire(numSireSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numDam(numDamSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type AFounder(AFounderSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numFounder(numFounderSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type IndivName(IndivNameSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_makeA(numSire, numDam, AFounder, numFounder, IndivName));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_makeA_lowMem
Rcpp::NumericMatrix rcpp_makeA_lowMem(const arma::ivec& numSire, const arma::ivec& numDam, const arma::mat& AFounder, const arma::ivec& numFounder, const Rcpp::CharacterVector IndivName, const arma::ivec& numKeep, const arma::ivec& ainKeep, const arma::ivec& anOff);
RcppExport SEXP optiSel_rcpp_makeA_lowMem(SEXP numSireSEXP, SEXP numDamSEXP, SEXP AFounderSEXP, SEXP numFounderSEXP, SEXP IndivNameSEXP, SEXP numKeepSEXP, SEXP ainKeepSEXP, SEXP anOffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type numSire(numSireSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numDam(numDamSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type AFounder(AFounderSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numFounder(numFounderSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type IndivName(IndivNameSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type numKeep(numKeepSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ainKeep(ainKeepSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type anOff(anOffSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_makeA_lowMem(numSire, numDam, AFounder, numFounder, IndivName, numKeep, ainKeep, anOff));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_nativecont
Rcpp::NumericVector rcpp_nativecont(std::string pathNative, int NFileN, int NC, const arma::ivec& ArmaIndexN, int M, const arma::vec& ArmaNkb);
RcppExport SEXP optiSel_rcpp_nativecont(SEXP pathNativeSEXP, SEXP NFileNSEXP, SEXP NCSEXP, SEXP ArmaIndexNSEXP, SEXP MSEXP, SEXP ArmaNkbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pathNative(pathNativeSEXP);
    Rcpp::traits::input_parameter< int >::type NFileN(NFileNSEXP);
    Rcpp::traits::input_parameter< int >::type NC(NCSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexN(ArmaIndexNSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmaNkb(ArmaNkbSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_nativecont(pathNative, NFileN, NC, ArmaIndexN, M, ArmaNkb));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_segBreedComp
Rcpp::NumericMatrix rcpp_segBreedComp(std::vector<std::string> pathNative, int Nfile, int N, const arma::ivec& ArmaIndexN, const arma::ivec& MatChr, const arma::vec& Armakb);
RcppExport SEXP optiSel_rcpp_segBreedComp(SEXP pathNativeSEXP, SEXP NfileSEXP, SEXP NSEXP, SEXP ArmaIndexNSEXP, SEXP MatChrSEXP, SEXP ArmakbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type pathNative(pathNativeSEXP);
    Rcpp::traits::input_parameter< int >::type Nfile(NfileSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexN(ArmaIndexNSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type MatChr(MatChrSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Armakb(ArmakbSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_segBreedComp(pathNative, Nfile, N, ArmaIndexN, MatChr, Armakb));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_segIBD
Rcpp::NumericMatrix rcpp_segIBD(std::string path1, std::string path2, int NFile1, int NFile2, const arma::ivec& ArmaIndex1, const arma::ivec& ArmaIndex2, int N1, int N2, int minSNP, double minL, const arma::vec& ArmacM, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip);
RcppExport SEXP optiSel_rcpp_segIBD(SEXP path1SEXP, SEXP path2SEXP, SEXP NFile1SEXP, SEXP NFile2SEXP, SEXP ArmaIndex1SEXP, SEXP ArmaIndex2SEXP, SEXP N1SEXP, SEXP N2SEXP, SEXP minSNPSEXP, SEXP minLSEXP, SEXP ArmacMSEXP, SEXP ArmakbSEXP, SEXP aSEXP, SEXP stdsymBSEXP, SEXP skipSEXP, SEXP cskipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path1(path1SEXP);
    Rcpp::traits::input_parameter< std::string >::type path2(path2SEXP);
    Rcpp::traits::input_parameter< int >::type NFile1(NFile1SEXP);
    Rcpp::traits::input_parameter< int >::type NFile2(NFile2SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndex1(ArmaIndex1SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndex2(ArmaIndex2SEXP);
    Rcpp::traits::input_parameter< int >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< int >::type N2(N2SEXP);
    Rcpp::traits::input_parameter< int >::type minSNP(minSNPSEXP);
    Rcpp::traits::input_parameter< double >::type minL(minLSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmacM(ArmacMSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Armakb(ArmakbSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type stdsymB(stdsymBSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type cskip(cskipSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_segIBD(path1, path2, NFile1, NFile2, ArmaIndex1, ArmaIndex2, N1, N2, minSNP, minL, ArmacM, Armakb, a, stdsymB, skip, cskip));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_segIBDandN
Rcpp::NumericMatrix rcpp_segIBDandN(std::string pathThisBreed, std::string pathNative, int NFileC, int NFileN, const arma::ivec& ArmaIndexC, const arma::ivec& ArmaIndexN, int NC, int minSNP, double minL, const arma::vec& ArmaPos, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip);
RcppExport SEXP optiSel_rcpp_segIBDandN(SEXP pathThisBreedSEXP, SEXP pathNativeSEXP, SEXP NFileCSEXP, SEXP NFileNSEXP, SEXP ArmaIndexCSEXP, SEXP ArmaIndexNSEXP, SEXP NCSEXP, SEXP minSNPSEXP, SEXP minLSEXP, SEXP ArmaPosSEXP, SEXP ArmakbSEXP, SEXP aSEXP, SEXP stdsymBSEXP, SEXP skipSEXP, SEXP cskipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pathThisBreed(pathThisBreedSEXP);
    Rcpp::traits::input_parameter< std::string >::type pathNative(pathNativeSEXP);
    Rcpp::traits::input_parameter< int >::type NFileC(NFileCSEXP);
    Rcpp::traits::input_parameter< int >::type NFileN(NFileNSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexC(ArmaIndexCSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexN(ArmaIndexNSEXP);
    Rcpp::traits::input_parameter< int >::type NC(NCSEXP);
    Rcpp::traits::input_parameter< int >::type minSNP(minSNPSEXP);
    Rcpp::traits::input_parameter< double >::type minL(minLSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmaPos(ArmaPosSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Armakb(ArmakbSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type stdsymB(stdsymBSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type cskip(cskipSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_segIBDandN(pathThisBreed, pathNative, NFileC, NFileN, ArmaIndexC, ArmaIndexN, NC, minSNP, minL, ArmaPos, Armakb, a, stdsymB, skip, cskip));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_segIBDandNVersion2
Rcpp::NumericMatrix rcpp_segIBDandNVersion2(std::string pathThisBreed, int NFileC, int NC, const arma::ivec& ArmaIndexC, const arma::mat& ArmaNat, int minSNP, double minL, const arma::vec& ArmaPos, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip);
RcppExport SEXP optiSel_rcpp_segIBDandNVersion2(SEXP pathThisBreedSEXP, SEXP NFileCSEXP, SEXP NCSEXP, SEXP ArmaIndexCSEXP, SEXP ArmaNatSEXP, SEXP minSNPSEXP, SEXP minLSEXP, SEXP ArmaPosSEXP, SEXP ArmakbSEXP, SEXP aSEXP, SEXP stdsymBSEXP, SEXP skipSEXP, SEXP cskipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pathThisBreed(pathThisBreedSEXP);
    Rcpp::traits::input_parameter< int >::type NFileC(NFileCSEXP);
    Rcpp::traits::input_parameter< int >::type NC(NCSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexC(ArmaIndexCSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ArmaNat(ArmaNatSEXP);
    Rcpp::traits::input_parameter< int >::type minSNP(minSNPSEXP);
    Rcpp::traits::input_parameter< double >::type minL(minLSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmaPos(ArmaPosSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Armakb(ArmakbSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type stdsymB(stdsymBSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type cskip(cskipSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_segIBDandNVersion2(pathThisBreed, NFileC, NC, ArmaIndexC, ArmaNat, minSNP, minL, ArmaPos, Armakb, a, stdsymB, skip, cskip));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_segInbreeding
Rcpp::NumericVector rcpp_segInbreeding(std::string path1, std::string path2, int NFile1, int NFile2, const arma::ivec& ArmaIndex1, const arma::ivec& ArmaIndex2, int N1, int N2, int M, int minSNP, double minL, const arma::vec& ArmacM, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip);
RcppExport SEXP optiSel_rcpp_segInbreeding(SEXP path1SEXP, SEXP path2SEXP, SEXP NFile1SEXP, SEXP NFile2SEXP, SEXP ArmaIndex1SEXP, SEXP ArmaIndex2SEXP, SEXP N1SEXP, SEXP N2SEXP, SEXP MSEXP, SEXP minSNPSEXP, SEXP minLSEXP, SEXP ArmacMSEXP, SEXP ArmakbSEXP, SEXP aSEXP, SEXP stdsymBSEXP, SEXP skipSEXP, SEXP cskipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path1(path1SEXP);
    Rcpp::traits::input_parameter< std::string >::type path2(path2SEXP);
    Rcpp::traits::input_parameter< int >::type NFile1(NFile1SEXP);
    Rcpp::traits::input_parameter< int >::type NFile2(NFile2SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndex1(ArmaIndex1SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndex2(ArmaIndex2SEXP);
    Rcpp::traits::input_parameter< int >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< int >::type N2(N2SEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type minSNP(minSNPSEXP);
    Rcpp::traits::input_parameter< double >::type minL(minLSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmacM(ArmacMSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Armakb(ArmakbSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type stdsymB(stdsymBSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type cskip(cskipSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_segInbreeding(path1, path2, NFile1, NFile2, ArmaIndex1, ArmaIndex2, N1, N2, M, minSNP, minL, ArmacM, Armakb, a, stdsymB, skip, cskip));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_segN
Rcpp::NumericMatrix rcpp_segN(std::string pathNative, int NFileN, int NC, const arma::ivec& ArmaIndexN, const arma::vec& ArmaNkb);
RcppExport SEXP optiSel_rcpp_segN(SEXP pathNativeSEXP, SEXP NFileNSEXP, SEXP NCSEXP, SEXP ArmaIndexNSEXP, SEXP ArmaNkbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pathNative(pathNativeSEXP);
    Rcpp::traits::input_parameter< int >::type NFileN(NFileNSEXP);
    Rcpp::traits::input_parameter< int >::type NC(NCSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ArmaIndexN(ArmaIndexNSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ArmaNkb(ArmaNkbSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_segN(pathNative, NFileN, NC, ArmaIndexN, ArmaNkb));
    return rcpp_result_gen;
END_RCPP
}
