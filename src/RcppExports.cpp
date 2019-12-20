// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dist_transform
arma::vec dist_transform(arma::mat myfield, int matl);
RcppExport SEXP _rdimexp_dist_transform(SEXP myfieldSEXP, SEXP matlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type myfield(myfieldSEXP);
    Rcpp::traits::input_parameter< int >::type matl(matlSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_transform(myfield, matl));
    return rcpp_result_gen;
END_RCPP
}
// coarsen
arma::mat coarsen(arma::mat myfield, int nk);
RcppExport SEXP _rdimexp_coarsen(SEXP myfieldSEXP, SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type myfield(myfieldSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen(myfield, nk));
    return rcpp_result_gen;
END_RCPP
}
// make_datamat
arma::mat make_datamat(arma::mat myfield, int snipsize, int nsmooth);
RcppExport SEXP _rdimexp_make_datamat(SEXP myfieldSEXP, SEXP snipsizeSEXP, SEXP nsmoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type myfield(myfieldSEXP);
    Rcpp::traits::input_parameter< int >::type snipsize(snipsizeSEXP);
    Rcpp::traits::input_parameter< int >::type nsmooth(nsmoothSEXP);
    rcpp_result_gen = Rcpp::wrap(make_datamat(myfield, snipsize, nsmooth));
    return rcpp_result_gen;
END_RCPP
}
// emp_vario
arma::mat emp_vario(arma::mat cov);
RcppExport SEXP _rdimexp_emp_vario(SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(emp_vario(cov));
    return rcpp_result_gen;
END_RCPP
}
// dist_euclid
arma::mat dist_euclid(arma::mat locs);
RcppExport SEXP _rdimexp_dist_euclid(SEXP locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type locs(locsSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_euclid(locs));
    return rcpp_result_gen;
END_RCPP
}
// fitkrig
arma::vec fitkrig(arma::vec data, arma::mat locs_new, arma::mat locs_obs, arma::mat sig_obs, arma::vec params);
RcppExport SEXP _rdimexp_fitkrig(SEXP dataSEXP, SEXP locs_newSEXP, SEXP locs_obsSEXP, SEXP sig_obsSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type locs_new(locs_newSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type locs_obs(locs_obsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_obs(sig_obsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(fitkrig(data, locs_new, locs_obs, sig_obs, params));
    return rcpp_result_gen;
END_RCPP
}
// lsq_vario_fit
double lsq_vario_fit(arma::vec params, arma::mat empvario, arma::mat ds);
RcppExport SEXP _rdimexp_lsq_vario_fit(SEXP paramsSEXP, SEXP empvarioSEXP, SEXP dsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type empvario(empvarioSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ds(dsSEXP);
    rcpp_result_gen = Rcpp::wrap(lsq_vario_fit(params, empvario, ds));
    return rcpp_result_gen;
END_RCPP
}
// prof_nll
double prof_nll(double phi, arma::mat datamat, arma::mat ds, double prange);
RcppExport SEXP _rdimexp_prof_nll(SEXP phiSEXP, SEXP datamatSEXP, SEXP dsSEXP, SEXP prangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type datamat(datamatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ds(dsSEXP);
    Rcpp::traits::input_parameter< double >::type prange(prangeSEXP);
    rcpp_result_gen = Rcpp::wrap(prof_nll(phi, datamat, ds, prange));
    return rcpp_result_gen;
END_RCPP
}
// optim_nll_rcpp
arma::vec optim_nll_rcpp(arma::vec& interval, arma::mat& data, arma::mat& ds, double& prange);
RcppExport SEXP _rdimexp_optim_nll_rcpp(SEXP intervalSEXP, SEXP dataSEXP, SEXP dsSEXP, SEXP prangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ds(dsSEXP);
    Rcpp::traits::input_parameter< double& >::type prange(prangeSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_nll_rcpp(interval, data, ds, prange));
    return rcpp_result_gen;
END_RCPP
}
// optim_z_rcpp
arma::vec optim_z_rcpp(arma::vec& init_val, arma::mat& locs_obs, arma::mat& emp_vario, double& lambda);
RcppExport SEXP _rdimexp_optim_z_rcpp(SEXP init_valSEXP, SEXP locs_obsSEXP, SEXP emp_varioSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type locs_obs(locs_obsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type emp_vario(emp_varioSEXP);
    Rcpp::traits::input_parameter< double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_z_rcpp(init_val, locs_obs, emp_vario, lambda));
    return rcpp_result_gen;
END_RCPP
}
// freqweight
arma::vec freqweight(arma::vec x, arma::vec w);
RcppExport SEXP _rdimexp_freqweight(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(freqweight(x, w));
    return rcpp_result_gen;
END_RCPP
}
// cov_specd
arma::mat cov_specd(arma::vec rvec, arma::vec rprof, arma::mat r);
RcppExport SEXP _rdimexp_cov_specd(SEXP rvecSEXP, SEXP rprofSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rvec(rvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rprof(rprofSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_specd(rvec, rprof, r));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rdimexp_dist_transform", (DL_FUNC) &_rdimexp_dist_transform, 2},
    {"_rdimexp_coarsen", (DL_FUNC) &_rdimexp_coarsen, 2},
    {"_rdimexp_make_datamat", (DL_FUNC) &_rdimexp_make_datamat, 3},
    {"_rdimexp_emp_vario", (DL_FUNC) &_rdimexp_emp_vario, 1},
    {"_rdimexp_dist_euclid", (DL_FUNC) &_rdimexp_dist_euclid, 1},
    {"_rdimexp_fitkrig", (DL_FUNC) &_rdimexp_fitkrig, 5},
    {"_rdimexp_lsq_vario_fit", (DL_FUNC) &_rdimexp_lsq_vario_fit, 3},
    {"_rdimexp_prof_nll", (DL_FUNC) &_rdimexp_prof_nll, 4},
    {"_rdimexp_optim_nll_rcpp", (DL_FUNC) &_rdimexp_optim_nll_rcpp, 4},
    {"_rdimexp_optim_z_rcpp", (DL_FUNC) &_rdimexp_optim_z_rcpp, 4},
    {"_rdimexp_freqweight", (DL_FUNC) &_rdimexp_freqweight, 2},
    {"_rdimexp_cov_specd", (DL_FUNC) &_rdimexp_cov_specd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_rdimexp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
