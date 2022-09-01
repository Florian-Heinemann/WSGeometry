// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// StabilizedScaling_Rcpp
Rcpp::List StabilizedScaling_Rcpp(const arma::mat& cost_matrix, const arma::vec& supply, const arma::vec& demand, int supply_div_type, int demand_div_type, const arma::vec supply_div_parameters, const arma::vec demand_div_parameters, int iter_max, const arma::vec& epsvec, double tol, int check_interval, double indicator_slack, double clamp);
RcppExport SEXP _WSGeometry_StabilizedScaling_Rcpp(SEXP cost_matrixSEXP, SEXP supplySEXP, SEXP demandSEXP, SEXP supply_div_typeSEXP, SEXP demand_div_typeSEXP, SEXP supply_div_parametersSEXP, SEXP demand_div_parametersSEXP, SEXP iter_maxSEXP, SEXP epsvecSEXP, SEXP tolSEXP, SEXP check_intervalSEXP, SEXP indicator_slackSEXP, SEXP clampSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type cost_matrix(cost_matrixSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< int >::type supply_div_type(supply_div_typeSEXP);
    Rcpp::traits::input_parameter< int >::type demand_div_type(demand_div_typeSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type supply_div_parameters(supply_div_parametersSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type demand_div_parameters(demand_div_parametersSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsvec(epsvecSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type check_interval(check_intervalSEXP);
    Rcpp::traits::input_parameter< double >::type indicator_slack(indicator_slackSEXP);
    Rcpp::traits::input_parameter< double >::type clamp(clampSEXP);
    rcpp_result_gen = Rcpp::wrap(StabilizedScaling_Rcpp(cost_matrix, supply, demand, supply_div_type, demand_div_type, supply_div_parameters, demand_div_parameters, iter_max, epsvec, tol, check_interval, indicator_slack, clamp));
    return rcpp_result_gen;
END_RCPP
}
// Sinkhorn_Rcpp
Rcpp::List Sinkhorn_Rcpp(const arma::mat& cost_matrix, const arma::vec& supply, const arma::vec& demand, int supply_div_type, int demand_div_type, const arma::vec supply_div_parameters, const arma::vec demand_div_parameters, int iter_max, const arma::vec& epsvec, double tol, int thread_cnt, int max_lines_per_work, double indicator_slack);
RcppExport SEXP _WSGeometry_Sinkhorn_Rcpp(SEXP cost_matrixSEXP, SEXP supplySEXP, SEXP demandSEXP, SEXP supply_div_typeSEXP, SEXP demand_div_typeSEXP, SEXP supply_div_parametersSEXP, SEXP demand_div_parametersSEXP, SEXP iter_maxSEXP, SEXP epsvecSEXP, SEXP tolSEXP, SEXP thread_cntSEXP, SEXP max_lines_per_workSEXP, SEXP indicator_slackSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type cost_matrix(cost_matrixSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< int >::type supply_div_type(supply_div_typeSEXP);
    Rcpp::traits::input_parameter< int >::type demand_div_type(demand_div_typeSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type supply_div_parameters(supply_div_parametersSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type demand_div_parameters(demand_div_parametersSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsvec(epsvecSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type thread_cnt(thread_cntSEXP);
    Rcpp::traits::input_parameter< int >::type max_lines_per_work(max_lines_per_workSEXP);
    Rcpp::traits::input_parameter< double >::type indicator_slack(indicator_slackSEXP);
    rcpp_result_gen = Rcpp::wrap(Sinkhorn_Rcpp(cost_matrix, supply, demand, supply_div_type, demand_div_type, supply_div_parameters, demand_div_parameters, iter_max, epsvec, tol, thread_cnt, max_lines_per_work, indicator_slack));
    return rcpp_result_gen;
END_RCPP
}
// bary_sinkhorn_arma
Rcpp::List bary_sinkhorn_arma(arma::mat weights, arma::mat frechet, int maxIter, double lambda, arma::mat C, double thresh, int threads);
RcppExport SEXP _WSGeometry_bary_sinkhorn_arma(SEXP weightsSEXP, SEXP frechetSEXP, SEXP maxIterSEXP, SEXP lambdaSEXP, SEXP CSEXP, SEXP threshSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type frechet(frechetSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bary_sinkhorn_arma(weights, frechet, maxIter, lambda, C, thresh, threads));
    return rcpp_result_gen;
END_RCPP
}
// rand_cxx
double rand_cxx();
RcppExport SEXP _WSGeometry_rand_cxx() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rand_cxx());
    return rcpp_result_gen;
END_RCPP
}
// index_cxx
int index_cxx(arma::mat x, double y);
RcppExport SEXP _WSGeometry_index_cxx(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(index_cxx(x, y));
    return rcpp_result_gen;
END_RCPP
}
// norm
double norm(arma::mat x, arma::mat y);
RcppExport SEXP _WSGeometry_norm(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(norm(x, y));
    return rcpp_result_gen;
END_RCPP
}
// expo
arma::mat expo(arma::mat a);
RcppExport SEXP _WSGeometry_expo(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(expo(a));
    return rcpp_result_gen;
END_RCPP
}
// gen_cost
arma::mat gen_cost(const arma::mat A, const arma::mat B);
RcppExport SEXP _WSGeometry_gen_cost(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_cost(A, B));
    return rcpp_result_gen;
END_RCPP
}
// wsbary_cxx_armaP
Rcpp::List wsbary_cxx_armaP(const Rcpp::List weightsR, arma::mat positions1, const Rcpp::List positionssetR, const arma::mat frechet_weights, const bool fixed_support, const int maxIter, const int weights_maxIter, const int pos_maxIter, const double stepsize, const double thresh, bool headstart, const int headstartlength, int threads);
RcppExport SEXP _WSGeometry_wsbary_cxx_armaP(SEXP weightsRSEXP, SEXP positions1SEXP, SEXP positionssetRSEXP, SEXP frechet_weightsSEXP, SEXP fixed_supportSEXP, SEXP maxIterSEXP, SEXP weights_maxIterSEXP, SEXP pos_maxIterSEXP, SEXP stepsizeSEXP, SEXP threshSEXP, SEXP headstartSEXP, SEXP headstartlengthSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type weightsR(weightsRSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type positions1(positions1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type positionssetR(positionssetRSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type frechet_weights(frechet_weightsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fixed_support(fixed_supportSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const int >::type weights_maxIter(weights_maxIterSEXP);
    Rcpp::traits::input_parameter< const int >::type pos_maxIter(pos_maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< bool >::type headstart(headstartSEXP);
    Rcpp::traits::input_parameter< const int >::type headstartlength(headstartlengthSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(wsbary_cxx_armaP(weightsR, positions1, positionssetR, frechet_weights, fixed_support, maxIter, weights_maxIter, pos_maxIter, stepsize, thresh, headstart, headstartlength, threads));
    return rcpp_result_gen;
END_RCPP
}
// krbary_subgrad_cxx
Rcpp::List krbary_subgrad_cxx(const Rcpp::List weightsR, const Rcpp::List costMatsR, const arma::mat frechet_weights, const int maxIter, const double stepsize, const double thresh, bool headstart, const int headstartlength, int threads);
RcppExport SEXP _WSGeometry_krbary_subgrad_cxx(SEXP weightsRSEXP, SEXP costMatsRSEXP, SEXP frechet_weightsSEXP, SEXP maxIterSEXP, SEXP stepsizeSEXP, SEXP threshSEXP, SEXP headstartSEXP, SEXP headstartlengthSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type weightsR(weightsRSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type costMatsR(costMatsRSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type frechet_weights(frechet_weightsSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< bool >::type headstart(headstartSEXP);
    Rcpp::traits::input_parameter< const int >::type headstartlength(headstartlengthSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(krbary_subgrad_cxx(weightsR, costMatsR, frechet_weights, maxIter, stepsize, thresh, headstart, headstartlength, threads));
    return rcpp_result_gen;
END_RCPP
}
// maaipm_fixed_cpp
Rcpp::List maaipm_fixed_cpp(arma::vec p, arma::vec s, arma::vec x, const arma::vec b, const arma::vec costvec, const arma::sp_mat constMat, const int N, const int m, const int M, const arma::vec sizes, const arma::vec sizescsum, const int nr, const int nc, const arma::sp_mat U, const int maxIter, const double thresh, const int threads);
RcppExport SEXP _WSGeometry_maaipm_fixed_cpp(SEXP pSEXP, SEXP sSEXP, SEXP xSEXP, SEXP bSEXP, SEXP costvecSEXP, SEXP constMatSEXP, SEXP NSEXP, SEXP mSEXP, SEXP MSEXP, SEXP sizesSEXP, SEXP sizescsumSEXP, SEXP nrSEXP, SEXP ncSEXP, SEXP USEXP, SEXP maxIterSEXP, SEXP threshSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type costvec(costvecSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizescsum(sizescsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type U(USEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(maaipm_fixed_cpp(p, s, x, b, costvec, constMat, N, m, M, sizes, sizescsum, nr, nc, U, maxIter, thresh, threads));
    return rcpp_result_gen;
END_RCPP
}
// maaipm_free_cpp
Rcpp::List maaipm_free_cpp(arma::vec p, arma::vec s, arma::vec x, const arma::vec b, arma::vec costvec, arma::mat support, arma::mat fullsupport, const arma::sp_mat constMat, const int N, const int m, const int M, const arma::vec sizes, const arma::vec sizescsum, const int nr, const int nc, const arma::sp_mat U, const int maxIter, const int outer_maxIter, const double thresh, const int threads);
RcppExport SEXP _WSGeometry_maaipm_free_cpp(SEXP pSEXP, SEXP sSEXP, SEXP xSEXP, SEXP bSEXP, SEXP costvecSEXP, SEXP supportSEXP, SEXP fullsupportSEXP, SEXP constMatSEXP, SEXP NSEXP, SEXP mSEXP, SEXP MSEXP, SEXP sizesSEXP, SEXP sizescsumSEXP, SEXP nrSEXP, SEXP ncSEXP, SEXP USEXP, SEXP maxIterSEXP, SEXP outer_maxIterSEXP, SEXP threshSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type costvec(costvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type support(supportSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fullsupport(fullsupportSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizescsum(sizescsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type U(USEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const int >::type outer_maxIter(outer_maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(maaipm_free_cpp(p, s, x, b, costvec, support, fullsupport, constMat, N, m, M, sizes, sizescsum, nr, nc, U, maxIter, outer_maxIter, thresh, threads));
    return rcpp_result_gen;
END_RCPP
}
// transport_network_dual_arma
arma::mat transport_network_dual_arma(arma::mat a, arma::mat b, arma::mat C);
RcppExport SEXP _WSGeometry_transport_network_dual_arma(SEXP aSEXP, SEXP bSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(transport_network_dual_arma(a, b, C));
    return rcpp_result_gen;
END_RCPP
}
// transport_network_primal_arma
arma::mat transport_network_primal_arma(arma::mat a, arma::mat b, arma::mat C);
RcppExport SEXP _WSGeometry_transport_network_primal_arma(SEXP aSEXP, SEXP bSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(transport_network_primal_arma(a, b, C));
    return rcpp_result_gen;
END_RCPP
}
// treegkr_Rcpp
Rcpp::List treegkr_Rcpp(Rcpp::List& tree, Rcpp::NumericVector& supply, Rcpp::NumericVector& demand, Rcpp::NumericVector& creation, Rcpp::NumericVector& destruction);
RcppExport SEXP _WSGeometry_treegkr_Rcpp(SEXP treeSEXP, SEXP supplySEXP, SEXP demandSEXP, SEXP creationSEXP, SEXP destructionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type creation(creationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type destruction(destructionSEXP);
    rcpp_result_gen = Rcpp::wrap(treegkr_Rcpp(tree, supply, demand, creation, destruction));
    return rcpp_result_gen;
END_RCPP
}
// umaaipm_free_cpp
Rcpp::List umaaipm_free_cpp(arma::vec p, arma::vec s, arma::vec x, const arma::vec b, arma::vec costvec, arma::mat support, arma::mat fullsupport, const arma::sp_mat constMat, const int N, const int m, const int M, const arma::vec sizes, const arma::vec sizescsum, const int nr, const int nc, const arma::sp_mat U, const double C, const int maxIter, const int outer_maxIter, const double thresh, const int threads);
RcppExport SEXP _WSGeometry_umaaipm_free_cpp(SEXP pSEXP, SEXP sSEXP, SEXP xSEXP, SEXP bSEXP, SEXP costvecSEXP, SEXP supportSEXP, SEXP fullsupportSEXP, SEXP constMatSEXP, SEXP NSEXP, SEXP mSEXP, SEXP MSEXP, SEXP sizesSEXP, SEXP sizescsumSEXP, SEXP nrSEXP, SEXP ncSEXP, SEXP USEXP, SEXP CSEXP, SEXP maxIterSEXP, SEXP outer_maxIterSEXP, SEXP threshSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type costvec(costvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type support(supportSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fullsupport(fullsupportSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizescsum(sizescsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type U(USEXP);
    Rcpp::traits::input_parameter< const double >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const int >::type outer_maxIter(outer_maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(umaaipm_free_cpp(p, s, x, b, costvec, support, fullsupport, constMat, N, m, M, sizes, sizescsum, nr, nc, U, C, maxIter, outer_maxIter, thresh, threads));
    return rcpp_result_gen;
END_RCPP
}
// umaaipm_fixed_cpp
Rcpp::List umaaipm_fixed_cpp(arma::vec p, arma::vec s, arma::vec x, const arma::vec b, arma::vec costvec, arma::mat support, arma::mat fullsupport, const arma::sp_mat constMat, const int N, const int m, const int M, const arma::vec sizes, const arma::vec sizescsum, const int nr, const int nc, const arma::sp_mat U, const double C, const int maxIter, const double thresh, const int threads);
RcppExport SEXP _WSGeometry_umaaipm_fixed_cpp(SEXP pSEXP, SEXP sSEXP, SEXP xSEXP, SEXP bSEXP, SEXP costvecSEXP, SEXP supportSEXP, SEXP fullsupportSEXP, SEXP constMatSEXP, SEXP NSEXP, SEXP mSEXP, SEXP MSEXP, SEXP sizesSEXP, SEXP sizescsumSEXP, SEXP nrSEXP, SEXP ncSEXP, SEXP USEXP, SEXP CSEXP, SEXP maxIterSEXP, SEXP threshSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type costvec(costvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type support(supportSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fullsupport(fullsupportSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sizescsum(sizescsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type U(USEXP);
    Rcpp::traits::input_parameter< const double >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(umaaipm_fixed_cpp(p, s, x, b, costvec, support, fullsupport, constMat, N, m, M, sizes, sizescsum, nr, nc, U, C, maxIter, thresh, threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_WSGeometry_StabilizedScaling_Rcpp", (DL_FUNC) &_WSGeometry_StabilizedScaling_Rcpp, 13},
    {"_WSGeometry_Sinkhorn_Rcpp", (DL_FUNC) &_WSGeometry_Sinkhorn_Rcpp, 13},
    {"_WSGeometry_bary_sinkhorn_arma", (DL_FUNC) &_WSGeometry_bary_sinkhorn_arma, 7},
    {"_WSGeometry_rand_cxx", (DL_FUNC) &_WSGeometry_rand_cxx, 0},
    {"_WSGeometry_index_cxx", (DL_FUNC) &_WSGeometry_index_cxx, 2},
    {"_WSGeometry_norm", (DL_FUNC) &_WSGeometry_norm, 2},
    {"_WSGeometry_expo", (DL_FUNC) &_WSGeometry_expo, 1},
    {"_WSGeometry_gen_cost", (DL_FUNC) &_WSGeometry_gen_cost, 2},
    {"_WSGeometry_wsbary_cxx_armaP", (DL_FUNC) &_WSGeometry_wsbary_cxx_armaP, 13},
    {"_WSGeometry_krbary_subgrad_cxx", (DL_FUNC) &_WSGeometry_krbary_subgrad_cxx, 9},
    {"_WSGeometry_maaipm_fixed_cpp", (DL_FUNC) &_WSGeometry_maaipm_fixed_cpp, 17},
    {"_WSGeometry_maaipm_free_cpp", (DL_FUNC) &_WSGeometry_maaipm_free_cpp, 20},
    {"_WSGeometry_transport_network_dual_arma", (DL_FUNC) &_WSGeometry_transport_network_dual_arma, 3},
    {"_WSGeometry_transport_network_primal_arma", (DL_FUNC) &_WSGeometry_transport_network_primal_arma, 3},
    {"_WSGeometry_treegkr_Rcpp", (DL_FUNC) &_WSGeometry_treegkr_Rcpp, 5},
    {"_WSGeometry_umaaipm_free_cpp", (DL_FUNC) &_WSGeometry_umaaipm_free_cpp, 21},
    {"_WSGeometry_umaaipm_fixed_cpp", (DL_FUNC) &_WSGeometry_umaaipm_fixed_cpp, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_WSGeometry(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
