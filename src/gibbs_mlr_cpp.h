#ifndef GIBBS_MLR_CPP_H
#define GIBBS_MLR_CPP_H

#include <RcppArmadillo.h>

Rcpp::S4 gibbs_mlr_cpp(uint32_t m, uint32_t burn, uint32_t thin,
                       const arma::colvec& y, const arma::mat& x,
                       const arma::colvec& mu0, const arma::mat& sigma0,
                       arma::colvec eta_start, float a0 = 0.001, float b0 = 0.001,
                       bool verbose = false, bool display_progress = false);

#endif
