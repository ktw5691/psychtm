#ifndef WAIC_H
#define WAIC_H

#include <RcppArmadillo.h>

double pwaic_d(const arma::colvec& like_pred);

double waic_d(const arma::colvec& like_pred, double p_effd);

Rcpp::NumericVector waic_all(uint32_t iter, const arma::mat& l_pred);

Rcpp::NumericVector waic_diff(const arma::mat& l_pred1, const arma::mat& l_pred2);

#endif
