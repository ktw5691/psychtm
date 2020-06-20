#ifndef POST_PRED_H
#define POST_PRED_H

#include <RcppArmadillo.h>

arma::rowvec post_pred_norm(const arma::colvec& y, const arma::mat& w,
                            const arma::colvec& eta, double sigma2);

arma::rowvec post_pred_logit(const arma::colvec& y, const arma::mat& w, const arma::colvec& eta);

#endif
