#include <RcppArmadillo.h>

#include "get_loglike.h"

//' @title Posterior predictive likelihood for SLDA/SLDAX/MLR
//'
//' @name post_pred_norm
//' @parm y A D x 1 vector of outcomes.
//' @param w A D x q matrix of additional predictors.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The residual variance.
//'
//' @return Posterior predictive likelihood.
//'
//' @noRd
arma::rowvec post_pred_norm(const arma::colvec& y, const arma::mat& w,
                            const arma::colvec& eta, double sigma2) {

  arma::colvec mu_hat = w * eta;
  // Log-likelihood of observed y given current posterior draw of eta, sigma2
  arma::colvec loglike_pred = -0.5 * (log(2.0 * M_PI) + log(sigma2) +
    arma::square(y - mu_hat) / sigma2);

  return exp(loglike_pred.t());
}

//' @title Posterior predictive likelihood for logistic SLDA/SLDAX/regression
//'
//' @name post_pred_logit
//' @parm y A D x 1 vector of outcomes.
//' @param w A D x q matrix of additional predictors.
//' @param eta A q x 1 vector of regression coefficients.
//' @return Predictive posterior likelihood of all D observations.
//'
//' @noRd
arma::rowvec post_pred_logit(const arma::colvec& y, const arma::mat& w, const arma::colvec& eta) {

  const uint16_t D = w.n_rows;
  arma::rowvec loglike_pred(D);
  arma::colvec mu_hat = w * eta;
  for (uint16_t d = 0; d < D; d++) {
    loglike_pred(d) = get_ll_logit_yd(y(d), arma::as_scalar(mu_hat(d)));
  }
  return exp(loglike_pred);
}
