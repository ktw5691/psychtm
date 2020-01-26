#include <RcppArmadillo.h>

#include "get_loglike.h"
#include "invlogit.h"
#include "rmvnorm_cpp.h"

//' @title Posterior predictive likelihood for sLDA/sLDAX/MLR
//'
//' @name post_pred_norm
//' @param w A D x q matrix of additional predictors.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The residual variance.
//'
//' @return Posterior predictive likelihood.
arma::rowvec post_pred_norm(const arma::mat& w,
                            const arma::colvec& eta, double sigma2) {

  const uint16_t D = w.n_rows;
  arma::colvec mu_hat = w * eta;

  arma::mat sigma(D, D, arma::fill::eye); // D x D identity matrix
  // Draw D new observations
  arma::colvec yhat = rmvnorm_cpp(1, mu_hat, sigma2 * sigma).t();
  arma::colvec loglike_pred = -0.5 / sigma2 * arma::square(yhat - mu_hat);

  return exp(loglike_pred.t());
}

//' @title Posterior predictive likelihood for logistic sLDA/sLDAX/regression
//'
//' @name post_pred_logit
//' @param w A D x q matrix of additional predictors.
//' @param eta A q x 1 vector of regression coefficients.
//' @return Predictive posterior likelihood of all D observations.
arma::rowvec post_pred_logit(const arma::mat& w, const arma::colvec& eta) {

  const uint16_t D = w.n_rows;
  arma::rowvec loglike_pred(D);
  arma::colvec mu_hat = w * eta;
  for (uint16_t d = 0; d < D; d++) {
    double phat = invlogit(arma::as_scalar(mu_hat(d)));
    uint16_t yhat = Rcpp::rbinom(1, 1, phat)(0);
    loglike_pred(d) = get_ll_logit_yd(yhat, arma::as_scalar(mu_hat(d)));
  }
  return exp(loglike_pred);
}
