#include <RcppArmadillo.h>
#include "rmvnorm_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Draw eta from full conditional posterior for sLDA/sLDAX/MLR models
//'
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param sigma2 The residual variance.
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//'
//' @return A q x 1 vector of draws of eta.
arma::colvec draw_eta_norm(const arma::mat& w, const arma::vec& y,
                           long double sigma2, const arma::vec& mu0,
                           const arma::mat& sigma0) {

  arma::mat wtw = w.t() * w;
  arma::mat sigma0_inv = arma::inv_sympd(sigma0);
  arma::mat sigma1 = arma::inv_sympd(sigma0_inv + wtw / sigma2);
  arma::colvec eta1 = sigma1 * (sigma0_inv * mu0 + w.t() * y / sigma2);

  return rmvnorm_cpp(1, eta1, sigma1).t();
}
