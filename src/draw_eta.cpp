#include <RcppArmadillo.h>

#include "get_logpost.h"
#include "rmvnorm_cpp.h"

//' @title Draw eta from full conditional posterior for SLDA/SLDAX/MLR models
//'
//' @name draw_eta_norm
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param sigma2 The residual variance.
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//'
//' @return A q x 1 vector of draws of eta.
//'
//' @noRd
arma::colvec draw_eta_norm(const arma::mat& w, const arma::vec& y,
                           long double sigma2, const arma::vec& mu0,
                           const arma::mat& sigma0) {

  arma::mat wtw = w.t() * w;
  arma::mat sigma0_inv = arma::inv_sympd(sigma0);
  arma::mat sigma1 = arma::inv_sympd(sigma0_inv + wtw / sigma2);
  arma::colvec eta1 = sigma1 * (sigma0_inv * mu0 + w.t() * y / sigma2);

  return rmvnorm_cpp(1, eta1, sigma1).t();
}

//' @title Draw eta from full conditional posterior for logistic SLDA/SLDAX/regression
//'   using Metropolis-Hastings (MH) algorithm
//'
//' @name draw_eta_logit
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param eta_prev A q x 1 vector of the previous draw of the regression
//'   coefficients.
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//' @param proposal_sd A q x 1 vector of proposal distribution standard deviations.
//' @param attempt A vector of the number of current attempted draws of eta.
//' @param accept A vector of the number of accepted draws of eta.
//'
//' @return A q x 1 vector of draws for eta.
//'
//' @noRd
arma::colvec draw_eta_logit(const arma::mat& w, const arma::colvec& y,
                            const arma::colvec& eta_prev,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            const arma::vec& proposal_sd,
                            arma::vec& attempt, arma::vec& accept) {

  const uint16_t q = w.n_cols;
  arma::colvec cand_eta = eta_prev; // Candidate draws of eta
  arma::colvec eta = eta_prev;
  double cur_logpost = eta_logpost_logit(w, y, eta_prev, mu0, sigma0);

  for (uint16_t j = 0; j < q; j++) {
    cand_eta(j) = Rcpp::rnorm(1, eta_prev(j), proposal_sd(j))(0);
    attempt(j)++; // Passed by reference to update outside function

    // Compute acceptance ratio
    double cand_logpost = eta_logpost_logit(w, y, cand_eta, mu0, sigma0);
    double log_r = cand_logpost - cur_logpost; // Symmetric proposals
    double log_u = log(Rcpp::runif(1)(0));
    if (log_r > log_u) {
      eta(j)      = cand_eta(j);
      cur_logpost = cand_logpost;
      accept(j)++; // Passed by reference to update outside function
    } else {
      eta(j)      = eta_prev(j);
      cand_eta(j) = eta_prev(j);
    }
  }
  return eta;
}
