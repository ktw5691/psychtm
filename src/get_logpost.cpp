#include <RcppArmadillo.h>

#include "get_loglike.h"

//' @title Log-posterior for regression coefficients with normal prior
//'
//' @name get_lpost_eta
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//'
//' @return The current log-posterior.
double get_lpost_eta(double ll, const arma::colvec& eta,
                     const arma::colvec& mu0, const arma::mat& sigma0) {

  double lp_temp = ll;
  // Add prior on eta
  double temp_prod = arma::as_scalar(
    (eta - mu0).t() * arma::inv_sympd(sigma0) * (eta - mu0)
  );
  lp_temp += (-0.5 * temp_prod);

  return lp_temp;
}

//' @title Log-posterior for normal outcome regression
//'
//' @name get_lpost_mlr
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param a0 The shape parameter for the prior on sigma2.
//' @param b0 The scale parameter for the prior on sigma2.
//'
//' @return The current log-posterior.
double get_lpost_mlr(double ll,
                     const arma::colvec& eta, double sigma2,
                     const arma::colvec& mu0, const arma::mat& sigma0,
                     double a0, double b0) {

  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);

  // Add prior on sigma2
  lp_temp += ((-0.5 * a0 - 1.0) * log(sigma2) - 0.5 * b0 / sigma2);

  return lp_temp;
}

//' @title Log-posterior for LDA model
//'
//' @name get_lpost_lda
//' @param ll A double of the current log-likelihood.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//' @param V The number of words in the vocabulary.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//'
//' @return The current log-posterior.
double get_lpost_lda(double ll, const arma::mat& theta, const arma::mat& beta,
                     double gamma_, double alpha_,
                     uint32_t V, const Rcpp::IntegerVector& docs_index) {

  // Add prior on beta matrix
  double temp_betapost = arma::accu(log(beta));
  double lp_temp = ll + ((gamma_ - 1.0) * temp_betapost);

  // Add prior on theta matrix
  double temp_thetapost = arma::accu(log(theta));
  lp_temp += ((alpha_ - 1) * temp_thetapost);

  return lp_temp;
}

//' @title Log-posterior for sLDA/sLDAX model
//'
//' @name get_lpost_slda_norm
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//' @param a0 The shape parameter for the prior on sigma2.
//' @param b0 The scale parameter for the prior on sigma2.
//' @param V The number of words in the vocabulary.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//'
//' @return The current log-posterior.
double get_lpost_slda_norm(double ll, const arma::colvec& eta, double sigma2,
                           const arma::mat& theta, const arma::mat& beta,
                           const arma::colvec& mu0, const arma::mat& sigma0,
                           double gamma_, double alpha_,
                           double a0, double b0, uint32_t V,
                           const Rcpp::IntegerVector& docs_index) {

  double lp_temp = get_lpost_mlr(ll, eta, sigma2, mu0, sigma0, a0, b0);
  lp_temp += get_lpost_lda(lp_temp, theta, beta, gamma_, alpha_, V, docs_index);
  return lp_temp;
}

//' @title Log-posterior for logistic sLDA/sLDAX model
//'
//' @name get_lpost_slda_logit
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//' @param V The number of words in the vocabulary.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//'
//' @return The current log-posterior.
double get_lpost_slda_logit(double ll, const arma::colvec& eta,
                            const arma::mat& theta, const arma::mat& beta,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            double gamma_, double alpha_,
                            uint32_t V, const Rcpp::IntegerVector& docs_index) {

  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);
  lp_temp += get_lpost_lda(lp_temp, theta, beta, gamma_, alpha_, V, docs_index);
  return lp_temp;
}

//' @title Compute full conditional log-posterior of eta for logistic sLDA/sLDAX/regression
//'
//' @name eta_logpost_logit
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param eta A q x 1 vector of regression coefficients
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//'
//' @return The full conditional log-posterior density of eta.
double eta_logpost_logit(const arma::mat& w, const arma::colvec& y,
                         const arma::vec& eta,
                         const arma::vec& mu0, const arma::mat& sigma0) {

  // Compute log-likelihood of y
  double loglike = get_ll_logit(y, w, eta);
  // Add log-prior on eta
  double logpost = get_lpost_eta(loglike, eta, mu0, sigma0);
  return logpost;
}
