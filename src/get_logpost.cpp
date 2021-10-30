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
//'
//' @noRd
double get_lpost_eta(double ll, const arma::colvec& eta,
                     const arma::colvec& mu0, const arma::mat& sigma0) {

  uint32_t q = eta.size();
  double lp_temp = ll;
  // Add log-prior on eta
  double temp_prod = arma::as_scalar(
    (eta - mu0).t() * arma::inv_sympd(sigma0) * (eta - mu0)
  );
  lp_temp += (-0.5 * ( static_cast<float>(q) * (log(2.0 * M_PI) + log(arma::det(sigma0))) + temp_prod ));

  return lp_temp;
}

//' @title Log-posterior contribution for normal outcome priors on coefs and sigma2 + log-likelihood
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
//' @return The current log-posterior contribution.
//'
//' @noRd
double get_lpost_mlr(double ll,
                     const arma::colvec& eta, double sigma2,
                     const arma::colvec& mu0, const arma::mat& sigma0,
                     double a0, double b0) {

  // Log-likelihood + log-prior on eta
  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);

  // Add log-prior on sigma2 ~ Inv-Gamma(a0 / 2, b0 / 2)
  lp_temp += (0.5 * a0) * log(0.5 * b0) - lgamma(0.5 * a0) +
    (-0.5 * a0 - 1.0) * log(sigma2) - 0.5 * b0 / sigma2;

  return lp_temp;
}

//' @title Log-posterior contribution from LDA model priors + log-likelihood
//'
//' @name get_lpost_lda
//' @param ll A double of the current log-likelihood.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//'
//' @return The current log-posterior contribution.
//'
//' @noRd
double get_lpost_lda(double ll, const arma::mat& theta, const arma::mat& beta,
                     double gamma_, double alpha_) {

  float D = static_cast<float>(theta.n_rows);
  float K = static_cast<float>(theta.n_cols);
  float V = static_cast<float>(beta.n_cols);

  // Add log-prior on beta matrix
  double sum_logbeta = arma::accu(log(beta));
  double beta_contrib = K * lgamma(V * gamma_) -
    V * K * lgamma(gamma_) +
    (gamma_ - 1.0) * sum_logbeta;

  // Add log-prior on theta matrix
  double sum_logtheta = arma::accu(log(theta)); // Sum log(theta_{di}) across docs and topics
  double theta_contrib = D * lgamma(K * alpha_) -
    K * D * lgamma(alpha_) +
    (alpha_ - 1) * sum_logtheta;

  return ll + beta_contrib + theta_contrib;
}

//' @title Log-posterior for SLDA/SLDAX model
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
//'
//' @return The current log-posterior.
//'
//' @noRd
double get_lpost_slda_norm(double ll, const arma::colvec& eta, double sigma2,
                           const arma::mat& theta, const arma::mat& beta,
                           const arma::colvec& mu0, const arma::mat& sigma0,
                           double gamma_, double alpha_,
                           double a0, double b0) {

  // Combine log-likelihood and log-priors for regression coefs and residual variance
  double lp_temp = get_lpost_mlr(ll, eta, sigma2, mu0, sigma0, a0, b0);
  // Add contribution from LDA log-priors
  lp_temp = get_lpost_lda(lp_temp, theta, beta, gamma_, alpha_);
  return lp_temp;
}

//' @title Log-posterior for logistic SLDA/SLDAX model
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
//'
//' @return The current log-posterior.
//'
//' @noRd
double get_lpost_slda_logit(double ll, const arma::colvec& eta,
                            const arma::mat& theta, const arma::mat& beta,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            double gamma_, double alpha_) {

  // Combine log-likelihood and log-prior for regression coefs
  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);
  // Add log-priors from LDA model
  lp_temp = get_lpost_lda(lp_temp, theta, beta, gamma_, alpha_);
  return lp_temp;
}

//' @title Compute full conditional log-posterior of eta for logistic SLDA/SLDAX/regression
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
//'
//' @noRd
double eta_logpost_logit(const arma::mat& w, const arma::colvec& y,
                         const arma::vec& eta,
                         const arma::vec& mu0, const arma::mat& sigma0) {

  // Compute log-likelihood of y
  double loglike = get_ll_logit(y, w, eta);
  // Add log-prior on eta
  double logpost = get_lpost_eta(loglike, eta, mu0, sigma0);
  return logpost;
}
