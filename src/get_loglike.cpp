#include <RcppArmadillo.h>

//' @title Log-likelihood for logistic regression for observation d
//'
//' @name get_ll_logit_yd
//' @param yd An integer 0/1 outcome to be predicted.
//' @param muhatd A double predicted outcome on logit scale.
//'
//' @return The current log-likelihood for observation d.
//'
//' @noRd
double get_ll_logit_yd(bool yd, double muhatd) {

  // Compute log-likelihood of y
  double ll_temp = yd * muhatd - log(1.0 + exp(muhatd));
  return ll_temp;
}

//' @title Log-likelihood for logistic regression
//'
//' @name get_ll_logit
//' @param y A D x 1 vector of 0/1 outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix.
//' @param eta A q x 1 vector of regression coefficients.
//'
//' @return The current log-likelihood.
//'
//' @noRd
double get_ll_logit(const arma::colvec& y, const arma::mat& w,
                    const arma::colvec& eta) {

  // Add likelihood of y
  uint32_t D = w.n_rows;
  arma::colvec muhat = w * eta;

  // Compute log-likelihood of y
  double ll_temp = 0.0;
  for (uint32_t d = 0; d < D; d++) {
    ll_temp += get_ll_logit_yd(y(d), arma::as_scalar(muhat(d)));
  }
  return ll_temp;
}

//' @title Log-likelihood for MLR
//'
//' @name get_ll_mlr
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//'
//' @return The current log-likelihood.
//'
//' @noRd
double get_ll_mlr(const arma::colvec& y, const arma::mat& w,
                  const arma::colvec& eta, double sigma2) {

  uint32_t D = w.n_rows;
  double temp_prod = arma::as_scalar(
    (y - w * eta).t() * (y - w * eta)
  );
  double ll_temp = -0.5 * (static_cast<float>(D) * (log(2.0 * M_PI) + log(sigma2)) + temp_prod / sigma2);

  return ll_temp;
}

//' @title Log-likelihood for LDA model
//'
//' @name get_ll_lda
//' @param zdocs A D x max(\eqn{N_d}) matrix of topic indicators for all documents.
//' @param docs A D x max(\eqn{N_d}) matrix of word indicators for all documents.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//' @param N A vector of length D containing the number of words in each document.
//'
//' @return The current log-likelihood.
//'
//' @noRd
double get_ll_lda(const arma::umat& zdocs, const arma::umat& docs,
                  const arma::mat& theta, const arma::mat& beta,
                  const Rcpp::IntegerVector& docs_index,
                  const arma::colvec& N) {

  double ll_temp = 0.0;
  // Add likelihood of documents
  for (uint32_t d : docs_index) {
    for (uint32_t n = 0; n < N(d); n++) {
      uint16_t zdn = zdocs(d, n) - 1; // Topic for word dn (0-indexed)
      uint32_t wdn = docs(d, n) - 1;  // Word for word dn (0-indexed)
      ll_temp += log(theta(d, zdn));  // f(z_{dn} | theta_d)
      ll_temp += log(beta(zdn, wdn)); // f(w_{dn} | z_{dn}, beta_{z_{dn}})
    }
  }
  return ll_temp;
}

//' @title Log-likelihood for SLDA/SLDAX model
//'
//' @name get_ll_slda_norm
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//' @param zdocs A D x max(\eqn{N_d}) matrix of topic indicators for all documents.
//' @param docs A D x max(\eqn{N_d}) matrix of word indicators for all documents.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//' @param N A vector of length D containing the number of words in each document.
//'
//' @return The current log-likelihood.
//'
//' @noRd
double get_ll_slda_norm(const arma::colvec& y, const arma::mat& w,
                        const arma::colvec& eta, double sigma2,
                        const arma::umat& zdocs, const arma::umat& docs,
                        const arma::mat& theta, const arma::mat& beta,
                        const Rcpp::IntegerVector& docs_index,
                        const arma::colvec& N) {

  double ll_temp = get_ll_mlr(y, w, eta, sigma2) +
    get_ll_lda(zdocs, docs, theta, beta, docs_index, N);
  return ll_temp;
}

//' @title Log-likelihood for logistic SLDA/SLDAX model
//'
//' @name get_ll_slda_logit
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param zdocs A D x max(\eqn{N_d}) matrix of topic indicators for all documents.
//' @param docs A D x max(\eqn{N_d}) matrix of word indicators for all documents.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//' @param N A vector of length D containing the number of words in each document.
//'
//' @return The current log-likelihood.
//'
//' @noRd
double get_ll_slda_logit(const arma::colvec& y, const arma::mat& w,
                         const arma::colvec& eta,
                         const arma::umat& zdocs, const arma::umat& docs,
                         const arma::mat& theta, const arma::mat& beta,
                         const Rcpp::IntegerVector& docs_index,
                         const arma::colvec& N) {

  double ll_temp = get_ll_logit(y, w, eta) +
    get_ll_lda(zdocs, docs, theta, beta, docs_index, N);
  return ll_temp;
}
