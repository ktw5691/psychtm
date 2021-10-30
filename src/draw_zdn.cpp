#include <RcppArmadilloExtensions/sample.h>

#include "get_loglike.h"
#include "get_log_numer_samplez.h"

//' @title Draw zdn from full conditional distribution for LDA/SLDA/SLDAX
//'
//' @name draw_zdn
//' @param log_num A K x 1 vector of the log-numerator for sampling from topics
//'   1, ..., K.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
//'
//' @noRd
uint16_t draw_zdn(arma::vec& log_num) {

  const uint16_t K = log_num.size(); // Number of topics
  long double denom = sum(exp(log_num));
  arma::vec pmf = exp(log_num - log(denom));

  bool good_pmf = true;
  for (uint16_t k = 0; k < K; k++) {
    if (std::isnan(pmf(k)) || std::isinf(pmf(k))) good_pmf = false;
  }

  if (!good_pmf) pmf = 1.0 / static_cast<float>(K);
  Rcpp::IntegerVector topics = Rcpp::seq_len(K);
  Rcpp::IntegerVector zdn = Rcpp::RcppArmadillo::sample(topics, 1, true, pmf);
  return zdn(0);
}

//' @title Draw zdn from full conditional distribution for LDA
//'
//' @name draw_zdn_lda
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
//'
//' @noRd
uint16_t draw_zdn_lda(uint32_t V,
                      const arma::vec& ndk_n, const arma::vec& nkm_n,
                      const arma::vec& nk_n, float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K
  arma::vec log_num = get_log_numer_samplez(V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_);
  uint16_t zdn = draw_zdn(log_num);
  return zdn;
}

//' @title Draw zdn from full conditional distribution for SLDA/SLDAX
//'
//' @name draw_zdn_slda_norm
//' @param yd The outcome variable for document \eqn{d}.
//' @param w_d A q x 1 vector containing row d of the predictor model matrix of
//'   assumed form (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The residual variance.
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
//'
//' @noRd
uint16_t draw_zdn_slda_norm(double yd, const arma::rowvec& w_d,
                            const arma::vec& eta, double sigma2,
                            uint32_t V, const arma::vec& ndk_n,
                            const arma::vec& nkm_n, const arma::vec& nk_n,
                            float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;
  //   this tends to occur if sigma^2 draw near 0
  if (sigma2 < 0.0) Rcpp::stop("sigma2 must be positive");
  arma::vec log_num = get_log_numer_samplez(
    V, ndk_n, nkm_n, nk_n, alpha_, gamma_) +
      R::dnorm(yd, dot(w_d, eta), sqrt(sigma2), true);
  uint16_t zdn = draw_zdn(log_num);
  return zdn;
}

//' @title Draw zdn from full conditional distribution for SLDA/SLDAX with binary outcome
//'
//' @name draw_zdn_slda_logit
//' @param yd A the outcome variable for document d.
//' @param w_d A q x 1 vector containing row d of the predictor model matrix of
//'   assumed form (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
//'
//' @noRd
uint16_t draw_zdn_slda_logit(double yd, const arma::rowvec& w_d,
                             const arma::vec& eta,
                             uint32_t V, const arma::vec& ndk_n,
                             const arma::vec& nkm_n, const arma::vec& nk_n,
                             float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;
  double muhat = dot(w_d, eta);
  // Compute log-likelihood of y_d
  double loglike = get_ll_logit_yd(yd, muhat);

  arma::vec log_num = get_log_numer_samplez(V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_) +
                                              loglike;
  uint16_t zdn = draw_zdn(log_num);
  return zdn;
}
