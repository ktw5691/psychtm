#include <RcppArmadillo.h>

//' @title Compute log-numerator vector for sampling zdn from full conditional distribution for LDA
//'
//' @name get_log_numer_samplez
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A K x 1 vector of the log-numerator from the LDA model to sample zdn.
//'
//' @noRd
arma::vec get_log_numer_samplez(uint32_t V, const arma::vec& ndk_n,
                                const arma::vec& nkm_n, const arma::vec& nk_n,
                                float alpha_, float gamma_) {

  const uint16_t K = nk_n.size(); // Number of topics
  if (K < 2) Rcpp::stop("number of topics must be at least 2");
  if (V < 2) Rcpp::stop("size of vocabulary V must be at least 2");
  if (ndk_n.size() != K) Rcpp::stop("ndk_n must be a vector of length K");
  if (nkm_n.size() != K) Rcpp::stop("nkm_n must be a vector of length K");
  if (nk_n.size() != K) Rcpp::stop("nk_n must be a vector of length K");
  if (alpha_ < 0.0) Rcpp::stop("alpha_ must be positive");
  if (gamma_ < 0.0) Rcpp::stop("gamma_ must be positive");

  arma::vec log_num =
    log(ndk_n + alpha_) +
    log(nkm_n + gamma_) -
    log(nk_n + static_cast<float>(V) * gamma_);

  return log_num;
}
