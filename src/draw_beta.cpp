#include <RcppArmadillo.h>

#include "draw_betak.h"

//' @title Sample \eqn{B} from full conditional distribution
//'
//' @name draw_beta
//' @param wz_co A K x V matrix of counts of word-topic co-occurrences
//'   (topics: columns; words: rows).
//' @param gamma_ The hyperparameter on the Dirichlet prior for \eqn{\beta_k}.
//'
//' @return A K x V matrix \eqn{B}.
//'
//' @noRd
// [[Rcpp::export(.draw_beta)]]

arma::mat draw_beta(const arma::mat& wz_co, float gamma_) {
  uint16_t K = wz_co.n_rows;
  uint32_t V = wz_co.n_cols;
  // Loop over topics
  arma::mat beta = arma::zeros(K, V);
  for (uint16_t k = 0; k < K; k++) {
    beta.row(k) = draw_betak(wz_co.row(k), gamma_);
  }
  return beta;
}
