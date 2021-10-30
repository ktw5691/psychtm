#include <RcppArmadillo.h>

#include "rdirichlet_cpp.h"

//' @title Sample \eqn{\beta_k} from full conditional distribution
//'
//' @name draw_betak
//' @param wz_co A V x 1 vector of counts of the draws of each word for topic
//'   k over all documents.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A V x 1 vector of estimates for \eqn{\beta_k}.
//'
//' @noRd
// [[Rcpp::export(.draw_betak)]]
arma::mat draw_betak(const arma::rowvec& wz_co, float gamma_) {

  const uint32_t V = wz_co.size(); // Number of topics
  if (gamma_ <= 0.0) Rcpp::stop("gamma_ must be positive");
  if (V < 2) Rcpp::stop("vocabulary size must be at least 2");

  arma::rowvec gam(V);
  gam.fill(gamma_);
  arma::mat betak = rdirichlet_cpp(1, wz_co + gam);

  return betak;
}
