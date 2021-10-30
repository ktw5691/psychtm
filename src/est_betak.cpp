#include <RcppArmadillo.h>

//' @title Estimate \eqn{\beta_k}
//'
//' @name est_betak
//' @param wz_co A V x 1 vector of counts of the draws of each word for topic
//'   k over all documents.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A V x 1 vector of estimates for \eqn{\beta_k}.
//'
//' @noRd
// [[Rcpp::export(.est_betak)]]
arma::rowvec est_betak(const arma::rowvec& wz_co, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic word to 0.0;
  if (gamma_ < 0.0) Rcpp::stop("gamma_ must be positive");
  const uint32_t V = wz_co.size(); // Vocabulary size
  if (V < 2) Rcpp::stop("vocabulary size V must be at least 2");

  arma::rowvec betak = exp(log(wz_co + gamma_) - log(sum(wz_co) + V * gamma_));

  // Check for impossible estimates and replace with 0
  betak.transform( [](double val) {
    if ((val > 1.0) | (std::isnan(val))) {
      return(0.0);
    } else {
      return val;
    }
  });

  return betak;
}
