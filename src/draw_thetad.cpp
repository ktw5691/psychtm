#include <RcppArmadillo.h>

#include "rdirichlet_cpp.h"

//' @title Sample \eqn{\theta_d} from full conditional distribution
//'
//' @name draw_thetad
//' @param z_count A K x 1 vector of counts of topic draw in document d.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//'
//' @return A K x 1 vector draw of \eqn{\theta_d}.
//'
//' @noRd
// [[Rcpp::export(.draw_thetad)]]
arma::mat draw_thetad(const arma::rowvec& z_count, float alpha_) {

  const uint16_t K = z_count.size(); // Number of topics
  if (alpha_ <= 0.0) Rcpp::stop("alpha_ must be positive");
  if (K < 2) Rcpp::stop("number of topics must be at least 2");

  arma::rowvec alp(K);
  alp.fill(alpha_);

  return rdirichlet_cpp(1, z_count + alp);;
}
