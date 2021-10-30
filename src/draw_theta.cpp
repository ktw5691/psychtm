#include <RcppArmadillo.h>

#include "draw_thetad.h"

//' @title Sample \eqn{\Theta} from full conditional distribution
//'
//' @name draw_theta
//' @param z_count A D x K matrix of counts of topic draws (columns) in
//'   documents (rows).
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//'
//' @return A D x K matrix \eqn{\Theta}.
//'
//' @noRd
// [[Rcpp::export(.draw_theta)]]
arma::mat draw_theta(const arma::mat& z_count, float alpha_) {
  uint32_t D = z_count.n_rows;
  uint16_t K = z_count.n_cols;
  arma::mat theta = arma::zeros(D, K);
  // Loop over documents
  for (uint32_t d = 0; d < D; d++) {
    theta.row(d) = draw_thetad(z_count.row(d), alpha_);
  }
  return theta;
}
