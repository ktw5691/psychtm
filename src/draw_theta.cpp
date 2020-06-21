#include <RcppArmadillo.h>

#include "draw_thetad.h"

//' @title Sample \eqn{\Theta} from full conditional distribution
//'
//' @name draw_theta
//' @param m Number of samples
//' @param z_count A D x K x m matrix of counts of topic draws (columns) in
//'   documents (rows) in m posterior samples.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//'
//' @return A D x K x m array of m draws of \eqn{\Theta}.
//' @export
// [[Rcpp::export(.draw_theta)]]
arma::cube draw_theta(uint32_t m, const arma::ucube& z_count, float alpha_) {
  uint32_t D = z_count.n_rows;
  uint16_t K = z_count.n_cols;
  arma::cube drawsm = arma::cube(D, K, m, arma::fill::zeros);
  // Loop over samples
  for (uint32_t i = 0; i < m; i++) {
    // Loop over documents
    arma::mat theta = arma::zeros(D, K);
    for (uint32_t d = 0; d < D; d++) {
      theta.row(d) = draw_thetad(z_count.slice(i).row(d), alpha_);
    }
    drawsm.slice(i) = theta;
  }
  return drawsm;
}
