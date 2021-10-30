#include <RcppArmadillo.h>

//' @title Sample from K-dimensional Dirichlet distribution Dir(\eqn{\alpha})
//'
//' @name rdirichlet_cpp
//' @param n The number of samples to draw.
//' @param alpha The K x 1 vector of concentration parameters (column vector).
//'
//' @return A n x K matrix of random draws.
//'
//' @noRd
arma::mat rdirichlet_cpp(uint32_t n, const arma::rowvec& alpha_) {

  // Draw Dirichlet rv using Gamma distribution
  // See https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation

  const uint16_t K = alpha_.size();
  arma::mat draws = arma::randn(n, K);

  for (uint32_t i = 0; i < n; i++) {
    long double norm = 0;
    // Draw each component from Gamma rv
    for (uint16_t k = 0; k < K; k++) {
      long double cur = R::rgamma(alpha_[k], 1.0); // Shape, Rate parameterization
      draws(i, k) = cur;
      norm += cur;
    }
    // Normalize draw over K components
    for (int k = 0; k < K; k++) {
      draws(i, k) = draws(i, k) / norm;
    }
  }

  return draws;
}
