#include <RcppArmadillo.h>

//' @title Sample from multivariate Gaussian N(\eqn{\mu}, \eqn{\Sigma})
//'
//' @name rmvnorm_cpp
//' @param n The number of samples to draw.
//' @param mu The q x 1 mean vector of the distribution (column vector).
//' @param sigma The q x q variance-covariance matrix of the distribution.
//'
//' @return A n x q matrix of random draws.
//'
//' @noRd
arma::mat rmvnorm_cpp(uint32_t n, const arma::colvec& mu,
                      const arma::mat& sigma) {

  // TODO: Check for positive-definite sigma
  const uint16_t ncols = sigma.n_cols;
  if ((mu.size() != sigma.n_rows) || (mu.size() != ncols)) {
    Rcpp::Rcerr <<
      "sigma must be a square matrix and mu must be a column vector with length equal to the number of rows and columns in sigma\n";
  }
  arma::mat y = arma::randn(n, ncols);

  return arma::repmat(mu, 1, n).t() + y * arma::chol(sigma);
}
