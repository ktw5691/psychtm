#include <RcppArmadillo.h>

//' @title Draw sigma2 from full conditional posterior for SLDA/SLDAX/MLR
//'
//' @name draw_sigma2
//' @param a0 The prior shape parameter for \eqn{\sigma^2}.
//' @param b0 The prior scale parameter for \eqn{\sigma^2}.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable.
//' @param eta A q x 1 vector of regression coefficients.
//'
//' @return A draw for \eqn{\sigma^2}.
//'
//' @noRd
long double draw_sigma2(float a0, float b0,
                        const arma::mat& w, const arma::colvec& y,
                        const arma::colvec& eta) {

  const uint32_t D = y.size(); // Number of observations
  if ((a0 <= 0.0) || (b0 <= 0.0)) Rcpp::stop("a0 and b0 must be positive");
  long double a = 0.5 * (static_cast<float>(D) + a0);

  if ((w.n_rows != D) || (w.n_cols != eta.size()))
    Rcpp::stop("zbar must be a D x q matrix and eta must be a q x 1 vector");

  arma::colvec resid(D);
  resid = y - w * eta;
  double b_update = arma::as_scalar(resid.t() * resid);

  long double scale = 0.5 * (b0 + b_update);
  // Parameterization of R::rgamma(a, b) is a = shape, b = 1 / scale
  long double b = 1.0 / scale;
  long double sigma2inv = R::rgamma(a, b);
  return 1.0 / sigma2inv;
}
