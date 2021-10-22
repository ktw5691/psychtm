#include <RcppArmadillo.h>

//' @title Estimate \eqn{\theta_d}
//'
//' @name est_thetad
//' @param z_count A K x 1 vector of counts of topic draw in document d.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//'
//' @return A K x 1 vector of estimate for \eqn{\theta_d}.
//'
//' @noRd
// [[Rcpp::export(.est_thetad)]]
arma::rowvec est_thetad(const arma::rowvec& z_count, float alpha_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 0.0;
  const uint16_t K = z_count.size(); // Number of topics
  if (alpha_ < 0.0) Rcpp::stop("alpha_ must be positive");
  if (K < 2) Rcpp::stop("number of topics must be at least 2");

  arma::rowvec thetad = exp(log(z_count + alpha_) - log(sum(z_count) +
    static_cast<float>(K) * alpha_));

  // Check for impossible estimates and replace with 0
  thetad.transform( [](double val) {
    if ((val > 1.0) | (std::isnan(val))) {
      return(0.0);
    } else {
      return val;
    }
  });

  return thetad;
}
