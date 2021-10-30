#include <RcppArmadillo.h>

//' @title Contribution to effective number of parameters for WAIC from observation y_d
//'
//' @name pwaic_d
//' @param like_pred A m x 1 vector of predictive likelihoods (NOT log-likelihoods).
//' @return The contribution of y_d (its predictive posterior likelihood variance)
//'   to the effective number of parameters.
//'
//' @noRd
double pwaic_d(const arma::colvec& like_pred) {

  // Get variance of log-predictive likelihood
  double mcmc_vary = arma::var(log(like_pred)); // Denominator of m - 1
  return mcmc_vary;
}

//' @title WAIC for observation y_d
//'
//' @name waic_d
//' @param like_pred A m x 1 vector of predictive likelihoods (NOT log-likelihoods) for y_d.
//' @param p_effd The contribution to the effective number of parameters from
//'   obs y_d.
//' @return WAIC contribution for observation d (on deviance scale).
double waic_d(const arma::colvec& like_pred, double p_effd) {

  // Log of mean posterior predictive density for y_d
  double lppd = log(arma::mean(like_pred));

  // Contribution to WAIC (on deviance scale, i.e., -2ll) for y_d
  return (-2.0 * (lppd - p_effd)); // See Gelman, Hwang, Vehtari (2014, p. 1003)
}

//' @title Compute WAIC for all outcomes.
//'
//' @name waic_all
//' @param iter The length of the sampled chain.
//' @param l_pred A `iter` x D matrix of predictive likelihoods (NOT log-likelihoods).
//'
//' @return Vector of (1) WAIC for model, (2) standard error for WAIC, and (3)
//'   the effective number of parameters.
//'
//' @examples
//' data(teacher_rate)
//' fit_mlr <- gibbs_mlr(rating ~ grade, data = teacher_rate, m = 5)
//' waic_all(iter = 5, t(lpd(fit_mlr)))
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector waic_all(uint32_t iter, const arma::mat& l_pred) {

  Rcpp::NumericVector full_waic_se(3);
  const uint16_t D = l_pred.n_cols; // Number of observations
  arma::colvec waic(D);
  arma::colvec peff(D);

  // Compute WAIC and p_eff for entire data set
  for (uint16_t d = 0; d < D; d++) {
    peff(d) = pwaic_d(l_pred.submat(0, d, iter - 1, d));
    waic(d) = waic_d(l_pred.submat(0, d, iter - 1, d), peff(d));
  }
  double peff_sum = arma::sum(peff);  // p_eff
  double waic_sum = arma::sum(waic);  // WAIC

  // Compute SE(WAIC)
  double se_waic = arma::var(waic);   // Denominator of D - 1
  se_waic *= static_cast<float>(D);   // Variance of WAIC times D
  se_waic = sqrt(se_waic);            // SE(WAIC)

  full_waic_se(0) = waic_sum;
  full_waic_se(1) = se_waic;
  full_waic_se(2) = peff_sum;

  return full_waic_se;
}

//' @title Compute difference (WAIC1 - WAIC2) in WAIC and its SE for two models.
//'
//' @name waic_diff
//' @param l_pred1 A m1 x D matrix of predictive likelihoods (NOT log-likelihoods) from model 1.
//' @param l_pred2 A m2 x D matrix of predictive likelihoods (NOT log-likelihoods) from model 2.
//'
//' @return A vector of (1) the difference in WAIC (on the deviance scale)
//'   between models and (2) the standard error of the difference in WAIC.
//'
//' @examples
//' data(teacher_rate)
//' fit_mlr <- gibbs_mlr(rating ~ grade, data = teacher_rate, m = 100)
//' fit_mlr2 <- gibbs_mlr(rating ~ grade + I(grade^2), data = teacher_rate, m = 100)
//' # Returns (1) D = WAIC(fit_mlr2) - WAIC(fit_mlr) and (2) SE(D)
//' #   Suggests that a linear relationship is preferable
//' waic_diff(t(lpd(fit_mlr2)), t(lpd(fit_mlr)))
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector waic_diff(const arma::mat& l_pred1, const arma::mat& l_pred2) {

  Rcpp::NumericVector diff_waic_se(2);
  const uint32_t m1 = l_pred1.n_rows; // Length of first chain
  const uint32_t m2 = l_pred2.n_rows; // Length of second chain
  const uint16_t D = l_pred1.n_cols;  // Number of observations
  arma::colvec waic1(D);
  arma::colvec waic2(D);
  arma::colvec peff1(D);
  arma::colvec peff2(D);
  arma::colvec waic_diff(D);

  // Compute WAIC and p_eff for entire data set
  for (uint16_t d = 0; d < D; d++) {
    peff1(d) = pwaic_d(l_pred1.submat(0, d, m1 - 1, d));
    peff2(d) = pwaic_d(l_pred2.submat(0, d, m2 - 1, d));
    waic1(d) = waic_d(l_pred1.submat(0, d, m1 - 1, d), peff1(d));
    waic2(d) = waic_d(l_pred2.submat(0, d, m2 - 1, d), peff2(d));
  }
  waic_diff = waic1 - waic2;
  double waic_diff_sum = arma::sum(waic_diff); // WAIC1 - WAIC2

  // Compute SE(WAIC1 - WAIC2)
  double se_waic_diff = arma::var(waic_diff);  // Denominator of D - 1
  se_waic_diff *= static_cast<float>(D);       // Multiply variance of diff by D
  se_waic_diff = sqrt(se_waic_diff);           // SE(WAIC1 - WAIC2)

  diff_waic_se(0) = waic_diff_sum; // Difference
  diff_waic_se(1) = se_waic_diff; // SE(WAIC1 - WAIC2)

  return diff_waic_se;
}
