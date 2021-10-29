#include "draw_eta.h"
#include "draw_sigma2.h"
#include "get_loglike.h"
#include "get_logpost.h"
#include "post_pred.h"
#include "waic.h"

// These header files need to be included AFTER draw_zdn.h (which then includes RcppArmadilloExtensions/sample.h)
#include <progress.hpp>
#include <progress_bar.hpp>

//' @title Collapsed Gibbs sampler for multiple linear regression
//'
//' @name gibbs_mlr_cpp
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param thin The period of iterations to keep after the burn-in period
//'   (default: `1`).
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param x A D x (p + 1) matrix of additional predictors.
//' @param mu0 A (p + 1) x 1 mean vector for the prior on the regression
//'   coefficients.
//' @param sigma0 A (p + 1) x (p + 1) variance-covariance matrix for the
//'   prior on the regression coefficients.
//' @param eta_start A (p + 1) x 1 vector of starting values for the
//'   regression coefficients.
//' @param a0 The shape parameter for the prior on sigma2 (default: `0.001`)
//' @param b0 The scale parameter for the prior on sigma2 (default: `0.001`)
//' @param verbose Should parameter draws be output during sampling? (default:
//'   `false`).
//' @param display_progress Show progress bar? (default: `false`). Do not use
//'   with `verbose = true`.
//'
//' @return An object of class `Mlr`.
//'
//' @noRd
// [[Rcpp::export(.gibbs_mlr_cpp)]]
Rcpp::S4 gibbs_mlr_cpp(uint32_t m, uint32_t burn, uint32_t thin,
                       const arma::colvec& y, const arma::mat& x,
                       const arma::colvec& mu0, const arma::mat& sigma0,
                       arma::colvec eta_start, float a0 = 0.001, float b0 = 0.001,
                       bool verbose = false, bool display_progress = false) {

  if (m <= burn) Rcpp::stop("Length of chain m not greater than burn-in period.");
  if ( (thin == 0) || (thin > (m - burn)) )
    Rcpp::stop("Thinning period thin must be at least 1 and must not exceed the length of chain after burn-in period (m - burn).");

  Rcpp::S4 results("Mlr"); // Create object slda of class Mlr

  const uint32_t D = x.n_rows;
  const uint16_t pp1 = x.n_cols;
  const uint32_t chain_outlength = (m - burn) / thin; // Truncates remainder

  arma::mat etam(chain_outlength, pp1);
  Rcpp::NumericVector sigma2m(chain_outlength);
  Rcpp::NumericVector loglike(chain_outlength); // Store log-likelihood (up to an additive constant)
  Rcpp::NumericVector logpost(chain_outlength); // Store log-posterior (up to an additive constant)
  arma::mat l_pred(chain_outlength, D);

  // Initialize sigma^2
  double sigma2 = var(y) / 2.0;

  if (verbose) {
    Rcpp::Rcout << 1 << " eta: " << eta_start.t() << " ~~~~ sigma2: " << sigma2 << "\n";
  }

  arma::colvec eta(pp1);

  Progress p(m, display_progress);
  for (uint32_t i = 1; i <= m; i++) {

    // Draw eta
    try {
      eta = draw_eta_norm(x, y, sigma2, mu0, sigma0);
    } catch (std::exception& e) {
      Rcpp::Rcerr << "Runtime Error: " << e.what() << " while drawing eta vector\n";
      forward_exception_to_r(e);
    }

    // Draw sigma2
    try {
      sigma2 = draw_sigma2(a0, b0, x, y, eta);
    } catch (std::exception& e) {
      Rcpp::Rcerr << "Runtime Error: " << e.what() << " while drawing sigma2\n";
      forward_exception_to_r(e);
    }

    if ( (i > burn) && (i % thin == 0) ) {
      // Likelihood
      uint32_t temp_pos = (i - burn - 1) / thin;
      if (temp_pos == chain_outlength) break; // Invalid index
      loglike(temp_pos) = get_ll_mlr(y, x, eta, sigma2);
      // Log-posterior
      logpost(temp_pos) = get_lpost_mlr(loglike(temp_pos), eta, sigma2,
              mu0, sigma0, a0, b0);

      l_pred.row(temp_pos) = post_pred_norm(y, x, eta, sigma2);

      etam.row(temp_pos) = eta.t();
      sigma2m(temp_pos) = sigma2;
    }
    if (i % 500 == 0) {
      if (verbose) {
        Rcpp::Rcout << i << " eta: " << eta.t() << " ~~~~ sigma2: " << sigma2 << "\n";
      }
    }
    if (display_progress) {
      // Check to see if user cancelled sampler
      if (Progress::check_abort()) {
        Rcpp::stop("");
      } else {
        p.increment();
      }
    }
  }

  // Compute WAIC and p_eff
  Rcpp::NumericVector waic_and_se = waic_all(chain_outlength, l_pred);

  results.slot("ndocs") = D;
  results.slot("nchain") = chain_outlength;
  results.slot("eta") = etam;
  results.slot("sigma2") = sigma2m;
  results.slot("mu0") = mu0;
  results.slot("sigma0") = sigma0;
  results.slot("a0") = a0;
  results.slot("b0") = b0;
  results.slot("eta_start") = eta_start;
  results.slot("loglike") = loglike;
  results.slot("logpost") = logpost;
  results.slot("waic") = waic_and_se(0);
  results.slot("se_waic") = waic_and_se(1);
  results.slot("p_eff") = waic_and_se(2);
  results.slot("lpd") = l_pred;

  return results;
}
