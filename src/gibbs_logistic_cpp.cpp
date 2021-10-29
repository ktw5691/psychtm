#include "draw_eta.h"
#include "get_loglike.h"
#include "get_logpost.h"
#include "post_pred.h"
#include "waic.h"

// These header files need to be included AFTER draw_zdn.h (which then includes RcppArmadilloExtensions/sample.h)
#include <progress.hpp>
#include <progress_bar.hpp>


//' @title Collapsed Gibbs sampler for logistic regression
//'
//' @name gibbs_logistic_cpp
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param thin The period of iterations to keep after the burn-in period
//'   (default: `1`).
//' @param y A D x 1 vector of binary outcomes (0/1) to be predicted.
//' @param x A D x p matrix of additional predictors (no column of 1s for
//'   intercept).
//' @param mu0 A (p + 1) x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A (p + 1) x (p + 1) variance-covariance matrix for the prior
//'   on the regression coefficients.
//' @param eta_start A (p + 1) x 1 vector of starting values for the
//'   regression coefficients.
//' @param proposal_sd The proposal standard deviation for drawing the
//'   regression coefficients, N(0, `proposal_sd`) (default: `2.38, ..., 2.38`).
//' @param verbose Should parameter draws be output during sampling? (default:
//'   `false`).
//' @param display_progress Show progress bar? (default: `false`). Do not use
//'   with `verbose = true`.
//'
//' @return An object of class Logistic.
//'
//' @noRd
// [[Rcpp::export(.gibbs_logistic_cpp)]]
Rcpp::S4 gibbs_logistic_cpp(uint32_t m, uint32_t burn, uint32_t thin,
                            const arma::colvec& y, const arma::mat& x,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            arma::colvec eta_start, arma::vec proposal_sd,
                            bool verbose = false, bool display_progress = false) {

  if (m <= burn) Rcpp::stop("Length of chain m not greater than burn-in period.");
  if ( (thin == 0) || (thin > (m - burn)) )
    Rcpp::stop("Thinning period thin must be at least 1 and must not exceed the length of chain after burn-in period (m - burn).");

  Rcpp::S4 results("Logistic"); // Create object slda of class Logistic

  const uint32_t D = y.size();
  const uint16_t pp1 = x.n_cols;
  const uint32_t chain_outlength = (m - burn) / thin; // Truncates remainder

  arma::mat etam = arma::mat(chain_outlength, pp1, arma::fill::zeros);
  Rcpp::NumericVector loglike(chain_outlength); // Store log-likelihood (up to an additive constant)
  Rcpp::NumericVector logpost(chain_outlength); // Store log-posterior (up to an additive constant)
  arma::mat l_pred(chain_outlength, D);

  arma::vec attempt = arma::vec(pp1, arma::fill::zeros);
  arma::vec accept  = arma::vec(pp1, arma::fill::zeros);
  arma::vec acc_rate(pp1);

  arma::colvec eta = eta_start;

  Progress prog(m, display_progress);
  for (uint32_t i = 1; i <= m; i++) {

    // Draw eta
    try {
      eta = draw_eta_logit(x, y, eta,
                           mu0, sigma0, proposal_sd,
                           attempt, accept);
    } catch (std::exception& e) {
      Rcpp::Rcerr << "Runtime Error: " << e.what() <<
        " while drawing eta vector\n";
      forward_exception_to_r(e);
    }

    // Update acceptance rates and tune proposal standard deviations
    acc_rate = accept / attempt;
    for (uint16_t j = 0; j < pp1; j++) {
      if (i < (static_cast<float>(burn) / 5.0) &&
          attempt(j) >= 50 &&
          (i % 50 == 0)) {
        if (acc_rate(j) < 0.159) {
          // too low acceptance, decrease jumping width
          proposal_sd(j) = proposal_sd(j) * 0.8;
        }
        if (acc_rate(j) > 0.309) {
          // too high acceptance, increase jumping width
          proposal_sd(j) = proposal_sd(j) * 1.2;
        }
        accept(j) = attempt(j) = 0;
      }
    }

    if ( (i > burn) && (i % thin == 0) ) {
      uint32_t temp_pos = (i - burn - 1) / thin;
      if (temp_pos == chain_outlength) break; // Invalid index
      // Likelihood
      loglike(temp_pos) = get_ll_logit(y, x, eta);
      // Log-posterior
      logpost(temp_pos) = get_lpost_eta(loglike(temp_pos), eta, mu0, sigma0);
      l_pred.row(temp_pos) = post_pred_logit(y, x, eta);
      etam.row(temp_pos) = eta.t();
    }
    if (i % 500 == 0) {
      if (verbose) {
        Rcpp::Rcout << i << " eta: " << eta.t() << "\n" <<
          "AR: " << acc_rate.t() << "\n" <<
            "Proposal SD: " << proposal_sd.t() << "\n";
      }
    }
    if (display_progress) {
      // Check to see if user cancelled sampler
      if (Progress::check_abort()) {
        Rcpp::stop("");
      } else {
        prog.increment();
      }
    }
  }

  // Compute WAIC and p_eff
  Rcpp::NumericVector waic_and_se = waic_all(chain_outlength, l_pred);

  results.slot("ndocs") = D;
  results.slot("nchain") = chain_outlength;
  results.slot("eta") = etam;
  results.slot("mu0") = mu0;
  results.slot("sigma0") = sigma0;
  results.slot("eta_start") = eta_start;
  results.slot("proposal_sd") = proposal_sd;
  results.slot("loglike") = loglike;
  results.slot("logpost") = logpost;
  results.slot("waic") = waic_and_se(0);
  results.slot("se_waic") = waic_and_se(1);
  results.slot("p_eff") = waic_and_se(2);
  results.slot("lpd") = l_pred;

  return results;
}
