#include "count.h"
#include "count_topic_word.h"
#include "draw_eta.h"
#include "draw_sigma2.h"
#include "draw_zdn.h"
#include "draw_beta.h"
#include "draw_theta.h"
#include "est_betak.h"
#include "est_thetad.h"
#include "get_loglike.h"
#include "get_logpost.h"
#include "post_pred.h"
#include "waic.h"

// These header files need to be included AFTER draw_zdn.h (which then includes RcppArmadilloExtensions/sample.h)
#include <progress.hpp>
#include <progress_bar.hpp>

//' @title Collapsed Gibbs sampler for the SLDAX model
//'
//' @name gibbs_sldax_cpp
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param thin The period of iterations to keep after the burn-in period
//'   (default: `1`).
//' @param y A D x 1 vector of binary outcomes (0/1) to be predicted.
//' @param x A D x p matrix of additional predictors (no column of 1s for
//'   intercept).
//' @param docs A D x max(\eqn{N_d}) matrix of word indices for all documents.
//' @param V The number of unique terms in the vocabulary.
//' @param K The number of topics.
//' @param model An integer denoting the type of model to fit.
//' @param mu0 A K x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A (K + p + 1) x (K + p + 1) variance-covariance matrix for the
//'   prior on the regression coefficients. The first p + 1 columns/rows
//'   correspond to predictors in X, while the last K columns/rows correspond to
//'   the K topic means.
//' @param a0 The shape parameter for the prior on sigma2.
//' @param b0 The scale parameter for the prior on sigma2.
//' @param eta_start A (K + p) x 1 vector of starting values for the
//'   regression coefficients. The first p elements correspond to predictors
//'   in X, while the last K elements correspond to the K topic means.
//' @param constrain_eta A logical (default = `false`): If `true`, the
//'   regression coefficients will be constrained so that they are in descending
//'   order; if `false`, no constraints will be applied.
//' @param sample_beta A logical (default = `true`): If `true`, the
//'   topic-vocabulary distributions are sampled from their full conditional
//'   distribution.
//' @param sample_theta A logical (default = `true`): If `true`, the
//'   topic proportions are sampled from their full conditional distribution.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions
//'   (default: `1.0`).
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities (default: `1.0`).
//' @param proposal_sd The proposal standard deviation for drawing the
//'   regression coefficients, N(0, `proposal_sd`) (default: `0.2`).
//' @param interaction_xcol The column number of the design matrix for the
//' additional predictors for which an interaction with the \eqn{K} topics is
//' desired (default: `-1L`, no interaction). Currently only supports a
//' single continuous predictor or a two-category categorical predictor
//' represented as a single dummy-coded column.
//' @param return_assignments A logical (default = `false`): If
//'   `true`, returns an N x \eqn{max N_d} x M array of topic assignments
//'   in slot `@topics`. CAUTION: this can be memory-intensive.
//' @param verbose Should parameter draws be output during sampling? (default:
//'   `false`).
//' @param display_progress Show progress bar? (default: `false`). Do not use
//'   with `verbose = true`.
//'
//' @noRd
// [[Rcpp::export(.gibbs_sldax_cpp)]]
Rcpp::S4 gibbs_sldax_cpp(const arma::umat& docs, uint32_t V,
                   uint32_t m, uint32_t burn, uint32_t thin,
                   uint16_t K, uint8_t model,
                   const arma::colvec& y,
                   const arma::mat& x,
                   const arma::colvec& mu0,
                   const arma::mat& sigma0,
                   float a0, float b0,
                   arma::colvec eta_start,
                   arma::vec proposal_sd,
                   int interaction_xcol = -1,
                   float alpha_ = 1.0, float gamma_ = 1.0,
                   bool constrain_eta = false,
                   bool sample_beta = true, bool sample_theta = true,
                   bool return_assignments = false,
                   bool verbose = false, bool display_progress = false) {

  if (m <= burn) Rcpp::stop("Length of chain m not greater than burn-in period.");
  if ( (thin == 0) || (thin > (m - burn)) )
    Rcpp::stop("Thinning period thin must be at least 1 and must not exceed the length of chain after burn-in period (m - burn).");

  const uint8_t lda = 1;
  const uint8_t slda = 2;
  const uint8_t sldax = 3;
  const uint8_t slda_logit = 4;
  const uint8_t sldax_logit = 5;

  Rcpp::S4 results("Sldax");

  const uint32_t chain_outlength = (m - burn) / thin; // Truncates remainder
  Rcpp::NumericVector sigma2m(chain_outlength);
  double sigma2 = var(y) / 2.0;

  const uint32_t D = docs.n_rows;
  arma::colvec N(D);
  const Rcpp::IntegerVector topics_index = Rcpp::seq_len(K);
  const Rcpp::IntegerVector docs_index   = Rcpp::seq_len(D) - 1;
  for (uint32_t d : docs_index) N(d) = sum(docs.row(d) > 0);
  const uint32_t maxNd = max(N);
  arma::mat theta(D, K);
  arma::mat beta(K, V);
  // Topic draw counts
  arma::mat ndk(D, K);
  // Counts of topic-word co-occurences in corpus (K x V)
  arma::mat nkm(K, V);
  // Topic draws for all words and docs
  arma::vec nk(K);
  arma::umat zdocs = arma::umat(D, maxNd, arma::fill::zeros);
  arma::mat zbar(D, K);

  // D x max(N_d) x m array to store topic draws
  uint32_t ntopics_slices;
  if (return_assignments) {
    ntopics_slices = chain_outlength;
  } else {
    ntopics_slices = 1;
  }
  arma::ucube topicsm = arma::ucube(D, maxNd, ntopics_slices, arma::fill::zeros);

  // D x K x m array to store topic proportion draws
  uint32_t ntheta_slices;
  if (sample_theta) {
    ntheta_slices = chain_outlength;
  } else {
    ntheta_slices = 1;
  }
  arma::cube thetam = arma::cube(D, K, ntheta_slices, arma::fill::zeros);

  // K x V x m array to store topic-vocabulary distribution draws
  uint32_t nbeta_slices;
  if (sample_beta) {
    nbeta_slices = chain_outlength;
  } else {
    nbeta_slices = 1;
  }
  arma::cube betam = arma::cube(K, V, nbeta_slices, arma::fill::zeros);

  Rcpp::NumericVector loglike(chain_outlength); // Store log-likelihood (up to an additive constant)
  Rcpp::NumericVector logpost(chain_outlength); // Store log-posterior (up to an additive constant)

  for (uint32_t d : docs_index) {
    // Randomly initialize topic draws from Uniform(1, K)
    zdocs.row(d).subvec(0, N(d) - 1) = arma::randi<arma::urowvec>( N(d), arma::distr_param(1, K) );
    for (uint16_t k = 0; k < K; k++) {
      // Count topic draws in each document
      ndk(d, k) = sum(zdocs.row(d) == (k + 1));
    }
  }

  if (model != lda) {
    // Get zbar matrix
    zbar = ndk.each_col() / N;
  }

  try {
    nkm = count_topic_word(K, V, zdocs, docs);
  } catch(std::exception& e) {
    Rcpp::Rcerr << "Runtime error: " << e.what() <<
      " while computing topic-word co-occurrences\n";
    forward_exception_to_r(e);
  }

  nk = arma::sum(ndk, 0).t(); // Column sums

  if (sample_theta) {
    try {
      theta = draw_theta(ndk, alpha_);
    } catch(std::exception& e) {
      Rcpp::Rcerr << "Runtime Error: " << e.what() << " when sampling Theta\n";
      forward_exception_to_r(e);
    }
  } else {
    // Estimate theta
    for (uint32_t d : docs_index) {
      try {
        theta.row(d) = est_thetad(ndk.row(d), alpha_);
      } catch(std::exception& e) {
        Rcpp::Rcerr << "Runtime Error: " << e.what() <<
          " when estimating theta vector for document " << d << "\n";
        forward_exception_to_r(e);
      }
    }
  }

  if (sample_beta) {
    try {
      beta = draw_beta(nkm, gamma_);
    } catch(std::exception& e) {
      Rcpp::Rcerr << "Runtime Error: " << e.what() << " when sampling Beta\n";
      forward_exception_to_r(e);
    }
  } else {
    // Estimate beta
    for (uint16_t k = 0; k < K; k++) {
      try {
        beta.row(k) = est_betak(nkm.row(k), gamma_);
      } catch(std::exception& e) {
        Rcpp::Rcerr << "Runtime Error: " << e.what() << " estimating row " << k <<
          "of beta matrix\n";
        forward_exception_to_r(e);
      }
    }
  }

  // For supervised models, set number of columns in design matrix
  uint16_t p = 0; // Num. cols. involving fixed predictors
  uint16_t q = 0; // Total num. cols. in design matrix
  if (model == sldax || model == sldax_logit) {
    if (interaction_xcol > 0) {
      p = x.n_cols + K - 1; // Need K - 1 interaction columns
    } else {
      p = x.n_cols;
    }
    q = p + K; // num. covariates + K + (K - 1) for interaction of one covariate w/ topics
  } else if (model == slda || model == slda_logit) {
    p = K; // Num. cols. involving fixed predictors
    q = p; // Total num. cols. in design matrix
  }
  arma::mat etam(chain_outlength, q);
  arma::colvec eta(q);
  arma::colvec etac(q);

  // For supervised models SLDA/SLDAX, set up regression coefficients
  if (constrain_eta) {
    if (model == sldax || model == sldax_logit) {
      // Constrain starting values of eta to be in descending order
      if (interaction_xcol > 0) {
        // Only sort topic "linear effects", not interactions
        std::partial_sort(eta_start.begin() + p - K + 1,
                          eta_start.end() - K + 1,
                          eta_start.end(), std::greater<float>());
      } else {
        std::sort(eta_start.begin() + p,
                  eta_start.end(),
                  std::greater<float>());
      }
    } else if (model == slda || model == slda_logit) {
      std::sort(eta_start.begin(), eta_start.end(), std::greater<float>());
    }
  }

  if (verbose && model != lda) Rcpp::Rcout << 1 << " eta: " << eta_start.t() << "\n";
  eta = eta_start;

  arma::vec attempt = arma::vec(q, arma::fill::zeros);
  arma::vec accept  = arma::vec(q, arma::fill::zeros);
  arma::vec acc_rate(q);

  arma::mat r = zbar;
  // Predictive posterior likelihood of y
  arma::mat l_pred(chain_outlength, D);
  arma::mat xz_int(D, K - 1); // Need K - 1 interaction columns
  if (model == sldax || model == sldax_logit) {
    r = join_rows(x, zbar); // Full predictor matrix `r`
    if (interaction_xcol > 0) {
      xz_int = zbar.cols(0, K - 2).each_col() % x.col(interaction_xcol - 1);
      r = join_rows(r, xz_int);
    }
  }

  Progress prog(m, display_progress);
  for (uint32_t i = 1; i <= m; i++) {

    // Draw z
    for (uint32_t d : docs_index) {
      for (uint32_t n = 0; n < N(d); n++) {
        uint32_t word = docs(d, n);
        uint16_t topic = zdocs(d, n);
        // Draw topic z_{dn} excluding current z_{dn} value
        update_zcounts(d, word, topic, d, ndk, nk, nkm);
        try {
          if (model == lda) {
            topic = draw_zdn_lda(V, ndk.row(d).t(), nkm.col(word - 1), nk,
                                 alpha_, gamma_);
          } else if (model == slda || model == sldax) {
            topic = draw_zdn_slda_norm(y(d), r.row(d), eta, sigma2, V,
                                       ndk.row(d).t(), nkm.col(word - 1), nk,
                                       alpha_, gamma_);
          } else if (model == slda_logit || model == sldax_logit) {
            topic = draw_zdn_slda_logit(y(d), r.row(d), eta, V,
                                        ndk.row(d).t(), nkm.col(word - 1), nk,
                                        alpha_, gamma_);
          } else {
            Rcpp::Rcerr << "Invalid model specified. Exiting.\n";
          }
        } catch(std::exception& e) {
          Rcpp::Rcerr << "Runtime Error: " << e.what() <<
            " occurred while drawing topic for word " << n << " in document "
            << d << "\n";
          forward_exception_to_r(e);
        }

        zdocs(d, n) = topic;
        // Update topic count in doc d
        ndk(d, topic - 1)++;
        // Update topic count in corpus
        nk(topic - 1)++;
        // Update topic-word counts in corpus
        nkm(topic - 1, word - 1)++;
      }
    }

    if (sample_theta) {
      // Sample theta
      try {
        theta = draw_theta(ndk, alpha_);
      } catch(std::exception& e) {
        Rcpp::Rcerr << "Runtime Error: " << e.what() << " when sampling Theta\n";
        forward_exception_to_r(e);
      }
    } else {
      // Estimate theta
      for (uint32_t d : docs_index) {
        try {
          theta.row(d) = est_thetad(ndk.row(d), alpha_);
        } catch(std::exception& e) {
          Rcpp::Rcerr << "Runtime Error: " << e.what() <<
            " when estimating theta vector for document " << d << "\n";
          forward_exception_to_r(e);
        }
      }
    }

    if (model != lda) {
      // Get zbar matrix
      zbar = ndk.each_col() / N;
    }

    if (model == sldax || model == sldax_logit) {
      r = join_rows(x, zbar); // Full predictor matrix `r`
      if (interaction_xcol > 0) {
        xz_int = zbar.cols(0, K - 2).each_col() % x.col(interaction_xcol - 1);
        r = join_rows(r, xz_int);
      }
    } else {
      r = zbar;
    }

    if (sample_beta) {
      try {
        beta = draw_beta(nkm, gamma_);
      } catch(std::exception& e) {
        Rcpp::Rcerr << "Runtime Error: " << e.what() << " when sampling Beta\n";
        forward_exception_to_r(e);
      }
    } else {
      // Estimate beta
      for (uint16_t k = 0; k < K; k++) {
        try {
          beta.row(k) = est_betak(nkm.row(k), gamma_);
        } catch(std::exception& e) {
          Rcpp::Rcerr << "Runtime Error: " << e.what() << " estimating row " << k <<
            "of beta matrix\n";
          forward_exception_to_r(e);
        }
      }
    }

    // Draw eta
    if (model != lda) {
      bool eta_order = false;
      uint16_t iter = 0;
      uint16_t max_iter = 1; // Only draw once
      if (constrain_eta) {
        max_iter = 10000;   // Max repeated draws if trying to constrain eta
      }
      while (!eta_order & (iter < max_iter)) { // Runs at least once
        iter++;
        if (model == slda || model == sldax) {
          try {
            etac = draw_eta_norm(r, y, sigma2, mu0, sigma0);
          } catch (std::exception& e) {
            Rcpp::Rcerr << "Runtime Error: " << e.what() <<
              " while drawing eta vector\n";
            forward_exception_to_r(e);
          }
        } else if (model == slda_logit || model == sldax_logit) {
          try {
            etac = draw_eta_logit(r, y, eta, mu0, sigma0,
                                  proposal_sd, attempt, accept);
          } catch (std::exception& e) {
            Rcpp::Rcerr << "Runtime Error: " << e.what() <<
              " while drawing eta vector\n";
            forward_exception_to_r(e);
          }
        }
        if (constrain_eta) {
          if (interaction_xcol > 0) {
            for (uint16_t k = p - K + 2; k < p + 1; k++) {
              // Force eta components to be in descending order (first is largest) to resolve label switching of topics
              eta_order = (etac(k - 1) > etac(k)) || (fabs(etac(k - 1) - etac(k)) < .000001);
              if (!eta_order) break;
            }
          } else {
            for (uint16_t k = p + 1; k < q; k++) {
              // Force eta components to be in descending order (first is largest)
              //   to resolve label switching of topics
              eta_order = (etac(k - 1) > etac(k)) || (fabs(etac(k - 1) - etac(k)) < .000001);
              if (!eta_order) break;
            }
          }
        }
      }
      if ( (iter >= max_iter - 1) & (max_iter > 1) ) Rcpp::Rcout << "Reached max_iter while drawing eta\n";
      eta = etac; // New draw

      if (model == slda || model == sldax) {
        // Draw sigma2
        try {
          sigma2 = draw_sigma2(a0, b0, r, y, eta);
        } catch (std::exception& e) {
          Rcpp::Rcerr << "Runtime Error: " << e.what() << " while drawing sigma2\n";
          forward_exception_to_r(e);
        }
      }
    }

    if ( (i > burn) && (i % thin == 0) ) {
      if (verbose && i == burn + 1 && burn > 0) Rcpp::Rcout << "Finished burn-in period\n";

      uint32_t temp_pos = (i - burn - 1) / thin;
      if (temp_pos == chain_outlength) break; // Invalid index

      if (return_assignments) topicsm.slice(temp_pos) = zdocs;
      if (sample_theta) thetam.slice(temp_pos) = theta;
      if (sample_beta) betam.slice(temp_pos) = beta;

      // Likelihood
      switch(model) {
        case lda:
          loglike(temp_pos) = get_ll_lda(zdocs, docs, theta, beta,
                                             docs_index, N);
          logpost(temp_pos) = get_lpost_lda(
            loglike(temp_pos), theta, beta, gamma_, alpha_);
          break;
        case slda: // Fall through
        case sldax:
          loglike(temp_pos) = get_ll_slda_norm(
            y, r, eta, sigma2, zdocs, docs, theta, beta, docs_index, N);
          logpost(temp_pos) = get_lpost_slda_norm(
            loglike(temp_pos), eta, sigma2, theta, beta, mu0, sigma0, gamma_,
            alpha_, a0, b0);
          l_pred.row(temp_pos) = post_pred_norm(y, r, eta, sigma2);
          sigma2m(temp_pos) = sigma2;
          break;
        case slda_logit: // Fall through
        case sldax_logit:
          loglike(temp_pos) = get_ll_slda_logit(
            y, r, eta, zdocs, docs, theta, beta, docs_index, N);
          logpost(temp_pos) = get_lpost_slda_logit(loglike(temp_pos),
            eta, theta, beta, mu0, sigma0, gamma_, alpha_);
          l_pred.row(temp_pos) = post_pred_logit(y, r, eta);
          break;
        default: Rcpp::Rcerr << "Invalid model specified. Exiting.\n";
      }
      if (model != lda) {
        etam.row(temp_pos) = eta.t();
      }
    }

    if (model == slda_logit || model == sldax_logit) {
      acc_rate = accept / attempt;
      for (uint16_t j = 0; j < q; j++) {
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
      if (verbose && i == burn / 5.0 && burn > 0) Rcpp::Rcout << "Finished tuning proposal distributions\n";
    }

    if (i % 500 == 0 && verbose) {
      switch(model) {
        case slda: // Fall through
        case sldax:
          Rcpp::Rcout << i << " eta: " << eta.t() <<
            " ~~~~ sigma2: " << sigma2 <<
            " ~~~~ zbar d1: " << zbar.row(0) << "\n";
          break;
        case slda_logit: // Fall through
        case sldax_logit:
          Rcpp::Rcout << i << " eta: " << eta.t() <<
            " ~~~~ zbar d1: " << zbar.row(0) <<
            " ~~~~ AR: " << acc_rate.t() << "\n" <<
            " ~~~~ Proposal SD: " << proposal_sd.t() << "\n";
          break;
        default:
          break;
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

  results.slot("ntopics") = K;
  results.slot("ndocs")   = D;
  results.slot("nvocab")  = V;
  results.slot("nchain")  = chain_outlength;
  if (return_assignments) {
    results.slot("topics") = topicsm; // Full chain
  } else {
    results.slot("topics") = zdocs; // Last draw
  }
  results.slot("alpha")   = alpha_;
  results.slot("gamma")   = gamma_;
  results.slot("loglike") = loglike;
  results.slot("logpost") = logpost;
  if (sample_theta) {
    results.slot("theta") = thetam; // Full chain of draws or estimates
  } else {
    results.slot("theta") = theta; // Last draw or estimate
  }
  if (sample_beta) {
    results.slot("beta") = betam; // Full chain of draws or estimates
  } else {
    results.slot("beta") = beta; // Last draw or estimate
  }

  if (model != lda) {
    results.slot("eta")       = etam;
    results.slot("mu0")       = mu0;
    results.slot("sigma0")    = sigma0;
    results.slot("eta_start") = eta_start;
    results.slot("waic")      = waic_and_se(0);
    results.slot("se_waic")   = waic_and_se(1);
    results.slot("p_eff")     = waic_and_se(2);
    results.slot("lpd")       = l_pred;

    if (model == slda || model == sldax) {
      results.slot("sigma2") = sigma2m;
      results.slot("a0")     = a0;
      results.slot("b0")     = b0;
    }

    if (model == slda_logit || model == sldax_logit)
      results.slot("proposal_sd") = proposal_sd;
  }

  return results;
}
