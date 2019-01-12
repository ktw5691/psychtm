#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <vector>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]

// Function to handle errors
void error(std::string s) {
  throw std::runtime_error(s);
}

//' Sample from multivariate Gaussian N(\eqn{\mu}, \eqn{\Sigma})
//'
//' @param n The number of samples to draw.
//' @param mu The mean vector of the distribution (column vector).
//' @param sigma The variance-covariance matrix of the distribution.
//' @export
// [[Rcpp::export]]
arma::mat rmvnorm_cpp(uint32_t n, const arma::colvec& mu,
                      const arma::mat& sigma) {

  // TODO: Check for positive-definite sigma
  if ((mu.size() != sigma.n_rows) || (mu.size() != sigma.n_cols)) {
    Rcerr <<
      "sigma must be a square matrix and mu must be a column vector with length equal to the number of rows and columns in sigma\n";
  }
  const uint16_t ncols = sigma.n_cols;
  arma::mat y = arma::randn(n, ncols);

  return arma::repmat(mu, 1, n).t() + y * arma::chol(sigma);
}

//' Draw eta from full conditional posterior
//'
//' @param zbar A D x K matrix with row d containing the mean number of draws of
//'   topics \eqn{z_1, \ldots, z_K} in document \eqn{d} where each row sums to
//'   1.
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param sigma2 The residual variance.
//' @param mu0 A K x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A K x K prior variance-covariance matrix for the regression
//'   coefficients.
arma::mat draw_eta_cpp(const arma::mat& zbar, const arma::vec& y,
                       long double sigma2, const arma::vec& mu0,
                       const arma::mat& sigma0) {

  const uint32_t D = zbar.n_rows;
  const uint16_t K = zbar.n_cols;
  if (K < 2) error("number of topics must be at least 2");
  if ((zbar.n_rows != D) || (zbar.n_cols != K))
    error("zbar must be a D x K matrix");
  if (y.size() != D) error("y must be of length D");
  if (mu0.size() != K) error("mu0 must be of length K");
  if ((sigma0.n_rows != K) || (sigma0.n_cols != K))
    error("sigma0 must be a K x K matrix");

  arma::mat ztz(K, K);
  ztz = zbar.t() * zbar;

  arma::colvec eta1(K);
  arma::mat sigma0_inv = sigma0.i();
  arma::mat sigma1 = (sigma0_inv + ztz / sigma2).i();
  eta1 = sigma1 * (sigma0_inv * mu0 + zbar.t() * y / sigma2);
  return rmvnorm_cpp(1, eta1, sigma1);
}

//' Draw sigma2 from full conditional posterior
//'
//' @param D The number of documents.
//' @param a0 The prior shape parameter for \eqn{\sigma^2}.
//' @param b0 The prior scale parameter for \eqn{\sigma^2}.
//' @param zbar A D x K matrix with row \eqn{d} containing the mean number of
//'   draws of topics \eqn{z_1, \ldots, z_K} in document \eqn{d} where each row
//'   sums to 1.
//' @param y A D x 1 vector of the outcome variable.
//' @param eta A K x 1 vector of regression coefficients.
//'
long double draw_sigma2_cpp(uint32_t D, float a0, float b0,
                            const arma::mat& zbar, const arma::colvec& y,
                            const arma::colvec& eta) {

  long double a = 0.5 * (static_cast<float>(D) + a0);

  if ((a0 < 0.0) || (b0 < 0.0)) error("a0 and b0 must be positive");

  if ((zbar.n_rows != D) || (zbar.n_cols != eta.size()))
    error("zbar must be a D x K matrix and eta must be a K x 1 vector");
  if (y.size() != D) error("y must be of length D");

  arma::colvec resid(D);
  resid = y - zbar * eta;
  double b_update = arma::as_scalar(resid.t() * resid);

  // Parameterization is 1 / rate
  long double b = 1.0 / (0.5 * (b0 + b_update));
  long double sigma2inv = R::rgamma(a, b);
  return 1.0 / sigma2inv;
}

//' Estimate beta_k (on log-scale)
//'
//' @param k The topic label (i.e., the row of the K x V beta matrix).
//' @param V The number of terms in the corpus vocabulary.
//' @param wz_co A V x 1 vector of counts of the draws of each word for topic
//'   \eqn{k} over all documents.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
arma::vec est_betak_cpp(uint16_t k, uint32_t V, const arma::vec& wz_co,
                        float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic word to 0.0;
  if (gamma_ < 0.0) error("gamma_ must be positive");
  if (V < 2) error("vocabulary size V must be at least 2");
  if (wz_co.size() != V) error("wz_co must be of length V");

  arma::vec betak = exp(log(wz_co + gamma_) - log(sum(wz_co) + V * gamma_));
  for (uint32_t v = 0; v < V; v++) {
    if ((betak[v] > 1.0) | (isnan(betak[v]))) betak[v] = 0.0;
  }
  return betak;
}

//' Estimate theta_d (on log scale)
//'
//' @param z_count A K x 1 vector of counts of topic draw in document \eqn{d}.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//' @param K The number of topics.
//'
arma::vec est_thetad_cpp(const arma::vec& z_count, float alpha_, uint16_t K) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 0.0;
  if (alpha_ < 0.0) error("alpha_ must be positive");
  if (K < 2) error("number of topics must be at least 2");
  if (z_count.size() != K) error("z_count must be of length K");

  arma::vec thetad = exp(log(z_count + alpha_) - log(sum(z_count) +
    static_cast<float>(K) * alpha_));
  for (uint32_t k = 0; k < K; k++) {
    if ((thetad[k] > 1.0) | (isnan(thetad[k]))) thetad[k] = 0.0;
  }

  return thetad;
}

//' Count topic-word co-occurences in corpus (ntopic x nvocab) (parallelizable)
//'
//' @param D The number of documents in the corpus.
//' @param K The number of topics.
//' @param V The number of terms in the corpus vocabulary.
//' @param doc_topic A D x max(\eqn{N_d}) matrix of topic assignments for
//'   the corpus.
//' @param doc_word A D x max(\eqn{N_d}) matrix of words for corpus.
//'
arma::mat count_topic_word_cpp(uint32_t D, uint16_t K, uint32_t V,
                               const arma::mat& doc_topic,
                               const arma::mat& doc_word) {

  if (K < 2) error("number of topics must be at least 2");
  if (V < 2) error("size of vocabulary V must be at least 2");
  if (doc_topic.n_rows != D) error("doc_topic must be a matrix with D rows");
  if (doc_word.n_rows != D) error("doc_word must be a matrix with D rows");
  if (doc_topic.n_cols != doc_word.n_cols)
    error("doc_topic and doc_word must have the same number of columns");

  const IntegerVector topics_index = seq_len(K);
  const IntegerVector docs_index = seq_len(D) - 1;
  // Matrix to store number of topic/word co-occurences
  arma::mat topic_word_freq = arma::zeros(K, V);
  // Loop through documents
  // for (uint32_t doc = 0; doc < D; doc++) {
  for (uint32_t doc : docs_index) {
    // Loop through words in document
    for (uint32_t pos = 0; pos < doc_word.n_cols; pos++) {
      // Loop through topics
      for (uint16_t topic = 1; topic <= K; topic++) {
        // Loop through vocabulary
        for (uint32_t v = 1; v <= V; v++) {
          topic_word_freq(topic - 1, v - 1) += ((doc_topic(doc, pos) == topic) *
            (doc_word(doc, pos) == v));
        }
      }
    }
  }
  return topic_word_freq;
}

//' Draw zdn from full conditional distribution
//'
//' @param yd A the outcome variable for document \eqn{d}.
//' @param zbar_d A K x 1 vector containing the empirical topic proportions in
//'   document \eqn{d} (should sum to 1).
//' @param eta A K x 1 vector of regression coefficients.
//' @param sigma2 The residual variance.
//' @param K The number of topics.
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document \eqn{d} excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word \eqn{m} in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
uint16_t draw_zdn_cpp(double yd, const arma::vec& zbar_d, const arma::vec& eta,
                      double sigma2, uint16_t K, uint32_t V,
                      const arma::vec& ndk_n, const arma::vec& nkm_n,
                      const arma::vec& nk_n, float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;
  //   this tends to occur if sigma^2 draw near 0

  if (K < 2) error("number of topics must be at least 2");
  if (V < 2) error("size of vocabulary V must be at least 2");
  if (zbar_d.size() != K) error("zbar_d must be a vector of length K");
  if (eta.size() != K) error("eta must be a vector of length K");
  if (sigma2 < 0.0) error("sigma2 must be positive");
  if (ndk_n.size() != K) error("ndk_n must be a vector of length K");
  if (nkm_n.size() != K) error("nkm_n must be a vector of length K");
  if (nk_n.size() != K) error("nk_n must be a vector of length K");
  if (alpha_ < 0.0) error("alpha_ must be positive");
  if (gamma_ < 0.0) error("gamma_ must be positive");

  arma::vec log_num =
    log(ndk_n + alpha_) +
    log(nkm_n + gamma_) -
    log(nk_n + static_cast<float>(V) * gamma_) +
    R::dnorm(yd, sum(zbar_d.t() * eta), sqrt(sigma2), true);

  long double denom = sum(exp(log_num));
  arma::vec pmf = exp(log_num - log(denom));

  bool good_pmf = true;
  for (uint16_t k = 0; k < K; k++) {
    if (std::isnan(pmf(k)) || std::isinf(pmf(k))) good_pmf = false;
  }

  if (!good_pmf) {
    for (uint16_t k = 0; k < K; k++) pmf(k) = 1.0 / static_cast<float>(K);
  }
  IntegerVector topics = seq_len(K);
  IntegerVector zdn = RcppArmadillo::sample(topics, 1, true, pmf);
  return zdn(0);
}

//' Collapsed Gibbs sampler for the sLDA model
//'
//' @include slda-class.R
//'
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param docs A D x max(\eqn{N_d}) matrix of word indices for all documents.
//' @param w A D x V matrix of counts for all documents and vocabulary terms.
//' @param K The number of topics.
//' @param mu0 A K x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A K x K variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param eta_start A K x 1 vector of starting values for the regression
//'   coefficients.
//' @param constrain_eta A logical (default = \code{FALSE}): If \code{TRUE}, the
//'   regression coefficients will be constrained so that they are in descending
//'   order; if \code{FALSE}, no constraints will be applied.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions
//'   (default: 0.1).
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities (default: 1.01).
//' @param a0 The shape parameter for the prior on sigma2 (default: 0.001)
//' @param b0 The scale parameter for the prior on sigma2 (default: 0.001)
//' @param verbose Should parameter draws be output during sampling? (default:
//'   \code{FALSE}).
//' @param display_progress Should percent progress of sampler be displayed
//'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
//'   \code{display_progress} be set to \code{TRUE} at any given time.
//' @export
// [[Rcpp::export]]
S4 cgibbs_slda_cpp(uint32_t m, uint16_t burn, const arma::colvec& y,
                     const arma::mat& docs, const arma::mat& w, uint16_t K,
                     const arma::colvec& mu0, const arma::mat& sigma0,
                     arma::colvec eta_start, bool constrain_eta = false,
                     float alpha_ = 0.1, float gamma_ = 1.01,
                     float a0 = 0.001, float b0 = 0.001,
                     bool verbose = false, bool display_progress = false) {

  S4 slda("Slda"); // Create object slda of class Slda

  const uint32_t D = w.n_rows;
  const uint32_t V = w.n_cols;
  NumericVector N(D);
  const IntegerVector topics_index = seq_len(K);
  const IntegerVector docs_index   = seq_len(D) - 1;
  const IntegerVector vocab_index  = seq_len(V);

  for (uint32_t d : docs_index) N(d) = sum(docs.row(d) > 0);
  const uint32_t maxNd = max(N);
  arma::mat etam(m, K);
  NumericVector sigma2m(m);
  NumericVector loglike(m); // Store log-likelihood (up to an additive constant)
  NumericVector logpost(m); // Store log-posterior (up to an additive constant)

  // D x K x m array
  arma::cube thetam = arma::zeros(D, K, m);
  // K x V x m array
  arma::cube betam = arma::zeros(K, V, m);
  // D x K x m array to store topic draw counts
  arma::cube ndk = arma::zeros(D, K, m);
  // Topic draws for all words and docs
  arma::mat zdocs = arma::zeros(D, maxNd);

  // Randomly assign topics
  NumericVector init_topic_probs(K);
  arma::mat zbar = arma::zeros(D, K);
  for (uint16_t k = 0; k < K; k++)
    init_topic_probs(k) = 1.0 / static_cast<float>(K);
  for (uint32_t d : docs_index) {
    for (uint32_t n = 0; n < N(d); n++) {
      zdocs(d, n) = RcppArmadillo::sample(
        topics_index, 1, true, init_topic_probs)(0);
    }
    for (uint16_t k = 0; k < K; k++) {
      // Count topic draws in each document
      ndk(d, k, 0) = sum(zdocs.row(d) == k);
      // Compute topic empirical proportions
      zbar(d, k) = static_cast<double>(ndk(d, k, 0)) / static_cast<double>(N(d));
    }
  }

  // Counts of topic-word co-occurences in corpus (K x V)
  arma::mat nkm = arma::zeros(K, V);
  try {
    nkm = count_topic_word_cpp(D, K, V, zdocs, docs);
  } catch(std::exception& e) {
    Rcerr << "Runtime error: " << e.what() <<
      " while computing topic-word co-occurrences\n";
  }
  NumericVector nk(K);
  for (uint16_t k = 0; k < K; k++) nk(k) = sum(ndk.slice(0).col(k));

  // Initialize sigma^2
  sigma2m(0) = var(y);

  if (constrain_eta) {
    // Constrain starting values of eta s.t. components are in descending order
    std::sort(eta_start.begin(), eta_start.end(), std::greater<int>());
  }
  etam.row(0) = eta_start.t();

  if (verbose) {
    Rcout << 1 << "eta: " << etam.row(0) << "~~~~ sigma2: " << sigma2m(0) <<
      "~~~~ zbar2" << zbar.row(1) << "\n";
  }

  // Estimate theta
  for (uint32_t d : docs_index) {
    try {
      thetam.slice(0).row(d) = est_thetad_cpp(
        ndk.slice(0).row(d).t(), alpha_, K).t();
    } catch(std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() <<
        " when estimating theta vector for document " << d << "\n";
    }
  }

  // Estimate beta
  for (uint16_t k = 0; k < K; k++) {
    try {
      betam.slice(0).row(k) = est_betak_cpp(k, V, nkm.row(k).t(),
                  gamma_ = gamma_).t();
    } catch(std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() << " estimating row " << k <<
        "of beta matrix\n";
    }

  }

  // Add likelihood of y
  double temp_prod = arma::as_scalar(
    (y - zbar * etam.row(0).t()).t() * (y - zbar * etam.row(0).t())
  );
  loglike(0) = -0.5 / sigma2m(0) * temp_prod;
  // Add likelihood of documents
  for (uint32_t d : docs_index) {
    for (uint32_t n = 0; n < N(d); n++) {
      uint16_t zdn = zdocs(d, n) - 1; // Topic for word dn (0-indexed)
      uint32_t wdn = docs(d, n) - 1;  // Word for word dn (0-indexed)
      // Not right?
      loglike(0) += log(thetam(d, zdn, 0));  // f(z_{dn} | theta_d)
      loglike(0) += log(betam(zdn, wdn, 0)); // f(w_{dn} | z_{dn}, beta_{z_{dn}})
    }
  }
  logpost(0) = loglike(0);
  // Add prior on eta
  temp_prod = arma::as_scalar(
    (etam.row(0).t() - mu0).t() * sigma0.i() * (etam.row(0).t() - mu0)
  );
  logpost(0) += (-0.5 * temp_prod);
  // Add prior on sigma2
  logpost(0) += ((-0.5 * a0 - 1.0) * log(sigma2m(0)) - 0.5 * b0 / sigma2m(0));
  // Add prior on beta matrix
  double temp_betapost = 0;
  for (uint16_t k = 0; k < K; k++) {
    for (uint32_t v = 0; v < V; v++) {
      temp_betapost += log(betam(k, v, 0));
    }
  }
  logpost(0) += ((gamma_ - 1.0) * temp_betapost);
  // Add prior on theta matrix
  double temp_thetapost = 0;
  for (uint32_t d : docs_index) {
    for (uint16_t k = 0; k < K; k++) {
      temp_thetapost = log(thetam(d, k, 0));
    }
  }
  logpost(0) += ((alpha_ - 1) * temp_thetapost);

  Progress p(m, display_progress);
  for (uint32_t i = 1; i < m; i++) {

    // Draw z
    for (uint32_t d : docs_index) {
      for (uint32_t n = 0; n < N(d); n++) {
        uint32_t word = docs(d, n);
        uint16_t topic = zdocs(d, n);
        // Exclude word n from topic counts in doc d
        arma::vec ndk_n = ndk.slice(i - 1).row(d).t();
        ndk_n(topic - 1)--;
        if (ndk_n(topic - 1) < 0) ndk_n(topic - 1) = 0;
        ndk.slice(i).row(d) = ndk_n.t();
        // Exclude word n from topic counts in corpus
        arma::vec nk_n = nk;
        nk_n(topic - 1)--;
        if (nk_n(topic - 1) < 0) nk_n(topic - 1) = 0;
        nk = nk_n;
        // Exclude word n from topic-word counts
        arma::mat nkm_n = nkm;
        nkm_n(topic - 1, word - 1)--;
        // Fix possible negative counts
        if (nkm_n(topic - 1, word - 1) < 0) nkm_n(topic - 1, word - 1) = 0;
        nkm_n = nkm_n.col(word - 1).t();
        nkm.col(word - 1) = nkm_n.t();

        try {
          topic = draw_zdn_cpp(y(d), zbar.row(d).t(), etam.row(i - 1).t(),
                               sigma2m(i - 1), K, V, ndk_n, nkm_n.t(), nk_n,
                               alpha_, gamma_);
        } catch(std::exception& e) {
          Rcerr << "Runtime Error: " << e.what() <<
            " occurred while drawing topic for word " << n << " in document " <<
            d << "\n";
        }
        zdocs(d, n) = topic;
        // Update topic count in doc d
        ndk(d, topic - 1, i)++;
        // Update topic count in corpus
        nk(topic - 1)++;
        // Update topic-word counts in corpus
        nkm(topic - 1, word - 1)++;
      }
    }

    for (uint32_t d: docs_index) {
      for (uint16_t k = 0; k < K; k++) {
        ndk(d, k, i) = sum(zdocs.row(d) == (k + 1));
        zbar(d, k) = static_cast<double>(ndk(d, k, i)) /
          static_cast<double>(N(d));
        // Estimate beta
        try {
          betam.slice(i).row(k) = est_betak_cpp(k, V, nkm.row(k).t(),
                                                gamma_).t();
        } catch(std::exception& e) {
          Rcerr << "Runtime Error: " << e.what() << " estimating row " << k <<
            "of beta matrix\n";
        }
      }
      // Estimate theta
      try {
        thetam.slice(i).row(d) = est_thetad_cpp(
          ndk.slice(i).row(d).t(), alpha_, K).t();
      } catch(std::exception& e) {
        Rcerr << "Runtime Error: " << e.what() <<
          " when estimating theta vector for document " << d << "\n";
      }
    }

    // Draw eta
    if (constrain_eta) {
      bool eta_order = false;
      uint16_t iter = 0;
      constexpr uint16_t max_iter = 1000;
      arma::vec etac(K);
      while (!eta_order & (iter < max_iter)) {
        iter++;
        try {
          etac = draw_eta_cpp(zbar, y, sigma2m(i - 1), mu0, sigma0).t();
        } catch (std::exception& e) {
          Rcerr << "Runtime Error: " << e.what() <<
            " while drawing eta vector\n";
        }
        for (uint16_t k = 1; k < K; k++) {
          // Force eta components to be in descending order (first is largest)
          //   to resolve label switching of topics
          eta_order = etac(k - 1) >= etac(k);
          if (!eta_order) break;
        }
      }
      etam.row(i) = etac.t();
    } else {
      try {
        etam.row(i) = draw_eta_cpp(zbar, y, sigma2m(i - 1), mu0, sigma0);
      } catch (std::exception& e) {
        Rcerr << "Runtime Error: " << e.what() <<
          " while drawing eta vector\n";
      }
    }

    // Draw sigma2
    try {
      sigma2m(i) = draw_sigma2_cpp(D, a0, b0, zbar, y, etam.row(i).t());
    } catch (std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() << " while drawing sigma2\n";
    }

    // Add likelihood of y
    double temp_prod = arma::as_scalar(
      (y - zbar * etam.row(i).t()).t() * (y - zbar * etam.row(i).t())
    );
    loglike(i) = -0.5 / sigma2m(i) * temp_prod;
    // Add likelihood of documents
    for (uint32_t d : docs_index) {
      for (uint32_t n = 0; n < N(d); n++) {
        uint16_t zdn = zdocs(d, n) - 1; // Topic for word dn (0-indexed)
        uint32_t wdn = docs(d, n) - 1;  // Word for word dn (0-indexed)
        loglike(i) += log(thetam(d, zdn, i));  // f(z_{dn} | theta_d)
        loglike(i) += log(betam(zdn, wdn, i)); // f(w_{dn} | z_{dn}, beta_{z_{dn}})
      }
    }
    logpost(i) = loglike(i);
    // Add prior on eta
    temp_prod = arma::as_scalar(
      (etam.row(i).t() - mu0).t() * sigma0.i() * (etam.row(i).t() - mu0)
    );
    logpost(i) += (-0.5 * temp_prod);
    // Add prior on sigma2
    logpost(i) += ((-0.5 * a0 - 1.0) * log(sigma2m(i)) - 0.5 * b0 / sigma2m(i));
    // Add prior on beta matrix
    double temp_betapost = 0;
    for (uint16_t k = 0; k < K; k++) {
      for (uint32_t v = 0; v < V; v++) {
        temp_betapost += log(betam(k, v, i));
      }
    }
    logpost(i) += ((gamma_ - 1.0) * temp_betapost);
    // Add prior on theta matrix
    double temp_thetapost = 0;
    for (uint32_t d : docs_index) {
      for (uint16_t k = 0; k < K; k++) {
        temp_thetapost = log(thetam(d, k, i));
      }
    }
    logpost(i) += ((alpha_ - 1) * temp_thetapost);

    if (i % 100 == 0) {
      if (verbose) {
        Rcout << i << "eta: " << etam.row(i) << "~~~~ sigma2: " << sigma2m(i) <<
          "~~~~ zbar2" << zbar.row(2) << "\n";
      }
    }
    if (display_progress) {
      p.increment();
    }
    Rcpp::checkUserInterrupt(); // Check to see if user cancelled sampler
  }
  const IntegerVector keep = seq(burn, m - 1);
  NumericVector keep_sigma2(m - burn);
  NumericVector keep_loglike(m - burn);
  NumericVector keep_logpost(m - burn);
  arma::cube keep_beta(K, V, m - burn);
  arma::cube keep_theta(D, K, m - burn);
  for (uint32_t t = 0; t < m - burn; t ++) {
    keep_sigma2(t) = sigma2m(t + burn);
    keep_loglike(t) = loglike(t + burn);
    keep_logpost(t) = logpost(t + burn);
    keep_beta.slice(t) = betam.slice(t + burn);
    keep_theta.slice(t) = thetam.slice(t + burn);
  }

  slda.slot("ntopics") = K;
  slda.slot("ndocs") = D;
  slda.slot("nvocab") = V;
  slda.slot("nchain") = m - burn;
  slda.slot("eta") = etam.rows(burn, m - 1);
  slda.slot("sigma2") = keep_sigma2;
  slda.slot("beta") = keep_beta;
  slda.slot("theta") = keep_theta;
  slda.slot("mu0") = mu0;
  slda.slot("sigma0") = sigma0;
  slda.slot("alpha") = alpha_;
  slda.slot("gamma") = gamma_;
  slda.slot("a0") = a0;
  slda.slot("b0") = b0;
  slda.slot("eta_start") = eta_start;
  slda.slot("loglike") = keep_loglike;
  slda.slot("logpost") = keep_logpost;

  return slda;
}

//' Simulate data from the sLDA model
//'
//' @param D The number of documents in the corpus.
//' @param V The number of terms in the corpus vocabulary.
//' @param N A D x 1 vector of the number of words in each document.
//' @param K The number of topics.
//' @param theta A D x K matrix of topic proportions (each row must sum to 1).
//' @param beta A K x V matrix of word probabilities per topic (each row must
//'   sum to 1).
//'
//' @export
// [[Rcpp::export]]
List sim_slda(uint32_t D, uint32_t V, arma::vec N, uint16_t K, arma::mat theta,
              arma::mat beta, arma::colvec eta, long double sigma2) {

  uint32_t maxN = max(N);
  arma::mat docs   = arma::zeros(D, maxN);
  arma::mat topics = arma::zeros(D, maxN);
  const IntegerVector topic_vec  = seq_len(K);
  const IntegerVector docs_index = seq_len(D) - 1;
  Rcout << "topic_vec" << topic_vec << std::endl;
  const IntegerVector vocab_vec = seq_len(V);
  arma::vec thetad = arma::zeros(K);
  uint32_t doc_len;
  arma::mat w = arma::zeros(D, V);
  arma::mat zbar = arma::zeros(D, K);
  arma::vec betan = arma::zeros(V);

  for (uint32_t d : docs_index) {
    doc_len = N(d);
    arma::vec topic(maxN);
    topic.fill(0);
    arma::vec words(maxN);
    words.fill(0);
    thetad = theta.row(d).t();
    Rcout << thetad << std::endl;
    IntegerVector topic_temp;
    IntegerVector word_temp;

    for (uint32_t n = 0; n < doc_len; n++) {
      betan = beta.row(topic(n)).t();
      topic_temp = RcppArmadillo::sample(topic_vec, 1, true, thetad);
      topic(n) = topic_temp(0);
      word_temp = RcppArmadillo::sample(vocab_vec, 1, true, betan);
      words(n) = word_temp(0);
    }
    docs.row(d) = words.t();
    topics.row(d) = topic.t();
    Rcout << topics.row(d) << std::endl;
  }

  for (uint32_t d : docs_index) {
    doc_len = N(d);
    for (uint32_t v = 0; v < V; v++) {
      for (uint32_t n = 0; n < N(d); n++) {
        if (docs(d, n) == (v + 1)) w(d, v)++;
        for (uint16_t k = 0; k < K; k++) {
          if (topics(d, n) == (k + 1)) zbar(d, k)++;
        }
      }
    }
    zbar.row(d) = zbar.row(d) / doc_len; // Normalize
    Rcout << zbar.row(d) << std::endl;
  }

  // for (uint32_t d : docs_index) {
  //   doc_len = N(d);
  //   // for (uint16_t k = 0; k < K; k++) {
  //   for (uint16_t k : topic_vec) {
  //     for (uint32_t n = 0; n < doc_len; n++) {
  //       if (topics(d, n) == (k + 1)) zbar(d, k)++;
  //     }
  //   }
  //   zbar.row(d) = zbar.row(d) / doc_len; // Normalize
  //   Rcout << zbar.row(d) << std::endl;
  // }

  arma::colvec y(D);
  y = rmvnorm_cpp(1, zbar * eta, sigma2 * arma::eye(D, D)).t();

  return List::create(Rcpp::Named("y")            = y,
                      Rcpp::Named("docs")         = docs,
                      Rcpp::Named("dtm")          = w,
                      Rcpp::Named("num_docs")     = D,
                      Rcpp::Named("num_topics")   = K,
                      Rcpp::Named("doc_lengths")  = N,
                      Rcpp::Named("vocab_length") = V,
                      Rcpp::Named("topics")       = topics,
                      Rcpp::Named("zbar")         = zbar,
                      Rcpp::Named("eta")          = eta,
                      Rcpp::Named("sigma2")       = sigma2,
                      Rcpp::Named("beta")         = beta,
                      Rcpp::Named("theta")        = theta);
}
