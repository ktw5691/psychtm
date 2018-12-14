#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <vector>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

//' Sample from multivariate Gaussian N(\eqn{\mu'}, \eqn{\Sigma})
//'
//' @param n The number of samples to draw.
//' @param mu The mean vector of the distribution (column vector).
//' @param sigma The variance-covariance matrix of the distribution.
//' @export
// [[Rcpp::export]]
arma::mat rmvnorm_cpp(int n, arma::colvec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat y = arma::randn(n, ncols);
  // TODO: Check for symmetric and positive-definite sigma
  // TODO: Check that dimensions of sigma match length of mu
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
//' @export
arma::mat draw_eta_cpp(arma::mat zbar, arma::vec y, double sigma2,
                       arma::vec mu0, arma::mat sigma0) {

  const int D = zbar.n_rows;
  const int K = zbar.n_cols;

  arma::mat ztz(K, K);
  ztz = zbar.t() * zbar;

  arma::colvec eta1(K);
  arma::mat sigma0_inv = sigma0.i();
  arma::mat sigma1 = (sigma0_inv + ztz / sigma2).i();
  eta1 = sigma1 * (sigma0_inv * mu0 + zbar.t() * y / sigma2);
  return rmvnorm_cpp(1, eta1, sigma1);
}

//' Draw sigma2 from full conditional posterior

//' @param D The number of documents.
//' @param a0 The prior shape parameter for \eqn{\sigma^2}.
//' @param b0 The prior scale parameter for \eqn{\sigma^2}.
//' @param zbar A D x K matrix with row \eqn{d} containing the mean number of
//'   draws of topics \eqn{z_1, \ldots, z_K} in document \eqn{d} where each row
//'   sums to 1.
//' @param y A D x 1 vector of the outcome variable.
//' @param eta A K x 1 vector of regression coefficients.
//'
//' @export
double draw_sigma2_cpp(int D, double a0, double b0, arma::mat zbar,
                       arma::colvec y, arma::colvec eta) {

  double a = 0.5 * (D + a0);

  arma::colvec resid(D);
  resid = y - zbar * eta;
  arma::mat b_update(1, 1);
  b_update = resid.t() * resid;

  // Parameterization is 1 / rate
  double b = 1.0 / (0.5 * (b0 + b_update(0, 0)));
  double sigma2inv = R::rgamma(a, b);
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
//' @export
arma::vec est_betak_cpp(int k, long V, arma::vec wz_co, double gamma_) {

  arma::vec betak = exp(log(wz_co + gamma_) - log(sum(wz_co) + V * gamma_));
  for (long v = 0; v < V; v++) {
    // Overflow occured if probability too small
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
//' @export
arma::vec est_thetad_cpp(arma::vec z_count, double alpha_, int K) {
  return exp(log(z_count + alpha_) - log(sum(z_count) + K * alpha_));
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
//' @export
arma::mat count_topic_word_cpp(
    int D, int K, int32_t V, arma::mat doc_topic, arma::mat doc_word) {

  // Matrix to store number of topic/word co-occurences
  arma::mat topic_word_freq = arma::zeros(K, V);
  // Loop through documents
  for (int32_t doc = 0; doc < D; doc++) {
    // Loop through words in document
    for (int32_t pos = 0; pos < doc_word.n_cols; pos++) {
      // Loop through topics
      for (int topic = 1; topic <= K; topic++) {
        // Loop through vocabulary
        for (long v = 1; v <= V; v++) {
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
//' @export
int draw_zdn_cpp(double yd, arma::vec zbar_d, arma::vec eta,
                 double sigma2, int K, long V, arma::vec ndk_n,
                 arma::vec nkm_n, arma::vec nk_n,
                 double alpha_, double gamma_) {

  arma::vec log_num =
    log(ndk_n + alpha_) +
    log(nkm_n + gamma_) -
    log(nk_n + V * gamma_) +
    R::dnorm(yd, sum(zbar_d.t() * eta), sqrt(sigma2), true);

  double denom = sum(exp(log_num));
  arma::vec pmf = exp(log_num - log(denom));

  // Possible to compute probabilities of 0 --> set pmf to {1/K, 1/K, ..., 1/K}
  //   Current state is poor, probabilities for all topics all near 0
  //   Tends to occur if sigma^2 draw near 0
  bool good_pmf = true;
  for (int k = 0; k < K; k++) {
    if (std::isnan(pmf(k)) || std::isinf(pmf(k))) {
      good_pmf = false;
    }
  }

  if (!good_pmf) {
    for (int k = 0; k < K; k++) pmf(k) = 1.0 / K;
  }
  IntegerVector topics = seq_len(K);
  IntegerVector zdn = RcppArmadillo::sample(topics, 1, true, pmf);
  return zdn(0);
}

//' Collapsed Gibbs sampler for the sLDA model
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
List cgibbs_slda_cpp(long m, long burn, arma::colvec y, arma::mat docs,
                     arma::mat w, int K, arma::colvec mu0, arma::mat sigma0,
                     arma::colvec eta_start, bool constrain_eta = false,
                     double alpha_ = 0.1, double gamma_ = 1.01,
                     double a0 = 0.001, double b0 = 0.001,
                     bool verbose = false, bool display_progress = false) {

  const long D = w.n_rows;
  const long V = w.n_cols;
  NumericVector N(D);
  const IntegerVector topics_index = seq_len(K);

  for (long d = 0; d < D; d++) N(d) = sum(docs.row(d) > 0);
  const long maxNd = max(N);
  arma::mat etam(m, K);
  NumericVector sigma2m(m);

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
  for (int k = 0; k < K; k++) init_topic_probs(k) = 1.0 / K;
  for (long d = 0; d < D; d++) {
    for (long n = 0; n < N(d); n++) {
      zdocs(d, n) = RcppArmadillo::sample(topics_index, 1, true, init_topic_probs)(0);
    }
  }

  // Count topic draws in each document
  for (long d = 0; d < D; d++) {
    for (int k = 0; k < K; k++) {
      ndk(d, k, 0) = sum(zdocs.row(d) == (k + 1));
    }
  }

  // Compute topic empirical proportions
  arma::mat zbar = arma::zeros(D, K);
  for (long d = 0; d < D; d++) {
    for (int k = 0; k < K; k++) {
      zbar(d, k) = ndk(d, k, 0) / N(d);
    }
  }

  // Counts of topic-word co-occurences in corpus (K x V)
  arma::mat nkm = count_topic_word_cpp(D, K, V, zdocs, docs);
  NumericVector nk(K);
  for (int k = 0; k < K; k++) nk(k) = sum(ndk.slice(0).col(k));

  // Initialize sigma^2
  sigma2m(0) = var(y);

  if (constrain_eta) {
    // Constrain starting values of eta s.t. components are in descending order
    bool eta_order = false;
    int iter = 0;
    arma::vec etac(K);
    const int max_iter = 1000;
    while (!eta_order & (iter < max_iter)) {
      iter += iter;
      etac = eta_start;
      for (int k = 1; k < K; k++) {
        // Force eta components to be in descending order (first is largest)
        //   to resolve label switching of topics
        eta_order = etac(k - 1) >= etac(k);
        if (!eta_order) {
          eta_start(k - 1) = etac(k);
          eta_start(k) = etac(k - 1);
        }
      }
    }
  }
  etam.row(0) = eta_start.t();

  if (verbose) {
    Rcout << 1 << "eta: " << etam.row(0) << "~~~~ sigma2: " << sigma2m(0) <<
      "~~~~ zbar2" << zbar.row(1) << "\n";
  }

  // Estimate theta
  for (long d = 0; d < D; d++) {
    thetam.slice(0).row(d) = est_thetad_cpp(ndk.slice(0).row(d).t(), alpha_, K).t();
  }

  // Estimate beta
  for (int k = 0; k < K; k++) {
    betam.slice(0).row(k) = est_betak_cpp(k, V, nkm.row(k).t(), gamma_ = gamma_).t();
  }

  Progress p(m, display_progress);
  for (long i = 1; i < m; i++) {

    // Draw z
    for (long d = 0; d < D; d++) {
      for (long n = 0; n < N(d); n++) {
        long word = docs(d, n);
        int topic = zdocs(d, n);
        // Exclude word n from topic counts in doc d
        arma::vec ndk_n = ndk.slice(i - 1).row(d).t();
        ndk_n(topic - 1) = ndk_n(topic - 1) - 1;
        if (ndk_n(topic - 1) < 0) ndk_n(topic - 1) = 0;
        ndk.slice(i).row(d) = ndk_n.t();
        // Exclude word n from topic counts in corpus
        arma::vec nk_n = nk;
        nk_n(topic - 1) = nk_n(topic - 1) - 1;
        if (nk_n(topic - 1) < 0) nk_n(topic - 1) = 0;
        nk = nk_n;
        // Exclude word n from topic-word counts
        arma::mat nkm_n = nkm;
        nkm_n(topic - 1, word - 1) = nkm_n(topic - 1, word - 1) - 1;
        if (nkm_n(topic - 1, word - 1) < 0) nkm_n(topic - 1, word - 1) = 0; // Fix possible negative counts
        nkm_n = nkm_n.col(word - 1).t();
        nkm.col(word - 1) = nkm_n.t();

        topic = draw_zdn_cpp(y(d), zbar.row(d).t(), etam.row(i - 1).t(), sigma2m(i - 1),
                             K, V, ndk_n, nkm_n.t(), nk_n, alpha_, gamma_);
        zdocs(d, n) = topic;
        // Update topic count in doc d
        ndk(d, topic - 1, i) = ndk(d, topic - 1, i) + 1;
        // Update topic count in corpus
        nk(topic - 1) = nk(topic - 1) + 1;
        // Update topic-word counts in corpus
        nkm(topic - 1, word - 1) = nkm(topic - 1, word - 1) + 1;
      }
    }

    for (long d = 0; d < D; d++) {
      for (int k = 0; k < K; k++) {
        ndk(d, k, i) = sum(zdocs.row(d) == (k + 1));
      }
    }

    for (long d = 0; d < D; d++) {
      for (int k = 0; k < K; k++) {
        zbar(d, k) = ndk(d, k, i) / N(d);
      }
    }

    // Estimate theta
    for (long d = 0; d < D; d++) thetam.slice(i).row(d) = est_thetad_cpp(
      ndk.slice(i).row(d).t(), alpha_, K).t();

    // Estimate beta
    for (int k = 0; k < K; k++) betam.slice(i).row(k) = est_betak_cpp(
      k, V, nkm.row(k).t(), gamma_).t();

    // Draw eta
    if (constrain_eta) {
      bool eta_order = false;
      int iter = 0;
      const int max_iter = 1000;
      arma::vec etac(K);
      while (!eta_order & (iter < max_iter)) {
        iter = iter + 1;
        etac = draw_eta_cpp(zbar, y, sigma2m(i - 1), mu0, sigma0).t();
        for (int k = 1; k < K; k++) {
          // Force eta components to be in descending order (first is largest)
          //   to resolve label switching of topics
          eta_order = etac(k - 1) >= etac(k);
          if (!eta_order) {
            break;
          }
        }
      }
      etam.row(i) = etac.t();
    } else {
      etam.row(i) = draw_eta_cpp(zbar, y, sigma2m(i - 1), mu0, sigma0);
    }

    // Draw sigma2
    sigma2m(i) = draw_sigma2_cpp(D, a0, b0, zbar, y, etam.row(i).t());

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
  arma::cube keep_beta(K, V, m - burn);
  arma::cube keep_theta(D, K, m - burn);
  for (long t = 0; t < m - burn; t ++) {
    keep_sigma2(t) = sigma2m(t + burn);
    keep_beta.slice(t) = betam.slice(t + burn);
    keep_theta.slice(t) = thetam.slice(t + burn);
  }
  return List::create(Rcpp::Named("eta")       = etam.rows(burn, m - 1),
                      Rcpp::Named("sigma2")    = keep_sigma2,
                      Rcpp::Named("beta")      = keep_beta,
                      Rcpp::Named("theta")     = keep_theta,
                      Rcpp::Named("mu0")       = mu0,
                      Rcpp::Named("sigma0")    = sigma0,
                      Rcpp::Named("alpha")     = alpha_,
                      Rcpp::Named("gamma")     = gamma_,
                      Rcpp::Named("a0")        = a0,
                      Rcpp::Named("b0")        = b0,
                      Rcpp::Named("eta_start") = eta_start);
}
