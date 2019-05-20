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
//' @param mu The q x 1 mean vector of the distribution (column vector).
//' @param sigma The q x q variance-covariance matrix of the distribution.
//'
//' @return A n x q matrix of random draws.
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

//' Draw eta from full conditional posterior for sLDA/sLDAX/MLR models
//'
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param sigma2 The residual variance.
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//'
//' @return A q x 1 vector of draws of eta.
arma::mat draw_eta_norm(const arma::mat& w, const arma::vec& y,
                        long double sigma2, const arma::vec& mu0,
                        const arma::mat& sigma0) {

  const uint16_t q = w.n_cols;

  arma::mat wtw(q, q);
  wtw = w.t() * w;

  arma::colvec eta1(q);
  arma::mat sigma0_inv = sigma0.i();
  arma::mat sigma1 = (sigma0_inv + wtw / sigma2).i();
  eta1 = sigma1 * (sigma0_inv * mu0 + w.t() * y / sigma2);
  return rmvnorm_cpp(1, eta1, sigma1);
}

//' Compute inverse logit
//' @param x A double
//'
//' @return Inverse-logit(x) on probability scale.
double invlogit(double x) {
  return exp(x) / (1.0 + exp(x));
}

//' Log-likelihood for logistic regression
//'
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix.
//' @param eta A q x 1 vector of regression coefficients.
//'
//' @return The current log-likelihood.
double get_ll_logit(const arma::colvec& y, const arma::mat& w,
                    const arma::colvec& eta) {

  // Add likelihood of y
  uint32_t D = w.n_rows;
  arma::colvec muhat(D);
  muhat = arma::as_scalar(w * eta);

  // Compute log-likelihood of y
  double ll_temp = 0.0;
  for (uint32_t d = 0; d < D; d++) {
    ll_temp += (y(d) * log(invlogit(muhat(d))) +
      (1.0 - y(d)) * log(1.0 / (1.0 + exp(muhat(d)))));
  }

  return ll_temp;
}

//' Log-posterior for regression coefficients with normal prior
//'
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//'
//' @return The current log-posterior.
double get_lpost_eta(double ll, const arma::colvec& eta,
                     const arma::colvec& mu0, const arma::mat& sigma0) {

  double lp_temp = ll;
  // Add prior on eta
  double temp_prod = arma::as_scalar(
    (eta - mu0).t() * sigma0.i() * (eta - mu0)
  );
  lp_temp += (-0.5 * temp_prod);

  return lp_temp;
}

//' Compute full conditional log-posterior of eta for logistic sLDA/sLDAX/regression
//'
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param eta A q x 1 vector of regression coefficients
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//'
//' @return The full conditional log-posterior density of eta.
double eta_logpost_logit(const arma::mat& w, const arma::vec& y,
                         const arma::vec& eta,
                         const arma::vec& mu0, const arma::mat& sigma0) {

  // Compute log-likelihood of y
  double loglike = get_ll_logit(y, w, eta);

  // Add log-prior on eta
  double logpost = get_lpost_eta(loglike, eta, mu0, sigma0);

  return logpost;
}

//' Draw eta from full conditional posterior for logistic sLDA/sLDAX/regression
//'   using Metropolis-Hastings (MH) algorithm
//'
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable for each document.
//' @param eta_prev A q x 1 vector of the previous draw of the regression
//'   coefficients.
//' @param mu0 A q x 1 vector of prior means for the regression coefficients.
//' @param sigma0 A q x q prior variance-covariance matrix for the regression
//'   coefficients.
//' @param proposal_sd A q x 1 vector of proposal distribution standard deviations.
//' @param attempt The number of current attempted draws of eta by MH.
//' @param accept The number of accepted draws of eta by MH.
//'
//' @return A q x 1 draw for eta.
arma::colvec draw_eta_logit(const arma::mat& w, const arma::colvec& y,
                         const arma::colvec& eta_prev,
                         const arma::colvec& mu0, const arma::mat& sigma0,
                         const arma::vec& proposal_sd,
                         arma::vec& attempt, arma::vec& accept) {

  const uint16_t q = w.n_cols;

  // Candidate draws of eta
  arma::colvec cand_eta(q);
  cand_eta = eta_prev;
  arma::vec eta(q);
  eta = eta_prev;
  double cur_logpost = eta_logpost_logit(w, y, eta_prev, mu0, sigma0);

  for (uint16_t j = 0; j < q; j++) {
    cand_eta(j) = rnorm(1, eta_prev(j), proposal_sd(j))(0);
    attempt(j)++; // Passed by reference to update outside function

    // Compute acceptance ratio
    double cand_logpost = eta_logpost_logit(w, y, cand_eta, mu0, sigma0);
    double log_r = cand_logpost - cur_logpost; // Symmetric proposals
    double log_u = log(runif(1)(0));
    if (log_r > log_u) {
      eta(j) = cand_eta(j);
      cur_logpost = cand_logpost;
      accept(j)++; // Passed by reference to update outside function
    } else {
      eta(j) = eta_prev(j);
      cand_eta(j) = eta_prev(j);
    }
  }
  return eta;
}

//' Draw sigma2 from full conditional posterior for sLDA/sLDAX/MLR
//'
//' @param D The number of documents.
//' @param a0 The prior shape parameter for \eqn{\sigma^2}.
//' @param b0 The prior scale parameter for \eqn{\sigma^2}.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable.
//' @param eta A q x 1 vector of regression coefficients.
//'
//' @return A draw for sigma2.
long double draw_sigma2(uint32_t D, float a0, float b0,
                        const arma::mat& w, const arma::colvec& y,
                        const arma::colvec& eta) {

  long double a = 0.5 * (static_cast<float>(D) + a0);

  if ((a0 < 0.0) || (b0 < 0.0)) error("a0 and b0 must be positive");

  if ((w.n_rows != D) || (w.n_cols != eta.size()))
    error("zbar must be a D x q matrix and eta must be a q x 1 vector");
  if (y.size() != D) error("y must be of length D");

  arma::colvec resid(D);
  resid = y - w * eta;
  double b_update = arma::as_scalar(resid.t() * resid);

  // Parameterization is 1 / rate
  long double b = 1.0 / (0.5 * (b0 + b_update));
  long double sigma2inv = R::rgamma(a, b);
  return 1.0 / sigma2inv;
}

//' Estimate beta_k
//'
//' @param k The topic label (i.e., the row of the K x V beta matrix).
//' @param V The number of terms in the corpus vocabulary.
//' @param wz_co A V x 1 vector of counts of the draws of each word for topic
//'   \eqn{k} over all documents.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A V x 1 vector of estimates for beta_k.
//'
//' @export
// [[Rcpp::export]]
arma::vec est_betak_cpp(uint16_t k, uint32_t V, const arma::vec& wz_co,
                        float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic word to 0.0;
  if (gamma_ < 0.0) error("gamma_ must be positive");
  if (V < 2) error("vocabulary size V must be at least 2");
  if (wz_co.size() != V) error("wz_co must be of length V");

  arma::vec betak = exp(log(wz_co + gamma_) - log(sum(wz_co) + V * gamma_));
  for (uint32_t v = 0; v < V; v++) {
    if ((betak[v] > 1.0) | (std::isnan(betak[v]))) betak[v] = 0.0;
  }
  return betak;
}

//' Estimate theta_d
//'
//' @param z_count A K x 1 vector of counts of topic draw in document \eqn{d}.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//' @param K The number of topics.
//'
//' @return A K x 1 vector of estimate for theta_d.
//'
//' @export
// [[Rcpp::export]]
arma::vec est_thetad_cpp(const arma::vec& z_count, float alpha_, uint16_t K) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 0.0;
  if (alpha_ < 0.0) error("alpha_ must be positive");
  if (K < 2) error("number of topics must be at least 2");
  if (z_count.size() != K) error("z_count must be of length K");

  arma::vec thetad = exp(log(z_count + alpha_) - log(sum(z_count) +
    static_cast<float>(K) * alpha_));
  for (uint32_t k = 0; k < K; k++) {
    if ((thetad[k] > 1.0) | (std::isnan(thetad[k]))) thetad[k] = 0.0;
  }

  return thetad;
}

//' Count topic-word co-occurences in corpus (ntopic x nvocab)
//'
//' @param D The number of documents in the corpus.
//' @param K The number of topics.
//' @param V The number of terms in the corpus vocabulary.
//' @param doc_topic A D x max(\eqn{N_d}) matrix of topic assignments for
//'   the corpus.
//' @param doc_word A D x max(\eqn{N_d}) matrix of words for corpus.
//'
//' @return A K x V matrix of topic-word co-occurence counts.
//'
//' @export
// [[Rcpp::export]]
arma::mat count_topic_word(uint32_t D, uint16_t K, uint32_t V,
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

//' Compute log-numerator vector for sampling zdn from full conditional distribution for LDA
//'
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
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A K x 1 vector of the log-numerator from the LDA model to sample zdn.
arma::vec get_log_numer_samplez(uint16_t K, uint32_t V, const arma::vec& ndk_n,
                                const arma::vec& nkm_n, const arma::vec& nk_n,
                                float alpha_, float gamma_) {

  if (K < 2) error("number of topics must be at least 2");
  if (V < 2) error("size of vocabulary V must be at least 2");
  if (ndk_n.size() != K) error("ndk_n must be a vector of length K");
  if (nkm_n.size() != K) error("nkm_n must be a vector of length K");
  if (nk_n.size() != K) error("nk_n must be a vector of length K");
  if (alpha_ < 0.0) error("alpha_ must be positive");
  if (gamma_ < 0.0) error("gamma_ must be positive");

  arma::vec log_num =
    log(ndk_n + alpha_) +
    log(nkm_n + gamma_) -
    log(nk_n + static_cast<float>(V) * gamma_);

  return log_num;
}

//' Draw zdn from full conditional distribution for LDA/sLDA/sLDAX
//'
//' @param log_num A K x 1 vector of the log-numerator for sampling from topics
//'   1, ..., K.
//' @param K The number of topics.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn(arma::vec& log_num, uint16_t K) {

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

//' Draw zdn from full conditional distribution for LDA
//'
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
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn_lda(uint16_t K, uint32_t V,
                      const arma::vec& ndk_n, const arma::vec& nkm_n,
                      const arma::vec& nk_n, float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K
  arma::vec log_num = get_log_numer_samplez(K, V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_);
  uint16_t zdn = draw_zdn(log_num, K);
  return zdn;
}

//' Draw zdn from full conditional distribution for sLDA/sLDAX
//'
//' @param yd A the outcome variable for document \eqn{d}.
//' @param w_d A q x 1 vector containing row d of the predictor model matrix of
//'   assumed form (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
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
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn_slda_norm(double yd, const arma::vec& w_d,
                            const arma::vec& eta, double sigma2, uint16_t K,
                            uint32_t V, const arma::vec& ndk_n,
                            const arma::vec& nkm_n, const arma::vec& nk_n,
                            float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;
  //   this tends to occur if sigma^2 draw near 0
  if (sigma2 < 0.0) error("sigma2 must be positive");

  arma::vec log_num = get_log_numer_samplez(K, V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_) +
    R::dnorm(yd, sum(w_d.t() * eta), sqrt(sigma2), true);
  uint16_t zdn = draw_zdn(log_num, K);
  return zdn;
}

//' Draw zdn from full conditional distribution for sLDA/sLDAX with binary outcome
//'
//' @param yd A the outcome variable for document \eqn{d}.
//' @param w_d A q x 1 vector containing row d of the predictor model matrix of
//'   assumed form (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
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
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn_slda_logit(double yd, const arma::vec& w_d,
                             const arma::vec& eta, uint16_t K,
                             uint32_t V, const arma::vec& ndk_n,
                             const arma::vec& nkm_n, const arma::vec& nk_n,
                             float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;

  double muhat;
  muhat = arma::as_scalar(w_d.t() * eta);

  // Compute log-likelihood of y_d
  double loglike = 0.0;
  loglike += (yd * log(invlogit(muhat) +
    (1.0 - yd) * log(1.0 / (1.0 + exp(muhat)))));

  arma::vec log_num = get_log_numer_samplez(K, V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_) +
    loglike;
    uint16_t zdn = draw_zdn(log_num, K);
    return zdn;
}

//' Posterior predictive likelihood for sLDA/sLDAX/MLR
//'
//' @param w A D x q matrix of additional predictors.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The residual variance.
//'
//' @return Posterior predictive likelihood.
arma::colvec post_pred_norm(const arma::mat& w, const arma::colvec& y,
                            const arma::colvec& eta, double sigma2) {

  const uint16_t D = w.n_rows;
  arma::colvec loglike_pred = arma::zeros(D);
  arma::colvec mu_hat(D);
  mu_hat = w * eta;

  for (uint16_t d = 0; d < D; d++) {
    double yhat;
    yhat = Rcpp::rnorm(1, arma::as_scalar(mu_hat(d)), sigma2)(0);
    double temp_prod = (yhat - mu_hat(d)) * (yhat - mu_hat(d));
    loglike_pred(d) = -0.5 / sigma2 * temp_prod;
  }

  return exp(loglike_pred);
}

//' Posterior predictive likelihood for logistic sLDA/sLDAX/regression
//'
//' @param w A D x q matrix of additional predictors.
//' @param eta A q x 1 vector of regression coefficients.
//' @return Predictive posterior likelihood of all D observations.
arma::colvec post_pred_logit(const arma::mat& w, const arma::colvec& eta) {

  const uint16_t D = w.n_rows;
  arma::colvec y_pred(D); // Store set of D predictions
  arma::colvec loglike_pred(D);
  arma::colvec mu_hat(D);
  mu_hat = w * eta;
  for (uint16_t d = 0; d < D; d++) {
    double phat = arma::as_scalar(invlogit(mu_hat(d)));
    uint16_t yhat;
    yhat = Rcpp::rbinom(1, 1, phat)(0);
    y_pred(d) = yhat;
    loglike_pred(d) = yhat * log(phat) +
      (1.0 - yhat) * log(1.0 / (1.0 + exp(arma::as_scalar(mu_hat(d)))));
  }

  return exp(loglike_pred);
}

//' Contribution to effective number of parameters for WAIC from observation y_d
//'
//' @param like_pred A m x 1 vector of predictive likelihoods (NOT log-likelihoods).
//' @export
//' @return The contribution of y_d (its predictive posterior likelihood variance)
//'   to the effective number of parameters.
// [[Rcpp::export]]
double pwaic_d(const arma::colvec& like_pred) {

  const uint32_t m = like_pred.size();
  double mcmc_meany = 0.0;
  double mcmc_vary = 0.0;
  // Get mean log-predictive likelihood
  for (uint32_t i = 0; i < m; i++) {
    mcmc_meany += log(like_pred(i));
  }
  mcmc_meany /= static_cast<float>(m);

  // Get variance of log-predictive likelihood
  for (uint32_t i = 0; i < m; i++) {
    mcmc_vary += ((log(like_pred(i)) - mcmc_meany) *
      (log(like_pred(i)) - mcmc_meany));
  }
  mcmc_vary /= static_cast<float>(m - 1);

  return mcmc_vary;
}

//' WAIC for observation y_d
//'
//' @param like_pred A m x 1 vector of predictive likelihoods (NOT log-likelihoods) for y_d.
//' @param p_eff The contribution to the effective number of parameters from
//'   obs y_d.
//' @return WAIC contribution for observation d (on deviance scale).
//' @export
// [[Rcpp::export]]
double waic_d(const arma::colvec& like_pred, double p_effd) {

  const uint32_t m = like_pred.n_rows;
  double likeyd = 0;

  // Get mean predictive likelihood for y_d
  for (uint32_t i = 0; i < m; i++) {
    likeyd += like_pred(i);
  }
  likeyd /= static_cast<float>(m);

  // Log of mean post. pred. density for y_d
  double lppd = log(likeyd);

  // Contribution to WAIC (on deviance scale, i.e., -2ll) for y_d
  return (-2.0 * (lppd - p_effd)); // See Gelman, Hwang, Vehtari (2014, p. 1003)
}

//' Compute WAIC for all outcomes.
//'
//' @param D The number of documents.
//' @param iter The current iteration of the chain.
//' @param l_pred A m x D matrix of predictive likelihoods (NOT log-likelihoods).
//' @return Vector of (1) WAIC for model, (2) standard error for WAIC, and (3)
//'   the effective number of parameters.
//' @export
// [[Rcpp::export]]
NumericVector waic_all(uint16_t D, uint32_t iter, const arma::mat& l_pred) {

  NumericVector full_waic_se(3);

  double peff_sum = 0.0;
  double waic_sum = 0.0;
  arma::colvec waic = arma::zeros(D);
  arma::colvec peff = arma::zeros(D);
  double se_waic = 0.0;

  // Compute WAIC and p_eff for entire data set
  for (uint16_t d = 0; d < D; d++) {
    peff(d) = pwaic_d(l_pred.submat(0, d, iter - 1, d));
    peff_sum += peff(d);
    waic(d) = waic_d(l_pred.submat(0, d, iter - 1, d), peff(d));
    waic_sum += waic(d);
  }

  // Compute SE(WAIC)
  double mean_waic = waic_sum / static_cast<float>(D); // Mean WAIC over all docs
  for (uint16_t d = 0; d < D; d++) {
    se_waic += ((waic(d) - mean_waic) * (waic(d) - mean_waic));
  }
  se_waic /= static_cast<float>(D - 1); // Variance of WAIC over all docs
  se_waic *= static_cast<float>(D); // Variance of WAIC times D
  se_waic = sqrt(se_waic); // SE(WAIC)

  full_waic_se(0) = waic_sum;
  full_waic_se(1) = se_waic;
  full_waic_se(2) = peff_sum;

  return full_waic_se;
}

//' Compute difference (WAIC1 - WAIC2) in WAIC and its SE for two models.
//'
//' @param D The number of documents.
//' @param m1 The length of the chain for model 1.
//' @param m2 The length of the chain for model 2.
//' @param l_pred1 A m x D matrix of predictive likelihoods (NOT log-likelihoods) from model 1.
//' @param l_pred2 A m x D matrix of predictive likelihoods (NOT log-likelihoods) from model 2.
//' @return A vector of (1) the difference in WAIC (on the deviance scale)
//'   between models and (2) the standard error of the difference in WAIC.
//' @export
// [[Rcpp::export]]
NumericVector waic_diff(uint16_t D, uint32_t m1, uint32_t m2,
                        const arma::mat& l_pred1, const arma::mat& l_pred2) {

  NumericVector diff_waic_se(2);

  double waic_diff_sum = 0.0;
  arma::colvec waic1 = arma::zeros(D);
  arma::colvec waic2 = arma::zeros(D);
  arma::colvec peff1 = arma::zeros(D);
  arma::colvec peff2 = arma::zeros(D);
  arma::colvec waic_diff = arma::zeros(D);
  double se_waic_diff = 0.0;

  // Compute WAIC and p_eff for entire data set
  for (uint16_t d = 0; d < D; d++) {
    peff1(d) = pwaic_d(l_pred1.submat(0, d, m1 - 1, d));
    peff2(d) = pwaic_d(l_pred2.submat(0, d, m2 - 1, d));
    waic1(d) = waic_d(l_pred1.submat(0, d, m1 - 1, d), peff1(d));
    waic2(d) = waic_d(l_pred2.submat(0, d, m2 - 1, d), peff2(d));
    waic_diff(d) = waic1(d) - waic2(d); // Difference for doc d
    waic_diff_sum += waic_diff(d); // Sum of differences
  }

  // Compute SE(WAIC1 - WAIC2)
  double mean_diff_waic = waic_diff_sum / static_cast<float>(D); // Mean difference
  for (uint16_t d = 0; d < D; d++) {
    se_waic_diff += ((waic_diff(d) - mean_diff_waic) * (waic_diff(d) - mean_diff_waic));
  }
  se_waic_diff /= static_cast<float>(D - 1); // Variance of difference
  se_waic_diff *= static_cast<float>(D); // Multiply variance of diff by D
  se_waic_diff = sqrt(se_waic_diff); // SE(WAIC1 - WAIC2)

  diff_waic_se(0) = waic_diff_sum; // Difference
  diff_waic_se(1) = se_waic_diff; // SE(WAIC1 - WAIC2)

  return diff_waic_se;
}

//' Log-likelihood for MLR
//'
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//'
//' @return The current log-likelihood.
double get_ll_mlr(const arma::colvec& y, const arma::mat& w,
                  const arma::colvec& eta, double sigma2) {

  double temp_prod = arma::as_scalar(
    (y - w * eta).t() * (y - w * eta)
  );
  double ll_temp = -0.5 / sigma2 * temp_prod;

  return ll_temp;
}

//' Log-likelihood for sLDA/sLDAX model
//'
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//' @param zdocs A D x max(\eqn{N_d}) matrix of topic indicators for all documents.
//' @param docs A D x max(\eqn{N_d}) matrix of word indicators for all documents.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//' @param N A vector of length D containing the number of words in each document.
//'
//' @return The current log-likelihood.
double get_ll_slda_norm(const arma::colvec& y, const arma::mat& w,
                        const arma::colvec& eta, double sigma2,
                        const arma::mat& zdocs, const arma::mat& docs,
                        const arma::mat& theta, const arma::mat& beta,
                        const IntegerVector& docs_index, const NumericVector& N) {

  double ll_temp = get_ll_mlr(y, w, eta, sigma2);
  // Add likelihood of documents
  for (uint32_t d : docs_index) {
    for (uint32_t n = 0; n < N(d); n++) {
      uint16_t zdn = zdocs(d, n) - 1; // Topic for word dn (0-indexed)
      uint32_t wdn = docs(d, n) - 1;  // Word for word dn (0-indexed)
      ll_temp += log(theta(d, zdn));  // f(z_{dn} | theta_d)
      ll_temp += log(beta(zdn, wdn)); // f(w_{dn} | z_{dn}, beta_{z_{dn}})
    }
  }
  return ll_temp;
}

//' Log-likelihood for logistic sLDA/sLDAX model
//'
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param zdocs A D x max(\eqn{N_d}) matrix of topic indicators for all documents.
//' @param docs A D x max(\eqn{N_d}) matrix of word indicators for all documents.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//' @param N A vector of length D containing the number of words in each document.
//'
//' @return The current log-likelihood.
// double get_ll_slda_logit(const arma::colvec& y, const arma::mat& w,
//                         const arma::colvec& eta,
//                         const arma::mat& zdocs, const arma::mat& docs,
//                         const arma::mat& theta, const arma::mat& beta,
//                         const IntegerVector& docs_index, const NumericVector& N) {
// }

//' Log-posterior for normal outcome regression
//'
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param a0 The shape parameter for the prior on sigma2.
//' @param b0 The scale parameter for the prior on sigma2.
//'
//' @return The current log-posterior.
double get_lpost_mlr(double ll,
                     const arma::colvec& eta, double sigma2,
                     const arma::colvec& mu0, const arma::mat& sigma0,
                     double a0, double b0) {

  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);

  // Add prior on sigma2
  lp_temp += ((-0.5 * a0 - 1.0) * log(sigma2) - 0.5 * b0 / sigma2);

  return lp_temp;
}

//' Log-posterior for sLDA/sLDAX model
//'
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The current draw of the residual variance of y.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//' @param a0 The shape parameter for the prior on sigma2.
//' @param b0 The scale parameter for the prior on sigma2.
//' @param V The number of words in the vocabulary.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//'
//' @return The current log-posterior.
double get_lpost_slda_norm(double ll, const arma::colvec& eta, double sigma2,
                           const arma::mat& theta, const arma::mat& beta,
                           const arma::colvec& mu0, const arma::mat& sigma0,
                           double gamma_, double alpha_,
                           double a0, double b0, uint32_t V,
                           const IntegerVector& docs_index) {

  uint16_t K = theta.n_cols;
  double lp_temp = get_lpost_mlr(ll, eta, sigma2, mu0, sigma0, a0, b0);

  // Add prior on beta matrix
  double temp_betapost = 0;
  for (uint16_t k = 0; k < K; k++) {
    for (uint32_t v = 0; v < V; v++) {
      temp_betapost += log(beta(k, v));
    }
  }
  lp_temp += ((gamma_ - 1.0) * temp_betapost);
  // Add prior on theta matrix
  double temp_thetapost = 0;
  for (uint32_t d : docs_index) {
    for (uint16_t k = 0; k < K; k++) {
      temp_thetapost = log(theta(d, k));
    }
  }
  lp_temp += ((alpha_ - 1) * temp_thetapost);

  return lp_temp;
}

//' Log-posterior for logistic sLDA/sLDAX model
//'
//' @param ll A double of the current log-likelihood.
//' @param eta A q x 1 vector of regression coefficients.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param mu0 A q x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A q x q variance-covariance matrix for the prior on the
//'   regression coefficients.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//' @param V The number of words in the vocabulary.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//'
//' @return The current log-posterior.
double get_lpost_slda_logit(double ll, const arma::colvec& eta,
                            const arma::mat& theta, const arma::mat& beta,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            double gamma_, double alpha_,
                            uint32_t V, const IntegerVector& docs_index) {

  uint16_t K = theta.n_cols;
  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);

  // Add prior on beta matrix
  double temp_betapost = 0;
  for (uint16_t k = 0; k < K; k++) {
    for (uint32_t v = 0; v < V; v++) {
      temp_betapost += log(beta(k, v));
    }
  }
  lp_temp += ((gamma_ - 1.0) * temp_betapost);
  // Add prior on theta matrix
  double temp_thetapost = 0;
  for (uint32_t d : docs_index) {
    for (uint16_t k = 0; k < K; k++) {
      temp_thetapost = log(theta(d, k));
    }
  }
  lp_temp += ((alpha_ - 1) * temp_thetapost);

  return lp_temp;
}

//' Update number of times topic was drawn in document excluding current word
//'
//' @param d The current document index.
//' @param n The current word index in document d.
//' @param topic The current topic index.
//' @param ndk A vector of the current number of draws of each topic in document d.
//'
//' @return A vector of the current number of draws of each topic in document d
//'   excluding word n.
arma::vec count_topicd(uint32_t d, uint32_t n, uint16_t topic,
                       const arma::vec& ndk) {
  // Exclude word n from topic counts in doc d
  arma::vec ndk_n = ndk;
  ndk_n(topic - 1)--;
  if (ndk_n(topic - 1) < 0) ndk_n(topic - 1) = 0;
  return ndk_n;
}

//' Update number of times topic was drawn in corpus excluding current word
//'
//' @param topic The current topic index.
//' @param nk A vector of the current number of draws of each topic in the corpus.
//'
//' @return A vector of the current number of draws of each topic in the corpus
//'   excluding the current word.
arma::vec count_topic_corpus(uint16_t topic, const arma::vec& nk) {
  // Exclude word n from topic counts in corpus
  arma::vec nk_n = nk;
  nk_n(topic - 1)--;
  if (nk_n(topic - 1) < 0) nk_n(topic - 1) = 0;
  return nk_n;
}

//' Update number of times word and topic co-occur in corpus
//'
//' @param word The current word index.
//' @param topic The current topic index.
//' @param nkm A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus.
//'
//' @return A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus excluding the current word.
arma::mat count_word_topic(uint32_t word, uint16_t topic, const arma::mat& nkm) {
  // Exclude word n from topic-word counts
  arma::mat nkm_n = nkm;
  nkm_n(topic - 1, word - 1)--;
  // Fix possible negative counts
  if (nkm_n(topic - 1, word - 1) < 0) nkm_n(topic - 1, word - 1) = 0;
  nkm_n = nkm_n.col(word - 1).t();
  return nkm_n;
}

//' Update all topic and topic-word counts in document and corpus
//'
//' @param d The current document index.
//' @param n The current word index in document d.
//' @param word The current word index.
//' @param topic The current topic index.
//' @param ndk A vector of the current number of draws of each topic in document d.
//' @param nk A vector of the current number of draws of each topic in the corpus.
//' @param nkm A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus.
void update_zcounts(uint32_t d, uint32_t n, uint32_t word, uint16_t topic,
                    arma::mat& ndk, NumericVector& nk, arma::mat& nkm) {

  arma::vec ndktemp = ndk.row(d).t();
  NumericVector nktemp = nk;
  arma::mat nkmtemp = nkm;

  ndk.row(d) = count_topicd(d, n, topic, ndktemp).t();
  nk = count_topic_corpus(topic, nktemp);
  nkm.col(word - 1) = count_word_topic(word, topic, nkmtemp).t();

}

//////////////////////////////// Gibbs Samplers ////////////////////////////////

//' Collapsed Gibbs sampler for multiple linear regression
//'
//' @include slda-class.R
//'
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param y A D x 1 vector of outcomes to be predicted.
//' @param x A D x (p + 1) matrix of additional predictors.
//' @param mu0 A (p + 1) x 1 mean vector for the prior on the regression
//'   coefficients.
//' @param sigma0 A (p + 1) x (p + 1) variance-covariance matrix for the
//'   prior on the regression coefficients.
//' @param eta_start A (p + 1) x 1 vector of starting values for the
//'   regression coefficients.
//' @param a0 The shape parameter for the prior on sigma2 (default: 0.001)
//' @param b0 The scale parameter for the prior on sigma2 (default: 0.001)
//' @param verbose Should parameter draws be output during sampling? (default:
//'   \code{FALSE}).
//' @param display_progress Should percent progress of sampler be displayed
//'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
//'   \code{display_progress} be set to \code{TRUE} at any given time.
//'
//' @return An object of class Mlr.
//' @export
// [[Rcpp::export]]
S4 gibbs_mlr(uint32_t m, uint32_t burn, const arma::colvec& y,
             const arma::mat& x,
             const arma::colvec& mu0, const arma::mat& sigma0,
             arma::colvec eta_start, float a0 = 0.001, float b0 = 0.001,
             bool verbose = false, bool display_progress = false) {

  if (m <= burn) {
    stop("Length of chain m not greater than burn-in period.");
  }

  S4 slda("Mlr"); // Create object slda of class Mlr

  const uint32_t D = x.n_rows;
  const uint16_t pp1 = x.n_cols;

  arma::mat etam(m - burn, pp1);
  NumericVector sigma2m(m - burn);
  NumericVector loglike(m - burn); // Store log-likelihood (up to an additive constant)
  NumericVector logpost(m - burn); // Store log-posterior (up to an additive constant)
  arma::mat l_pred(m - burn, D);

  // Initialize sigma^2
  double sigma2 = var(y) / 2.0;

  if (verbose) {
    Rcout << 1 << "eta: " << eta_start.t() << " ~~~~ sigma2: " << sigma2 << "\n";
  }

  arma::colvec eta(pp1);

  Progress p(m, display_progress);
  for (uint32_t i = 1; i <= m; i++) {

    // Draw eta
    try {
      // etam.row(i) = draw_eta_norm(x, y, sigma2m(i - 1), mu0, sigma0);
      eta = draw_eta_norm(x, y, sigma2, mu0, sigma0).t();
    } catch (std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() <<
        " while drawing eta vector\n";
    }

    // Draw sigma2
    try {
      // sigma2m(i) = draw_sigma2(D, a0, b0, x, y, etam.row(i).t());
      sigma2 = draw_sigma2(D, a0, b0, x, y, eta);
    } catch (std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() << " while drawing sigma2\n";
    }

    if (i > burn) {
      // Likelihood
      loglike(i - burn - 1) = get_ll_mlr(y, x, eta, sigma2);
      // Log-posterior
      logpost(i - burn - 1) = get_lpost_mlr(loglike(i - burn - 1), eta, sigma2,
                                            mu0, sigma0, a0, b0);

      l_pred.row(i - burn - 1) = post_pred_norm(x, y, eta, sigma2).t();

      etam.row(i - burn - 1) = eta.t();
      sigma2m(i - burn - 1) = sigma2;
    }
    if (i % 500 == 0) {
      if (verbose) {
        Rcout << i << "eta: " << eta.t() << " ~~~~ sigma2: " << sigma2 << "\n";
      }
    }
    if (display_progress) {
      p.increment();
    }
    Rcpp::checkUserInterrupt(); // Check to see if user cancelled sampler
  }

  // Compute WAIC and p_eff
  NumericVector waic_and_se(3);
  waic_and_se = waic_all(D, m - burn, l_pred);

  slda.slot("ndocs") = D;
  slda.slot("nchain") = m - burn;
  slda.slot("eta") = etam;
  slda.slot("sigma2") = sigma2m;
  slda.slot("mu0") = mu0;
  slda.slot("sigma0") = sigma0;
  slda.slot("a0") = a0;
  slda.slot("b0") = b0;
  slda.slot("eta_start") = eta_start;
  slda.slot("loglike") = loglike;
  slda.slot("logpost") = logpost;
  slda.slot("p_eff") = waic_and_se(2);
  slda.slot("waic") = waic_and_se(0);
  slda.slot("se_waic") = waic_and_se(1);
  slda.slot("lpd") = l_pred;

  return slda;
}

//' Collapsed Gibbs sampler for logistic regression
//'
//' @include slda-class.R
//'
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param y A D x 1 vector of binary outcomes (0/1) to be predicted.
//' @param x A D x p matrix of additional predictors (no column of 1s for
//'   intercept).
//' @param mu0 A (p + 1) x 1 mean vector for the prior on the regression coefficients.
//' @param sigma0 A (p + 1) x (p + 1) variance-covariance matrix for the prior
//'   on the regression coefficients.
//' @param eta_start A (p + 1) x 1 vector of starting values for the
//'   regression coefficients.
//' @param proposal_sd The proposal standard deviation for drawing the
//'   regression coefficients, N(0, proposal_sd) (default: 0.2).
//' @param verbose Should parameter draws be output during sampling? (default:
//'   \code{FALSE}).
//' @param display_progress Should percent progress of sampler be displayed
//'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
//'   \code{display_progress} be set to \code{TRUE} at any given time.
//' @export
// [[Rcpp::export]]
S4 gibbs_logistic(uint32_t m, uint32_t burn, const arma::colvec& y,
                  const arma::mat& x,
                  const arma::colvec& mu0, const arma::mat& sigma0,
                  arma::colvec eta_start, arma::vec proposal_sd,
                  bool verbose = false, bool display_progress = false) {

  if (m <= burn) {
    stop("Length of chain m not greater than burn-in period.");
  }

  S4 slda("Logistic"); // Create object slda of class Logistic

  const uint32_t D = y.size();
  const uint16_t pp1 = x.n_cols;

  arma::mat etam(m - burn, pp1);
  NumericVector loglike(m - burn); // Store log-likelihood (up to an additive constant)
  NumericVector logpost(m - burn); // Store log-posterior (up to an additive constant)
  arma::mat l_pred(m - burn, D);

  arma::vec attempt = arma::zeros(pp1);
  arma::vec accept = arma::zeros(pp1);

  arma::colvec eta(pp1);

  Progress prog(m, display_progress);
  for (uint32_t i = 1; i <= m; i++) {

    // Draw eta
    try {
      eta = draw_eta_logit(x, y, eta,
                           mu0, sigma0, proposal_sd,
                           attempt, accept);
    } catch (std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() <<
        " while drawing eta vector\n";
    }

    // Update acceptance rates and tune proposal standard deviations
    arma::vec acc_rate(pp1);
    for (uint16_t j = 0; j < pp1; j++) {
      acc_rate(j) = static_cast<float>(accept(j)) / static_cast<float>(attempt(j));
    }
    for (uint16_t j = 0; j < pp1; j++) {
      if (i < (static_cast<float>(burn) / 5.0) && attempt(j) >= 50 && (i % 50 == 0)) {
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

    if (i > burn) {
      // Likelihood
      loglike(i - burn - 1) = get_ll_logit(y, x, eta);
      // Log-posterior
      logpost(i - burn - 1) = get_lpost_eta(loglike(i - burn - 1), eta,
                                            mu0, sigma0);

      l_pred.row(i - burn - 1) = post_pred_logit(x, eta).t();

      etam.row(i - burn - 1) = eta.t();
    }
    if (i % 500 == 0) {
      if (verbose) {
        Rcout << i << "eta: " << eta.t() << "\n" <<
          "accept rate: " << acc_rate.t() << "\n" <<
          "prop_sd: " << proposal_sd.t() << "\n";
      }
    }
    if (display_progress) {
      prog.increment();
    }
    Rcpp::checkUserInterrupt(); // Check to see if user cancelled sampler
  }

  // Compute WAIC and p_eff
  NumericVector waic_and_se(3);
  waic_and_se = waic_all(D, m - burn, l_pred);

  slda.slot("ndocs") = D;
  slda.slot("nchain") = m - burn;
  slda.slot("eta") = etam;
  slda.slot("mu0") = mu0;
  slda.slot("sigma0") = sigma0;
  slda.slot("eta_start") = eta_start;
  slda.slot("proposal_sd") = proposal_sd;
  slda.slot("loglike") = loglike;
  slda.slot("logpost") = logpost;
  slda.slot("p_eff") = waic_and_se(2);
  slda.slot("waic") = waic_and_se(0);
  slda.slot("se_waic") = waic_and_se(1);
  slda.slot("lpd") = l_pred;

  return slda;
}

//' gibbs_lda()/gibbs_slda()/gibbs_sldax()/gibbs_slda_logit()/gibbs_sldax_logit()
