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
  const uint16_t ncols = sigma.n_cols;
  if ((mu.size() != sigma.n_rows) || (mu.size() != ncols)) {
    Rcerr <<
      "sigma must be a square matrix and mu must be a column vector with length equal to the number of rows and columns in sigma\n";
  }
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

  arma::mat wtw = w.t() * w;
  arma::mat sigma0_inv = sigma0.i();
  arma::mat sigma1 = (sigma0_inv + wtw / sigma2).i();
  arma::colvec eta1 = sigma1 * (sigma0_inv * mu0 + w.t() * y / sigma2);

  return rmvnorm_cpp(1, eta1, sigma1);
}

//' Compute inverse logit
//' @param x A double
//'
//' @return Inverse-logit(x) on probability scale.
double invlogit(double x) {
  double temp = exp(x) / (1.0 + exp(x));
  return temp;
}

//' Log-likelihood for logistic regression for observation d
//'
//' @param yd An integer 0/1 outcome to be predicted.
//' @param muhatd A double predicted outcome on logit scale.
//'
//' @return The current log-likelihood for observation d.
double get_ll_logit_yd(int yd, double muhatd) {

  // Compute log-likelihood of y
  double ll_temp = yd * muhatd - log(1.0 + exp(muhatd));
  return ll_temp;
}

//' Log-likelihood for logistic regression
//'
//' @param y A D x 1 vector of 0/1 outcomes to be predicted.
//' @param w A D x q matrix containing a predictor model matrix.
//' @param eta A q x 1 vector of regression coefficients.
//'
//' @return The current log-likelihood.
double get_ll_logit(const arma::colvec& y, const arma::mat& w,
                    const arma::colvec& eta) {

  // Add likelihood of y
  uint32_t D = w.n_rows;
  arma::colvec muhat = w * eta;

  // Compute log-likelihood of y
  double ll_temp = 0.0;
  for (uint32_t d = 0; d < D; d++) {
    ll_temp += get_ll_logit_yd(y(d), arma::as_scalar(muhat(d)));
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
//' @param attempt A vector of the number of current attempted draws of eta.
//' @param accept A vector of the number of accepted draws of eta.
//'
//' @return A q x 1 vector of draws for eta.
arma::colvec draw_eta_logit(const arma::mat& w, const arma::colvec& y,
                         const arma::colvec& eta_prev,
                         const arma::colvec& mu0, const arma::mat& sigma0,
                         const arma::vec& proposal_sd,
                         arma::vec& attempt, arma::vec& accept) {

  const uint16_t q = w.n_cols;
  arma::colvec cand_eta = eta_prev; // Candidate draws of eta
  arma::vec eta = eta_prev;
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
//' @param a0 The prior shape parameter for \eqn{\sigma^2}.
//' @param b0 The prior scale parameter for \eqn{\sigma^2}.
//' @param w A D x q matrix containing a predictor model matrix of assumed form
//'   (X, Zbar, XZbarInteractions).
//' @param y A D x 1 vector of the outcome variable.
//' @param eta A q x 1 vector of regression coefficients.
//'
//' @return A draw for \eqn{\sigma^2}.
long double draw_sigma2(float a0, float b0,
                        const arma::mat& w, const arma::colvec& y,
                        const arma::colvec& eta) {

  const uint32_t D = y.size(); // Number of observations
  if ((a0 < 0.0) || (b0 < 0.0)) error("a0 and b0 must be positive");
  long double a = 0.5 * (static_cast<float>(D) + a0);

  if ((w.n_rows != D) || (w.n_cols != eta.size()))
    error("zbar must be a D x q matrix and eta must be a q x 1 vector");

  arma::colvec resid(D);
  resid = y - w * eta;
  double b_update = arma::as_scalar(resid.t() * resid);

  // Parameterization is 1 / rate
  long double b = 1.0 / (0.5 * (b0 + b_update));
  long double sigma2inv = R::rgamma(a, b);
  return 1.0 / sigma2inv;
}

//' Estimate \eqn{\beta_k}
//'
//' @param wz_co A V x 1 vector of counts of the draws of each word for topic
//'   k over all documents.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A V x 1 vector of estimates for \eqn{\beta_k}.
//'
//' @export
// [[Rcpp::export]]
arma::vec est_betak(const arma::vec& wz_co, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic word to 0.0;
  if (gamma_ < 0.0) error("gamma_ must be positive");
  const uint32_t V = wz_co.size(); // Vocabulary size
  if (V < 2) error("vocabulary size V must be at least 2");

  arma::vec betak = exp(log(wz_co + gamma_) - log(sum(wz_co) + V * gamma_));
  for (uint32_t v = 0; v < V; v++) {
    if ((betak[v] > 1.0) | (std::isnan(betak[v]))) betak[v] = 0.0;
  }
  return betak;
}

//' Estimate \eqn{\theta_d}
//'
//' @param z_count A K x 1 vector of counts of topic draw in document d.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//'
//' @return A K x 1 vector of estimate for \eqn{\theta_d}.
//'
//' @export
// [[Rcpp::export]]
arma::vec est_thetad(const arma::vec& z_count, float alpha_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 0.0;
  const uint16_t K = z_count.size(); // Number of topics
  if (alpha_ < 0.0) error("alpha_ must be positive");
  if (K < 2) error("number of topics must be at least 2");

  arma::vec thetad = exp(log(z_count + alpha_) - log(sum(z_count) +
    static_cast<float>(K) * alpha_));
  for (uint32_t k = 0; k < K; k++) {
    if ((thetad[k] > 1.0) | (std::isnan(thetad[k]))) thetad[k] = 0.0;
  }

  return thetad;
}

//' Count topic-word co-occurences in corpus.
//'
//' Computes topic-word co-occurence matrix for a corpus of \eqn{D} documents
//' with the maximum length of a document in the corpus equal to max(\eqn{N_d})
//' and a vocabulary of \eqn{V} unique terms in the corpus.
//'
//' Indices in \code{doc_topic} and \code{doc_word} where no word exists in the
//' document must be set to 0.
//'
//' @param K The number of topics.
//' @param V The number of terms in the corpus vocabulary.
//' @param doc_topic A \eqn{D} x max(\eqn{N_d}) matrix of topic assignments for
//'   the corpus.
//' @param doc_word A \eqn{D} x max(\eqn{N_d}) matrix of words for corpus.
//'
//' @return A \eqn{K} x \eqn{V} matrix of topic-word co-occurence counts.
//'
//' @export
// [[Rcpp::export]]
arma::mat count_topic_word(uint16_t K, uint32_t V,
                           const arma::mat& doc_topic,
                           const arma::mat& doc_word) {
  const uint32_t D = doc_topic.n_rows; // Number of documents
  if (K < 2) error("number of topics must be at least 2");
  if (V < 2) error("size of vocabulary V must be at least 2");
  const uint32_t maxnd = doc_word.n_cols; // Maximum document length
  if (doc_word.n_rows != D) error("'doc_word' and 'doc_topic' must have the same number of rows");
  if (doc_topic.n_cols != maxnd)
    error("'doc_topic' and 'doc_word' must have the same number of columns");

  const IntegerVector topics_index = seq_len(K);
  const IntegerVector docs_index = seq_len(D) - 1;
  // Matrix to store number of topic/word co-occurences
  arma::mat topic_word_freq = arma::zeros(K, V);
  // Loop through documents
  for (uint32_t doc : docs_index) {
    // Loop through words in document
    for (uint32_t pos = 0; pos < maxnd; pos++) {
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
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//' @param alpha_ The hyperparameter on the Dirichlet prior for \eqn{\theta_d}.
//' @param gamma_ The hyperparameter for the Dirichlet priors on \eqn{\beta_k}.
//'
//' @return A K x 1 vector of the log-numerator from the LDA model to sample zdn.
arma::vec get_log_numer_samplez(uint32_t V, const arma::vec& ndk_n,
                                const arma::vec& nkm_n, const arma::vec& nk_n,
                                float alpha_, float gamma_) {

  const uint16_t K = nk_n.size(); // Number of topics
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
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn(arma::vec& log_num) {

  const uint16_t K = log_num.size(); // Number of topics
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
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn_lda(uint32_t V,
                      const arma::vec& ndk_n, const arma::vec& nkm_n,
                      const arma::vec& nk_n, float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K
  arma::vec log_num = get_log_numer_samplez(V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_);
  uint16_t zdn = draw_zdn(log_num);
  return zdn;
}

//' Draw zdn from full conditional distribution for sLDA/sLDAX
//'
//' @param yd The outcome variable for document \eqn{d}.
//' @param w_d A q x 1 vector containing row d of the predictor model matrix of
//'   assumed form (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param sigma2 The residual variance.
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn_slda_norm(double yd, const arma::vec& w_d,
                            const arma::vec& eta, double sigma2,
                            uint32_t V, const arma::vec& ndk_n,
                            const arma::vec& nkm_n, const arma::vec& nk_n,
                            float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;
  //   this tends to occur if sigma^2 draw near 0
  if (sigma2 < 0.0) error("sigma2 must be positive");
  arma::vec log_num = get_log_numer_samplez(V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_) +
    R::dnorm(yd, sum(w_d.t() * eta), sqrt(sigma2), true);
  uint16_t zdn = draw_zdn(log_num);
  return zdn;
}

//' Draw zdn from full conditional distribution for sLDA/sLDAX with binary outcome
//'
//' @param yd A the outcome variable for document d.
//' @param w_d A q x 1 vector containing row d of the predictor model matrix of
//'   assumed form (X, Zbar, XZbarInteractions).
//' @param eta A q x 1 vector of regression coefficients.
//' @param V The number of terms in the corpus vocabulary.
//' @param ndk_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} in
//'   document d excluding the current word \eqn{w_n} from the counts.
//' @param nkm_n A K x 1 vector of counts of topic \eqn{k = 1, \ldots, K} and
//'   word m in the corpus excluding the current word \eqn{w_n} from the
//'   counts.
//' @param nk_n A K x 1 vector of counts of draws of topic
//'   \eqn{k = 1, \ldots, K} in the corpus excluding the current word \eqn{w_n}
//'   from the counts.
//'
//' @return Indicator for the topic draw from {1, 2, ..., K}.
uint16_t draw_zdn_slda_logit(double yd, const arma::vec& w_d,
                             const arma::vec& eta,
                             uint32_t V, const arma::vec& ndk_n,
                             const arma::vec& nkm_n, const arma::vec& nk_n,
                             float alpha_, float gamma_) {

  // Warning: Overflow caused by probabilities near 0 handled by setting
  //   probability for the problematic topic to 1 / K;
  double muhat = arma::as_scalar(w_d.t() * eta);
  // Compute log-likelihood of y_d
  double loglike = get_ll_logit_yd(yd, muhat);

  arma::vec log_num = get_log_numer_samplez(V, ndk_n, nkm_n, nk_n,
                                            alpha_, gamma_) +
    loglike;
  uint16_t zdn = draw_zdn(log_num);
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
  arma::colvec mu_hat = w * eta;

  for (uint16_t d = 0; d < D; d++) {
    double yhat = Rcpp::rnorm(1, arma::as_scalar(mu_hat(d)), sigma2)(0);
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
  arma::colvec loglike_pred(D);
  arma::colvec mu_hat = w * eta;
  for (uint16_t d = 0; d < D; d++) {
    double phat = invlogit(arma::as_scalar(mu_hat(d)));
    uint16_t yhat = Rcpp::rbinom(1, 1, phat)(0);
    loglike_pred(d) = get_ll_logit_yd(yhat, arma::as_scalar(mu_hat(d)));
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
//' @param p_effd The contribution to the effective number of parameters from
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
//' @param iter The current iteration of the chain.
//' @param l_pred A m x D matrix of predictive likelihoods (NOT log-likelihoods).
//'
//' @return Vector of (1) WAIC for model, (2) standard error for WAIC, and (3)
//'   the effective number of parameters.
//' @export
// [[Rcpp::export]]
NumericVector waic_all(uint32_t iter, const arma::mat& l_pred) {

  NumericVector full_waic_se(3);
  double peff_sum = 0.0;
  double waic_sum = 0.0;
  const uint16_t D = l_pred.n_cols; // Number of observations
  arma::colvec waic(D);
  arma::colvec peff(D);
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
//' @param l_pred1 A m1 x D matrix of predictive likelihoods (NOT log-likelihoods) from model 1.
//' @param l_pred2 A m2 x D matrix of predictive likelihoods (NOT log-likelihoods) from model 2.
//'
//' @return A vector of (1) the difference in WAIC (on the deviance scale)
//'   between models and (2) the standard error of the difference in WAIC.
//' @export
// [[Rcpp::export]]
NumericVector waic_diff(const arma::mat& l_pred1, const arma::mat& l_pred2) {

  NumericVector diff_waic_se(2);
  double waic_diff_sum = 0.0;
  const uint32_t m1 = l_pred1.n_rows; // Length of first chain
  const uint32_t m2 = l_pred2.n_rows; // Length of second chain
  const uint16_t D = l_pred1.n_cols;  // Number of observations
  arma::colvec waic1(D);
  arma::colvec waic2(D);
  arma::colvec peff1(D);
  arma::colvec peff2(D);
  arma::colvec waic_diff(D);
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
    se_waic_diff += ((waic_diff(d) - mean_diff_waic) *
      (waic_diff(d) - mean_diff_waic));
  }
  se_waic_diff /= static_cast<float>(D - 1); // Variance of difference
  se_waic_diff *= static_cast<float>(D);     // Multiply variance of diff by D
  se_waic_diff = sqrt(se_waic_diff);         // SE(WAIC1 - WAIC2)

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

//' Log-likelihood for LDA model
//'
//' @param zdocs A D x max(\eqn{N_d}) matrix of topic indicators for all documents.
//' @param docs A D x max(\eqn{N_d}) matrix of word indicators for all documents.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//' @param N A vector of length D containing the number of words in each document.
//'
//' @return The current log-likelihood.
double get_ll_lda(const arma::mat& zdocs, const arma::mat& docs,
                  const arma::mat& theta, const arma::mat& beta,
                  const IntegerVector& docs_index, const NumericVector& N) {

  double ll_temp = 0.0;
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

  double ll_temp = get_ll_mlr(y, w, eta, sigma2) +
    get_ll_lda(zdocs, docs, theta, beta, docs_index, N);
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
double get_ll_slda_logit(const arma::colvec& y, const arma::mat& w,
                         const arma::colvec& eta,
                         const arma::mat& zdocs, const arma::mat& docs,
                         const arma::mat& theta, const arma::mat& beta,
                         const IntegerVector& docs_index, const NumericVector& N) {

 double ll_temp = get_ll_logit(y, w, eta) +
   get_ll_lda(zdocs, docs, theta, beta, docs_index, N);
 return ll_temp;
}

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

//' Log-posterior for LDA model
//'
//' @param ll A double of the current log-likelihood.
//' @param theta A D x K matrix of the current estimates of the document topic proportions.
//' @param beta a K x V matrix of the current estimates of the word-topic probabilities.
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions.
//' @param V The number of words in the vocabulary.
//' @param docs_index A vector of length D containing elements 1, 2, ..., D.
//'
//' @return The current log-posterior.
double get_lpost_lda(double ll, const arma::mat& theta, const arma::mat& beta,
                     double gamma_, double alpha_,
                     uint32_t V, const IntegerVector& docs_index) {

  uint16_t K = theta.n_cols;

  // Add prior on beta matrix
  double temp_betapost = 0;
  for (uint16_t k = 0; k < K; k++) {
    for (uint32_t v = 0; v < V; v++) {
      temp_betapost += log(beta(k, v));
    }
  }
  double lp_temp = ll + ((gamma_ - 1.0) * temp_betapost);

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

  double lp_temp = get_lpost_mlr(ll, eta, sigma2, mu0, sigma0, a0, b0);
  lp_temp += get_lpost_lda(lp_temp, theta, beta, gamma_, alpha_, V, docs_index);
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

  double lp_temp = get_lpost_eta(ll, eta, mu0, sigma0);
  lp_temp += get_lpost_lda(lp_temp, theta, beta, gamma_, alpha_, V, docs_index);
  return lp_temp;
}

//' Update number of times topic was drawn in document excluding current word
//'
//' @param topic The current topic index.
//' @param ndk A vector of the current number of draws of each topic in document d.
//'
//' @return A vector of the current number of draws of each topic in document d
//'   excluding word n.
arma::vec count_topicd(uint16_t topic, const arma::vec& ndk) {
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
//' @param word The current word index.
//' @param topic The current topic index.
//' @param ndk A vector of the current number of draws of each topic in document d.
//' @param nk A vector of the current number of draws of each topic in the corpus.
//' @param nkm A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus.
void update_zcounts(uint32_t d, uint32_t word, uint16_t topic,
                    arma::mat& ndk, NumericVector& nk, arma::mat& nkm) {

  arma::vec ndktemp = ndk.row(d).t();
  NumericVector nktemp = nk;
  arma::mat nkmtemp = nkm;

  ndk.row(d) = count_topicd(topic, ndktemp).t();
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
//' @return An object of class \code{Mlr}.
//' @export
//' @family Gibbs sampler
// [[Rcpp::export]]
S4 gibbs_mlr_cpp(uint32_t m, uint32_t burn, const arma::colvec& y,
                 const arma::mat& x,
                 const arma::colvec& mu0, const arma::mat& sigma0,
                 arma::colvec eta_start, float a0 = 0.001, float b0 = 0.001,
                 bool verbose = false, bool display_progress = false) {

  if (m <= burn) {
    stop("Length of chain m not greater than burn-in period.");
  }

  S4 results("Mlr"); // Create object slda of class Mlr

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
    Rcout << 1 << " eta: " << eta_start.t() << " ~~~~ sigma2: " << sigma2 << "\n";
  }

  arma::colvec eta(pp1);

  Progress p(m, display_progress);
  for (uint32_t i = 1; i <= m; i++) {

    // Draw eta
    try {
      eta = draw_eta_norm(x, y, sigma2, mu0, sigma0).t();
    } catch (std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() <<
        " while drawing eta vector\n";
    }

    // Draw sigma2
    try {
      sigma2 = draw_sigma2(a0, b0, x, y, eta);
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
        Rcout << i << " eta: " << eta.t() << " ~~~~ sigma2: " << sigma2 << "\n";
      }
    }
    if (display_progress) {
      p.increment();
    }
    Rcpp::checkUserInterrupt(); // Check to see if user cancelled sampler
  }

  // Compute WAIC and p_eff
  NumericVector waic_and_se = waic_all(m - burn, l_pred);

  results.slot("ndocs") = D;
  results.slot("nchain") = m - burn;
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
//'   regression coefficients, N(0, proposal_sd) (default: 2.38, ..., 2.38).
//' @param verbose Should parameter draws be output during sampling? (default:
//'   \code{FALSE}).
//' @param display_progress Should percent progress of sampler be displayed
//'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
//'   \code{display_progress} be set to \code{TRUE} at any given time.
//'
//' @return An object of class \code{Logistic}.
//' @export
//' @family Gibbs sampler
// [[Rcpp::export]]
S4 gibbs_logistic_cpp(uint32_t m, uint32_t burn, const arma::colvec& y,
                      const arma::mat& x,
                      const arma::colvec& mu0, const arma::mat& sigma0,
                      arma::colvec eta_start, arma::vec proposal_sd,
                      bool verbose = false, bool display_progress = false) {

  if (m <= burn) {
    stop("Length of chain m not greater than burn-in period.");
  }

  S4 results("Logistic"); // Create object slda of class Logistic

  const uint32_t D = y.size();
  const uint16_t pp1 = x.n_cols;

  arma::mat etam(m - burn, pp1);
  NumericVector loglike(m - burn); // Store log-likelihood (up to an additive constant)
  NumericVector logpost(m - burn); // Store log-posterior (up to an additive constant)
  arma::mat l_pred(m - burn, D);

  arma::vec attempt = arma::zeros(pp1);
  arma::vec accept = arma::zeros(pp1);
  arma::vec acc_rate(pp1);

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
        Rcout << i << " eta: " << eta.t() << "\n" <<
          "AR: " << acc_rate.t() << "\n" <<
          "Proposal SD: " << proposal_sd.t() << "\n";
      }
    }
    if (display_progress) {
      prog.increment();
    }
    Rcpp::checkUserInterrupt(); // Check to see if user cancelled sampler
  }

  // Compute WAIC and p_eff
  NumericVector waic_and_se = waic_all(m - burn, l_pred);

  results.slot("ndocs") = D;
  results.slot("nchain") = m - burn;
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

//' Collapsed Gibbs sampler for the sLDA-X model
//'
//' In general, don't use this directly. Instead, use \code{\link{gibbs_sldax}}.
//'
//' @include slda-class.R
//'
//' @param m The number of iterations to run the Gibbs sampler.
//' @param burn The number of iterations to discard as the burn-in period.
//' @param y A D x 1 vector of binary outcomes (0/1) to be predicted.
//' @param x A D x p matrix of additional predictors (no column of 1s for
//'   intercept).
//' @param docs A D x max(\eqn{N_d}) matrix of word indices for all documents.
//' @param w A D x V matrix of counts for all documents and vocabulary terms.
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
//' @param constrain_eta A logical (default = \code{TRUE}): If \code{TRUE}, the
//'   regression coefficients will be constrained so that they are in descending
//'   order; if \code{FALSE}, no constraints will be applied.
//' @param alpha_ The hyper-parameter for the prior on the topic proportions
//'   (default: 0.1).
//' @param gamma_ The hyper-parameter for the prior on the topic-specific
//'   vocabulary probabilities (default: 1.01).
//' @param proposal_sd The proposal standard deviation for drawing the
//'   regression coefficients, N(0, proposal_sd) (default: 0.2).
//' @param interaction_xcol The column number of the design matrix for the
//' additional predictors for which an interaction with the \eqn{K} topics is
//' desired (default: \eqn{-1L}, no interaction). Currently only supports a
//' single continuous predictor or a two-category categorical predictor
//' represented as a single dummy-coded column.
//' @param verbose Should parameter draws be output during sampling? (default:
//'   \code{FALSE}).
//' @param display_progress Should percent progress of sampler be displayed
//'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
//'   \code{display_progress} be set to \code{TRUE} at any given time.
//' @export
// [[Rcpp::export]]
S4 gibbs_sldax_cpp(const arma::mat& docs,
                   const arma::mat& w,
                   uint32_t m, uint32_t burn,
                   uint16_t K, uint8_t model,
                   const arma::colvec& y,
                   const arma::mat& x,
                   const arma::colvec& mu0,
                   const arma::mat& sigma0,
                   float a0, float b0,
                   arma::colvec eta_start,
                   arma::vec proposal_sd,
                   int interaction_xcol = -1,
                   float alpha_ = 0.1, float gamma_ = 1.01,
                   bool constrain_eta = true,
                   bool verbose = false, bool display_progress = false) {

  if (m <= burn) {
    stop("Length of chain m not greater than burn-in period.");
  }

  const uint8_t lda = 1;
  const uint8_t slda = 2;
  const uint8_t sldax = 3;
  const uint8_t slda_logit = 4;
  const uint8_t sldax_logit = 5;

  S4 results("Sldax");

  NumericVector sigma2m(m - burn);
  double sigma2 = var(y) / 2.0;

  const uint32_t D = docs.n_rows;
  const uint32_t V = w.n_cols;
  NumericVector N(D);
  const IntegerVector topics_index = seq_len(K);
  const IntegerVector docs_index   = seq_len(D) - 1;
  for (uint32_t d : docs_index) N(d) = sum(docs.row(d) > 0);
  const uint32_t maxNd = max(N);
  arma::mat theta(D, K);
  arma::mat beta(K, V);
  // Topic draw counts
  arma::mat ndk(D, K);
  // Counts of topic-word co-occurences in corpus (K x V)
  arma::mat nkm(K, V);
  // Topic draws for all words and docs
  NumericVector nk(K);
  arma::mat zdocs = arma::zeros(D, maxNd);
  arma::mat zbar(D, K);

  // D x max(N_d) x m array to store topic draws
  arma::cube topicsm = arma::zeros(D, maxNd, m - burn);
  NumericVector loglike(m - burn); // Store log-likelihood (up to an additive constant)
  NumericVector logpost(m - burn); // Store log-posterior (up to an additive constant)

  // Randomly assign topics
  NumericVector init_topic_probs(K);
  for (uint16_t k = 0; k < K; k++)
    init_topic_probs(k) = 1.0 / static_cast<float>(K);

  for (uint32_t d : docs_index) {
    for (uint32_t n = 0; n < N(d); n++) {
      zdocs(d, n) = RcppArmadillo::sample(
        topics_index, 1, true, init_topic_probs)(0);
    }
    for (uint16_t k = 0; k < K; k++) {
      // Count topic draws in each document
      ndk(d, k) = sum(zdocs.row(d) == (k + 1));
      // Compute topic empirical proportions
      if (model != lda)
        zbar(d, k) = static_cast<double>(ndk(d, k)) / static_cast<double>(N(d));
    }
  }
  try {
    nkm = count_topic_word(K, V, zdocs, docs);
  } catch(std::exception& e) {
    Rcerr << "Runtime error: " << e.what() <<
      " while computing topic-word co-occurrences\n";
  }
  for (uint16_t k = 0; k < K; k++) nk(k) = sum(ndk.col(k));

  // Estimate theta
  for (uint32_t d : docs_index) {
    try {
      theta.row(d) = est_thetad(ndk.row(d).t(), alpha_).t();
    } catch(std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() <<
        " when estimating theta vector for document " << d << "\n";
    }
  }

  // Estimate beta
  for (uint16_t k = 0; k < K; k++) {
    try {
      beta.row(k) = est_betak(nkm.row(k).t(), gamma_).t();
    } catch(std::exception& e) {
      Rcerr << "Runtime Error: " << e.what() << " estimating row " << k <<
        "of beta matrix\n";
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
  arma::mat etam(m - burn, q);
  arma::colvec eta(q);
  arma::colvec etac(q);

  // For supervised models sLDA/sLDAX, set up regression coefficients
  if (constrain_eta) {
    if (model == sldax || model == sldax_logit) {
      // Constrain starting values of eta to be in descending order
      if (interaction_xcol > 0) {
        // Only sort topic "linear effects", not interactions
        std::partial_sort(eta_start.begin() + p - K + 1,
                          eta_start.end() - K + 1,
                          eta_start.end(), std::greater<float>());
      } else {
        std::sort(eta_start.begin() + p, eta_start.end(),
                  std::greater<float>());
      }
    } else if (model == slda || model == slda_logit) {
      std::sort(eta_start.begin(), eta_start.end(), std::greater<float>());
    }
  }

  if (verbose) Rcout << 1 << " eta: " << eta_start.t() << "\n";
  eta = eta_start;

  arma::vec attempt = arma::zeros(q);
  arma::vec accept = arma::zeros(q);
  arma::vec acc_rate(q);

  arma::mat r = zbar;
  // Predictive posterior likelihood of y
  arma::mat l_pred(m - burn, D);
  arma::mat xz_int(D, K - 1); // Need K - 1 interaction columns
  if (model == sldax || model == sldax_logit) {
    r = join_rows(x, zbar); // Full predictor matrix `r`
    if (interaction_xcol > 0) {
      for (int k = 0; k < (K - 1); k++) { // Need K - 1 interaction columns
        xz_int.col(k) = x.col(interaction_xcol - 1) % zbar.col(k);
      }
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
        update_zcounts(d, word, topic, ndk, nk, nkm);
        try {
          if (model == lda) {
            topic = draw_zdn_lda(V, ndk.row(d).t(), nkm.col(word - 1), nk,
                                 alpha_, gamma_);
          } else if (model == slda || model == sldax) {
            topic = draw_zdn_slda_norm(y(d), r.row(d).t(), eta, sigma2, V,
                                       ndk.row(d).t(), nkm.col(word - 1), nk,
                                       alpha_, gamma_);
          } else if (model == slda_logit || model == sldax_logit) {
            topic = draw_zdn_slda_logit(y(d), r.row(d).t(), eta, V,
                                        ndk.row(d).t(), nkm.col(word - 1), nk,
                                        alpha_, gamma_);
          } else {
            Rcerr << "Invalid model specified. Exiting.\n";
          }
        } catch(std::exception&  e) {
          Rcerr << "Runtime Error: " << e.what() <<
            " occurred while drawing topic for word " << n << " in document "
            << d << "\n";
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

    for (uint32_t d: docs_index) {
      if (model != lda) {
        // Get zbar matrix
        for (uint16_t k = 0; k < K; k++) {
          zbar(d, k) = static_cast<double>(ndk(d, k)) / static_cast<double>(N(d));
        }
      }
      // Estimate theta for doc d
      try {
        theta.row(d) = est_thetad(ndk.row(d).t(), alpha_).t();
      } catch(std::exception& e) {
        Rcerr << "Runtime Error: " << e.what() <<
          " when estimating theta vector for document " << d << "\n";
      }
    }

    if (model == sldax || model == sldax_logit) {
      r = join_rows(x, zbar); // Full predictor matrix `r`
      if (interaction_xcol > 0) {
        for (int k = 0; k < (K - 1); k++) { // Need K - 1 interaction columns
          xz_int.col(k) = x.col(interaction_xcol - 1) % zbar.col(k);
        }
        r = join_rows(r, xz_int);
      }
    } else {
      r = zbar;
    }

    // Estimate beta
    for (uint16_t k = 0; k < K; k++) {
      try {
        beta.row(k) = est_betak(nkm.row(k).t(), gamma_).t();
      } catch(std::exception& e) {
        Rcerr << "Runtime Error: " << e.what() << " estimating row " << k <<
          "of beta matrix\n";
      }
    }

    // Draw eta
    if (model != lda) {
      bool eta_order = false;
      uint16_t iter = 0;
      uint16_t max_iter = 1; // Only draw once
      if (constrain_eta) {
        max_iter = 1000;   // Max repeated draws if trying to constrain eta
      }
      while (!eta_order & (iter < max_iter)) { // Runs at least once
        iter++;
        if (model == slda || model == sldax) {
          try {
            etac = draw_eta_norm(r, y, sigma2, mu0, sigma0).t();
          } catch (std::exception& e) {
            Rcerr << "Runtime Error: " << e.what() <<
              " while drawing eta vector\n";
          }
        } else if (model == slda_logit || model == sldax_logit) {
          try {
            etac = draw_eta_logit(r, y, eta, mu0, sigma0,
                                  proposal_sd, attempt, accept);
          } catch (std::exception& e) {
            Rcerr << "Runtime Error: " << e.what() <<
              " while drawing eta vector\n";
          }
        }
        if (constrain_eta) {
          if (interaction_xcol > 0) {
            for (uint16_t k = p - K + 2; k < p + 1; k++) {
              // Force eta components to be in descending order (first is largest) to resolve label switching of topics
              eta_order = etac(k - 1) >= etac(k);
              if (!eta_order) break;
            }
          } else {
            for (uint16_t k = p + 1; k < q; k++) {
              // Force eta components to be in descending order (first is largest)
              //   to resolve label switching of topics
              eta_order = etac(k - 1) >= etac(k);
              if (!eta_order) break;
            }
          }
        }
      }
      eta = etac; // New draw

      if (model == slda || model == sldax) {
        // Draw sigma2
        try {
          sigma2 = draw_sigma2(a0, b0, r, y, eta);
        } catch (std::exception& e) {
          Rcerr << "Runtime Error: " << e.what() << " while drawing sigma2\n";
        }
      }
    }

    if (i > burn) {
      if (i == burn + 1 && burn > 0) Rcout << "Finished burn-in period\n";
      topicsm.slice(i - burn - 1) = zdocs;
      // Likelihood
      switch(model) {
        case lda:
          loglike(i - burn - 1) = get_ll_lda(zdocs, docs, theta, beta,
                                             docs_index, N);
          logpost(i - burn - 1) = get_lpost_lda(loglike(i - burn - 1),
            theta, beta, gamma_, alpha_, V, docs_index);
          break;
        case slda: // Fall through
        case sldax:
          loglike(i - burn - 1) = get_ll_slda_norm(
            y, r, eta, sigma2, zdocs, docs, theta, beta, docs_index, N);
          logpost(i - burn - 1) = get_lpost_slda_norm(loglike(i - burn - 1),
            eta, sigma2, theta, beta, mu0, sigma0, gamma_, alpha_, a0, b0, V,
            docs_index);
          l_pred.row(i - burn - 1) = post_pred_norm(r, y, eta, sigma2).t();
          sigma2m(i - burn - 1) = sigma2;
          break;
        case slda_logit: // Fall through
        case sldax_logit:
          loglike(i - burn - 1) = get_ll_slda_logit(
            y, r, eta, zdocs, docs, theta, beta, docs_index, N);
          logpost(i - burn - 1) = get_lpost_slda_logit(loglike(i - burn - 1),
            eta, theta, beta, mu0, sigma0, gamma_, alpha_, V, docs_index);
          l_pred.row(i - burn - 1) = post_pred_logit(r, eta).t();
          break;
        default: Rcerr << "Invalid model specified. Exiting.\n";
      }
      if (model != lda) {
        etam.row(i - burn - 1) = eta.t();
      }
    }

    if (model == slda_logit || model == sldax_logit) {
      for (uint16_t j = 0; j < q; j++) {
        acc_rate(j) = static_cast<float>(accept(j)) / static_cast<float>(attempt(j));
      }
      for (uint16_t j = 0; j < q; j++) {
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
      if (i == burn / 5.0 && burn > 0) Rcout << "Finished tuning proposal distributions\n";
    }

    if (i % 500 == 0 && verbose) {
      switch(model) {
        case slda: // Fall through
        case sldax:
          Rcout << i << " eta: " << eta.t() <<
            " ~~~~ sigma2: " << sigma2 <<
            " ~~~~ zbar d1: " << zbar.row(0) << "\n";
          break;
        case slda_logit: // Fall through
        case sldax_logit:
          Rcout << i << " eta: " << eta.t() <<
            " ~~~~ zbar d1: " << zbar.row(0) <<
            " ~~~~ AR: " << acc_rate.t() << "\n" <<
            " ~~~~ Proposal SD: " << proposal_sd.t() << "\n";
          break;
        default:
          Rcout << i << " zbar d1: " << zbar.row(0) << "\n";
          break;
      }
    }
    if (display_progress) {
      prog.increment();
    }
    Rcpp::checkUserInterrupt(); // Check to see if user cancelled sampler
  }

  // Compute WAIC and p_eff
  NumericVector waic_and_se = waic_all(m - burn, l_pred);

  results.slot("ntopics") = K;
  results.slot("ndocs") = D;
  results.slot("nvocab") = V;
  results.slot("nchain") = m - burn;
  results.slot("topics") = topicsm;
  results.slot("alpha") = alpha_;
  results.slot("gamma") = gamma_;
  results.slot("loglike") = loglike;
  results.slot("logpost") = logpost;

  if (model != lda) {
    results.slot("eta") = etam;
    results.slot("mu0") = mu0;
    results.slot("sigma0") = sigma0;
    results.slot("eta_start") = eta_start;
    results.slot("waic") = waic_and_se(0);
    results.slot("se_waic") = waic_and_se(1);
    results.slot("p_eff") = waic_and_se(2);
    results.slot("lpd") = l_pred;

    if (model == slda || model == sldax) {
      results.slot("sigma2") = sigma2m;
      results.slot("a0") = a0;
      results.slot("b0") = b0;
    }

    if (model == slda_logit || model == sldax_logit)
      results.slot("proposal_sd") = proposal_sd;
  }

  return results;
}
