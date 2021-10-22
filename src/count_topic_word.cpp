#include <RcppArmadillo.h>

#include "table_cpp.h"

//' @title Count topic-word co-occurences in corpus.
//'
//' Computes topic-word co-occurence matrix for a corpus of \eqn{D} documents
//' with the maximum length of a document in the corpus equal to max(\eqn{N_d})
//' and a vocabulary of \eqn{V} unique terms in the corpus.
//'
//' Indices in `doc_topic` and `doc_word` where no word exists in the
//' document must be set to 0.
//'
//' @name count_topic_word
//' @param K The number of topics.
//' @param V The number of terms in the corpus vocabulary.
//' @param doc_topic A \eqn{D} x max(\eqn{N_d}) matrix of topic assignments for
//'   the corpus.
//' @param doc_word A \eqn{D} x max(\eqn{N_d}) matrix of words for corpus.
//'
//' @return A \eqn{K} x \eqn{V} matrix of topic-word co-occurence counts.
//'
//' @noRd
// [[Rcpp::export(.count_topic_word)]]
arma::mat count_topic_word(uint16_t K, uint32_t V,
                           const arma::umat& doc_topic,
                           const arma::umat& doc_word) {

  const uint32_t D = doc_topic.n_rows; // Number of documents
  if (K < 2) Rcpp::stop("number of topics must be at least 2");
  if (V < 2) Rcpp::stop("size of vocabulary V must be at least 2");
  const uint32_t maxnd = doc_word.n_cols; // Maximum document length
  if (doc_word.n_rows != D) Rcpp::stop("'doc_word' and 'doc_topic' must have the same number of rows");
  if (doc_topic.n_cols != maxnd)
    Rcpp::stop("'doc_topic' and 'doc_word' must have the same number of columns");

  const arma::ucolvec topics_index = arma::linspace<arma::ucolvec>(1, K, K);
  // Matrix to store number of topic/word co-occurences
  arma::mat topic_word_freq = arma::mat(K, V, arma::fill::zeros);

  for (uint16_t topic : topics_index) {
    arma::uvec ids = arma::find(doc_topic == topic);
    arma::urowvec word_ids = doc_word.elem(ids).t();
    // Update word frequencies for topic
    std::map<uint32_t, uint32_t> counts = table_cpp(word_ids);
    for (const auto &p : counts) {
      topic_word_freq(topic - 1, p.first - 1) += p.second;
    }
  } // Topics
  return topic_word_freq;
}
