#include <RcppArmadillo.h>

//' @title Get the top number of observations
//'
//' @name table_cpp
//' @param v A row vector of values to be tabulated.
//'
//' @return A map of (values, frequencies)
//'
//' @noRd
std::map<uint32_t, uint32_t> table_cpp(const arma::urowvec& v) {

  // Create a map to store frequencies
  std::map<uint32_t, uint32_t> Elt;
  Elt.clear();

  // Fill the map with occurrences per number.
  for (arma::uword i = 0; i < v.size(); ++i) Elt[ v(i) ] += 1;

  return Elt;
}
