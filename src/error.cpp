#include <Rcpp.h>

// Function to handle errors
void error(std::string s) {
  throw std::runtime_error(s);
}
