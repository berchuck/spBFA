#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
bool Play(arma::vec X) {
  arma::vec Y = X.elem(arma::find_finite(X));
  return Y.size() == 0;
}


/*** R
Play(c(-Inf, -Inf))
*/
