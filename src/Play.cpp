#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Play(int f, arma::cube Cube) {
  return Cube(arma::span::all, arma::span(f - 1, f - 1), arma::span::all);
}


/*** R
Trials <- array(dim = c(5, 2, 3))
Trials[ , , ] <- 1:30
Play(1, Trials)
*/
