#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//Function to get Lambda matrix-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetLambda(arma::mat const& Theta, arma::umat const& Xi, int K, int M) {
  arma::mat Lambda(M, K);
  for (arma::uword i = 0; i < M; i++) {
    for (arma::uword j = 0; j < K; j++) {
      arma::uword XiLoc = Xi(i, j);
      Lambda(i, j) = Theta(XiLoc, j);
    }
  }
  return Lambda;
}



//Function to get weights, w_jl(s_i)-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube GetWeights(arma::cube const& Alpha, int K, int M, int L) {
  arma::cube Weights(L, M, K);
  arma::cube UpperPhiAlpha(L, M, K);
  for (arma::uword j = 0; j < K; j++) {
    for (arma::uword i = 0; i < M; i++) {
      for (arma::uword l = 0; l < L; l ++) {
        UpperPhiAlpha(l, i, j) = UpperpnormRcpp(Alpha(l, i, j));
      }
    }
  }
  for (arma::uword j = 0; j < K; j++) {
    for (arma::uword i = 0; i < M; i++) {
      for (arma::uword l = 0; l < L; l ++) {
        if (l == 0) Weights(l, i, j) = pnormRcpp(Alpha(l, i, j));
        if (l > 0) Weights(l, i, j) = pnormRcpp(Alpha(l, i, j)) * arma::prod(UpperPhiAlpha.slice(j)(arma::span(0, l - 1), i));
      }
    }
  }
  return Weights;
}

  
  
//Function to get log weights, log(w_jl(s_i))-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube GetlogWeights(arma::cube const& Alpha, int K, int M, int L) {
  arma::cube logWeights(L, M, K);
  arma::cube logUpperPhiAlpha(L, M, K);
  for (arma::uword j = 0; j < K; j++) {
    for (arma::uword i = 0; i < M; i++) {
      for (arma::uword l = 0; l < L; l ++) {
        double AlphaLIJ = Alpha(l, i, j);
        logUpperPhiAlpha(l, i, j) = lUpperpnormRcpp(AlphaLIJ);
      }
    }
  }
  for (arma::uword j = 0; j < K; j++) {
    for (arma::uword i = 0; i < M; i++) {
      for (arma::uword l = 0; l < L; l ++) {
        double AlphaLIJ = Alpha(l, i, j);
        if (l == 0) logWeights(l, i, j) = lpnormRcpp(AlphaLIJ);
        if (l > 0) logWeights(l, i, j) = lpnormRcpp(AlphaLIJ) + arma::sum(logUpperPhiAlpha.slice(j)(arma::span(0, l - 1), i));
      }
    }
  }
  return logWeights;
}

  
  
//Matrix inverse using cholesky decomoposition for covariances-------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CholInv(arma::mat const& Cov) {
  return arma::inv_sympd(Cov);
}



//Matrix inverse of for 3x3 matrix-----------------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Inv3(arma::mat const& A) {
  arma::mat result = arma::mat(3, 3);
  double determinant = A(0, 0) * ( A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2) ) - A(0, 1) * ( A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0) ) + A(0, 2) * ( A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0) );
  double invdet = 1 / determinant;
  result(0, 0) =  ( A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2) ) * invdet;
  result(1, 0) = -( A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1) ) * invdet;
  result(2, 0) =  ( A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1) ) * invdet;
  result(0, 1) = -( A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0) ) * invdet;
  result(1, 1) =  ( A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) ) * invdet;
  result(2, 1) = -( A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2) ) * invdet;
  result(0, 2) =  ( A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1) ) * invdet;
  result(1, 2) = -( A(0, 0) * A(2, 1) - A(2, 0) * A(0, 1) ) * invdet;
  result(2, 2) =  ( A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1) ) * invdet;
  return result;
}



//Matrix inverse of for 2x2 matrix-----------------------------------------------------------------------
arma::mat Inv2(arma::mat const& A) {
  arma::mat result = arma::mat(2, 2);
  double determinant = ( A(0, 0) * A(1, 1) ) - ( A(0, 1) * A(0, 1) );
  double invdet = 1 / determinant;
  result(0, 0) =  A(1, 1) * invdet;
  result(0, 1) = -A(1, 0) * invdet;
  result(1, 0) = -A(0, 1) * invdet;
  result(1, 1) =  A(0, 0) * invdet;
  return result;
}



//Function for making an upper diagonal matrix symmetric-------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat makeSymm(arma::mat const& A) {
  return arma::symmatu(A);
}



//Function that checks numerical equality of two objects against a tolerance-----------------------------
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}
