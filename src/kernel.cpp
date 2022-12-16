// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef MSpMat::InnerIterator InIterMat;

// the sparse matrix is of class "dgCMatrix"
// number of rows = number of points of interest
// number of columns = number of cells
// each element of the input sparse matrix is the ratio of distance to bandwidth

// [[Rcpp::export("GaussianKernelDensity")]]
NumericVector GaussianKernelDensity(Eigen::MappedSparseMatrix<double> X){
  int ncol = X.cols();
  NumericVector ret(ncol, 0.0); //vector of size ncol filled with 0.0
  for (int j = 0; j < ncol; j++) {
    // for each non-zero element in the j column
    double Ai = 0.0; // Ai has the initial value 0
    for (InIterMat i_(X, j); i_; ++i_) {
      // sum of each column
      Ai += exp(-0.5*pow(i_.value(), 2)); // apply the Gaussian Kernel
    };
  ret[j] = Ai;
  };
  return ret;
}

// [[Rcpp::export("UniformKernelDensity")]]
NumericVector UniformKernelDensity(Eigen::MappedSparseMatrix<double> X){
  int ncol = X.cols();
  NumericVector ret(ncol, 0.0); //vector of size ncol filled with 0.0
  for (int j = 0; j < ncol; j++) {
    // for each non-zero element in the j column
    double Ai = 0.0; // Ai has the initial value 0
    for (InIterMat i_(X, j); i_; ++i_) {
      // sum of each column
      Ai += i_.value(); // apply the Uniform Kernel
    };
    ret[j] = Ai;
  };
  return ret;
}

// [[Rcpp::export("QuarticKernelDensity")]]
NumericVector QuarticKernelDensity(Eigen::MappedSparseMatrix<double> X){
  int ncol = X.cols();
  NumericVector ret(ncol, 0.0); //vector of size ncol filled with 0.0
  for (int j = 0; j < ncol; j++) {
    // for each non-zero element in the j column
    double Ai = 0.0; // Ai has the initial value 0
    for (InIterMat i_(X, j); i_; ++i_) {
      // sum of each column
      Ai += pow(1.-pow(i_.value(), 2), 2); // apply the Quartic Kernel
    };
    ret[j] = Ai;
  };
  return ret;
}
// [[Rcpp::export("TriweightKernelDensity")]]
NumericVector TriweightKernelDensity(Eigen::MappedSparseMatrix<double> X){
  int ncol = X.cols();
  NumericVector ret(ncol, 0.0); //vector of size ncol filled with 0.0
  for (int j = 0; j < ncol; j++) {
    // for each non-zero element in the j column
    double Ai = 0.0; // Ai has the initial value 0
    for (InIterMat i_(X, j); i_; ++i_) {
      // sum of each column
      Ai += pow(1.-pow(i_.value(), 2), 3); // apply the triweight Kernel
    };
    ret[j] = Ai;
  };
  return ret;
}

NumericVector EpanechnikovKernelDensity(Eigen::MappedSparseMatrix<double> X){
  int ncol = X.cols();
  NumericVector ret(ncol, 0.0); //vector of size ncol filled with 0.0
  for (int j = 0; j < ncol; j++) {
    // for each non-zero element in the j column
    double Ai = 0.0; // Ai has the initial value 0
    for (InIterMat i_(X, j); i_; ++i_) {
      // sum of each column
      Ai += (1.-pow(i_.value(), 2)); // apply the Epanechnikov Kernel
    };
    ret[j] = Ai;
  };
  return ret;
}
