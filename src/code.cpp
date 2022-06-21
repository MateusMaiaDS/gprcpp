#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace RcppEigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::Upper;
using Eigen::VectorXd;


// [[Rcpp::depends(RcppEigen)]]


// Type defining as standard template
typedef Map<MatrixXd> MapMatd; // Mapping a double matrix


// [[Rcpp::export]]
inline MatrixXd MtM(const MapMatd& M){
  int n(M.cols()); // Getting the number of columns from x

  return MatrixXd(n,n).setZero().selfadjointView<Lower>().
         rankUpdate(M.adjoint());
}


// [[Rcpp::export]]
inline MatrixXd AtB(const MapMatd& A,const MapMatd& B){
    return A.adjoint() * B;
}

// [[Rcpp::export]]
inline MatrixXd A_solve(const MapMatd& M){
  int n(M.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  // return M.inverse(); // Solve the system Mx=I;
  return M.llt().solve(I);
}

// [[Rcpp::export]]
inline MatrixXd A_solve_B(const MapMatd& A, const MapMatd& B){
  int n(A.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  return A.llt().solve(I) * B; // Solve the system Mx=I;
}

// [[Rcpp::export]]
inline MatrixXd A_solve_B_simple(const MapMatd& A, const MapMatd& B){
  int n(A.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  return A.llt().solve(B) ; // Solve the system Mx=I;
}







