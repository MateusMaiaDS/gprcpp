#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
// [[Rcpp::depends(RcppEigen)]]
// The line above (depends) it will make all the dependcies be included on the file
#include "gprcpp_types.h"
using namespace Rcpp;
using namespace RcppEigen;
using namespace std;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;



// [[Rcpp::export]]
MatrixXd MtM(const MapMatd& M){
  int n(M.cols()); // Getting the number of columns from x

  return MatrixXd(n,n).setZero().selfadjointView<Lower>().
         rankUpdate(M.adjoint());
}




// [[Rcpp::export]]
MatrixXd A_solve(const MapMatd& M){
  int n(M.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  // return M.inverse(); // Solve the system Mx=I;
  return M.llt().solve(I);
}

// [[Rcpp::export]]
MatrixXd A_solve_B(const MapMatd& A, const MapMatd& B){
  int n(A.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  return A.llt().solve(I) * B; // Solve the system Mx=I;
}

// [[Rcpp::export]]
MatrixXd A_solve_B_simple(const MapMatd& A, const MapMatd& B){
  int n(A.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  return A.llt().solve(B) ; // Solve the system Mx=I;
}

// [[Rcpp::export]]
MatrixXd A_solve_B_simple_matrixXd(const MatrixXd& A, const MapMatd& B){
  int n(A.cols());
  MatrixXd I =  MatrixXd::Identity(n,n); // Creating a Identity matrix
  return A.llt().solve(B) ; // Solve the system Mx=I;
}


// Creating a function to create the squared sum matrix
// [[Rcpp::export]]
MatrixXd symm_distance_matrix(const MapMatd &A){
  int nrow(A.rows());
  double squared_norm_aux;
  MatrixXd D = MatrixXd(nrow,nrow).setZero();
  for(int i = 0; i < nrow; i++){
    for(int c = i+1; c<nrow; c++){
    squared_norm_aux = (A.row(i)-A.row(c)).squaredNorm();// Calculate the squared DISTANCE matrix
    D(i,c) =  squared_norm_aux ;
    D(c,i) = squared_norm_aux;
    }
  }
  return D;
}

// Creating a function to create the squared sum matrix
// [[Rcpp::export]]
MatrixXd distance_matrix(const MapMatd &A, const MapMatd &B){
  int nrow1(A.rows());
  int nrow2(B.rows());
  int ncol1(A.cols());
  int ncol2(B.cols());

  if (ncol1 != ncol2) {
    throw std::runtime_error("Incompatible number of dimensions");
  }

  MatrixXd D = MatrixXd(nrow1,nrow2).setZero();
  for(int i = 0; i < nrow1; i++){
    for(int c = 0; c < nrow2; c++){
      D(i,c) = (A.row(i) - B.row(c)).squaredNorm();
    }
  }
  return D;
}




// [[Rcpp::export]]
NumericMatrix symm_distance_matrix_old(NumericMatrix m1) {
  int nrow = m1.nrow();
  int ncol = m1.ncol();

  NumericMatrix out(nrow, nrow);

  for (int r = 0; r < nrow; r++) {
    out(r,r) = 0;
    for (int c = (r+1); c < nrow; c++) {
      double total = 0;
      for(int c_aux = 0; c_aux < ncol; c_aux++){
        total += pow( (m1(r,c_aux)-m1(c,c_aux)) ,2);
      }
      out(r, c) = total;
      out(c, r) = total;
    }
  }

  return out;
};

// Building the kernel matrix

// [[Rcpp::export]]
MatrixXd k_y_nugget(const MapMatd &A,
                                     const double phi,
                                     const double nu,
                                     const double nugget){
  int nrow(A.rows());
  double kernel_input;
  MatrixXd D = MatrixXd(nrow,nrow).setOnes();
  for(int i = 0; i < nrow; i++){
    D(i,i) = nu+nugget;
    for(int c = i+1; c<nrow; c++){
      kernel_input = nu*exp(-((A.row(i)-A.row(c)).squaredNorm())/(2*phi*phi));// Calculate the squared DISTANCE matrix
      D(i,c) =  kernel_input;
      D(c,i) = kernel_input;
    }
  }
  return D;
}

// [[Rcpp::export]]
MatrixXd k_A_B(const MapMatd &A,
                      const MapMatd &B,
                      const double phi,
                      const double nu,
                      const double nugget){
  int nrow(A.rows());
  int ncol(B.rows());

  double kernel_input;
  MatrixXd D = MatrixXd(nrow,ncol).setOnes();
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
        kernel_input = nu*exp(-((A.row(i)-B.row(j)).squaredNorm())/(2*phi*phi));// Calculate the squared DISTANCE matrix
        D(i,j) = kernel_input;
      }
  }

  return D;
}

// Get the GP-mean
// [[Rcpp::export]]
MatrixXd gp_mean(const MapMatd K_y_nug,
                const MapMatd K_A_B,
                const MapMatd y){

    return K_A_B.adjoint()*(A_solve_B_simple(K_y_nug,y));
    // return (A_solve_B_simple(K_y_nug,y));

}

// Get the GP-cov
// [[Rcpp::export]]
MatrixXd gp_cov( const MapMatd K_y_nug,
                        const MapMatd K_new,
                        const MapMatd K_A_B){

  return K_new-K_A_B.adjoint()*(A_solve_B_simple(K_y_nug,K_A_B));
}

// [[Rcpp::export]]
double get_log_D(MatrixXd X){
      VectorXd Dvec(X.ldlt().vectorD());
      return Dvec.array().log().sum();
}


//[[Rcpp::export]]
double phi_log_post(const MapMatd& X,
                    const MapMatd& y,
                    const double phi,
                    const double nu,
                    const double nugget){

  // Getting the covariance matrix
  MatrixXd K_y_value = k_y_nugget(X,phi,nu,nugget);
  return -0.5*get_log_D(2*3.1415926*K_y_value)-
    0.5*(y.adjoint()*A_solve_B_simple_matrixXd(K_y_value,y))(1,1);
}

//[[Rcpp::export]]
NumericVector phi_post_sample( const MapMatd X,
                 const MapMatd y,
                 const int n_mcmc,
                 const int n_burn,
                 const double nu,
                 const double nugget,
                 double phi_init = 0.1) {

  // Getting the n_post
  NumericVector phi_post_samples;
  double l_old, l_new, phi_new,acceptance; // Calculating the log_likelihoods

  for(int i=0;i<n_mcmc;i++){

    phi_new = R::runif(0,10);
    l_old = phi_log_post(X,y,phi_init,nu,nugget);
    l_new = phi_log_post(X,y,phi_new,nu,nugget);

    // Checking the acceptance
    acceptance = exp(l_new-l_old);

    if(R::runif(0,1)<=acceptance){
      phi_init = phi_new;
    }

    // Saving the post-burn samples
    if(i >=n_burn){
      phi_post_samples.push_back(phi_init);
    }

  }

  return phi_post_samples;

}



