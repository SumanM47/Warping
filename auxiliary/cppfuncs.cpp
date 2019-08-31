#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix submat(NumericMatrix X_, NumericVector ind_) {
  
  int n = X_.nrow(), k = X_.ncol();
  
  arma::mat X(X_.begin(), n, k, false);
  arma::uvec ind = as<arma::uvec>(ind_);
  arma::mat submat = X.rows(ind - 1);
  
  return wrap(submat);
}

//[[Rcpp::export]]
arma::cube subarr(arma::cube A, NumericVector ind_){
  arma::uvec ind = as<arma::uvec>(ind_);
  int nr = ind.size();
  int nc = A.n_cols;
  int ns = A.n_slices;
  arma::cube X(nr,nc,ns);
  for(int i=0;i<ns;i++){
    X.slice(i) = A.slice(i).rows(ind-1);
  }
  return(X);
}

//[[Rcpp::export]]
arma::vec get_sig2(arma::mat Xf, arma::mat C){
            
            int np = Xf.n_rows;
            
            arma::rowvec temp(np);
            arma::vec out(np);
            for(int i=0;i<np;i++){
            temp = Xf.row(i);
            out(i) = arma::sum((temp*C)%temp);
            }
            return(out);
            }