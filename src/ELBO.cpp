//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double lower_bound_2(arma::cube epsilon_expect, arma::mat sigma_log_S_expect, arma::mat sigma_expect,arma::cube S_2_expect, arma::mat S_expect, arma::mat mu_S_expect, 
                     arma::mat mu_S_2_expect){
  double pi = 3.141592653589793238462643383280;
  double res = 0;
  int K = epsilon_expect.n_rows;
  int M = epsilon_expect.n_cols;
  int N = epsilon_expect.n_slices;
  
  for (int i = 0; i < K; ++i){
    for (int m = 0; m < M; ++m){
      for (int n = 0; n < N; ++n){
          res = res+epsilon_expect(i,m,n)*(-0.5*log(2*pi) + 0.5*sigma_log_S_expect(i,m)-0.5*sigma_expect(i,m)*(S_2_expect(i,i,n)
                    -2*S_expect(i,n)*mu_S_expect(i,m)+mu_S_2_expect(i,m) ));
      }//for n
    }//for m
  }//for i
  return(res);
}//lower_bound_2


// [[Rcpp::export]]
double lower_bound_3(arma::mat lambda_A_0, arma::cube var_A_expect, arma::mat eta_A_0){
  int P = var_A_expect.n_slices;
  int K = var_A_expect.n_rows;
  double res = P/2*accu(log(lambda_A_0));
  double val = 0;
  NumericVector prob(1);
  
  for (int j = 0; j < K; ++j){
    for (int i = 0; i < K; ++i){
      val = -eta_A_0(j,i)*sqrt(lambda_A_0(0,i));
      prob = R::pnorm(val,0,1,1,0);
      res = res-0.5*lambda_A_0(0,i)*(var_A_expect(i,i,j)-pow(eta_A_0(j,i),2))-log(1-prob(0));     
    }//for i
  }//for j
  return(res);
}//lower_bound_3


// [[Rcpp::export]]
double lower_bound_9(arma::mat lambda_A, arma::cube var_A_expect, arma::mat eta_A){
  int P = var_A_expect.n_slices;
  int K = var_A_expect.n_rows;
  double res = P/2*accu(log(lambda_A));
  double val = 0;
  NumericVector prob(1);
  
  for (int j = 0; j < K; ++j){
    for (int i = 0; i < K; ++i){
      val = -eta_A(j,i)*sqrt(lambda_A(0,i));
      prob = R::pnorm(val,0,1,1,0);
      res = res-0.5*lambda_A(0,i)*(var_A_expect(i,i,j)-pow(eta_A(j,i),2))-log(1-prob(0));     
    }//for i
  }//for j
  return(res);
}//lower_bound_3




