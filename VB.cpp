// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat eta_A_update_fast(NumericVector dim_eta_A, mat lambda_A_0, mat eta_A_0, 
                         mat lambda_A, double beta_expect, mat S_expect, 
                         mat mean_A_expect, mat X){
  
  mat res(dim_eta_A[0], dim_eta_A[1],fill::zeros);
  mat out_prod(1,S_expect.n_cols, fill::zeros);
  int P = dim_eta_A[0];
  int K = dim_eta_A[1];
  int N = S_expect.n_cols;
  mat prod(P,N,fill::zeros);
  double num = 0;
  double tmp = 0;
  
  prod = mean_A_expect.t()*S_expect;
  
  for (int i = 0; i < P; ++i){
    for (int j = 0; j < K; ++j){
      num = 0;
      out_prod = mean_A_expect(j,i) * S_expect.row(j);
        
      for (int n = 0; n < N; ++n){
        num = num+(prod(i,n)-out_prod(0,n)-X(i,n))*S_expect(j,n);
      }//for n
      tmp = lambda_A_0(0,j)*eta_A_0(i,j)-beta_expect*num;
      res(i,j) = tmp/lambda_A(0,j);
    }//for j
  }//for i
    
  return res;
}//eta_A_update_fast


// [[Rcpp::export]]
mat eta_A_update_cycle_fast(NumericVector dim_eta_A, mat lambda_A_0, mat eta_A_0, 
                      mat lambda_A, double beta_expect, mat S_expect, 
                      mat mean_A_expect, mat X){
  
  mat res(dim_eta_A[0], dim_eta_A[1],fill::zeros);
  mat out_prod(1,S_expect.n_cols, fill::zeros);
  int P = dim_eta_A[0];
  int K = dim_eta_A[1];
  int N = S_expect.n_cols;
  mat prod(P,N,fill::zeros);
  double num = 0;
  double tmp = 0;
  double mu_A_cycle = 0;
  double Z = 0;
  double temp = 0;
  
  mat mean_A_expect_cycle = mean_A_expect;
  mat sigma = sqrt(1/lambda_A);
    
  prod = mean_A_expect.t()*S_expect;
  
  for (int i = 0; i < P; ++i){
    for (int j = 0; j < K; ++j){
      num = 0;
      out_prod = mean_A_expect(j,i) * S_expect.row(j);
      
      for (int n = 0; n < N; ++n){
        num = num+(prod(i,n)-out_prod(0,n)-X(i,n))*S_expect(j,n);
      }//for n
      tmp = lambda_A_0(0,j)*eta_A_0(i,j)-beta_expect*num;
      res(i,j) = tmp/lambda_A(0,j);
      
      //update mu_A
      Z = 1-R::pnorm(-res(i,j)/sigma(0,j), 0, 1, true, false);
      
      if (Z==0){
        mu_A_cycle = 0;
      }else{
        temp = R::dnorm(-res(i,j)/sigma(0,j),0, 1, false);
        mu_A_cycle = res(i,j)+temp/Z*sigma(0,j);
      }//else
      
      //update mean_A_expect_cycle
      mean_A_expect(j,i) = mu_A_cycle;
      
    }//for j
  }//for i
  
  return res;
}//eta_A_update_cycle_fast



// [[Rcpp::export]]
double b_noise_update_fast(double b_noise, mat X, mat mean_A_expect, cube var_A_expect, cube S_2_expect, mat S_expect){
  double res = 0;
  int P = X.n_rows;
  int N = X.n_cols;
  mat mean_A_expect_t(mean_A_expect.n_cols, mean_A_expect.n_rows,fill::zeros);
  //mat S_expect_t(S_expect.n_cols, S_expect.n_rows, fill::zeros);
  mat tmp(1,1,fill::zeros);
    
  mean_A_expect_t = mean_A_expect.t();
  //S_expect_t = S_expect.t();
  
  for (int j = 0; j < P; ++j){
    for (int n = 0; n < N; ++n){
      //tmp = -2*X(j,n)*mean_A_expect_t.row(j)*S_expect.col(n)+S_expect_t.row(n)*var_A_expect.slice(j)*S_expect.col(n);
      tmp = -2*X(j,n)*mean_A_expect_t.row(j)*S_expect.col(n)+trace( var_A_expect.slice(j)* S_2_expect.slice(n));
      res = res+0.5*( pow(X(j,n),2)+as_scalar(tmp));
    }//for n
  }//for j
  return(b_noise+res);
}//b_noise_update_fast


// [[Rcpp::export]]
mat mu_S_update_fast(NumericVector dim_mu_S, cube sigma_S, mat X, mat mu_S_expect, cube epsilon_expect,
                     mat sigma_expect, double beta_expect, mat mean_A_expect){
  mat res(dim_mu_S[0], dim_mu_S[1], fill::zeros);
  int K = mu_S_expect.n_rows;
  mat left(K,1,fill::zeros);
  mat right(K,1,fill::zeros);
  mat temp(sigma_expect.n_rows, sigma_expect.n_cols, fill::zeros);
  mat temp_t(sigma_expect.n_cols, sigma_expect.n_rows, fill::zeros);
  mat epsilon_expect_n(epsilon_expect.n_rows, epsilon_expect.n_cols, fill::zeros);
    
  temp = sigma_expect%mu_S_expect;
  temp_t = temp.t();
  
  for (int n = 0; n < dim_mu_S[1]; ++n){
    left = beta_expect*mean_A_expect*X.col(n);
    epsilon_expect_n = epsilon_expect.slice(n);
    
    for (int i = 0; i < dim_mu_S[0]; ++i){
      right(i,0) = as_scalar(epsilon_expect_n.row(i)*temp_t.col(i));
    }//for i
    res.col(n) = sigma_S.slice(n)*(left+right);
  }//for n
  return(res);
}//mu_S_update_fast



// [[Rcpp::export]]
cube lambda_S_update_fast(NumericVector dim_lambda_S, mat sigma_log_S_expect, mat pi_S, mat sigma_expect
                            ,mat S_expect, cube S_2_expect, mat mu_S_expect, mat mu_S_2_expect){
  //https://stackoverflow.com/questions/14947431/how-to-use-pi-in-rcppeigen
  double pi = 3.141592653589793238462643383280;
  double inf = 9999999999999999;
  double constant = 0;//-0.5*log(2*pi);
  cube res(dim_lambda_S[0], dim_lambda_S[1], dim_lambda_S[2], fill::zeros);
  cube res_tmp = res;
  
  for (int i = 0; i < dim_lambda_S[0]; ++i){
    for (int m = 0; m < dim_lambda_S[1]; ++m){
      for (int n = 0; n < dim_lambda_S[2]; ++n){
        res_tmp(i,m,n) = 0.5*sigma_log_S_expect(i,m)+log(pi_S(i,m))+constant-0.5*sigma_expect(i,m)*(S_2_expect(i,i,n)-2*S_expect(i,n)*mu_S_expect(i,m)+mu_S_2_expect(i,m));
      }//for n
    }//for j
  }//for i

  double sum_along_m = 0;
  for (int i = 0; i < dim_lambda_S[0]; ++i){
    for (int n = 0; n < dim_lambda_S[2]; ++n){
      sum_along_m = accu(exp(res_tmp.slice(n).row(i)));
      
      //if ( (sum_along_m==0) || (sum_along_m==R_PosInf)){
      if ( (sum_along_m==0) || (sum_along_m>inf)){
        for (int m = 0; m < dim_lambda_S[1]; ++m){
          res(i,m,n) = 1/accu(exp(res_tmp.slice(n).row(i)-res_tmp(i,m,n)));
        }//for m
      }else{
         for (int m = 0; m < dim_lambda_S[1]; ++m){
           res(i,m,n) = exp(res_tmp(i,m,n))/sum_along_m;
         }//for m
      }//else
    }//for m
  }//for i
  return(res);
}//lambda_S_update_fast




