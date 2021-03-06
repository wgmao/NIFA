# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

lower_bound_2 <- function(epsilon_expect, sigma_log_S_expect, sigma_expect, S_2_expect, S_expect, mu_S_expect, mu_S_2_expect) {
    .Call(`_NIFA_lower_bound_2`, epsilon_expect, sigma_log_S_expect, sigma_expect, S_2_expect, S_expect, mu_S_expect, mu_S_2_expect)
}

lower_bound_3 <- function(lambda_A_0, var_A_expect, eta_A_0) {
    .Call(`_NIFA_lower_bound_3`, lambda_A_0, var_A_expect, eta_A_0)
}

lower_bound_9 <- function(lambda_A, var_A_expect, eta_A) {
    .Call(`_NIFA_lower_bound_9`, lambda_A, var_A_expect, eta_A)
}

eta_A_update_fast <- function(dim_eta_A, lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X) {
    .Call(`_NIFA_eta_A_update_fast`, dim_eta_A, lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X)
}

eta_A_update_cycle_fast <- function(dim_eta_A, lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X) {
    .Call(`_NIFA_eta_A_update_cycle_fast`, dim_eta_A, lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X)
}

b_noise_update_fast <- function(b_noise, X, mean_A_expect, var_A_expect, S_2_expect, S_expect) {
    .Call(`_NIFA_b_noise_update_fast`, b_noise, X, mean_A_expect, var_A_expect, S_2_expect, S_expect)
}

mu_S_update_fast <- function(dim_mu_S, sigma_S, X, mu_S_expect, epsilon_expect, sigma_expect, beta_expect, mean_A_expect) {
    .Call(`_NIFA_mu_S_update_fast`, dim_mu_S, sigma_S, X, mu_S_expect, epsilon_expect, sigma_expect, beta_expect, mean_A_expect)
}

lambda_S_update_fast <- function(dim_lambda_S, sigma_log_S_expect, pi_S, sigma_expect, S_expect, S_2_expect, mu_S_expect, mu_S_2_expect) {
    .Call(`_NIFA_lambda_S_update_fast`, dim_lambda_S, sigma_log_S_expect, pi_S, sigma_expect, S_expect, S_2_expect, mu_S_expect, mu_S_2_expect)
}

