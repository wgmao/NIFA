\name{NIFA}
\alias{NIFA}
\title{Main NIFA function}
\usage{
NIFA(X, K = 6, M = 4, init="sd", S.init=NULL, A.init=NULL, L1.sd = NULL, L2.sd = NULL , verbose=F, beta_threshold = 2e-5, 
     S_threshold = 6e-5, max.iter = 1000, ref=NULL, rho_S_diff=1.5,phi_S_prior=1, a_S_prior=1,b_S_prior=1, a_noise_prior=1, 
     b_noise_prior=1,lambda_A_prior=1, eta_A_prior=0, ELBO=F, A.immediate =F, beta_expect_flag = NULL)
}
\arguments{
\item{X}{the data to be processed with genes in rows and samples in columns.}
\item{K}{number of latent factors (default: 6)}
\item{M}{maximum number of mixtures associated with each latent factor (default: 4)}
\item{init}{initialization method and there are three options: "random", "ICA" and "sd" (default: "sd"-simple decomposition)}
\item{S.init}{initialization values for latent factors S (default: NULL)}
\item{A.init}{initialization values for loadings A (default: NULL)}
\item{L1.sd}{regularization parameter L1 in simple decomposition (default: NULL)}
\item{L2.sd}{regularization parameter L2 in simple decomposition (default: NULL)}
\item{verbose}{print out intermediate messages (default: F)}
\item{beta_threshold}{threshold for relateive change of beta (default: 2e-5)}
\item{S_threshold}{threshold for relative change of S (default: 6e-5)}
\item{max.iter}{maximum number of iterations (default: 1000)}
\item{ref}{independent measurement that can be associated with S, e.g.,  one_hot_encode(celltype) (default: NULL)}
\item{rho_S_diff}{hyperparameter to initialize rho (default: 1.5)}
\item{phi_S_prior}{hyperparameter to initialize phi (default: 1)}
\item{a_S_prior}{prior parameter for S (default: 1)}
\item{b_S_prior}{prior parameter for S (default: 1)}
\item{a_noise_prior}{prior paramter for beta (default: 1)}
\item{b_noise_prior}{prior parameter for beta (default: 1, recommended: prod(dim(X))*5)}
\item{lambda_A_prior}{prior parameter for A (default: 1)}
\item{eta_A_prior}{prior parameter for A (default: 0)}
\item{ELBO}{whether calculate ELBO or not (default: F)}
\item{A.immediate}{where update A till convergence within each iteration (default: F, recommended: F)}
\item{beta_expect_flag}{you can specify a fixed beta_expect by assigning a value, e.g., 1e-2 (default: NULL)}

}
\value{
\item{X}{same as input}
\item{K}{same as input}
\item{M}{same as input}
\item{beta_threshold}{same as input}
\item{S_threshold}{same as input}
\item{rho_S_diff}{same as input}
\item{phi_S_prior}{same as input}
\item{a_S_prior}{same as input}
\item{b_S_prior}{same as input}
\item{a_noise_prior}{same as input}
\item{b_noise_prior}{same as input}
\item{lambda_A_prior}{same as input}
\item{eta_A_prior}{same as input}
\item{loop}{number of iterations before converge}
\item{S_expect}{expectation of latent factors S, see the matching table}
\item{S_2_expect}{see the matching table}
\item{mu_S_expect}{see the matching table}
\item{epsilon_expect}{see the matching table}
\item{sigma_log_S_expect}{see the matching table}
\item{beta_expect}{see the matching table}
\item{sigma_expect}{see the matching table}
\item{mu_S_2_expect}{see the matching table}
\item{mean_A_expect}{expectation of loading A, see the matching table}
\item{var_A_expect}{see the matching table}
\item{sigma_S}{see the matching table}
\item{mu_S}{see the matching table}
\item{lambda_S}{see the matching table}
\item{pi_S}{see the matching table}
\item{rho_S}{see the matching table}
\item{phi_S}{see the matching table}
\item{a_S}{see the matching table}
\item{b_S}{see the matching table}
\item{a_noise}{see the matching table}
\item{b_noise}{see the matching table}
\item{lambda_A}{see the matching table}
\item{eta_A}{see the matching table}
\item{mu_A}{see the matching table}
\item{sigma_A}{see the matching table}
\item{rho_S_0}{initial values of rho_S}
\item{phi_S_0}{initial values of phi_S}
\item{a_S_0}{initial values of a_S}
\item{b_S_0}{initial values of b_S}
\item{lambda_S_0}{initial values of lambda_S}
\item{pi_S_0}{initial values of pi_S}
\item{a_noise_0}{initial values of a_noise}
\item{b_noise_0}{initial values of b_noise}
\item{lambda_A_0}{initial values of lambda_A}
\item{eta_A_0}{initial values of eta_A}
\item{mu_S_change}{relative change of mu_S, normF(mu_S-mu_S_old)/normF(mu_S_old)}
\item{beta_change}{relative change of beta}
\item{ELBO_accu}{ELBO value}
}
\examples{
NIFA.res <- NIFA(tscale(rnaseq.norm), K = 25, M = 4, max.iter = 500, S_threshold = 6e-5, init = "sd", 
                 ref = one_hot_encode(celltype), b_noise_prior = prod(dim(rnaseq.norm))*5)
}
\description{
Main NIFA function
}
