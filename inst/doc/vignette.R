library(NIFA)
library(fastICA)
library(NMF)



P <- 2000
N <- 500
K <- 6
M <- 2


set.seed(1)
simu <- simulateData(ngenes = P, nsamples = N, k = K, M = M, sd = 0.1, colinear.sd = 2, shape = 5, scale = 1)
X <- tscale(simu$X)



K <- 8
M <- 3


svdres <- svd(X, nu = K, nv = K, LINPACK = T)
svd.res <- apply(abs(cor(svdres$v, t(simu$S))),2,max)
print(svd.res)


ICAres <- fastICA(X,n.comp = K)
ICA.res <- apply(abs(cor(t(ICAres$A), t(simu$S))),2,max)
print(ICA.res)


ICAres.def <- fastICA(X, n.comp = K,  alg.typ = "deflation")
ICA.res.def <- apply(abs(cor(t(ICAres.def$A), t(simu$S))),2,max)
print(ICA.res.def)


NMFres <- nmf(simu$X+abs(min(simu$X)), rank = K)
NMF.res <- apply(abs(cor(t(NMFres@fit@H), t(simu$S))),2,max)
print(NMF.res)

NIFAres <- NIFA(tscale(X), K = K, M = M, max.iter = 500, S_threshold = 6e-5, init = "sd", A.init = NULL, S.init = NULL, verbose = T, ref = t(simu$S), beta_expect_flag = NULL, L1.sd = NULL, L2.sd = NULL, b_noise_prior = 1)
NIFA.res <- apply(abs(cor(t(NIFAres$mu_S), t(simu$S))),2,max)
print(NIFA.res)



