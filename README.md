# NIFA

R library dependencies
```
MASS
mclust
Rcpp
RcppArmadillo
```

To install the R package
```
library(devtools)
install_github("wgmao/NIFA")
```

The main function is `NIFA()` and there is a short [vignette](blob/master/inst/doc/vignette.pdf) based on a test dataset (simulated scRNA-seq) called `SimKumar4easy` which is publicly available via the bioconductor package `DuoClustering2018`.

# Parameters that have major effects on the result
There are five parameters that are more sensitive than others: `K`, `S_threshold`/`max.iter` and `b_noise_prior`/`beta_expect_flag`.
- `K`  number of latent factors. 
- `S_threshold` and `max.iter` control the number of iterations.
- `b_noise_prior` Based on experience, the recommendation is to set it as `prod(dim(X))*5)`. If the result doesn't look good, you can manually set up a fixed value for the noise parameter using `beta_expect_flag` This is highly data-dependent. 

# Matching Table
![Matching Table](https://github.com/wgmao/NIFA/blob/master/match.png)
