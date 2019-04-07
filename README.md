# tfHuber

Tuning-Free Huber Regression

## Goal of the package

This package implements Huber mean estimation, adaptive Huber regression and regularized Huber regression (Huber-Lasso) efficiently. For all the above methods, the robustness parameter tau is calibrated by a tuning-free principal.

Specifically, for Huber regression, assume that the observed data (Y, X) follow a linear model Y = beta_0 + X * beta + epsilon, where Y is an n-dimensional response vector, X is an n by d design matrix, and epsilon is an n-vector of noise variables whose distributions can be asymmetric and/or heavy-tailed. The package can compute the standard Huber's M-estimator if d < n and Huber-Lasso regression estimator, in particular, the intercept term beta_0 is estimated via a two-step procedure. See the reference paper Wang et al. (2018) for more details of the tuning-free framework.

## Installation

Install `tfHuber` from github:

```{r gh-installation, eval = FALSE}
install.packages("devtools")
library(devtools)
devtools::install_github("XiaoouPan/tfHuber")
library(tfHuber)
```

## Getting help

Help on the functions can be accessed by typing `?`, followed by function name at the R command prompt. 

For example, `?huberReg` will present a detailed documentation with inputs, outputs and examples of the function `huberReg`.

## Common error messages

The package `tfHuber` is implemented in `Rcpp` and `RcppArmadillo`, so the following error messages might appear when you first install it (we'll keep updating common error messages with feedback from users):

* Error: "...could not find build tools necessary to build tfHuber": For Windows you need Rtools, for Mac OS X you need to install Command Line Tools for XCode. See (https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites). 

* Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out of date. This is a common environment specific issue. 

    1. In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions here: http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/. 

    2. For >= R 3.4.* : download the installer from the here: https://gcc.gnu.org/wiki/GFortranBinaries#MacOS. Then run the installer.


## Functions

There are three functions in this package: 

* `huberMean`: Huber mean estimation. 
* `huberReg`: Adaptive Huber regression.
* `cvHuberLasso`: K-fold cross-validation for Huber-Lasso regression.

## Simple examples 

Fist we show an example for Huber mean estimation. We generate data from log-normal distribution, which is asymmetric and heavy-tailed. Then we estimate its mean by adaptive Huber loss.

```{r}
library(tfHuber)
n = 1000
X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
meanList = huberMean(X)
hMean = meanList$mu
```

Then we present an example for adaptive Huber regression. Here we generate data from linear model Y = X * theta + epsilon, where epsilon is from log-normal distribution, and estimate intercept and coefficients by tuning-free Huber regression.

```{r}
library(tfHuber)
n = 500
d = 5
thetaStar = rep(3, d + 1)
X = matrix(rnorm(n * d), n, d)
error = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
Y = as.numeric(cbind(rep(1, n), X) %*% thetaStar + error)
listHuber = huberReg(X, Y)
thetaHuber = listHuber$theta
```

Finally we illustrate regularized Huber regression via an example. Again, we generate data from linear model Y = X * theta + epsilon, where X is a high-dimensional design matrix, and epsilon is from log-normal distribution. Then we estimate intercept and coefficients by Huber-Lasso regression, where the sparsity parameter lambda is calibrated by k-folds cross-validation, and robustness parameter tau is determined by a tuning-free procedure.

```{r}
library(tfHuber)
n = 200
d = 500
s = 5
thetaStar = c(rep(3, s + 1), rep(0, d - s))
X = matrix(rnorm(n * d), n, d)
error = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
Y = as.numeric(cbind(rep(1, n), X) %*% thetaStar + error)
listHuberLasso = cvHuberLasso(X, Y)
thetaHuberLasso = listHuberLasso$theta
```

## License

GPL (>= 2)

## Authors

Xiaoou Pan <xip024@ucsd.edu>, Wen-Xin Zhou <wez243@ucsd.edu> 

## Reference

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ Integration. J. Stat. Softw. 40(8) 1-18. (http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Comp. Stat. Dat. Ana. 71 1054-1063. (http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Liu, H., Sun, Q. and Zhang, T. (2018). I-LAMM for sparse learning: Simultaneous control of algorithmic complexity and statistical error. Ann. Stat. 46 814–841. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6007998/)

Pan, X., Sun, Q. and Zhou, W.-X. (2019). Nonconvex regularized robust regression with oracle properties in polynomial time. Preprint. (https://www.math.ucsd.edu/~wez243/NH.pdf).

Sanderson, C. and Curtin, R. (2016). Armadillo: a template-based C++ library for linear algebra. J. Open. Src. Softw. 1 26. (http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Sun, Q., Zhou, W.-X. and Fan, J. (2018). Adaptive Huber regression, J. Amer. Stat. Assoc, to appear. (https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1543124)

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. J. R. Stat. Soc. Ser. B. Stat. Methodol. 58 267–288. (https://www.jstor.org/stable/2346178?seq=1#metadata_info_tab_contents)

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint. (https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf)
