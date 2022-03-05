# adaHuber

**Ada**ptive **Huber** Estimation and Regression

## Description

This package implements the Huber mean estimator, Huber covariance matrix estimation, adaptive Huber regression and *l<sub>1</sub>*-regularized Huber regression (Huber-Lasso) estimators efficiently. For all these methods, the robustification parameter *&tau;* is calibrated via a tuning-free principle.

Specifically, for Huber regression, assume the observed data vectors (*Y*, *X*) follow a linear model *Y = &theta;<sub>0</sub> + X &theta; + &epsilon;*, where *Y* is an *n*-dimensional response vector, *X* is an *n* &times; *d* design matrix, and *&epsilon;* is an *n*-vector of noise variables whose distributions can be asymmetric and/or heavy-tailed. The package computes the standard Huber's *M*-estimator when *d < n* and the Huber-Lasso estimator when *d > n*. The vector of coefficients *&theta;* and the intercept term *&theta;<sub>0</sub>* are estimated successively via a two-step procedure. See [Wang et al., 2020](https://www.math.ucsd.edu/~wez243/tfHuber.pdf) for more details of the two-step tuning-free framework.

## Recent update

2022-02-28

We are updating the package and will submit a new version to CRAN soon.

The most efficient implementation of three functions `huberMean`, `huberCov`, `huberReg` in this package have been merged into another `R` library `FarmTest`, which has a [CRAN binary release](https://CRAN.R-project.org/package=FarmTest). To avoid the annoying compiling issues caused by this source package, and experience faster and more stable computation, we recommend installing `FarmTest`.

## Installation

Install `tfHuber` from GitHub:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("XiaoouPan/tfHuber")
library(tfHuber)
```

## Common error messages

First of all, to avoid most unexpected error messages, it is **strongly** recommended to update `R` to version >= 3.6.1.

Besides, since the library `tfHuber` is coded in `Rcpp` and `RcppArmadillo`, when you first install it, the following two build tools are required:

1. Rtools for Windows OS or XCode Command Line Tools for Mac OS. See [this link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for details.

2. gfortran binaries: see [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS) for instructions.

`tfHuber` should be working well after these steps. Some common error messages along with their solutions are collected below, and we'll keep updating them based on users' feedback:

* Error: "...could not find build tools necessary to build FarmTest": Please see step 1 above.

* Error: "library not found for -lgfortran/..": Please see step 2 above.

## Functions

There are four functions in this package: 

* `huberMean`: Huber mean estimation.
* `huberCov`: Huber covariance matrix estimation.
* `huberReg`: Adaptive Huber regression.
* `cvHuberLasso`: *K*-fold cross-validated Huber-Lasso regression.

## Getting help

Help on the functions can be accessed by typing `?`, followed by function name at the R command prompt. 

For example, `?huberReg` will present a detailed documentation with inputs, outputs and examples of the function `huberReg`.

## Examples 

First, we present an example of Huber mean estimation. We generate data from a log-normal distribution, which is asymmetric and heavy-tailed. We estimate its mean by the tuning-free Huber mean estimator.

```r
library(tfHuber)
n = 1000
X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
meanList = huberMean(X)
hMean = meanList$mu
```

Then we present an example of Huber covariance matrix estimation. We generate data from *t* distribution with df = 3, which is heavy-tailed. We estimate its covariance matrix by the method proposed in [Ke et al., 2019](https://arxiv.org/abs/1811.01520).

```r
library(tfHuber)
n = 100
d = 50
X = matrix(rt(n * d, df = 3), n, d) / sqrt(3)
hubCov = huberCov(X)
```

Next, we present an example of adaptive Huber regression. Here we generate data from a linear model *Y = X &theta; + &epsilon;*, where *&epsilon;* follows a log-normal distribution, and estimate the intercept and coefficients by tuning-free Huber regression.

```r
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

Finally, we illustrate the use of *l<sub>1</sub>*-regularized Huber regression. Again, we generate data from a linear model *Y = X &theta; + &epsilon;*, where *&theta;* is a high-dimensional vector, and *&epsilon;* is from a log-normal distribution. We estimate the intercept and coefficients by Huber-Lasso regression, where the regularization parameter *&lambda;* is calibrated by *K*-fold cross-validation, and the robustification parameter *&tau;* is chosen by a tuning-free procedure.

```r
library(tfHuber)
n = 100
d = 200
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

## Author(s)

Xiaoou Pan <xip024@ucsd.edu>, Wen-Xin Zhou <wez243@ucsd.edu> 

## References

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ integration. *J. Stat. Softw.* **40** 1-18. [Paper](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. *Comput. Statist. Data Anal.* **71** 1054-1063. [Paper](http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Liu, H., Sun, Q. and Zhang, T. (2018). I-LAMM for sparse learning: Simultaneous control of algorithmic complexity and statistical error. *Ann. Statist.* **46** 814–841. [Paper](https://projecteuclid.org/euclid.aos/1522742437)

Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions. *Statis. Sci.* **34** 454-471. [Paper](https://projecteuclid.org/euclid.ss/1570780979)

Pan, X., Sun, Q. and Zhou, W.-X. (2021). Iteratively reweighted l1-penalized robust regression. *Electron. J. Stat.* **15** 3287-3348. [Paper](https://doi.org/10.1214/21-EJS1862)

Sanderson, C. and Curtin, R. (2016). Armadillo: A template-based C++ library for linear algebra. *J. Open Source Softw.* **1** 26. [Paper](http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. *J. Amer. Stat. Assoc.* **115** 254-265. [Paper](https://doi.org/10.1080/01621459.2018.1543124)

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *J. R. Stat. Soc. Ser. B. Stat. Methodol.* **58** 267–288. [Paper](https://www.jstor.org/stable/2346178?seq=1#metadata_info_tab_contents)

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2020). A new principle for tuning-free Huber regression. *Stat. Sinica* to appear. [Paper](https://www.math.ucsd.edu/~wez243/tfHuber.pdf)
