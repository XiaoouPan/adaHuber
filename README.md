# tfHuber

Tuning-Free Huber Estimation and Regression

## Description

This package implements the Huber mean estimator, Huber covariance matrix estimation, adaptive Huber regression and l<sub>1</sub>-regularized Huber regression (Huber-Lasso) estimators efficiently. For all these methods, the robustification parameter tau is calibrated by a tuning-free principle.

Specifically, for Huber regression, assume the observed data vectors (Y, X) follow a linear model Y = &theta;<sub>0</sub> + X &theta; + &epsilon;, where Y is an n-dimensional response vector, X is an n by d design matrix, and &epsilon; is an n-vector of noise variables whose distributions can be asymmetric and/or heavy-tailed. The package computes the standard Huber's M-estimator when d < n and the Huber-Lasso estimator when d > n. The vector of coefficients beta and the intercept term &beta;<sub>0</sub> are estimated successively via a two-step procedure. See [Wang et al., 2018](https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf) for more details of the two-step tuning-free framework.

## Installation

Install `tfHuber` from GitHub:

```r
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

* Error: "...could not find build tools necessary to build tfHuber": For Windows you need Rtools, for Mac OS X you need to install Command Line Tools for XCode. See [this link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for details. 

* Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out of date. This is a common environment specific issue. 

    1. In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

    2. For >= R 3.4.* : download the installer from [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS). Then run the installer.

## Functions

There are four functions in this package: 

* `huberMean`: Huber mean estimation.
* `huberCov`: Huber covariance matrix estimation.
* `huberReg`: Adaptive Huber regression.
* `cvHuberLasso`: K-fold cross-validated Huber-Lasso regression.

## Examples 

First, we present an example of Huber mean estimation. We generate data from a log-normal distribution, which is asymmetric and heavy-tailed. We estimate its mean by the tuning-free Huber mean estimator.

```r
library(tfHuber)
n = 1000
X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
meanList = huberMean(X)
hMean = meanList$mu
```

Then we present an example of Huber covariance matrix estimation. We generate data from t distribution with df = 3, which is heavy-tailed. We estimate its covariance matrix by the method proposed in [Ke et al., 2019](https://arxiv.org/abs/1811.01520).

```r
library(tfHuber)
n = 100
d = 50
X = matrix(rt(n * d, df = 3), n, d) / sqrt(3)
hubCov = huberCov(X)
```

Next, we present an example of adaptive Huber regression. Here we generate data from a linear model Y = X &theta; + &epsilon;, where &epsilon; follows a log-normal distribution, and estimate the intercept and coefficients by tuning-free Huber regression.

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

Finally, we illustrate the use of l<sub>1</sub>-regularized Huber regression. Again, we generate data from a linear model Y = X &theta; + &epsilon;, where &theta; is a high-dimensional vector, and &epsilon; is from a log-normal distribution. We estimate the intercept and coefficients by Huber-Lasso regression, where the regularization parameter &lambda; is calibrated by K-fold cross-validation, and the robustification parameter &tau; is chosen by a tuning-free procedure.

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

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ Integration. J. Stat. Softw. 40(8) 1-18. [Paper](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Comp. Stat. Dat. Ana. 71 1054-1063. [Paper](http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Liu, H., Sun, Q. and Zhang, T. (2018). I-LAMM for sparse learning: Simultaneous control of algorithmic complexity and statistical error. Ann. Statist. 46 814–841. [Paper](https://projecteuclid.org/euclid.aos/1522742437)

Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-Friendly Covariance Estimation for Heavy-Tailed Distributions: A Survey and Recent Results. Statis. Sci. To appear. [Paper](https://arxiv.org/abs/1811.01520)

Pan, X., Sun, Q. and Zhou, W.-X. (2019). Nonconvex regularized robust regression with oracle properties in polynomial time. Preprint. [Paper](https://arxiv.org/abs/1907.04027).

Sanderson, C. and Curtin, R. (2016). Armadillo: a template-based C++ library for linear algebra. J. Open. Src. Softw. 1(2) 26. [Paper](http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Sun, Q., Zhou, W.-X. and Fan, J. (2019). Adaptive Huber regression. J. Amer. Stat. Assoc. 0 1-12. [Paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1543124)

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. J. R. Stat. Soc. Ser. B. Stat. Methodol. 58 267–288. [Paper](https://www.jstor.org/stable/2346178?seq=1#metadata_info_tab_contents)

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A new principle for tuning-free Huber regression. Preprint. [Paper](https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf)
