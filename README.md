# adaHuber

**Ada**ptive **Huber** Estimation and Regression

## Description

This package implements the Huber-type estimator for mean, covariance matrix, regression and *l<sub>1</sub>*-regularized Huber regression (Huber-Lasso). For all these methods, the robustification parameter *&tau;* is calibrated via a tuning-free principle.

Specifically, for Huber regression, assume the observed data vectors (*Y*, *X*) follow a linear model *Y = &theta;<sub>0</sub> + X &theta; + &epsilon;*, where *Y* is an *n*-dimensional response vector, *X* is an *n* &times; *d* design matrix, and *&epsilon;* is an *n*-vector of noise variables whose distributions can be asymmetric and/or heavy-tailed. The package computes the standard Huber's *M*-estimator when *d < n* and the Huber-Lasso estimator when *d > n*. The vector of coefficients *&theta;* and the intercept term *&theta;<sub>0</sub>* are estimated successively via a two-step procedure. See [Wang et al., 2021](https://doi.org/10.5705/ss.202019.0045) for more details.

## Recent update

**2022-03-04**

Version 1.1 is submitted to CRAN.


## Installation

Install `adaHuber` from CRAN

```r
install.packages("adaHuber")
```

## Common error messages

* Error: Compilation failed (with messages involving lgfortran, clang, etc.). **Solution**: This is a compilation error of Rcpp-based source packages. It happens when we recently submit a new version to CRAN, but it usually takes 3-5 days to build the binary package. Please use an older version or patiently wait for 3-5 days and then install the updated version.

* Error: unable to load shared object.. Symbol not found: _EXTPTR_PTR. **Solution**: This issue is common in some specific versions of `R` when we load Rcpp-based libraries. It is an error in R caused by a minor change about `EXTPTR_PTR`. Upgrading R to 4.0.2 will solve the problem.


## Functions

There are five functions in this package: 

* `adaHuber.mean`: Adaptive Huber mean estimation.
* `adaHuber.cov`: Adaptive Huber covariance estimation.
* `adaHuber.reg`: Adaptive Huber regression.
* `adaHuber.lasso`: Adaptive Huber-Lasso regression.
* `adaHuber.cv.lasso`: Cross-validated adaptive Huber-Lasso regression.

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
