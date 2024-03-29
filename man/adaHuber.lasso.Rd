% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaHuber.R
\name{adaHuber.lasso}
\alias{adaHuber.lasso}
\title{Regularized Adaptive Huber Regression}
\usage{
adaHuber.lasso(
  X,
  Y,
  lambda = 0.5,
  tau = 0,
  phi0 = 0.01,
  gamma = 1.2,
  epsilon = 0.001,
  iteMax = 500
)
}
\arguments{
\item{X}{A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates.}

\item{Y}{An \eqn{n}-dimensional response vector.}

\item{lambda}{(\strong{optional}) Regularization parameter. Must be positive. Default is 0.5.}

\item{tau}{(\strong{optional}) The robustness parameter. If not specified or the input value is non-positive, a tuning-free principle is applied. Default is 0 (hence, tuning-free).}

\item{phi0}{(\strong{optional}) The initial quadratic coefficient parameter in the local adaptive majorize-minimize algorithm. Default is 0.01.}

\item{gamma}{(\strong{optional}) The adaptive search parameter (greater than 1) in the local adaptive majorize-minimize algorithm. Default is 1.2.}

\item{epsilon}{(\strong{optional}) Tolerance level of the gradient-based algorithm. The iteration will stop when the maximum magnitude of all the elements of the gradient is less than \code{tol}. Default is 1e-03.}

\item{iteMax}{(\strong{optional}) Maximum number of iterations. Default is 500.}
}
\value{
An object containing the following items will be returned:
\describe{
\item{\code{coef}}{A \eqn{(p + 1)} vector of estimated sparse regression coefficients, including the intercept.}
\item{\code{tau}}{The robustification parameter calibrated by the tuning-free principle (if the input is non-positive).}
\item{\code{iteration}}{Number of iterations until convergence.}
\item{\code{phi}}{The quadratic coefficient parameter in the local adaptive majorize-minimize algorithm.}
}
}
\description{
Sparse regularized Huber regression models in high dimensions with \eqn{\ell_1} (lasso) penalty. The function implements a localized majorize-minimize algorithm with a gradient-based method.
}
\examples{
n = 200; p = 500; s = 10
beta = c(rep(1.5, s + 1), rep(0, p - s))
X = matrix(rnorm(n * p), n, p)
err = rt(n, 2)
Y = cbind(rep(1, n), X) \%*\% beta + err 

fit.lasso = adaHuber.lasso(X, Y, lambda = 0.5)
beta.lasso = fit.lasso$coef
}
\references{
Pan, X., Sun, Q. and Zhou, W.-X. (2021). Iteratively reweighted l1-penalized robust regression. Electron. J. Stat., 15, 3287-3348.

Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115 254-265.

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
}
\seealso{
See \code{\link{adaHuber.cv.lasso}} for regularized adaptive Huber regression with cross-validation.
}
