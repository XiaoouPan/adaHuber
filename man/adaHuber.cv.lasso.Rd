% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaHuber.R
\name{adaHuber.cv.lasso}
\alias{adaHuber.cv.lasso}
\title{Cross-Validated Regularized Adaptive Huber Regression.}
\usage{
adaHuber.cv.lasso(
  X,
  Y,
  lambdaSeq = NULL,
  kfolds = 5,
  numLambda = 50,
  phi0 = 0.01,
  gamma = 1.2,
  epsilon = 0.001,
  iteMax = 500
)
}
\arguments{
\item{X}{A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates.}

\item{Y}{An \eqn{n}-dimensional response vector.}

\item{lambdaSeq}{(\strong{optional}) A sequence of candidate regularization parameters. If unspecified, a reasonable sequence will be generated.}

\item{kfolds}{(\strong{optional}) Number of folds for cross-validation. Default is 5.}

\item{numLambda}{(\strong{optional}) Number of \eqn{\lambda} values for cross-validation if \code{lambdaSeq} is unspeficied. Default is 50.}

\item{phi0}{(\strong{optional}) The initial quadratic coefficient parameter in the local adaptive majorize-minimize algorithm. Default is 0.01.}

\item{gamma}{(\strong{optional}) The adaptive search parameter (greater than 1) in the local adaptive majorize-minimize algorithm. Default is 1.2.}

\item{epsilon}{(\strong{optional}) A tolerance level for the stopping rule. The iteration will stop when the maximum magnitude of the change of coefficient updates is less than \code{epsilon}. Default is 0.001.}

\item{iteMax}{(\strong{optional}) Maximum number of iterations. Default is 500.}
}
\value{
An object containing the following items will be returned:
\describe{
\item{\code{coef}}{A \eqn{(p + 1)} vector of estimated sparse regression coefficients, including the intercept.}
\item{\code{lambdaSeq}}{The sequence of candidate regularization parameters.}
\item{\code{lambda}}{Regularization parameter selected by cross-validation.}
\item{\code{tau}}{The robustification parameter calibrated by the tuning-free principle.}
\item{\code{iteration}}{Number of iterations until convergence.}
\item{\code{phi}}{The quadratic coefficient parameter in the local adaptive majorize-minimize algorithm.}
}
}
\description{
Sparse regularized adaptive Huber regressionwith "lasso" penalty. The function implements a localized majorize-minimize algorithm with a gradient-based method. The regularization parameter \eqn{\lambda} is selected by cross-validation, and the robustification parameter \eqn{\tau} is determined by a tuning-free principle.
}
\examples{
n = 100; p = 200; s = 5
beta = c(rep(1.5, s + 1), rep(0, p - s))
X = matrix(rnorm(n * p), n, p)
err = rt(n, 2)
Y = cbind(rep(1, n), X) \%*\% beta + err 

fit.lasso = adaHuber.cv.lasso(X, Y)
beta.lasso = fit.lasso$coef
}
\references{
Pan, X., Sun, Q. and Zhou, W.-X. (2021). Iteratively reweighted l1-penalized robust regression. Electron. J. Stat., 15, 3287-3348.

Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115 254-265.

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
}
\seealso{
See \code{\link{adaHuber.lasso}} for regularized adaptive Huber regression with a specified \eqn{lambda}.
}
