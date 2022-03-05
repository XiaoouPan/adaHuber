#' @title Adaptive Huber Mean Estimation
#' @description Adaptive Huber mean estimator from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' @param X An \eqn{n}-dimensional data vector.
#' @param epsilon (\strong{optional}) The tolerance level in the iterative estimation procedure, iteration will stop when \eqn{|\mu_new - \mu_old| < \epsilon}. The defalut value is 1e-4.
#' @param iteMax (\strong{optional}) Maximum number of iterations. Default is 500.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{mu}}{The Huber mean estimator.}
#' \item{\code{tau}}{The robustness parameter determined by the tuning-free principle.}
#' \item{\code{iteration}}{The number of iterations in the estimation procedure.}
#' }
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
#' @examples 
#' n = 1000
#' mu = 2
#' X = rt(n, 2) + mu
#' fit.mean = adaHuber.mean(X)
#' fit.mean$mu
#' @export 
adaHuber.mean = function(X, epsilon = 0.0001, iteMax = 500) {
  return (huberMeanList(X, epsilon, iteMax))
}

#' @title Adaptive Huber Covariance Estimation
#' @description Adaptive Huber covariance estimator from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' @details The observed data \eqn{X} is an \eqn{n} by \eqn{p} matrix. The distribution of each entry can be asymmetrix and/or heavy-tailed. The function outputs a robust estimator for the covariance matrix of \eqn{X}. For the input matrix \code{X}, both low-dimension (\eqn{p < n}) and high-dimension (\eqn{p > n}) are allowed.
#' @param X An \eqn{n} by \eqn{p} data matrix.
#' @param epsilon (\strong{optional}) The tolerance level in the iterative estimation procedure. The problem is converted to mean estimation, and the stopping rule is the same as \code{adaHuber.mean}. The defalut value is 1e-4.
#' @param iteMax (\strong{optional}) Maximum number of iterations. Default is 500.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{means}}{The Huber estimators for column means. A \eqn{p}-dimensional vector.}
#' \item{\code{cov}}{The Huber estimator for covariance matrix. A \eqn{p} by \eqn{p} matrix.}
#' }
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions. Statis. Sci., 34, 454-471.
#' @seealso \code{\link{adaHuber.mean}} for adaptive Huber mean estimation.
#' @examples 
#' n = 100
#' p = 5
#' X = matrix(rt(n * p, 3), n, p)
#' fit.cov = adaHuber.cov(X)
#' fit.cov$means
#' fit.cov$cov
#' @export 
adaHuber.cov = function(X, epsilon = 0.0001, iteMax = 500) {
  fit = huberCov(X, epsilon, iteMax)
  return (list(means = as.numeric(fit$means), cov = fit$cov))
}

#' @title Adaptive Huber Regression
#' @description Adaptive Huber regression from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param Y An \eqn{n}-dimensional response vector.
#' @param method (\strong{optional}) A character string specifying the method to calibrate the robustification parameter \eqn{\tau}. Two choices are "standard"(default) and "adaptive". See Wang et al.(2021) for details.
#' @param epsilon (\strong{optional}) Tolerance level of the gradient descent algorithm. The iteration will stop when the maximum magnitude of all the elements of the gradient is less than \code{tol}. Default is 1e-04.
#' @param iteMax (\strong{optional}) Maximum number of iterations. Default is 500.
#' @return An object containing the following items will be returned:
#' \describe{
#' \item{\code{coef}}{A \eqn{(p + 1)}-vector of estimated regression coefficients, including the intercept.}
#' \item{\code{tau}}{The robustification parameter calibrated by the tuning-free principle.}
#' \item{\code{iteration}}{Number of iterations until convergence.}
#' }
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115, 254-265.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
#' @examples
#' n = 200
#' p = 10
#' beta = rep(1.5, p + 1)
#' X = matrix(rnorm(n * p), n, p)
#' err = rt(n, 2)
#' Y = cbind(1, X) %*% beta + err
#' 
#' fit.huber = adaHuber.reg(X, Y, method = "standard")
#' beta.huber = fit.huber$coef
#' 
#' fit.adahuber = adaHuber.reg(X, Y, method = "adaptive")
#' beta.adahuber = fit.adahuber$coef
#' @export 
adaHuber.reg = function(X, Y, method = c("standard", "adaptive"), epsilon = 0.0001, iteMax = 500) {
  if (nrow(X) != length(Y)) {
    stop("Error: the length of Y must be the same as the number of rows of X.")
  }
  if (ncol(X) >= nrow(X)) {
    stop("Error: the number of columns of X cannot exceed the number of rows of X.")
  }
  method = match.arg(method)
  fit = NULL
  if (method == "standard") {
    fit = huberReg(X, Y, epsilon, constTau = 1.345, iteMax)
  } else {
    fit = adaHuberReg(X, Y, epsilon, iteMax)
  }
  return (list(coef = as.numeric(fit$coef), tau = fit$tau, iteration = fit$iteration))
}

#' @title Regularized Adaptive Huber Regression
#' @description Sparse regularized Huber regression models in high dimensions with \eqn{\ell_1} (lasso) penalty. The function implements a localized majorize-minimize algorithm with a gradient-based method.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. 
#' @param Y An \eqn{n}-dimensional response vector.
#' @param lambda (\strong{optional}) Regularization parameter. Must be positive. Default is 0.5.
#' @param tau (\strong{optional}) The robustness parameter. If not specified or the input value is non-positive, a tuning-free principle is applied. Default is 0 (hence, tuning-free).
#' @param phi0 (\strong{optional}) The initial quadratic coefficient parameter in the local adaptive majorize-minimize algorithm. Default is 0.01.
#' @param gamma (\strong{optional}) The adaptive search parameter (greater than 1) in the local adaptive majorize-minimize algorithm. Default is 1.2.
#' @param epsilon (\strong{optional}) Tolerance level of the gradient-based algorithm. The iteration will stop when the maximum magnitude of all the elements of the gradient is less than \code{tol}. Default is 1e-03.
#' @param iteMax (\strong{optional}) Maximum number of iterations. Default is 500.
#' @return An object containing the following items will be returned:
#' \describe{
#' \item{\code{coef}}{A \eqn{(p + 1)} vector of estimated sparse regression coefficients, including the intercept.}
#' \item{\code{tau}}{The robustification parameter calibrated by the tuning-free principle (if the input is non-positive).}
#' \item{\code{iteration}}{Number of iterations until convergence.}
#' \item{\code{phi}}{The quadratic coefficient parameter in the local adaptive majorize-minimize algorithm.}
#' }
#' @references Pan, X., Sun, Q. and Zhou, W.-X. (2021). Iteratively reweighted l1-penalized robust regression. Electron. J. Stat., 15, 3287-3348.
#' @references Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115 254-265.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
#' @seealso See \code{\link{adaHuber.cv.lasso}} for regularized adaptive Huber regression with cross-validation.
#' @examples 
#' n = 200; p = 500; s = 10
#' beta = c(rep(1.5, s + 1), rep(0, p - s))
#' X = matrix(rnorm(n * p), n, p)
#' err = rt(n, 2)
#' Y = cbind(rep(1, n), X) %*% beta + err 
#' 
#' fit.lasso = adaHuber.lasso(X, Y, lambda = 0.5)
#' beta.lasso = fit.lasso$coef
#' @export 
adaHuber.lasso = function(X, Y, lambda = 0.5, tau = 0, phi0 = 0.01, gamma = 1.2, epsilon = 0.001, iteMax = 500) {
  if (nrow(X) != length(Y)) {
    stop("Error: the length of Y must be the same as the number of rows of X.")
  }
  if (lambda <= 0) {
    stop("Error: lambda must be positive.")
  }
  fit = NULL
  if (tau <= 0) {
    fit = adaHuberLassoList(X, Y, lambda, phi0, gamma, epsilon, iteMax)
    tau = fit$tau
  } else {
    fit = huberLassoList(X, Y, lambda, tau, phi0, gamma, epsilon, iteMax)
  }
  return (list(coef = as.numeric(fit$coef), tau = tau, iteration = fit$iteration, phi = fit$phi))
}

#' @title Cross-Validated Regularized Adaptive Huber Regression.
#' @description Sparse regularized adaptive Huber regressionwith "lasso" penalty. The function implements a localized majorize-minimize algorithm with a gradient-based method. The regularization parameter \eqn{\lambda} is selected by cross-validation, and the robustification parameter \eqn{\tau} is determined by a tuning-free principle.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. 
#' @param Y An \eqn{n}-dimensional response vector.
#' @param lambdaSeq (\strong{optional}) A sequence of candidate regularization parameters. If unspecified, a reasonable sequence will be generated.
#' @param kfolds (\strong{optional}) Number of folds for cross-validation. Default is 5.
#' @param numLambda (\strong{optional}) Number of \eqn{\lambda} values for cross-validation if \code{lambdaSeq} is unspeficied. Default is 50.
#' @param phi0 (\strong{optional}) The initial quadratic coefficient parameter in the local adaptive majorize-minimize algorithm. Default is 0.01.
#' @param gamma (\strong{optional}) The adaptive search parameter (greater than 1) in the local adaptive majorize-minimize algorithm. Default is 1.2.
#' @param epsilon (\strong{optional}) A tolerance level for the stopping rule. The iteration will stop when the maximum magnitude of the change of coefficient updates is less than \code{epsilon}. Default is 0.001.
#' @param iteMax (\strong{optional}) Maximum number of iterations. Default is 500.
#' @return An object containing the following items will be returned:
#' \describe{
#' \item{\code{coef}}{A \eqn{(p + 1)} vector of estimated sparse regression coefficients, including the intercept.}
#' \item{\code{lambdaSeq}}{The sequence of candidate regularization parameters.}
#' \item{\code{lambda}}{Regularization parameter selected by cross-validation.}
#' \item{\code{tau}}{The robustification parameter calibrated by the tuning-free principle.}
#' \item{\code{iteration}}{Number of iterations until convergence.}
#' \item{\code{phi}}{The quadratic coefficient parameter in the local adaptive majorize-minimize algorithm.}
#' }
#' @references Pan, X., Sun, Q. and Zhou, W.-X. (2021). Iteratively reweighted l1-penalized robust regression. Electron. J. Stat., 15, 3287-3348.
#' @references Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115 254-265.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
#' @seealso See \code{\link{adaHuber.lasso}} for regularized adaptive Huber regression with a specified \eqn{lambda}.
#' @examples 
#' n = 100; p = 200; s = 5
#' beta = c(rep(1.5, s + 1), rep(0, p - s))
#' X = matrix(rnorm(n * p), n, p)
#' err = rt(n, 2)
#' Y = cbind(rep(1, n), X) %*% beta + err 
#' 
#' fit.lasso = adaHuber.cv.lasso(X, Y)
#' beta.lasso = fit.lasso$coef
#' @export 
adaHuber.cv.lasso = function(X, Y, lambdaSeq = NULL, kfolds = 5, numLambda = 50, phi0 = 0.01, gamma = 1.2, epsilon = 0.001, iteMax = 500) {
  if (nrow(X) != length(Y)) {
    stop("Error: the length of Y must be the same as the number of rows of X.")
  }
  if (!is.null(lambdaSeq) && min(lambdaSeq) <= 0) {
    stop("Error: all lambda's must be positive.")
  }
  n = nrow(X)
  if (is.null(lambdaSeq)) {
    lambdaMax = max(abs(t(X) %*% Y)) / n
    lambdaMin = 0.01 * lambdaMax
    lambdaSeq = exp(seq(log(lambdaMin), log(lambdaMax), length.out = numLambda))
  }
  folds = sample(rep(1:kfolds, ceiling(n / kfolds)), n)
  fit = cvAdaHuberLasso(X, Y, lambdaSeq, folds, kfolds, phi0, gamma, epsilon, iteMax)
  return (list(coef = as.numeric(fit$coef), lambdaSeq = lambdaSeq, lambda = fit$lambda, tau = fit$tau, iteration = fit$iteration, phi = fit$phi))
}
