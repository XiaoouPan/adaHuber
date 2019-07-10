# include <RcppArmadillo.h>
# include <cmath>
# include <string>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double f1(const double x, const arma::vec& resSq, const int n) {
  return arma::sum(arma::min(resSq, x * arma::ones(n))) / (n * x) - std::log(n) / n;
}

// [[Rcpp::export]]
double rootf1(const arma::vec& resSq, const int n, double low, double up, 
              const double tol = 0.00001, const int maxIte = 500) {
  int ite = 0;
  while (ite <= maxIte && up - low > tol) {
    double mid = (up + low) / 2;
    double val = f1(mid, resSq, n);
    if (val == 0) {
      return mid;
    } else if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return (low + up) / 2;
}

//' The function calculates adaptive Huber mean estimator from a data sample, with \eqn{\tau} determined by a tuning-free principle.
//'
//' The observed data are \eqn{X}, which is an \eqn{n}-dimensional vector whose distribution can be asymmetrix and/or heavy-tailed. The function outputs a robust estimator for the mean of \eqn{X}.
//'
//' @title Tuning-free Huber mean estimation
//' @param X An \eqn{n}-dimensional data vector.
//' @param epsilon The tolerance level in the iterative estimation procedure, iteration will stop when \eqn{|\mu_new - \mu_old| < \epsilon} or \eqn{|\tau_new - \tau_old| < \epsilon}. The defalut value is 1e-5.
//' @param iteMax The maximal number of iteration in the iterative estimation procedure, iteration will stop when this number is reached. The defalut value is 500.
//' @return A list including the following terms will be returned:
//' \itemize{
//' \item \code{mu} The Huber mean estimator.
//' \item \code{tau} The robustness parameter determined by the tuning-free principle.
//' \item \code{iteration} The number of iterations in the estimation procedure.
//' }
//' @author Xiaoou Pan, Wen-Xin Zhou
//' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint.
//' @examples
//' n = 1000
//' X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
//' meanList = huberMean(X)
//' hMean = meanList$mu
//' @export
// [[Rcpp::export]]
Rcpp::List huberMean(const arma::vec& X, const double epsilon = 0.00001, const int iteMax = 500) {
  int n = X.size();
  double muOld = 0;
  double muNew = arma::mean(X);
  double tauOld = 0;
  double tauNew = arma::stddev(X) * std::sqrt((long double)n / std::log(n));
  int iteNum = 0;
  while (((std::abs(muNew - muOld) > epsilon) || (std::abs(tauNew - tauOld) > epsilon)) && iteNum < iteMax) {
    muOld = muNew;
    tauOld = tauNew;
    arma::vec res = X - muOld * arma::ones(n);
    arma::vec resSq = arma::square(res);
    tauNew = std::sqrt((long double)rootf1(resSq, n, arma::min(resSq), arma::sum(resSq)));
    arma::vec w = arma::min(tauNew * arma::ones(n) / arma::abs(res), arma::ones(n));
    muNew = arma::as_scalar(X.t() * w) / arma::sum(w);
    iteNum++;
  }
  return Rcpp::List::create(Rcpp::Named("mu") = muNew, Rcpp::Named("tau") = tauNew, 
                            Rcpp::Named("iteration") = iteNum);
}

//' The function fits adaptive Huber regression via iterative weighted least square, with \eqn{\tau} determined by a tuning-free principle and the intercept term \eqn{\beta_0} estimated via a two-step procedure.
//'
//' The observed data are \eqn{(Y, X)}, where \eqn{Y} is an \eqn{n}-dimensional response vector and \eqn{X} is an \eqn{n} by \eqn{d} design matrix with \eqn{d < n}. We assume that \eqn{Y} depends on \eqn{X} through a linear model \eqn{Y = X \beta + \epsilon}, where \eqn{\epsilon} is an \eqn{n}-dimensional noise vector whose distribution can be asymmetrix and/or heavy-tailed. All the arguments except for \eqn{X} and \eqn{Y} have default settings.
//'
//' @title Tuning-free Huber regression
//' @param X An \eqn{n} by \eqn{d} design matrix with each row being a sample and each column being a variable and \eqn{d < n}.
//' @param Y A continuous response vector with length \eqn{n}.
//' @param epsilon The tolerance level for the iterative weighted least square, the iteration will stop when \eqn{||\theta_new - \theta_old||_inf < \epsilon} or \eqn{|\tau_new - \tau_old| < \epsilon}. The defalut value is 1e-5.
//' @param constTau The constant term used to update \eqn{\tau} in the tuning-free procedure. In each round of iteration, \eqn{\tau} is updated to be \code{constTau} \eqn{* \sigma_MAD}, where \eqn{\sigma_MAD = median(|R - median(R)|) / \Phi^(-1)(3/4)} is the median absolute deviation estimator, and \eqn{R} is the residual from last round of iteration. The defalut value is 1.345.
//' @param iteMax The maximal number of iteration in iterative weighted least square, the iteration stops if this number is reached. The defalut value is 500.
//' @return A list including the following terms will be returned:
//' \itemize{
//' \item \code{theta} The estimated \eqn{\theta}, a vector with length \eqn{d + 1}, with the first one being the value of intercept.
//' \item \code{tauCoef} The robustness parameter \eqn{\tau} determined by the tuning-free principle to estimate coefficients except for the intercept.
//' \item \code{tauItcp} The robustness parameter \eqn{\tau} determined by the tuning-free principle to estimate the intercept.
//' \item \code{iteCoef} The number of iterations in the iterative weighted least square procedure.
//' \item \code{iteItcp} The number of iterations to estimate the intercept.
//' }
//' @author Xiaoou Pan, Wen-Xin Zhou
//' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint.
//' @seealso \code{\link{cvHuberLasso}}
//' @examples
//' n = 500
//' d = 5
//' thetaStar = rep(3, d + 1)
//' X = matrix(rnorm(n * d), n, d)
//' error = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
//' Y = as.numeric(cbind(rep(1, n), X) %*% thetaStar + error)
//' listHuber = huberReg(X, Y)
//' thetaHuber = listHuber$theta
//' @export
// [[Rcpp::export]]
Rcpp::List huberReg(const arma::mat& X, const arma::vec& Y, const double epsilon = 0.00001, 
                        const double constTau = 1.345, const int iteMax = 500) {
  int n = X.n_rows;
  int d = X.n_cols;
  arma::mat Z = arma::ones(n, d + 1);
  Z.cols(1, d) = X;
  arma::vec thetaOld = arma::zeros(d + 1);
  arma::vec thetaNew = arma::solve(Z.t() * Z, Z.t() * Y);
  double tauOld = 0;
  double tauNew = std::sqrt((long double)arma::sum(arma::square(Y - Z * thetaNew)) / (n - d)) *
    std::sqrt((long double)n / std::log((long double)(d + std::log(n * d))));
  double mad = 0;
  int iteNum = 0;
  while ((arma::norm(thetaNew - thetaOld, "inf") > epsilon || std::abs(tauNew - tauOld) > epsilon) 
           && iteNum < iteMax) {
    thetaOld = thetaNew;
    tauOld = tauNew;
    arma::vec res = Y - Z * thetaOld;
    mad = arma::median(arma::abs(res - arma::median(res))) / 0.6744898;
    tauNew = constTau * mad;
    arma::mat WZ = Z;
    arma::vec WY = Y;
    for (int i = 0; i < n; i++) {
      double w = tauNew / std::abs(res(i));
      if (w < 1) {
        WZ.row(i) *= w;
        WY(i) *= w;
      }
    }
    thetaNew = arma::solve(Z.t() * WZ, Z.t() * WY);
    iteNum++;
  }
  Rcpp::List listMean = huberMean(Y - X * thetaNew.rows(1, d));
  thetaNew(0) = listMean["mu"];
  return Rcpp::List::create(Rcpp::Named("theta") = thetaNew, Rcpp::Named("tauCoef") = tauNew, 
                            Rcpp::Named("tauItcp") = listMean["tau"], Rcpp::Named("iteCoef") = iteNum,
                            Rcpp::Named("iteItcp") = listMean["iteration"]);
}

// [[Rcpp::export]]
int sgn(const double x) {
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
arma::vec softThresh(const arma::vec& x, const arma::vec& lambda) {
  return arma::sign(x) % arma::max(arma::abs(x) - lambda, arma::zeros(x.size()));
}

// [[Rcpp::export]]
arma::vec cmptLambda(const arma::vec& beta, const double lambda) {
  arma::vec rst = lambda * arma::ones(beta.size());
  rst(0) = 0;
  return rst;
}

// [[Rcpp::export]]
double loss(const arma::vec& Y, const arma::vec& Ynew, const std::string lossType,
            const double tau) {
  double rst = 0;
  if (lossType == "l2") {
    rst = arma::mean(arma::square(Y - Ynew)) / 2;
  } else if (lossType == "Huber") {
    arma::vec res = Y - Ynew;
    for (int i = 0; i < Y.size(); i++) {
      if (std::abs(res(i)) <= tau) {
        rst += res(i) * res(i) / 2;
      } else {
        rst += tau * std::abs(res(i)) - tau * tau / 2;
      }
    }
    rst /= Y.size();
  }
  return rst;
}

// [[Rcpp::export]]
arma::vec gradLoss(const arma::mat& X, const arma::vec& Y, const arma::vec& beta,
                   const std::string lossType, const double tau) {
  arma::vec res = Y - X * beta;
  arma::vec rst = arma::zeros(beta.size());
  if (lossType == "l2") {
    rst = -1 * (res.t() * X).t();
  } else if (lossType == "Huber") {
    for (int i = 0; i < Y.size(); i++) {
      if (std::abs(res(i)) <= tau) {
        rst -= res(i) * X.row(i).t();
      } else {
        rst -= tau * sgn(res(i)) * X.row(i).t();
      }
    }
  }
  return rst / Y.size();
}

// [[Rcpp::export]]
arma::vec updateBeta(const arma::mat& X, const arma::vec& Y, arma::vec beta, const double phi,
                     const arma::vec& Lambda, const std::string lossType, const double tau) {
  arma::vec first = beta - gradLoss(X, Y, beta, lossType, tau) / phi;
  arma::vec second = Lambda / phi;
  return softThresh(first, second);
}

// [[Rcpp::export]]
double cmptPsi(const arma::mat& X, const arma::vec& Y, const arma::vec& betaNew,
               const arma::vec& beta, const double phi, const std::string lossType,
               const double tau) {
  arma::vec diff = betaNew - beta;
  double rst = loss(Y, X * beta, lossType, tau)
    + arma::as_scalar((gradLoss(X, Y, beta, lossType, tau)).t() * diff)
    + phi * arma::as_scalar(diff.t() * diff) / 2;
    return rst;
}

// [[Rcpp::export]]
Rcpp::List LAMM(const arma::mat& X, const arma::vec& Y, const arma::vec& Lambda, arma::vec beta,
                const double phi, const std::string lossType, const double tau, 
                const double gamma) {
  double phiNew = phi;
  arma::vec betaNew = arma::vec();
  while (true) {
    betaNew = updateBeta(X, Y, beta, phiNew, Lambda, lossType, tau);
    double FVal = loss(Y, X * betaNew, lossType, tau);
    double PsiVal = cmptPsi(X, Y, betaNew, beta, phiNew, lossType, tau);
    if (FVal <= PsiVal) {
      break;
    }
    phiNew *= gamma;
  }
  return Rcpp::List::create(Rcpp::Named("beta") = betaNew, Rcpp::Named("phi") = phiNew);
}

// [[Rcpp::export]]
arma::vec lasso(const arma::mat& X, const arma::vec& Y, const double lambda,
                const double phi0 = 0.001, const double gamma = 1.5, 
                const double epsilon_c = 0.001, const int iteMax = 500) {
  int d = X.n_cols - 1;
  arma::vec beta = arma::zeros(d + 1);
  arma::vec betaNew = arma::zeros(d + 1);
  arma::vec Lambda = cmptLambda(beta, lambda);
  double phi = phi0;
  int ite = 0;
  while (ite < iteMax) {
    ite++;
    Rcpp::List listLAMM = LAMM(X, Y, Lambda, beta, phi, "l2", 1, gamma);
    betaNew = Rcpp::as<arma::vec>(listLAMM["beta"]);
    phi = listLAMM["phi"];
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon_c) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// [[Rcpp::export]]
Rcpp::List huberLasso(const arma::mat& X, const arma::vec& Y, const double lambda,
                      double tau = -1, const double constTau = 1.345, const double phi0 = 0.001, 
                      const double gamma = 1.5, const double epsilon_c = 0.001, 
                      const int iteMax = 500) {
  int d = X.n_cols - 1;
  arma::vec beta = arma::zeros(d + 1);
  arma::vec betaNew = arma::zeros(d + 1);
  arma::vec Lambda = cmptLambda(beta, lambda);
  if (tau < 0) {
    arma::vec betaLasso = lasso(X, Y, lambda, phi0, gamma, epsilon_c, iteMax);
    arma::vec res = Y - X * betaLasso;
    double mad = arma::median(arma::abs(res - arma::median(res))) / 0.6744898;
    tau = constTau * mad;
  }
  double phi = phi0;
  int ite = 0;
  while (ite < iteMax) {
    ite++;
    Rcpp::List listLAMM = LAMM(X, Y, Lambda, beta, phi, "Huber", tau, gamma);
    betaNew = Rcpp::as<arma::vec>(listLAMM["beta"]);
    phi = listLAMM["phi"];
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon_c) {
      break;
    }
    beta = betaNew;
    arma::vec res = Y - X * beta;
    double mad = arma::median(arma::abs(res - arma::median(res))) / 0.6744898;
    tau = constTau * mad;
  }
  return Rcpp::List::create(Rcpp::Named("beta") = betaNew, Rcpp::Named("tau") = tau,
                            Rcpp::Named("iteration") = ite);
}

// [[Rcpp::export]]
arma::uvec getIndex(const int n, const int low, const int up) {
  arma::vec seq = arma::regspace(0, n - 1);
  return arma::find(seq >= low && seq <= up);
}

// [[Rcpp::export]]
arma::uvec getIndexComp(const int n, const int low, const int up) {
  arma::vec seq = arma::regspace(0, n - 1);
  return arma::find(seq < low || seq > up);
}

// [[Rcpp::export]]
double pairPred(const arma::mat& X, const arma::vec& Y, const arma::vec& beta) {
  int n = X.n_rows;
  int d = X.n_cols - 1;
  int m = n * (n - 1) >> 1;
  arma::mat pairX(m, d + 1);
  arma::vec pairY(m);
  int k = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      pairX.row(k) = X.row(i) - X.row(j);
      pairY(k++) = Y(i) - Y(j);
    }
  }
  arma::vec predY = pairX * beta;
  return arma::sum(arma::square(pairY - predY));
}

//' The function fits Huber-Lasso regression via I-LAMM algorithm, with \eqn{\tau} determined by a tuning-free principle, \eqn{\lambda} calibrated by k-folds cross-validation, and the intercept term \eqn{\beta_0} estimated via a two-step procedure.
//'
//' The observed data are \eqn{(Y, X)}, where \eqn{Y} is an \eqn{n}-dimensional response vector and \eqn{X} is an \eqn{n} by \eqn{d} design matrix. We assume that \eqn{Y} depends on \eqn{X} through a linear model \eqn{Y = X \beta + \epsilon}, where \eqn{\beta} is a sparse vector and \eqn{\epsilon} is an \eqn{n}-dimensional noise vector whose distribution can be asymmetrix and/or heavy-tailed. All the arguments except for \eqn{X} and \eqn{Y} have default settings.
//'
//' @title Tuning-free Huber-Lasso regression
//' @param X An \eqn{n} by \eqn{d} design matrix with each row being a sample and each column being a variable, either low-dimensional data (\eqn{d \le n}) or high-dimensional data (\eqn{d > n}) are allowed..
//' @param Y A continuous response vector with length \eqn{n}.
//' @param lSeq Sequence of tuning parameter of regularized regression \eqn{\lambda}, every element should be positive. If it's not specified, the default sequence is generated in this way: define \eqn{\lambda_max = max(|Y^T X|) / n}, and \eqn{\lambda_min = 0.01 * \lambda_max}, then \code{lseq} is a sequence from \eqn{\lambda_max} to \eqn{\lambda_min} that decreases uniformly on log scale.
//' @param nlambda Number of \eqn{\lambda} to generate the default sequence \code{lSeq}. It's not necessary if \code{lSeq} is specified. The default value is 30.
//' @param constTau The constant term used to update \eqn{\tau} in the tuning-free procedure. In each round of iteration, \eqn{\tau} is updated to be \code{constTau} \eqn{* \sigma_MAD}, where \eqn{\sigma_MAD = median(|R - median(R)|) / \Phi^(-1)(3/4)} is the median absolute deviation estimator, and \eqn{R} is the residual from last round of iteration. The defalut value is 1.345.
//' @param phi0 The initial value of the isotropic parameter \eqn{\phi} in I-LAMM algorithm. The defalut value is 0.001.
//' @param gamma The inflation parameter in I-LAMM algorithm, in each iteration of I-LAMM, we will inflate \eqn{\phi} by \eqn{\gamma}. The defalut value is 1.5.
//' @param epsilon_c The tolerance level for I-LAMM algorithm, iteration will stop when \eqn{||\theta_new - \theta_old||_inf < \epsilon_c}. The defalut value is 1e-3.
//' @param iteMax The maximal number of iteration in I-LAMM algorithm, the iteration stops if this number is reached. The defalut value is 500.
//' @param nfolds The number of folds to conduct cross validation for \eqn{\lambda}, values that are greater than 10 are not recommended, and it'll be modified to 10 if the input is greater than 10. The default value is 3.
//' @return A list including the following terms will be returned:
//' \itemize{
//' \item \code{theta} The estimated \eqn{\theta}, a vector with length \eqn{d + 1}, with the first one being the value of intercept.
//' \item \code{lambdaSeq} The sequence of \eqn{\lambda}'s for cross validation.
//' \item \code{lambdaMin} The value of \eqn{\lambda} in \code{lSeq} that minimized mse in k-fold cross-validation.
//' \item \code{tauCoef} The robustness parameter \eqn{\tau} determined by the tuning-free principle to estimate coefficients except for the intercept.
//' \item \code{tauItcp} The robustness parameter \eqn{\tau} determined by the tuning-free principle to estimate the intercept.
//' \item \code{iteCoef} The number of iterations in I-LAMM algirithm to estimate coefficients.
//' \item \code{iteItcp} The number of iterations to estimate the intercept.
//' }
//' @author Xiaoou Pan, Wen-Xin Zhou
//' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint.
//' @references Fan, J., Liu, H., Sun, Q. and Zhang, T. (2018). I-LAMM for sparse learning: Simultaneous control of algorithmic complexity and statistical error. Ann. Statist. 46 814–841.
//' @seealso \code{\link{huberReg}}
//' @examples
//' n = 200
//' d = 500
//' s = 5
//' thetaStar = c(rep(3, s + 1), rep(0, d - s))
//' X = matrix(rnorm(n * d), n, d)
//' error = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
//' Y = as.numeric(cbind(rep(1, n), X) %*% thetaStar + error)
//' listHuberLasso = cvHuberLasso(X, Y)
//' thetaHuberLasso = listHuberLasso$theta
//' @export
// [[Rcpp::export]]
Rcpp::List cvHuberLasso(const arma::mat& X, const arma::vec& Y,
                        Rcpp::Nullable<Rcpp::NumericVector> lSeq = R_NilValue, int nlambda = 30,
                        const double constTau = 1.345, const double phi0 = 0.001, 
                        const double gamma = 1.5, const double epsilon_c = 0.001, 
                        const int iteMax = 500, int nfolds = 3) {
  int n = X.n_rows;
  int d = X.n_cols;
  arma::mat Z = arma::ones(n, d + 1);
  Z.cols(1, d) = X;
  arma::vec lambdaSeq = arma::vec();
  if (lSeq.isNotNull()) {
    lambdaSeq = Rcpp::as<arma::vec>(lSeq);
    nlambda = lambdaSeq.size();
  } else {
    double lambdaMax = arma::max(arma::abs(Y.t() * Z)) / n;
    double lambdaMin = 0.01 * lambdaMax;
    lambdaSeq = arma::exp(arma::linspace(std::log((long double)lambdaMin),
                                   std::log((long double)lambdaMax), nlambda));
  }
  if (nfolds > 10 || nfolds > n) {
    nfolds = n < 10 ? n : 10;
  }
  int size = n / nfolds;
  arma::vec mse = arma::zeros(nlambda);
  for (int i = 0; i < nlambda; i++) {
    for (int j = 0; j < nfolds; j++) {
      int low = j * size;
      int up = (j == (nfolds - 1)) ? (n - 1) : ((j + 1) * size - 1);
      arma::uvec idx = getIndex(n, low, up);
      arma::uvec idxComp = getIndexComp(n, low, up);
      Rcpp::List hLassoList = huberLasso(Z.rows(idxComp), Y.rows(idxComp), lambdaSeq(i), -1, 
                                         constTau, phi0, gamma, epsilon_c, iteMax);
      arma::vec thetaHat = Rcpp::as<arma::vec>(hLassoList["beta"]);
      mse(i) += pairPred(Z.rows(idx), Y.rows(idx), thetaHat);
    }
  }
  arma::uword cvIdx = mse.index_min();
  Rcpp::List hLassoList = huberLasso(Z, Y, lambdaSeq(cvIdx), -1, constTau, phi0, gamma, 
                                     epsilon_c, iteMax);
  arma::vec theta = Rcpp::as<arma::vec>(hLassoList["beta"]);
  Rcpp::List listMean = huberMean(Y - Z.cols(1, d) * theta.rows(1, d));
  theta(0) = listMean["mu"];
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("lambdaSeq") = lambdaSeq,
                            Rcpp::Named("lambdaMin") = lambdaSeq(cvIdx), 
                            Rcpp::Named("tauCoef") = hLassoList["tau"], 
                            Rcpp::Named("tauItcp") = listMean["tau"], 
                            Rcpp::Named("iteCoef") = hLassoList["iteration"],
                            Rcpp::Named("iteItcp") = listMean["iteration"]);
}
