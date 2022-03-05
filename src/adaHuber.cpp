# include <RcppArmadillo.h>
# include <algorithm>
# include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int sgn(const double x) {
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
double f1(const double x, const arma::vec& resSq, const int n, const double rhs) {
  return arma::mean(arma::min(resSq / x, arma::ones(n))) - rhs;
}

// [[Rcpp::export]]
double rootf1(const arma::vec& resSq, const int n, const double rhs, double low, double up, const double tol = 0.001, const int maxIte = 500) {
  int ite = 1;
  while (ite <= maxIte && up - low > tol) {
    double mid = 0.5 * (up + low);
    double val = f1(mid, resSq, n, rhs);
    if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

// Functions for high-dim tunning-free calibration
/*// [[Rcpp::export]]
double f2(const double x, const arma::vec& resSq, const int n, const int s, const double rhs) {
  return arma::accu(arma::min(resSq / x, arma::ones(n))) / (n - s) - rhs;
}

// [[Rcpp::export]]
double rootf2(const arma::vec& resSq, const int n, const int s, const double rhs, double low, double up, const double tol = 0.001, 
              const int maxIte = 500) {
  int ite = 1;
  while (ite <= maxIte && up - low > tol) {
    double mid = 0.5 * (up + low);
    double val = f2(mid, resSq, n, s, rhs);
    if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}*/

// [[Rcpp::export]]
double huberDer(const arma::vec& res, const double tau, const int n) {
  double rst = 0.0;
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    rst -= std::abs(cur) <= tau ? cur : tau * sgn(cur);
  }
  return rst / n;
}

// [[Rcpp::export]]
double huberMean(arma::vec X, const int n, const double tol = 0.0001, const int iteMax = 500) {
  double rhs = std::log(n) / n;
  double mx = arma::mean(X);
  X -= mx;
  double tau = arma::stddev(X) * std::sqrt((long double)n / std::log(n));
  double derOld = huberDer(X, tau, n);
  double mu = -derOld, muDiff = -derOld;
  arma::vec res = X - mu;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  double derNew = huberDer(res, tau, n);
  double derDiff = derNew - derOld;
  int ite = 1;
  while (std::abs(derNew) > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = muDiff * derDiff;
    if (cross > 0) {
      double a1 = cross / derDiff * derDiff;
      double a2 = muDiff * muDiff / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    derOld = derNew;
    muDiff = -alpha * derNew;
    mu += muDiff;
    res = X - mu;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    derNew = huberDer(res, tau, n);
    derDiff = derNew - derOld;
    ite++;
  }
  return mu + mx;
}

// The function that is called in R
// [[Rcpp::export]]
Rcpp::List huberMeanList(arma::vec X, const double tol = 0.0001, const int iteMax = 500) {
  int n = X.size();
  double rhs = std::log(n) / n;
  double mx = arma::mean(X);
  X -= mx;
  double tau = arma::stddev(X) * std::sqrt((long double)n / std::log(n));
  double derOld = huberDer(X, tau, n);
  double mu = -derOld, muDiff = -derOld;
  arma::vec res = X - mu;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  double derNew = huberDer(res, tau, n);
  double derDiff = derNew - derOld;
  int ite = 1;
  while (std::abs(derNew) > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = muDiff * derDiff;
    if (cross > 0) {
      double a1 = cross / derDiff * derDiff;
      double a2 = muDiff * muDiff / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    derOld = derNew;
    muDiff = -alpha * derNew;
    mu += muDiff;
    res = X - mu;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    derNew = huberDer(res, tau, n);
    derDiff = derNew - derOld;
    ite++;
  }
  mu += mx;
  return Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named("tau") = tau, Rcpp::Named("iteration") = ite);
}

// [[Rcpp::export]]
arma::vec huberMeanVec(const arma::mat& X, const int n, const int p, const double epsilon = 0.0001, const int iteMax = 500) {
  arma::vec rst(p);
  for (int i = 0; i < p; i++) {
    rst(i) = huberMean(X.col(i), n, epsilon, iteMax);
  }
  return rst;
}

// [[Rcpp::export]]
double hMeanCov(const arma::vec& Z, const int n, const int d, const int N, double rhs, const double epsilon = 0.0001, const int iteMax = 500) {
  double muOld = 0;
  double muNew = arma::mean(Z);
  double tau = arma::stddev(Z) * std::sqrt((long double)n / (2 * std::log(d) + std::log(n)));
  int iteNum = 0;
  arma::vec res(n), resSq(n), w(n);
  while ((std::abs(muNew - muOld) > epsilon) && iteNum < iteMax) {
    muOld = muNew;
    res = Z - muOld;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, N, rhs, arma::min(resSq), arma::accu(resSq)));
    w = arma::min(tau / arma::abs(res), arma::ones(N));
    muNew = arma::as_scalar(Z.t() * w) / arma::accu(w);
    iteNum++;
  }
  return muNew;
}

// The function for covariance estimation called in R.
// [[Rcpp::export]]
Rcpp::List huberCov(const arma::mat& X, const double epsilon = 0.0001, const int iteMax = 500) {
  int n = X.n_rows;
  int p = X.n_cols;
  double rhs2 = (2 * std::log(p) + std::log(n)) / n;
  arma::vec mu(p);
  arma::mat sigmaHat(p, p);
  for (int j = 0; j < p; j++) {
    mu(j) = huberMean(X.col(j), n, epsilon, iteMax);
    double theta = huberMean(arma::square(X.col(j)), n, epsilon, iteMax);
    double temp = mu(j) * mu(j);
    if (theta > temp) {
      theta -= temp;
    }
    sigmaHat(j, j) = theta;
  }
  int N = n * (n - 1) >> 1;
  arma::mat Y(N, p);
  for (int i = 0, k = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      Y.row(k++) = X.row(i) - X.row(j);
    }
  }
  for (int i = 0; i < p - 1; i++) {
    for (int j = i + 1; j < p; j++) {
      sigmaHat(i, j) = sigmaHat(j, i) = hMeanCov(0.5 * Y.col(i) % Y.col(j), n, p, N, rhs2);
    }
  }
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("cov") = sigmaHat);
}

// [[Rcpp::export]]
double mad(const arma::vec& x) {
  return 1.482602 * arma::median(arma::abs(x - arma::median(x)));
}

// [[Rcpp::export]]
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx1, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) * sx1(i);
  }
  return X;
}

// [[Rcpp::export]]
void updateHuber(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double n1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (std::abs(cur) <= tau) {
      der(i) = -cur;
    } else {
      der(i) = -tau * sgn(cur);
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
Rcpp::List adaHuberReg(const arma::mat& X, arma::vec Y, const double epsilon = 0.0001, const int iteMax = 500) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const double n1 = 1.0 / n;
  double rhs = n1 * (p + std::log(n * p));
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  Y -= my;
  double tau = 1.345 * mad(Y);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > epsilon && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) %= sx1;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p), n, epsilon, iteMax);
  return Rcpp::List::create(Rcpp::Named("coef") = beta, Rcpp::Named("tau") = tau, Rcpp::Named("iteration") = ite);
}

// [[Rcpp::export]]
Rcpp::List huberReg(const arma::mat& X, arma::vec Y, const double epsilon = 0.0001, const double constTau = 1.345, const int iteMax = 500) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  Y -= my;
  double tau = constTau * mad(Y);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  tau = constTau * mad(res);
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > epsilon && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    tau = constTau * mad(res);
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) %= sx1;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p), n, epsilon, iteMax);
  return Rcpp::List::create(Rcpp::Named("coef") = beta, Rcpp::Named("tau") = tau, Rcpp::Named("iteration") = ite);
}

// [[Rcpp::export]]
arma::vec softThresh(const arma::vec& x, const arma::vec& lambda, const int p) {
  return arma::sign(x) % arma::max(arma::abs(x) - lambda, arma::zeros(p + 1));
}

// [[Rcpp::export]]
arma::vec cmptLambdaLasso(const double lambda, const int p) {
  arma::vec rst = lambda * arma::ones(p + 1);
  rst(0) = 0;
  return rst;
}

// [[Rcpp::export]]
double lossHuber(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta, const double n1, const double tau) {
  arma::vec res = Y - Z * beta;
  double rst = 0.0;
  for (int i = 0; i < Y.size(); i++) {
    double temp = std::abs(res(i));
    if (temp <= tau) {
      rst += 0.5 * temp * temp;
    } else {
      rst += tau * temp - 0.5 * tau * tau;
    }
  }
  return n1 * rst;
}

// [[Rcpp::export]]
double updateHuberHd(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta, arma::vec& grad, const double n1, const double tau) {
  arma::vec res = Y - Z * beta;
  double rst = 0.0;
  grad = arma::zeros(grad.size());
  for (int i = 0; i < Y.size(); i++) {
    double cur = res(i);
    if (std::abs(cur) <= tau) {
      grad -= cur * Z.row(i).t();
      rst += 0.5 * cur * cur;
    } else {
      grad -= tau * sgn(cur) * Z.row(i).t();
      rst += tau * std::abs(cur) - 0.5 * tau * tau;
    }
  }
  grad *= n1;
  return n1 * rst;
}

// [[Rcpp::export]]
double lamm(const arma::mat& Z, const arma::vec& Y, const arma::vec& Lambda, arma::vec& beta, const double tau, const double phi, const double gamma, 
            const int p, const double n1) {
  double phiNew = phi;
  arma::vec betaNew(p + 1);
  arma::vec grad(p + 1);
  double loss = updateHuberHd(Z, Y, beta, grad, n1, tau);
  while (true) {
    arma::vec first = beta - grad / phiNew;
    arma::vec second = Lambda / phiNew;
    betaNew = softThresh(first, second, p);
    double fVal = lossHuber(Z, Y, betaNew, n1, tau);
    arma::vec diff = betaNew - beta;
    double psiVal = loss + arma::as_scalar(grad.t() * diff) + 0.5 * phiNew * arma::as_scalar(diff.t() * diff);
    if (fVal <= psiVal) {
      break;
    }
    phiNew *= gamma;
  }
  beta = betaNew;
  return phiNew;
}

// Huber-lasso with a specified lambda and tau
// [[Rcpp::export]]
arma::vec huberLasso(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const int p, const double n1, const double phi0 = 0.1, 
                     const double gamma = 1.2, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lamm(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// The function that is called in R
// [[Rcpp::export]]
Rcpp::List huberLassoList(const arma::mat& X, arma::vec& Y, const double lambda, const double tau, const double phi0 = 0.1, const double gamma = 1.2, 
                          const double epsilon = 0.001, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lamm(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  betaNew.rows(1, p) %= sx1;
  betaNew(0) += my - arma::as_scalar(mx * betaNew.rows(1, p));
  return Rcpp::List::create(Rcpp::Named("coef") = betaNew, Rcpp::Named("iteration") = ite, Rcpp::Named("phi") = phi);
}

// [[Rcpp::export]]
int sparsity(const arma::vec& x) {
  arma::vec temp = arma::nonzeros(x);
  return temp.size();
}

// [[Rcpp::export]]
arma::vec adaHuberLasso(const arma::mat& Z, const arma::vec& Y, const double lambda, const int p, const double phi0 = 0.1, const double gamma = 1.2, 
                        const double epsilon = 0.001, const int iteMax = 500) {
  const int n = Z.n_rows;
  const double n1 = 1.0 / n;
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double rhs = n1 * (std::log(n * p));
  double tau = 1.345 * mad(Y);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lamm(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
    //int sparse = sparsity(beta);
    arma::vec res = Y - Z * beta;
    arma::vec resSq = arma::square(Y);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  }
  return betaNew;
}

// Huber-lasso with a specified lambda and a tuning-free tau. It is called in R.
// [[Rcpp::export]]
Rcpp::List adaHuberLassoList(const arma::mat& X, arma::vec& Y, const double lambda, const double phi0 = 0.1, const double gamma = 1.2, 
                             const double epsilon = 0.001, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double rhs = n1 * (std::log(n * p));
  double tau = 1.345 * mad(Y);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lamm(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
    arma::vec res = Y - Z * beta;
    arma::vec resSq = arma::square(Y);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  }
  betaNew.rows(1, p) %= sx1;
  betaNew(0) += my - arma::as_scalar(mx * betaNew.rows(1, p));
  return Rcpp::List::create(Rcpp::Named("coef") = betaNew, Rcpp::Named("iteration") = ite, Rcpp::Named("tau") = tau, Rcpp::Named("phi") = phi);
}

// [[Rcpp::export]]
double lossL2(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta) {
  arma::vec res = Y - Z * beta;
  return arma::accu(arma::square(res));
}

// cross-validated l1-penalized Huber regression
// [[Rcpp::export]]
Rcpp::List cvAdaHuberLasso(const arma::mat& X, arma::vec& Y, const arma::vec& lambdaSeq, const arma::vec& folds, const int kfolds, 
                           const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.001, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec mse = arma::zeros(nlambda);
  arma::vec betaHat(p + 1);
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      betaHat = adaHuberLasso(trainZ, trainY, lambdaSeq(i), p, phi0, gamma, epsilon, iteMax);
      mse(i) += lossL2(testZ, testY, betaHat);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambdaSeq(cvIdx), p);
  double rhs = n1 * (std::log(n * p));
  double tau = 1.345 * mad(Y);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lamm(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
    arma::vec res = Y - Z * beta;
    arma::vec resSq = arma::square(Y);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  }
  betaNew.rows(1, p) %= sx1;
  betaNew(0) += my - arma::as_scalar(mx * betaNew.rows(1, p));
  return Rcpp::List::create(Rcpp::Named("coef") = betaNew, Rcpp::Named("tau") = tau, Rcpp::Named("lambda") = lambdaSeq(cvIdx), 
                            Rcpp::Named("iteration") = ite, Rcpp::Named("phi") = phi);
}
