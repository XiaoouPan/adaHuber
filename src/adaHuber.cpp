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

// [[Rcpp::export]]
double f2(const double x, const arma::vec& resSq, const int N, const double rhs) {
  return arma::mean(arma::min(resSq / x, arma::ones(N))) - rhs;
}

// [[Rcpp::export]]
double rootf2(const arma::vec& resSq, const int n, const int d, const int N, const double rhs, double low, double up, const double tol = 0.001, 
              const int maxIte = 500) {
  int ite = 0;
  while (ite <= maxIte && up - low > tol) {
    double mid = 0.5 * (up + low);
    double val = f2(mid, resSq, N, rhs);
    if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

// [[Rcpp::export]]
double g1(const double x, const arma::vec& resSq, const int n, const double rhs) {
  return arma::mean(arma::min(resSq / x, arma::ones(n))) - rhs;
}

// [[Rcpp::export]]
double rootg1(const arma::vec& resSq, const int n, const double rhs, double low, double up, const double tol = 0.001, const int maxIte = 500) {
  int ite = 0;
  while (ite <= maxIte && up - low > tol) {
    double mid = 0.5 * (up + low);
    double val = g1(mid, resSq, n, rhs);
    if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

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
Rcpp::List huberMean(arma::vec X, const int n, const double tol = 0.0001, const int iteMax = 500) {
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
arma::vec huberMeanVec(const arma::mat& X, const int n, const int p, const double epsilon = 0.001, const int iteMax = 500) {
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
    tau = std::sqrt((long double)rootf2(resSq, n, d, N, rhs, arma::min(resSq), arma::accu(resSq)));
    w = arma::min(tau / arma::abs(res), arma::ones(N));
    muNew = arma::as_scalar(Z.t() * w) / arma::accu(w);
    iteNum++;
  }
  return muNew;
}

// [[Rcpp::export]]
Rcpp::List huberCov(const arma::mat& X, const int n, const int p) {
  double rhs2 = (2 * std::log(p) + std::log(n)) / n;
  arma::vec mu(p);
  arma::mat sigmaHat(p, p);
  for (int j = 0; j < p; j++) {
    mu(j) = huberMean(X.col(j), n);
    double theta = huberMean(arma::square(X.col(j)), n);
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
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) / sx(i);
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
arma::vec adaHuberReg(const arma::mat& X, arma::vec Y, const int n, const int p, const double tol = 0.0001, const int iteMax = 5000) {
  const double n1 = 1.0 / n;
  double rhs = n1 * (p + std::log(n * p));
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  Y -= my;
  double tau = 1.345 * mad(Y);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootg1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
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
    tau = std::sqrt((long double)rootg1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) /= sx;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p), n);
  return beta;
}

// [[Rcpp::export]]
arma::vec huberReg(const arma::mat& X, arma::vec Y, const int n, const int p, const double tol = 0.0001, const double constTau = 1.345, 
                   const int iteMax = 5000) {
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
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
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
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
  beta.rows(1, p) /= sx;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p), n);
  return beta;
}
