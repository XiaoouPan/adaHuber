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
  return Rcpp::List::create(Rcpp::Named("theta") = thetaNew, Rcpp::Named("tau1") = tauNew, 
                            Rcpp::Named("tau2") = listMean["tau"], Rcpp::Named("iteration1") = iteNum,
                            Rcpp::Named("iteration2") = listMean["iteration"]);
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

//' @export
// [[Rcpp::export]]
Rcpp::List cvHuberLasso(const arma::mat& X, const arma::vec& Y,
                        Rcpp::Nullable<Rcpp::NumericVector> lSeq = R_NilValue, int nlambda = 30,
                        const double constTau = 1.345, const double phi0 = 0.001, 
                        const double gamma = 1.5, const double epsilon_c = 0.001, 
                        const int iteMax = 50, int nfolds = 3) {
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
