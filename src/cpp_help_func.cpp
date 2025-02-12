#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title C++ implementation of the photometric measurement uncertainty.
//' 
//' @description Compute the photometric measurement uncertainty based on \deqn{\sigma(M) = \beta_0\exp(\beta_1(M - m_0)).}
//' 
//' @param M Numeric with length 1. The true GC magnitude
//' @param m0 Numeric with length 1. The offset magnitude. Default to `m0 = 25.5`. 
//' @param b0 Numeric with length 1. The uncertainty when M = m0. Default to `b0 = 0.08836`.
//' @param b1 Numeric with length 1. Rate at which the uncertainty increases with M. Default to `b1 = 0.645`.
//' @export
//' @return The measurement uncertainty at the true magnitude `M`.
// [[Rcpp::export]]
double err_M_cpp(double M, double m0 = 25.5, double b0 = 0.08836, double b1 = 0.645){
  return b0*exp(b1*(M - m0));
}

//' @title C++ implementation of the completeness fraction.
//' 
//' @description Compute the completeness fraction based on \deqn{f(M)=\frac{1}{\exp(\alpha(M- M_{\text{Lim}}))}.}
//' 
//' @param M A numeric vector. The true GC magnitude
//' @param Lim Numeric with length 1. 50% completeness fraction.
//' @param a Numeric with length 1. Rate at which the completeness fraction drops with M. Default to `a = 1.5`.
//' @export
//' @return The completeness fraction evaluated at the true magnitude `M`.
//[[Rcpp::export]]
NumericVector f_cpp(NumericVector M, double Lim, double a = 1.5){
  return 1/(1 + exp(a*(M - Lim)));
}


//[[Rcpp::export]]
NumericVector f_cpp_trunc(NumericVector M, double Lim, double a = 1.5){
  NumericVector ind = ifelse(M < 26.3, 1.0, 0.0);
  return 1/(1 + exp(a*(M - Lim)))*ind;
}

// [[Rcpp::export]]
double dnorm_cpp(double x, double mu, double sigma){
  if(sigma == 0 & x == mu){
    return R_PosInf;
  }
  else if(sigma == 0){
    return 0;
  }
  else{
    return 1/sqrt(2*M_PI)/sigma*exp(-(x-mu)*(x-mu)/2/sigma/sigma);
  }
}

//' @title C++ implementation to calculate the value of the noisy GCLF given the observed noisy magnitude `M`.
//'
//' @description Compute the noisy GCLF evaluated at `M` based on \deqn{\varphi(m; \mu_{\text{TO}}, \sigma^2) = \int_{-\infty}^\infty\phi(m; m_t, \sigma_M^2(m_t))\phi(m_t; \mu_{\text{TO}}, \sigma^2)dm_t,}
//' where \eqn{\phi(m; m_t, \sigma_M^2(m_t))} is a Gaussian density that gives the probability density of the noisy magnitude \eqn{m} given the true magnitude \eqn{m_t}. \eqn{\phi(m_t; \mu_{\text{TO}}, \sigma^2)} is the true GCLF with TO point \eqn{\mu_{\text{TO}}} and dispersion \eqn{\sigma}, evaluated at the true magnitude \eqn{m_t}.
//'
//' @param M A numeric vector. The observed GC magnitude
//' @param mu Numeric with length 1. True GCLF TO point
//' @param sigma Numeric with length 1. True GCLF dispersion
//' @param m0 Numeric with length 1. The offset magnitude
//' @param b0 Numeric with length 1. The uncertainty when M = m0
//' @param b1 Numeric with length 1. Rate at which the uncertainty increases with M
//' @param n Integer. Number of integration points used to marginalize out the true magnitude. Default to 20.
//' @export
//' @return A numeric vector that gives the value of the noisy GCLF given the observed noisy magnitude `M`.
// [[Rcpp::export]]
NumericVector phi_eM_cpp(NumericVector M, double mu, double sigma, double m0, double b0, double b1, int n = 20){
  NumericVector s (M.size());
  double dx;
  double sig;
  NumericVector x;
  for(int j = 0; j <= M.size() - 1; j++){
    sig = err_M_cpp(M(j));
    if(sig == R_PosInf){
      sig = 100;
    }
    else{
      sig = sig*(sig > 0.0001 & sig < 100) + 0.0001*(sig < 0.0001) + 100*(sig > 100);
    }
    x = wrap(arma::linspace(M(j) - 6*sig, M(j) + 6*sig, (int)n*std::max(1.0, 5*sig)));
    dx = x(1) - x(0);
    for(int i = 0; i <= x.size() - 1; i++){
      s(j) = s(j) + dx*dnorm_cpp(x(i), mu, sigma)*dnorm_cpp(x(i), M(j), err_M_cpp(x(i), m0, b0, b1));
    }
  }
  return s;
}

//[[Rcpp::export]]
double Phi_f_cpp(double Lim, double mu, double sigma, int n = 50, double a = 1.5){
  NumericVector x = wrap(arma::linspace(mu - 6*sigma, mu + 6*sigma, n));
  double dx = x(1) - x(0);
  return sum(dx*Rcpp::dnorm(x, mu, sigma)*f_cpp(x, Lim, a));
}

//' @title C++ implementation to calculate the proportion of remaining GCs after removal of faint GCs.
//' 
//' @description Compute the proportion of the remaining GCs after removing the faint GCs based on their noisy magnitudes and the completeness fraction \eqn{f(M)}. 
//' It computes the following integral \deqn{\int_{-\infty}^\infty\varphi(m; \mu_{\text{TO}}, \sigma^2)f(m)dm,} where the noisy GCLF \eqn{\varphi(m; \mu_{\text{TO}}, \sigma^2)} is obtained from [phi_eM_cpp()].
//' 
//' @param Lim Numeric with length 1. 50% completeness fraction
//' @param mu Numeric with length 1. True GCLF TO point
//' @param sigma Numeric with length 1. True GCLF dispersion
//' @param m0 Numeric with length 1. The offset magnitude
//' @param b0 Numeric with length 1. The uncertainty when M = m0
//' @param b1 Numeric with length 1. Rate at which the uncertainty increases with M
//' @param a Numeric with length 1. Rate at which the completeness fraction drops with M
//' @param n Integer. Number of integration points used to marginalize out the true magnitude. Default to 50.
//' @export
//' @return A numeric value between \eqn{(0,1)}: the proportion of the GCs remaining after removing the faint GCs.
// [[Rcpp::export]]
double p_eM_cpp(double Lim, double mu, double sigma, double m0, double b0, double b1, double a, int n = 50){
  NumericVector x = wrap(arma::linspace(mu - 6*sigma, mu + 6*sigma, n));
  double dx = x(1) - x(0);
  double res = 0.0;
  NumericVector y = phi_eM_cpp(x, mu, sigma, m0, b0, b1)*f_cpp(x, Lim, a);
  for(int j = 0; j <= y.size() - 1; j++){
    res = res + dx*y(j);
  }
 return res;
}

//[[Rcpp::export]]
double Phi_f_trun_cpp(double Lim, double mu, double sigma, int n = 50, double a = 1.5){
  NumericVector x = wrap(arma::linspace(mu - 6*sigma, 26.3, n));
  x = 26.3 - pow(26.3 - x, 3.0)/40.0;
  double res = 0.0;
  NumericVector diffx = diff(x);
  NumericVector y = Rcpp::dnorm(x, mu, sigma)*f_cpp_trunc(x, Lim, a);
  double dx; 
  for(int j = 1; j <= diffx.size() - 1; j++){
    dx = diffx(j-1);
    res = res + 0.5*dx*(y(j - 1) + y(j));
  }
  return res;
}

// [[Rcpp::export]]
double p_eM_trun_cpp(double Lim, double mu, double sigma, double m0, double b0, double b1, double a, double TO = 26.3, int n = 100){
  NumericVector x = wrap(arma::linspace(mu - 6*sigma, TO, n));
  x = TO - pow(TO - x, 3.0)/40.0;
  double res = 0.0;
  NumericVector diffx = diff(x);
  NumericVector y = phi_eM_cpp(x, mu, sigma, m0, b0, b1)*f_cpp_trunc(x, Lim, a);
  double dx; 
  for(int j = 1; j <= diffx.size() - 1; j++){
    dx = diffx(j-1);
    res = res + 0.5*dx*(y(j - 1) + y(j));
  }
  return res;
}


