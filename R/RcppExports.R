# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title C++ implementation of the photometric measurement uncertainty.
#' 
#' @description Compute the photometric measurement uncertainty based on \deqn{\sigma(M) = \beta_0\exp(\beta_1(M - m_0)).}
#' 
#' @param M Numeric with length 1. The true GC magnitude
#' @param m0 Numeric with length 1. The offset magnitude. Default to `m0 = 25.5`. 
#' @param b0 Numeric with length 1. The uncertainty when M = m0. Default to `b0 = 0.08836`.
#' @param b1 Numeric with length 1. Rate at which the uncertainty increases with M. Default to `b1 = 0.645`.
#' @export
#' @return The measurement uncertainty at the true magnitude `M`.
err_M_cpp <- function(M, m0 = 25.5, b0 = 0.08836, b1 = 0.645) {
    .Call(`_MATHPOP_err_M_cpp`, M, m0, b0, b1)
}

#' @title C++ implementation of the completeness fraction.
#' 
#' @description Compute the completeness fraction based on \deqn{f(M)=\frac{1}{\exp(\alpha(M- M_{\text{Lim}}))}.}
#' 
#' @param M A numeric vector. The true GC magnitude
#' @param Lim Numeric with length 1. 50% completeness fraction.
#' @param a Numeric with length 1. Rate at which the completeness fraction drops with M. Default to `a = 1.5`.
#' @export
#' @return The completeness fraction evaluated at the true magnitude `M`.
f_cpp <- function(M, Lim, a = 1.5) {
    .Call(`_MATHPOP_f_cpp`, M, Lim, a)
}

f_cpp_trunc <- function(M, Lim, a = 1.5) {
    .Call(`_MATHPOP_f_cpp_trunc`, M, Lim, a)
}

dnorm_cpp <- function(x, mu, sigma) {
    .Call(`_MATHPOP_dnorm_cpp`, x, mu, sigma)
}

#' @title C++ implementation to calculate the value of the noisy GCLF given the observed noisy magnitude `M`.
#'
#' @description Compute the noisy GCLF evaluated at `M` based on \deqn{\varphi(m; \mu_{\text{TO}}, \sigma^2) = \int_{-\infty}^\infty\phi(m; m_t, \sigma_M^2(m_t))\phi(m_t; \mu_{\text{TO}}, \sigma^2)dm_t,}
#' where \eqn{\phi(m; m_t, \sigma_M^2(m_t))} is a Gaussian density that gives the probability density of the noisy magnitude \eqn{m} given the true magnitude \eqn{m_t}. \eqn{\phi(m_t; \mu_{\text{TO}}, \sigma^2)} is the true GCLF with TO point \eqn{\mu_{\text{TO}}} and dispersion \eqn{\sigma}, evaluated at the true magnitude \eqn{m_t}.
#'
#' @param M A numeric vector. The observed GC magnitude
#' @param mu Numeric with length 1. True GCLF TO point
#' @param sigma Numeric with length 1. True GCLF dispersion
#' @param m0 Numeric with length 1. The offset magnitude
#' @param b0 Numeric with length 1. The uncertainty when M = m0
#' @param b1 Numeric with length 1. Rate at which the uncertainty increases with M
#' @param n Integer. Number of integration points used to marginalize out the true magnitude. Default to 20.
#' @export
#' @return A numeric vector that gives the value of the noisy GCLF given the observed noisy magnitude `M`.
phi_eM_cpp <- function(M, mu, sigma, m0, b0, b1, n = 20L) {
    .Call(`_MATHPOP_phi_eM_cpp`, M, mu, sigma, m0, b0, b1, n)
}

Phi_f_cpp <- function(Lim, mu, sigma, n = 50L, a = 1.5) {
    .Call(`_MATHPOP_Phi_f_cpp`, Lim, mu, sigma, n, a)
}

#' @title C++ implementation to calculate the proportion of remaining GCs after removal of faint GCs.
#' 
#' @description Compute the proportion of the remaining GCs after removing the faint GCs based on their noisy magnitudes and the completeness fraction \eqn{f(M)}. 
#' It computes the following integral \deqn{\int_{-\infty}^\infty\varphi(m; \mu_{\text{TO}}, \sigma^2)f(m)dm,} where the noisy GCLF \eqn{\varphi(m; \mu_{\text{TO}}, \sigma^2)} is obtained from [phi_eM_cpp()].
#' 
#' @param Lim Numeric with length 1. 50% completeness fraction
#' @param mu Numeric with length 1. True GCLF TO point
#' @param sigma Numeric with length 1. True GCLF dispersion
#' @param m0 Numeric with length 1. The offset magnitude
#' @param b0 Numeric with length 1. The uncertainty when M = m0
#' @param b1 Numeric with length 1. Rate at which the uncertainty increases with M
#' @param a Numeric with length 1. Rate at which the completeness fraction drops with M
#' @param n Integer. Number of integration points used to marginalize out the true magnitude. Default to 50.
#' @export
#' @return A numeric value between \eqn{(0,1)}: the proportion of the GCs remaining after removing the faint GCs.
p_eM_cpp <- function(Lim, mu, sigma, m0, b0, b1, a, n = 50L) {
    .Call(`_MATHPOP_p_eM_cpp`, Lim, mu, sigma, m0, b0, b1, a, n)
}

Phi_f_trun_cpp <- function(Lim, mu, sigma, n = 50L, a = 1.5) {
    .Call(`_MATHPOP_Phi_f_trun_cpp`, Lim, mu, sigma, n, a)
}

p_eM_trun_cpp <- function(Lim, mu, sigma, m0, b0, b1, a, TO = 26.3, n = 100L) {
    .Call(`_MATHPOP_p_eM_trun_cpp`, Lim, mu, sigma, m0, b0, b1, a, TO, n)
}

