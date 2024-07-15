#' @title Function to compute the (unnormalized) truncated noisy GCLF
#' 
#' @description Compute the (unnormalized) truncated noisy GCLF: \deqn{\varphi(m; \mu_{\text{TO}}, \sigma^2)f(m),}
#' where \eqn{\varphi(m; \mu_{\text{TO}}, \sigma^2)} is obtained from [phi_eM_cpp()].
#' 
#' 
#' @param M A numeric vector. Observed GC magnitudes.
#' @param mu A numeric value. True GCLF TO point.
#' @param sigma A numeric value. True GCLF dispersion.
#' @param Lim A numeric value. 50% completeness limit.
#' @param m0 A numeric value. Offset magnitude for measurement uncertainty.
#' @param b0 A numeric value. Measurement uncertainty when `M = m0`.
#' @param b1 A numeric value. Rate at which the measurement uncertainty increases as `M` increases.
#' @param a A numeric value. Rate at which the completeness fraction decreases as `M` increases.
#' @export
#' @return A numeric vector giving (unnormalized) truncated noisy GCLF evaluated at `M`.
psi_f <- function(M, mu, sigma, Lim, m0, b0, b1, a){
  phi_eM_cpp(M, mu, sigma, m0, b0, b1)*f(M, Lim, a)
}
