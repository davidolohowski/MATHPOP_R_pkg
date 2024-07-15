#' @title Simulate point pattern from a given Sersic profile 
#' 
#' @description Simulate point pattern from a given Sersic profile given by \deqn{\Lambda(s; N, R_h, n) = \frac{N b_n^{2n}}{2\pi R_h^2 n \Gamma(2n)e}\exp\left(-b_n\left(\frac{r(s)}{R_h}\right)^{1/n}\right),} 
#' where \deqn{r^2(s) = ((s_x - c_x)\cos(\vartheta) - (s_y - c_y)\sin(\vartheta))^2 + ((s_x- c_x)\sin(\vartheta)+ (s_y-c_y)\cos(\vartheta))^2/e^2.}.
#' 
#' @param N A non-negative numeric value. Number of GCs to simulate.
#' @param c A vector of length two for the central location of a galaxy. The first entry is the x coordinate of the galactic center, the second entry is the y coordinate.
#' @param R_eff A non-negative numeric value. Half-number radius of the GC system of a galaxy.
#' @param e A non-negative numeric value. Aspect ratio of the GC system of a galaxy.
#' @param n A non-negative numeric value. Sersic index of the GC system of a galaxy.
#' @param theta A numeric value in \eqn{(0, 2\pi)}. Orientation angle (in radian) of the GC system of a galaxy.
#' @keywords Sersic Intensity
#' @export
#' @return A data frame with `N` number of row and two columns that gives the locations of the simulated point pattern.
simulate_Sersic <- function(N, c, R_eff, e, n, theta){
  bn <- zipfR::Rgamma.inv(2*n, 0.5)
  z <- runif(N)
  x <- zipfR::Rgamma.inv(2*n, z)
  R <- (x/bn)^n*R_eff
  t <- runif(N, 0, 2*pi)
  A <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)
  X <- solve(A, t(cbind(R*cos(t), R*sin(t)*e)))
  return(data.frame(x = c[1] + X[1,], y = c[2] + X[2,]))
}
