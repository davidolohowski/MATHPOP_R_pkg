#' Simulate point pattern from a given Sersic profile
#' 
#' @param N A non-negative numeric value. Number of GCs to simulate.
#' @param c A vector of length two for the central location of a galaxy. The first entry is the x coordinate of the galactic center, the second entry is the y coordinate.
#' @param R_eff A non-negative numeric value. Half-number radius of the GC system of a galaxy.
#' @param e A non-negative numeric value. Aspect ratio of the GC system of a galaxy.
#' @param n A non-negative numeric value. Sersic index of the GC system of a galaxy.
#' @param theta A numeric value in (0, 2*pi). Orientation angle (in radian) of the GC system of a galaxy.
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
