#' Function to compute the intensity at the locations of a point pattern X given a Sersic intensity profile
#'
#' @param X A two-column matrix of the GC point pattern locations. The first column is the x coordinate, the second column is the y coordinate.
#' @param c A vector of length two for the central location of a galaxy. The first entry is the x coordinate of the galactic center, the second entry is the y coordinate.
#' @param N A non-negative numeric value. Mean number of GC in a galaxy.
#' @param R_eff A non-negative numeric value. Half-number radius of the GC system of a galaxy.
#' @param e A non-negative numeric value. Aspect ratio of the GC system of a galaxy.
#' @param n A non-negative numeric value. Sersic index of the GC system of a galaxy.
#' @param theta A numeric value in (0, 2*pi). Orientation angle (in radian) of the GC system of a galaxy.
#' @keywords Sersic Intensity
#' @export
#' @return A numeric vector with the length equal to the number of row of X that gives the Sersic intensity at the locations specified by X. 
Sersic_ints <- function(X, # Point pattern
                        c, # Galactic Center
                        N, # Mean number of GCs
                        R_eff, # Half-number radius
                        e = 1, # Ellipticity
                        n = 0.5, # Sersic index
                        theta = 0 # Orientation angle
){
  bn <- zipfR::Rgamma.inv(2*n, 0.5)
  r <- (((X[,1] - c[1])*cos(theta) - (X[,2]-c[2])*sin(theta))^2/R_eff^2 + ((X[,1] - c[1])*sin(theta) + (X[,2]-c[2])*cos(theta))^2/R_eff^2/e^2)^0.5
  return(exp(log(N) + (2*n)*log(bn) - bn*(r)^(1/n) - log(n) - lgamma(2*n) - log(2*pi) - 2*log(R_eff)- log(e)))
}
