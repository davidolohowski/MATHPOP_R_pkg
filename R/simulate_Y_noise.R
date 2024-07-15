#' Simulate GC locations and magnitudes with noisy magnitudes based on a Sersic profile and Gaussian GCLF.
#' 
#' @param S A `SpatialPolygonsDataFrame` that gives the spatial domain on which the point pattern resides.
#' @param l0 A positive numeric value. The intensity (/kpc^2) of GCs in the IGM.
#' @param c A matrix with two columns for the central locations of galaxies in `S`. The first column is the x coordinate of the galactic centers, the second column is the y coordinate.
#' @param N A (non-negative) numeric vector with length `nrow(c)`. Each element is the mean number of GCs in a galaxy.
#' @param R_eff A (non-negative) numeric vector with length `nrow(c)`. Each element is the half-number radius of the GC system of a galaxy.
#' @param e A (non-negative) numeric vector with length `nrow(c)`. Each element is the aspect ratio of the GC system of a galaxy.
#' @param n A (non-negative) numeric vector with length `nrow(c)`. Each element is the Sersic index of the GC system of a galaxy.
#' @param theta A numeric vector with length `nrow(c)`. Each element is within (0, 2*pi) and is the orientation angle of the GC system of a galaxy.
#' @param mu A numeric vector with length `nrow(c) + 1`. Each element is the GCLF TO point of each of the `nrow(c) + 1` GC sub-populations.
#' @param sigma A numeric vector with length `nrow(c) + 1`. Each element is the GCLF dispersion of each of the `nrow(c) + 1` GC sub-populations.
#' @param m0 A numeric value. The offset magnitude of the photometric measurement uncertainty.
#' @param b0 A numeric value. The uncertainty when `M = m0`.
#' @param b1 Rate at which the uncertainty increases with `M`.
#' @keywords Simulate Point Patterns Sersic
#' @export
simulate_Y_noisy <- function(S, l0, c, N, R_eff, e, n, theta, mu, sigma, m0 = 25.5, b0 = 0.0884, b1 = 0.645){
  A <- sf::st_area(st_as_sf(S))
  K <- nrow(c)
  N_b <- floor(l0*A)
  if(N_b == 0){
    Y <- data.frame()
  }
  else{
    Y0 <- as.data.frame(sp::spsample(S, N_b, type = 'random'))
    M <- rnorm(N_b, mu[1], sigma[1])
    M_err <- rnorm(N_b, M, err_M(M, m0, b0, b1))
    Y <- data.frame(x = Y0[,1], y = Y0[,2], M = M_err, Mt = M, id = rep(0,N_b), e_M = err_M(M, m0, b0, b1))
  }
  for (i in 1:K) {
    if(N[i] == 0)
    {
      Y <- Y
    }
    else{
      M <- rnorm(N[i], mu[i+1], sigma[i+1])
      M_err <- rnorm(N[i], M, err_M(M, m0, b0, b1))
      Y <- bind_rows(Y, data.frame(as.data.frame(simulate_Sersic(N[i], c[i,], R_eff[i], e[i], n[i], theta[i])), M = M_err, Mt = M, id = rep(i, N[i]), e_M = err_M(M, m0, b0, b1)))
    }
  }
  Y <- Y[unlist(sf::st_intersects(sf::st_as_sf(S), sf::st_as_sf(Y, coords = c('x','y')))),]
  return(Y)
}
