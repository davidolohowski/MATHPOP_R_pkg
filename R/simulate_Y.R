#' Simulate GC locations and magnitudes (without measurement uncertainty) based on a Sersic profile and Gaussian GCLF.
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
#' @keywords Simulate Point Patterns Sersic
#' @export
#' @return A data frame with four columns. `x`, `y` columns give the locations of the simulated GCs. `M` gives the GC magnitudes. `id` is the identifier for GC sub-population; `id = 0` means GCs from the IGM, while others are from galaxies.
simulate_Y <- function(S, l0, c, N, R_eff, e, n, theta, mu, sigma){
  A <- sf::st_area(st_as_sf(S))
  K <- nrow(c)
  N_y <- rpois(K, N)
  N_b <- rpois(1, l0*A)
  if(N_b == 0){
    Y <- data.frame()
  }
  else{
    Y0 <- as.data.frame(sp::spsample(S, N_b, type = 'random'))
    Y <- data.frame(x = Y0[,1], y = Y0[,2], M = rnorm(N_b, mu[1], sigma[1]), id = rep(0,N_b))
  }
  for (i in 1:K) {
    if(N_y[i] == 0)
    {
      Y <- Y
    }
    else{
      Y <- bind_rows(Y, data.frame(as.data.frame(simulate_Sersic(N_y[i], c[i,], R_eff[i], e[i], n[i], theta[i])), M = rnorm(N_y[i], mu[i+1], sigma[i+1]), id = rep(i, N_y[i])))
    }
  }
  Y <- Y[unlist(sf::st_intersects(st_as_sf(S), st_as_sf(Y, coords = c('x','y')))),]
  return(Y)
}
