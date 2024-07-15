# simulate GC data
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
