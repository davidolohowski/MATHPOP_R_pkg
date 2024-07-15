# log likelihood function of MATHPOP model
log_lik_MATHPOP <- function(S, grid, l0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error){

  alpha <- cf_error$alpha
  m50 <- cf_error$m50
  beta0 <- cf_error$beta0
  beta1 <- cf_error$beta1
  m1 <- cf_error$m1

  if(l0 <= 0 || any(N < 0) || any(R_eff < 0) || any(sigma < 0)){
    return(-Inf)
  }
  else{
    A <- sf::st_area(sf::st_as_sf(S))
    K <- nrow(c)
    Theta <- list()
    for (i in 1:K) {
      Theta[[i]] <- list(c = c[i,], N = N[i], R_eff = R_eff[i], e = e[i], n = n[i], theta = theta[i], mu = mu[i+1], sigma = sigma[i+1])
    }

    norm_const <- sum(unlist(lapply(Theta, function(x){integrate_Sersic(grid, x$c, x$N, x$R_eff, x$e, x$n, x$theta)*p*p_eM_cpp(Lim = m50, x$mu, x$sigma, m0 = m1, b0 = beta0, b1 = beta1, a = alpha)})))
    L <- - norm_const - l0*A*sum(p)*p_eM_cpp(Lim = m50, mu[1], sigma[1], m0 = m1, b0 = beta0, b1 = beta1, a = alpha)
    loc <- l0*p*psi_f(Y_obs$M, mu[1], sigma[1], Lim = m50, m0 = m1, b0 = beta0, b1 = beta1, a = alpha)
    loc <- loc + vapply(data.table::transpose(lapply(Theta, function(x){Sersic_ints(Y_obs[,c('x','y')], x$c, x$N, x$R_eff, x$e, x$n, x$theta)*p*psi_f(Y_obs$M, x$mu, x$sigma, Lim = m50, m0 = m1, b0 = beta0, b1 = beta1, a = alpha)})), sum, 0)
    L <- L + sum(log(loc))
    return(L)
  }
}
