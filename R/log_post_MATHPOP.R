# log posterior of MATHPOP model
log_post_MATHPOP <- function(S, grid, l0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, prior, Ng, cf_error){
  if(Ng == 0)
  {
    lp <- log_lik_MATHPOP(S, grid, l0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error) +
      dlnorm(l0, prior$IGM$l0[1], prior$IGM$l0[2], log = T) +
      dnorm(mu[1], prior$IGM$mu[1], prior$IGM$mu[2], log = T) +
      dlnorm(sigma[1], prior$IGM$sigma[1], prior$IGM$sigma[2], log = T) +
      sum(dnorm(mu[-1], prior$UDG$mu[,1], prior$UDG$mu[,2], log = T)) +
      sum(dlnorm(sigma[-1], prior$UDG$sigma[,1], prior$UDG$sigma[,2], log = T)) + sum(log(sigma)) +
      sum(VGAM::dfoldnorm(N, mean = prior$UDG$N[,1], sd = prior$UDG$N[,2], log = T)) + sum(log(N)) +
      sum(dlnorm(R_eff, prior$UDG$R_eff[,1], prior$UDG$R_eff[,2], log = T)) + sum(log(R_eff)) +
      sum(dlnorm(n, prior$UDG$n[,1], prior$UDG$n[,2], log = T)) + sum(log(n))
    return(lp)
  }
  else{
    K <- nrow(c)
    lp <- log_lik_MATHPOP(S, grid, l0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error) +
      dlnorm(l0, prior$IGM$l0[1], prior$IGM$l0[2], log = T) +
      dnorm(mu[1], prior$IGM$mu[1], prior$IGM$mu[2], log = T) +
      dlnorm(sigma[1], prior$IGM$sigma[1], prior$IGM$sigma[2], log = T) +
      sum(dnorm(mu[-1], c(prior$gal$mu[,1], prior$UDG$mu[,1]), c(prior$gal$mu[,2], prior$UDG$mu[,2]), log = T)) +
      sum(dlnorm(sigma[-1], c(prior$gal$sigma[,1], prior$UDG$sigma[,1]), c(prior$gal$sigma[,2], prior$UDG$sigma[,2]), log = T)) + sum(log(sigma)) +
      sum(dlnorm(N[1:Ng], prior$gal$N[,1], prior$gal$N[,2], log = T)) +
      sum(VGAM::dfoldnorm(N[(Ng+1):K], prior$UDG$N[,1], prior$UDG$N[,2], log = T)) + sum(log(N)) +
      sum(dlnorm(R_eff, c(prior$gal$R_eff[,1], prior$UDG$R_eff[,1]), c(prior$gal$R_eff[,2], prior$UDG$R_eff[,2]), log = T)) + sum(log(R_eff)) +
      sum(dlnorm(n, c(prior$gal$n[,1], prior$UDG$n[,1]), c(prior$gal$n[,2], prior$UDG$n[,2]), log = T)) + sum(log(n))
    return(lp)
  }
}
