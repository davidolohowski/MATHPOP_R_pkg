# log-likelihood of parametric finite mixture model
mix_func <- function(par, dat){
  w <- par[1]
  wr <- par[2]
  GCLF_TO <- 26.3
  GCLF_sig <- par[3]
  GC_color_mu_r <- par[4]
  GC_color_sig_r <- par[5]
  GC_color_mu_b <- par[6]
  GC_color_sig_b <- par[7]
  mu <- par[8]
  sigma <- par[9]

  M <- dat$F814W
  C <- dat$C

  nll <- -sum(log(w*dnorm(M, GCLF_TO, GCLF_sig)*f_cpp(M, 26.69, 6.56)/Phi_f_cpp(26.69, GCLF_TO, GCLF_sig, a = 6.56)*(wr*dnorm(C, GC_color_mu_r, GC_color_sig_r) + (1-wr)*dnorm(C, GC_color_mu_b, GC_color_sig_b)) +
                    (1-w)*dnorm(M, mu, sigma)*dunif(C, 0.8, 2.4)))

  return(nll)
}
