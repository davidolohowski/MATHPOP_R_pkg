#' @title Fitting the mixture model with noisy point source data
#' 
#' @description Function to obtain probabilistic GC catalog based on the parametric finite-mixture model in Li et al. (2024).
#' 
#' @param dat A data frame that contains the point source data used to obtain the probabilistic GC catalog. Should contain the columns `C` for GC color; `F814W` for GC magnitudes; `C_err` for measurement uncertainty in color; `M_err` for measurement uncertainty in magnitudes.
#' @param n_iter An integer for number of iteration to jitter the sources with their measurement uncertainty.
#' @param seed An integer for the random seed. Default to 12345.
#' @export
#' @return A list with three objects `prob`, `par`, and `sim`. `prob` contains the computed probability that a source is a GC. `par` contains the inferred parameter values for the finite-mixture model. `sim` contains the jittered point source data.
meas_uncertain_mix_func <- function(dat, n_iter, seed = 12345){
  set.seed(seed)
  p_mat <- matrix(0, ncol = n_iter, nrow = nrow(dat))
  par_mat <- matrix(0, ncol = 9, nrow = n_iter)
  sim_dat <- data.frame()
  for(i in 1:n_iter){
    sim_CM <- as.data.frame(t(apply(dat[,c('C', 'F814W', 'M_err', 'C_err')], 1,
                                    function(x){MASS::mvrnorm(n = 1, x[c(1,2)], Sigma = matrix(c(x[4]^2, x[3]^2, x[3]^2, x[3]^2), 2))})))

    data <- filter(sim_CM, C > 0.8 & C < 2.4)[,c('C', 'F814W')]
    sim_dat <- bind_rows(sim_dat, data)
    M <- data$F814W
    C <- data$C
    idx <- which(sim_CM$C > 0.8 & sim_CM$C < 2.4)
    res <- optim(c(0.6, 0.5, 1.2, 1.5, 0.16, 1.6, 0.1, 26.3, 0.3), mix_func, dat = data,
                 method = "L-BFGS-B", lower = c(0, 0, 0.9, 1.25, 0.1, 1.5, 0.1, 26, 0.1), upper = c(1, 1, 1.5, 1.5, 0.5, 1.8, 0.5, 27, 0.6))$par

    par_mat[i,] <- res

    w <- res[1]
    wr <- res[2]
    GCLF_TO <- 26.3
    GCLF_sig <- res[3]
    GC_color_mu_r <- res[4]
    GC_color_sig_r <- res[5]
    GC_color_mu_b <- res[6]
    GC_color_sig_b <- res[7]
    mu <- res[8]
    sigma <- res[9]

    GCLF_comp <- w*dnorm(M, GCLF_TO, GCLF_sig)*f_cpp(M, 26.69, 6.56)/Phi_f_cpp(26.69, GCLF_TO, GCLF_sig, a = 6.56)*
      (wr*dnorm(C, GC_color_mu_r, GC_color_sig_r) + (1-wr)*dnorm(C, GC_color_mu_b, GC_color_sig_b))
    cont_comp <- (1-w)*dnorm(M, mu, sigma)*dunif(C, 0.8, 2.4)

    prob <- GCLF_comp/(GCLF_comp + cont_comp)

    p_mat[idx,i] <- prob
  }
  return(list(prob = p_mat, par = par_mat, sim = sim_dat))
}
