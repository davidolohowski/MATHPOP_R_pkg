# Fit MATHPOP model with adaptive MCMC.
#' 
#' Function that fits the MATHPOP model
#' @import tidyverse
#' @import sf
#' @import sp
#' @import raster
#' @import spatstat
#' @import VGAM
#' @param Data A data.frame: Each row is an observed GC or point source in the data. Requires at least the columns named `x`, `y` for the spatial coordinates of GCs (in physical coordinates), and the magnitudes named `M`. If there are more columns, they need to be the probabilities a point source is a GC.
#' @param spat_dom A List. A list object containing a list called `vertices` that gives the vertices of the spatial domain, and a list of the number of integration grid called `n_grid`.
#' @param fixed_Theta A List. A list that specifies the known parameters of GC system, which contains two list objects `gal` and `UDG` that specify the respective known parameters of normal galaxies and UDGs. If there are no normal galaxies in the data, `gal` does not need to be specified.
#' @param prior A List. A list that specifies the parameter values of the prior distributions of the model parameters, specified in a similar fashion as `fixed_Theta`.
#' @param p A numeric value or vector. Crowding effect. Either a numeric value between \eqn{(0,1)}, or a numeric vector whose entries are all in \eqn{(0,1)} and length equal to `n_grid` in `spat_dom` with each numeric element being the crowding effect at the location of the spatial grid. In the current implementation, it is default to 1 (no crowding).
#' @param cf_error A List. List of parameters for completeness fraction and measurement uncertainties.
#' @param M An integer. Total number of iteration to run the MCMC algorithm.
#' @param Theta A List. Starting values of the MCMC chain. Default to `NULL`, and specified internally.
#' @param tune A List. Tuning parameters for initial MCMC pilot run.
#' @param n An integer. Initial MCMC pilot run iteration. Default to 1000.
#' @param prob_model Logical. Whether the GC data used is a probabilistic catalog or a binary catalog. Default to TRUE.
#' @param seed An integer. Random seed value. Default to 12345.
#' @param burnin A numeric value. A real number between 0 and 1. Percentage of the sample to be discarded as burn-in for MCMC. Default to 0.1
#' @keywords MATHPOP adaptive MCMC
#' @export
#' @return A data frame with `floor(burn_in*M)` number of rows that gives the posterior sample of the fitted MATHPOP model.
fit_MATHPOP <- function(Data, spat_dom, fixed_Theta, prior, p = 1, cf_error, M, Theta = NULL, tune, n = 1000, prob_model = TRUE, seed = 12345, burnin = 0.1){
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta || Current MCMC iteration: :current]",
                         total = M,
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 200)

  S <- sp::Polygon(spat_dom$vertices)
  S <- sp::SpatialPolygons(list(Polygons(list(S),'region')))
  S <- sp::SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

  grid <- sp::makegrid(S, n = spat_dom$n_grid)
  grid <- sp::SpatialPoints(grid, proj4string = CRS(proj4string(S)))
  grid <- raster::crop(grid, S)
  gridded(grid) <- T

  dx1 <- unname(grid@grid@cellsize[1])
  dx2 <- unname(grid@grid@cellsize[2])

  grid <- as.data.frame(grid)
  names(grid) <- c('x', 'y')
  grid <- as.matrix(grid)

  grid <- list(grid, c(dx1, dx2))

  if(is.null(fixed_Theta$gal)){
    center <- fixed_Theta$UDG$center
    e <- fixed_Theta$UDG$e
    theta <- fixed_Theta$UDG$theta
    UDG_ID <- fixed_Theta$UDG$UDG_ID
    Ng <- 0
    if(is.null(Theta)){
      set.seed(seed)
      Theta <- list(l0 = runif(1, exp(prior$IGM$l0) - 0.01, exp(prior$IGM$l0) + 0.01),
                    N = matrix(rnorm(nrow(center), log(10), 0.5), nrow = 1),
                    R_eff = matrix(rnorm(nrow(center), log(2), 0.5), nrow = 1),
                    n = matrix(rnorm(nrow(center), log(1), 0.75), nrow = 1),
                    mu = matrix(rnorm(nrow(center) + 1, c(prior$IGM$mu[1], prior$UDG$mu[,1]), 0.5), nrow = 1),
                    sigma = matrix(rnorm(nrow(center) + 1, c(prior$IGM$sigma[1], prior$UDG$sigma[,1]), 0.1), nrow = 1))
    }
  }
  else{
    center <- rbind(fixed_Theta$gal$center, fixed_Theta$UDG$center)
    e <- c(fixed_Theta$gal$e, fixed_Theta$UDG$e)
    theta <- c(fixed_Theta$gal$theta, fixed_Theta$UDG$theta)
    UDG_ID <- fixed_Theta$UDG$UDG_ID
    Ng <- nrow(fixed_Theta$gal$center)
    if(is.null(Theta)){
      set.seed(seed)
      Theta <- list(l0 = runif(1, exp(prior$IGM$l0) - 0.005, exp(prior$IGM$l0) + 0.005),
                    N = matrix(c(rnorm(Ng, prior$gal$N[,1], 0.25), rnorm(nrow(center) - Ng, log(10), 0.5)), nrow = 1),
                    R_eff = matrix(c(rnorm(Ng, prior$gal$R_eff[,1], 0.25), rnorm(nrow(center) - Ng, log(2), 0.5)), nrow = 1),
                    n = matrix(c(rnorm(Ng, log(0.5), 0.5), rnorm(nrow(center) - Ng, log(1), 0.75)), nrow = 1),
                    mu = matrix(rnorm(nrow(center) + 1, c(prior$IGM$mu[1], prior$gal$mu[,1], prior$UDG$mu[,1]), 0.5), nrow = 1),
                    sigma = matrix(rnorm(nrow(center) + 1, c(prior$IGM$sigma[1], prior$gal$sigma[,1], prior$UDG$sigma[,1]), 0.1), nrow = 1))
    }
  }

  K <- nrow(center)
  res <- matrix(0, nrow = M+1, ncol = 5*K+3)
  if(Ng > 0){
    col_names <- c('l0', paste0('N_gal_', 1:Ng), paste0('N_', UDG_ID), paste0('R_gal_', 1:Ng), paste0('R_', UDG_ID), paste0('n_gal_', 1:Ng), paste0('n_', UDG_ID), 'mu_0', paste0('mu_gal_', 1:Ng), paste0('mu_', UDG_ID), 'sigma_0', paste0('sigma_gal_', 1:Ng), paste0('sigma_', UDG_ID))
  }
  else{
    col_names <- c('l0', paste0('N_', UDG_ID), paste0('R_', UDG_ID), paste0('n_', UDG_ID), 'mu_0', paste0('mu_', UDG_ID), 'sigma_0', paste0('sigma_', UDG_ID))
  }
  colnames(res) <- col_names
  res[1,] <- c(Theta$l0, Theta$N, Theta$R_eff, Theta$n, Theta$mu, Theta$sigma)
  t_b <- tune$l0
  t_Ng <- tune$Ng
  t_Nu <- tune$Nu
  t_R <- tune$R
  t_n <- tune$n
  t_m <- tune$mu
  t_s <- tune$sigma
  gamma <- 2.38^2/(5*K+3)
  if(prob_model){
    Y_dat <- Data[, c('x', 'y', 'M')]
    Y_prob <- as.matrix(subset(Data, select = -c(x, y, M)))
    idx <- sample(1:ncol(Y_prob), M, replace = T)
  }
  for (i in 1:M) {
    if(i <= n){
      if(prob_model){
        Y_obs <- Y_dat[rbernoulli(nrow(Y_dat), p = Y_prob[,idx[i]]),]
      }
      else{
        Y_obs <- Data
      }
      if(Ng == 0)
      {
        l0_old <- res[i, 1]
        N_old <- exp(res[i, 2:(K+1)])
        R_eff_old <- exp(res[i,(K+2):(2*K+1)])
        n_old <- exp(res[i,(2*K+2):(3*K+1)])
        mu_old <- res[i, (3*K+2):(4*K+2)]
        sigma_old <- exp(res[i, (4*K+3):(5*K+3)])

        l0_new <- rnorm(1, l0_old, t_b)
        N_new <- exp(rnorm(K, log(N_old), t_Nu))
        R_eff_new <- exp(rnorm(K, log(R_eff_old), t_R))
        n_new <- exp(rnorm(K, log(n_old), t_n))
        mu_new <- rnorm(K+1, mu_old, t_m)
        sigma_new <- exp(rnorm(K+1, log(sigma_old), t_s))
      }
      else{
        l0_old <- res[i, 1]
        N_old_g <- exp(res[i, 2:(Ng+1)])
        N_old_u <- exp(res[i, (Ng+2):(K+1)])
        N_old <- c(N_old_g, N_old_u)
        R_eff_old <- exp(res[i,(K+2):(2*K+1)])
        n_old <- exp(res[i,(2*K+2):(3*K+1)])
        mu_old <- res[i, (3*K+2):(4*K+2)]
        sigma_old <- exp(res[i, (4*K+3):(5*K+3)])

        l0_new <- rnorm(1, l0_old, t_b)
        N_new_g <- exp(rnorm(Ng, log(N_old_g), t_Ng))
        N_new_u <- exp(rnorm(K-Ng, log(N_old_u), t_Nu))
        N_new <- c(N_new_g, N_new_u)
        R_eff_new <- exp(rnorm(K, log(R_eff_old), t_R))
        n_new <- exp(rnorm(K, log(n_old), t_n))
        mu_new <- rnorm(K+1, mu_old, t_m)
        sigma_new <- exp(rnorm(K+1, log(sigma_old), t_s))
      }

      r <- log_post_MATHPOP(S, grid, l0_new, center, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, Ng, cf_error = cf_error) -
        log_post_MATHPOP(S, grid, l0_old, center, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, Ng, cf_error = cf_error)

      C <- cov(res[c(1:n),])

      Mu <- colMeans(res[c(1:n),])

      if(log(runif(1)) < r){
        res[i+1, ] <- c(l0_new, log(N_new),
                        log(R_eff_new), log(n_new), mu_new, log(sigma_new))
      }
      else{
        res[i+1,] <- res[i,]
      }
      pb$tick()
    }
    else{
      if(prob_model){
        Y_obs <- Y_dat[rbernoulli(nrow(Y_dat), p = Y_prob[,idx[i]]),]
      }
      else{
        Y_obs <- Data
      }
      if(Ng == 0){
        l0_old <- res[i, 1]
        N_old <- exp(res[i, 2:(K+1)])
        R_eff_old <- exp(res[i,(K+2):(2*K+1)])
        n_old <- exp(res[i,(2*K+2):(3*K+1)])
        mu_old <- res[i, (3*K+2):(4*K+2)]
        sigma_old <- exp(res[i, (4*K+3):(5*K+3)])

        X <- res[i,]
        C <- (i-2)/(i-1)*C + (X-Mu)%*%t(X-Mu)/i
        Mu <- (Mu*(i-1) + X)/i

        Theta_new <- MASS::mvrnorm(1, mu = X,
                                   Sigma = gamma*(C + diag(1e-7, nrow = 5*K+3)))

        l0_new <- Theta_new[1]
        N_new <- exp(Theta_new[2:(K+1)])
        R_eff_new <- exp(Theta_new[(K+2):(2*K+1)])
        n_new <- exp(Theta_new[(2*K+2):(3*K+1)])
        mu_new <- Theta_new[(3*K+2):(4*K+2)]
        sigma_new <- exp(Theta_new[(4*K+3):(5*K+3)])
      }
      else{
        l0_old <- res[i, 1]
        N_old_g <- exp(res[i, 2:(Ng+1)])
        N_old_u <- exp(res[i, (Ng+2):(K+1)])
        N_old <- c(N_old_g, N_old_u)
        R_eff_old <- exp(res[i,(K+2):(2*K+1)])
        n_old <- exp(res[i,(2*K+2):(3*K+1)])
        mu_old <- res[i, (3*K+2):(4*K+2)]
        sigma_old <- exp(res[i, (4*K+3):(5*K+3)])

        X <- res[i,]
        C <- (i-2)/(i-1)*C + (X-Mu)%*%t(X-Mu)/i
        Mu <- (Mu*(i-1) + X)/i

        Theta_new <- MASS::mvrnorm(1, mu = X,
                                   Sigma = gamma*(C + 1e-7*diag(1, nrow = 5*K+3)))

        l0_new <- Theta_new[1]
        N_new_g <- exp(Theta_new[2:(Ng+1)])
        N_new_u <- exp(Theta_new[(Ng+2):(K+1)])
        N_new <- c(N_new_g, N_new_u)
        R_eff_new <- exp(Theta_new[(K+2):(2*K+1)])
        n_new <- exp(Theta_new[(2*K+2):(3*K+1)])
        mu_new <- Theta_new[(3*K+2):(4*K+2)]
        sigma_new <- exp(Theta_new[(4*K+3):(5*K+3)])
      }

      r <- log_post_MATHPOP(S, grid, l0_new, center, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, Ng, cf_error = cf_error) -
        log_post_MATHPOP(S, grid, l0_old, center, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, Ng, cf_error = cf_error)

      if(log(runif(1)) < r){
        res[i+1, ] <- c(l0_new, log(N_new),
                        log(R_eff_new), log(n_new), mu_new, log(sigma_new))
      }
      else{
        res[i+1,] <- res[i,]
      }
      pb$tick()
    }
  }
  res[,2:(K+1)] <- exp(res[,2:(K+1)])
  res[,(K+2):(2*K+1)] <- exp(res[,(K+2):(2*K+1)])
  res[,(2*K+2):(3*K+1)] <- exp(res[,(2*K+2):(3*K+1)])
  res[,(4*K+3):(5*K+3)] <- exp(res[,(4*K+3):(5*K+3)])
  res <- res[-c(1:(floor(0.1*M)+1)), ]
  return(res)
}

