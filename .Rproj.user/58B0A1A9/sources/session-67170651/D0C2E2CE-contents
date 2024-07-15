# integrate Sersic intensity over a domain
integrate_Sersic <- function(grid, c, N, R_eff, e = 1, n = 0.5, theta = 0){
  dx1 <- grid[[2]][1]
  dx2 <- grid[[2]][2]

  return(unname(sum(Sersic_ints(grid[[1]], c, N, R_eff, e, n, theta))*dx1*dx2))
}
