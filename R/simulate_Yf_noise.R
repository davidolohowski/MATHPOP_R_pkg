# simulate GC data after removing faint GCs
simulate_Yf_noise <- function(Y, Lim, alpha = 1.5){
  return(Y[rbernoulli(nrow(Y), p = f(Y$M, Lim, alpha)),])
}
