# completeness fraction
f <- function(M, Lim, alpha = 1.5){
  1/(1 + exp(alpha*(M-Lim)))
}
