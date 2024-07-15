# logit function
logit <- function(u, a, b){
  p <- pmin(1, pmax(0, (u-a)/(b-a)))
  return(pmin(999, pmax(-999,log(p/(1-p)))))
}
