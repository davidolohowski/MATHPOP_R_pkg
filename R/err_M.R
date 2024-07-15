# Photometric measurement uncertainty
 
err_M <- function(M, m0 = 25.5, b0 = 0.0884, b1 = 0.645){
  b0*exp(b1*(M - m0))
}
