
lclz_matrix <- function(n ,c, Rem){
  
  #---------------------------------------
  # Create a lclz mx using Gaspari-Cohn lclz function, see
  # Gaspari and Cohn (1999, Eq.(4.10))
  # 
  # Args
  # 
  # n - number of grid points on the circular domain
  # c - lclz length scale (NB: the lclz functions vaishes at distances > 2*c)
  # Rem - Earth radius, m
  #---------------------------------------
  
  C <- matrix(0, ncol = n, nrow = n)

  h = 2*pi*Rem/n
  for(i in 1:n){
    for(j in 1:n){
      
      if(abs(i-j)<n/2){
        z <- abs(i-j)             # greate-circle distance
      }else{
        z <- n - abs(i-j)         # greate-circle distance
      }
      z <- 2*sin(z*h/(Rem*2))*Rem # chordal distance
      
      if(z>=0 & z<=c)   C[i,j] <- -1/4*(z/c)^5+1/2*(z/c)^4+5/8*(z/c)^3-5/3*(z/c)^2+1
      if(z>=c & z<=2*c) C[i,j] <- 1/12*(z/c)^5-1/2*(z/c)^4+5/8*(z/c)^3+5/3*(z/c)^2-5*(z/c)+4-2/3*c/z
      # if(z>=2*c) C[i,j] <- 0
      
      # C[i,j] <- exp(-0.5*(z/c)^2)
    }
  }
  return(C)
}
