fft_Rspe_rearrange = function(Rspe_vector) {
  
  # For plotting purposes, rearrange the R fft spectral cplx n-vector
  # f0, f1, ..., f(n/2),f(-n/2 +1),..., f(-1)
  # to the wvn-symmetric form (an (n+1) cplx vector)
  # f(n/2),f(-n/2 +1), ..., f(-1), f(0), f(1),..., f(n/2)
  #
  # NB: n needs to be even.
  #
  # M Tsy 3 Jul 2017
  
  n = length(Rspe_vector)
  
  if(n %% 2 == 1){
    message("Odd n")
    print(n)  
    stop("Make n an even number")
  }
  
  nd2=n/2
  np1=n+1
  
  symm_spe_vector=c(1:np1)
  
  symm_spe_vector[1]  =Rspe_vector[nd2 +1]  # m=-n/2
  symm_spe_vector[np1]=Rspe_vector[nd2 +1]  # m=n/2
  
  symm_spe_vector[2:nd2]=Rspe_vector[(nd2 +2):n]      # m<0
  
  symm_spe_vector[(nd2+1):np1]=Rspe_vector[1:(nd2+1)] # m >= 0
  
  return(symm_spe_vector)
}
  
  
  
  