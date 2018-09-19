ext2mod_params4transform = function(pi_psi, kappa_psi) {
  
  #--------------------------------------------------------------------------
  # From the probability pi_psi=P(psi(t,s) <0), compute the g-function parameter epsilon_psi
  # computed as follows:
  # 
  # 1. Compute 
  # G = g [ Q(pi_psi) * log(kappa_psi)],
  # where g is the g-function (gfunction)
  # 2. Compute
  # epsilon_psi = G /(1-G)
  # 
  # Return:
  # 
  # epsilon_psi
  # 
  # M Tsy 2018 June
  #--------------------------------------------------------------------------
  
  
  G = gfunction( qnorm(pi_psi) * log(kappa_psi) )
  epsilon_psi = G /(1-G)
  
  return(epsilon_psi)
}