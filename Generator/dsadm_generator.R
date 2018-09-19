dsadm_generator <- function(x_start, n, N, ntime, dt, U, rho, nu, sigma, Rem, forcing = TRUE){
  
  # Integrate DSADM FOR MULTIPLE TIME STEPs.
  # An ensemble of N members is generated.
  # Can be used as the FCST model (ie unforced) or as 
  # the generating model of TRUTH (forced), or as 
  # the generating model for an ensemble member (possible truth, forced).
  # In the two latter (forced) cases, the forcing is generated within this routine
  #  as the white noise multiplied by sigma(t,s).
  # By default, forcing is ON.
  #-------------------------------------------
  # Arguments:
  #
  # x_start[1:n, 1:N] - ensm of initial fields
  # n - dim-ty of the state vector x
  # N - ensm size
  # ntime - number of model time steps to be generated including the starting time instant
  #   (so that ntime should normally be >=2)
  # dt - time step, s
  # U[1:n], rho[1:n], nu[1:n], sigma[1:n] - coefficient fields of the DSADM
  # Rem - Earth radius, m
  # forcing - logical switch: if TRUE, then forcing is computed here and added to the solution,
  #                           if FALSE, the pure FCST is computed.
  # 
  # Return: the generated space-time field: field[space, time].
  #
  # M Tsyrulnikov
  # June 2018
  #***************************************************************
  
  field = matrix(0, nrow=n, ncol=ntime)
  field3D = array(0, dim=c(n, N, ntime))
  field3D[,,1] = x_start
  
  if(ntime > 1){
    for (t in 2:ntime){
      field3D[,,t] = dsadm_step(field3D[,,t-1], n, N, dt, U[,t], rho[,t], nu[,t], sigma[,t], Rem)
    }
  }
  
  # If N=1, return 2D array, otherwise 3D array
  
  if(N == 1){
    field = field3D[,1,]
  }else{
    field = field3D
  }
  return(field)
}
