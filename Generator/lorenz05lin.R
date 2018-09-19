lorenz05lin <- function(X_ref_start, n, U_start, M, ntime, dt, F_Lorenz, J, noise_spatim){
  
  # Linearized Lorenz-2005 Model-II system
  # Integrate using RK2.
  # Simultaneously, compute (and return!) the reference trajectory.
  # Recursively call lorenz05lin_step.
  #==============
  # Arguments
  # 
  # X_ref_start[1:n] - the initial condition at the current time step for the REFERENCE trajectory,
  #   around which the nlin mdl is linearized.
  # n is the dim-ty of the state vector X
  # U_start[1:n, 1:M] - M perturbation initial conditions (X_m=X_ref + u[,m])
  # M - number of ini perturbations 
  # ntime is the number of time steps (inclu the starting step 1)
  # dt is the time step (unitless, "Lorenz time", dt=0.05 corresponds to 6h in the atmosphere)
  # F_Lorenz is the constant
  # J (or, equiv-ly, K=2*J+1) - integer, defines the spatial length scale (in the range 0...3)
  # noise_spatim[1:n, 1:M, 1:ntime] - white noise to be added to the model field AFTER each time step 
  # 
  # NB: noise_spatim is added ONLY to UUU, not the reference solution!
  # 
  # Return the forecasts - for both the reference solution Xref and the M perturbation solutions U.
  # Return the whole history X[:,:]
  #
  # M Tsyrulnikov
  # May 2018
  #***************************************************************

  XX_ref = matrix(0, nrow=n, ncol=ntime)
  UUU    = array(0, dim=c(n,M,ntime))
  
  XX_ref[,1] = X_ref_start   # reference ini cond
  UUU[,,1]   = U_start[,, drop=FALSE]
  
  K=2*J+1

  # recursively call the one-step solver
 
  for(k in 2:ntime){
    LIN = lorenz05lin_step(XX_ref[,k-1], n, as.matrix(UUU[,,k-1, drop=FALSE]), M, dt, F_Lorenz, J, 
                           as.matrix(noise_spatim[1:n, 1:M, k, drop=FALSE]))
    XX_ref[,k] = LIN$Xref_fcst
    UUU  [,,k] = LIN$U_fcst
  }
  
  return(list(XX_ref=XX_ref, UUU=UUU))
}