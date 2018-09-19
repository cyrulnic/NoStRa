lorenz05 <- function(X_start, n, ntime, dt, F_Lorenz, J, noise_spatim){
  
  # Lorenz-2005 Model-II system
  # Integrate using 2nd-order Runge-Kutta
  # Call lorenz05_step
  #
  # Arguments:
  #
  # X_start[1:n] is the initial condition
  # n is the dim-ty of the state vector X
  # ntime is the number of time steps (inclu the starting step 1)
  # dt is the time step (unitless, dt=0.05 corresponds to 6h in the atmosphere)
  # F_Lorenz is the constant
  # J (or, equiv-ly, K) - integer, defines the spatial length scale (in the range 0...3)
  # noise_spatim[1:n, 1:ntime] - white noise to be added to the model field AFTER each time step 
  # 
  # Return the whole history X[:,:]
  #  
  #***************************************************************
  # Specification of initial conditions.
  # NB: In the steady-state regime, 
  #     the mean X value for Lorenz96 (and about so in Lorenz05) is about 
  # 
  # assumed_mean = 1.2 * F_Lorenz^(1/3) # mean, from Lorenz-05, p.1577, bef Eq(5)
  # assumed_ms = assumed_mean * F_Lorenz # mean square, from Lorenz-05, p.1577, bef Eq(5)
  # assumed_var=assumed_ms - assumed_mean^2
  # assumed_sd = sqrt(assumed_var)
  # X1=rnorm(n, mean=assumed_mean, sd=assumed_sd) # ini condition
  #
  # M Tsyrulnikov
  # May 2018
  #***************************************************************

  source('lorenz05_step.R')
  
  X=matrix(0, nrow=n, ncol=ntime)
  X[,1] = X_start   # ini cond
  
  K=2*J+1

  # recursively call the one-step solver
 
  for(k in 2:ntime){
    X[,k] = lorenz05_step(X[,k-1], n, dt, F_Lorenz, J, noise_spatim[,k])
  }
  
  return(X)
}
