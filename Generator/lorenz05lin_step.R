lorenz05lin_step = function(X_ref_start, n, U_start, M, dt, F_Lorenz, J, noise){
  
  # Linearized Lorenz-2005 Model-II.
  # Integrate FOR ONE TIME STEP (using 2nd-order Runge-Kutta).
  # Simultaneously, compute (and return!) the tangent linear oprt applied
  # to the initial perturbations U_start.
  # 
  # Both the reference trajectory and the forecast perturbations are computed here.
  #----------
  # Lorenz-05:
  # 
  # Set K<<n. Let K be Odd. Then, with J=(K-1)/2, the model is:
  # 
  # dx/dt=f(x), where x and f are the n-vectors
  #***********************************************************************
  # f(x)_i = - W(i-2*K) * W(i-K)  +  1/K * sum_{j=-J}^J  W(i-K+j)*X(i+K+j) - X(i) + F,  
  # where
  # W(i) = 1/K * sum_{j=-J}^J  X(i-j)
  #*********************************************************************** 
  # (i is the space index, all spatial indices are  modulo n)
  # 
  # Time scale: dt=0.05 corresponds to 6h in the atmosphere (Lorenz&Emanuel 1998)
  #---------
  # RK2 (2nd-order Runge-Kutta, the half-step version):
  # 
  # (1) x[k+1/2]=x[k] + dt/2 * f(x[k])
  # (2) x[k+1] = x[k] + dt * f(x[k+1/2])
  # (k is the time index)
  #=============================
  # Linearization.
  # 
  #  X = Xref + u
  #  W = Wref + v
  #   where
  # v(i) = 1/K * sum_{j=-J}^J  u(i-j)
  #  
  #             du/dt = Jacobian(Xref) * u = F(Xref, u),
  #   where
  #  
  #  F = -( Wref(i-2*K) * v(i-K)  + Wref(i-K) * v(i-2*K) ) +
  #      1/K * ( Wref(i+j-K) * u(i+j+K) + Xref(i+j+K) * w(i+j-K) ) - u(i)
  # 
  # Runge-Kutta for u.
  #=============================
  # Arguments:
  #
  # X_ref_start[1:n] - the initial condition at the current time step for the REFERENCE trajectory,
  #   it's also the reference point in phase space around which the nlin mdl is linearized.
  # n is the dim-ty of the state vector X
  # U_start[1:n, 1:M] - M perturbation initial conditions (X_m=X_ref + u[,m])
  # M - number of ini perturbations 
  # dt is the time step (unitless, "Lorenz time", dt=0.05 corresponds to 6h in the atmosphere)
  # F_Lorenz is the constant
  # J (or, equiv-ly, K=2*J+1) - integer, defines the spatial length scale (in the range 0...3)
  # noise[1:n, 1:M] -  noises to be added to the model field AFTER each time step (0 for deterministic fcst)
  #
  # NB: noise is added ONLY to the pertbns U, not the reference solution!
  # 
  # Return the forecasts - for both the reference solution Xref and the M perturbation solutions U.
  #
  # M Tsyrulnikov
  # May 2018
  #***************************************************************
  
  K=2*J+1
 
  # classical Runge-Kutta 2nd order, the half-step version
 
  #==================
  # half-step
  
  RHS_Lorenz1 = RHS_Lorenz05_lin(X_ref_start, n, U_start, M, F_Lorenz, J,K)
  
  # 1. Reference solution
  Xref_half = X_ref_start + (dt/2)*RHS_Lorenz1$RHS
  
  # 2. Perturbation solutions
  U_half    = U_start     + (dt/2)*RHS_Lorenz1$RHS_lin
  
  #==================
  # full-step
  
  RHS_Lorenz2 = RHS_Lorenz05_lin(Xref_half, n, U_half, M, F_Lorenz, J,K)
    
  # 1. Reference solution
  Xref_fcst = X_ref_start + dt * RHS_Lorenz2$RHS
 
  # 2. Perturbation solutions, with the noise added
  U_fcst    = U_start     + dt * RHS_Lorenz2$RHS_lin    + noise

  
  return(list(Xref_fcst=Xref_fcst, U_fcst=U_fcst))
}



RHS_Lorenz05_lin = function(x, n, U, M, F_Lorenz, J,K)  {
  
  # The rhs operator of the Lorenz-05 Model-II system
  # x is the n-vector of the current state of the system (the Reference solution)
  # U[1:n, 1:M] - M perturbation initial conditions (x_m=X + u[,m])
  # F_Lorenz is the constant.
  # J,K are two integer constants
  
  # First, compute the spatially smoothed X_i and U[i,], i.e. W_i and V[i,]:
  # W[i]    = 1/K * sum_{j=-J}^J  X[i-j]
  # V[i, m] = 1/K * sum_{j=-J}^J  U[i-j, m]
  # (V is the space-smoothed U)
  
  W=c(1:n)
  RHS=c(1:n)
  
  V       = matrix(0, nrow=n, ncol=M)
  RHS_lin = matrix(0, nrow=n, ncol=M)

  for(i in 1:n){
    ind=((i-J-1):(i+J-1)) %%n +1
    W[i]=mean(x[ind])
    V[i,1:M] = apply( U, 2, function(vec) mean(vec[ind]) )
  }
  
  # Compute the rhs:
  # RHS[i]= - W(i-2*K) * W(i-K)  +  1/K * sum_{j=-J}^J  W(i-K+j)*X(i+K+j) - X(i) + F
  # RHS_lin[i,m]    = - ( W[i-2*K] * V[i-K, m]   + W[i-K]  * V[i-2*K, m] ) +  
  #   1/K * sum_{j=-J}^J ( W[i+j-K] * U[i+j+K, m] + X[i+j+K]* V[i+j-K, m] ) - U[i,1:M]

  for(i in 1:n){
    #ind=((i-J-1):(i+J-1)) %%n +1 # averaging support for current i
    
    ind1=(i-(2*K)-1) %%n +1  # i-2K modulo n
    ind2=(i-   K -1) %%n +1  # i-K  modulo n
    
    iind3=((i-J-K-1):(i+J-K-1)) %%n +1  # i+j-K,  j \in [-J,J]
    iind4=((i-J+K-1):(i+J+K-1)) %%n +1  # i+j+K,  j \in [-J,J]
    
    RHS[i]= - W[ind1]*W[ind2] +  mean( W[iind3] * x[iind4]) - x[i] + F_Lorenz
    
    RHS_lin[i,1:M]= - W[ind1]*V[ind2,1:M] - W[ind2]*V[ind1,1:M] 
                    + apply( U, 2, function(vec) mean(W[iind3] * vec[iind4]) )
                    + apply( V, 2, function(vec) mean(x[iind4] * vec[iind3]) )
                    - U[i,1:M]
  }
  
  return(list(RHS=RHS, RHS_lin=RHS_lin))
}

