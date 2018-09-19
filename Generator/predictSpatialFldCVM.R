predictSpatialFldCVM = function(ntime_filter, n, dt, stride, Rem, 
                                UU, rrho, nnu, ssigma,  
                                Bx_start){
  #-----------------------------------------------------------------------------------
  # Predict the spatial field (signal) CVM, Bx, by time stepping using the eqn
  # (see the implicit numerical scheme used to solve the time-continuous DSADM eqn)
  #............................... 
  # x_k = F_k * (x_k-1 + epsilon_k)                                         (1)
  #............................... 
  # Here
  # x is the truth (field, signal), everything is zero-mean, 
  # F_k = G_k^{-1} ...
  # and
  # 
  # epsilon_k[] = dt/h * sigma_k[] * N(0,I) == dt/h * Sigma_k * N(0,I)
  # 
  # Sigma_k = diag(sigma_k[]) (a matrix)  ==>
  #......................................................................
  # CVM(epsilon_k) = (dt/h)^2 * Sigma_k^2 == (dt/h)^2 * diag(sigma_k[]^2)    (2)
  #......................................................................
  # Then, from (1) and (2),
  #..............................................
  # Bx_k = F_k * [Bx_k-1 + CVM(epsilon_k)] * F_k              (*)
  #..............................................
  #  
  # Note that the first time step starts  from a timestep-0 CVM Bx_start
  # 
  # Args:
  # 
  # ntime_filter - nu of FILTER (ANLS) time steps 
  #               (the intvl between the consecutive analyses is
  #               a multiple (stride) of the model time step)
  # n - dim-ty of the state vector x
  # dt - MODEL time step (atmospheric time), sec.
  # stride - nu of model time steps between consecutive time instants
  #          when the returned array BBx is to be stored
  # Rem - Earth radius, m               
  # UU, rrho, nnu,  ssigma - scnd DSADM flds [space, time]
  # Bx_start - the imposed field CVM at the virtual time step 0
  # 
  # return: BBx[,,ind_time_store] -- at the anls time step only.
  # 
  # M Tsyrulnikov 
  # Jul 2018
  #-----------------------------------------------------------------------------------
  
  h           = 2*pi*Rem/n
  ntime_model = ntime_filter * stride
  
  BBx = array(0, dim=c(n, n, ntime_filter)) # every stride model time step starting from 1
  
  Bx_prev      = Bx_start
  Q_DSADM_diag = ssigma^2  /h *dt # model-error variances per model time step
  
  #---------------------------------------------------
  # main loop over model time steps

  i_filter=0
  
  for(i in 1:ntime_model){    
    
      # In the implicit DSADM time integration scheme, 
      # model error is added Before the mdl oprt is applied ==>
      
      CVM_forcing=diag( Q_DSADM_diag[,i] ) 
      BpQ = Bx_prev + CVM_forcing  
      
      # Bx = F * BpQ * F^T 
      # (1) PHI := F * BpQ
      # (2) B = PHI * F^T = (F * PHI^T)^T == F * PHI^T
      
      PHI = dsadm_step(BpQ,    n, n, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = FALSE)
      Bx  = dsadm_step(t(PHI), n, n, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = FALSE)
      
      # Store Bx
      
      if((i-1) %% stride == 0){
        i_filter=i_filter +1
        
        BBx[,,i_filter] = Bx
      } 

      Bx_prev = Bx

  }  # end time loop
  
  return(BBx)
}
