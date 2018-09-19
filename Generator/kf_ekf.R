KF_EKF = function(ntime_filter, n, dt, stride, ind_obs_space, ind_time_anls, Rem, 
                  UU, rrho, nnu, ssigma,  
                  F_Lorenz, J_Lorenz, sd_noise,
                  R_diag, m, OBS, 
                  X_flt_start, A_start, 
                  model_type, filter_type,
                  predict_BB_KF){
  #-----------------------------------------------------------------------------------
  # KF/EKF
  # 
  # Note that the "anls-fcst" cycle starts here from an timestep-0 field X_flt_start
  # and one cycle is {(i) fcst, (ii) anls}.
  # 
  # Args:
  # 
  # ntime_filter - nu of ANLS (filter) time steps 
  #                (onefilter time step is
  #               a multiple (stride) of the model time step)
  # n - dim-ty of the state vector x
  # dt - MODEL time step (atmospheric time), sec.
  # stride - nu of model time steps between consecutive analyses
  # ind_obs_space - vector of indices of the state vector, where (in space) OBS are present
  # ind_time_anls - model time steps at which the anls is to be performed 
  #                (= 1, 1 + stride, 1 + 2*stride, ...)
  # Rem - Earth radius, m
  # UU, rrho, nnu,  ssigma - scnd flds for DSADM
  # F_Lorenz, J_Lorenz, sd_noise - Lorenz-2005 params
  # 
  # R_diag - diagonal (a vector) of the obs-err CVM
  # m - distance between adjacent obs in space (in grid meshes, integer)
  # OBS - obs at ALL MODEL time steps at the obs locations defined by ind_obs_space
  # X_flt_start - at time step 1, the fcst is to be started from X_flt_start
  # A_start - the imposed "anls-err CVM" at the virtual time step 0
  # model_type = "DSADM" or "Lorenz05" or "Lorenz05lin"
  # filter_type - "KF" or "EKF" 
  #     NB: with model_type="Lorenz05", only filter_type ="EKF" is acceptable
  # predict_BB_KF - compute & return prior filtering cvms (TRUE/FALSE)?
  # 
  # return: XXf, XXa, BB_KF (at the anls times only) and B_mean, A_mean.
  #         NB: BB_KF are computed only if predict_BB_KF=TRUE.
  # 
  # A Rakitko, 
  # M Tsyrulnikov (current code owner)
  # June 2018
  #-----------------------------------------------------------------------------------
  
  h = 2*pi*Rem/n
  
  ntime_model = ntime_filter * stride
  
  XXf    = matrix(NA, nrow = n, ncol = ntime_model)
  XXa    = matrix(NA, nrow = n, ncol = ntime_model)
  B_mean = matrix(0,  nrow = n, ncol = n)
  A_mean = matrix(0,  nrow = n, ncol = n)
  
  if(predict_BB_KF) {
    BB  = array(NA, c(n, n, ntime_filter))
  }else{                                   
    BB  = array(NA, c(n, n, 1))  # select a dummy array to save space
  }
  
  #AA         <- array(NA,c(n, n, ntime_model))
  
  n_obs=length(ind_obs_space)        # number of obs
  H = matrix(0, nrow=n_obs, ncol=n)  # obs oprt
  for (i in 1:n_obs){
    H[i, ind_obs_space[i]] = 1
  }
  
  R=diag(R_diag)                     # obs-err CVM
  
  Xa = X_flt_start # the 1st fcst starts at time step 0 from this field
  A  = A_start

  eps=1e-9 # EKF: Jacobian assessment through finite differences:
           #      scale down anls-CVM columns to reach the linear regime 1e-7..1e-9 are ok
 
  #---------------------------------------------------
  # Checks
  
  if(model_type != "DSADM" & model_type != "Lorenz05" & model_type != "Lorenz05lin"){
    print(model_type)  
    stop("KF_EKF: wrong model_type")
  }
  
  if(filter_type != "KF" & filter_type != "EKF"){
    print(filter_type)  
    stop("KF_EKF: wrong filter_type")
  }
  
  if(model_type == "Lorenz05" & filter_type == "KF"){
    print(model_type)
    print(filter_type)  
    stop("KF_EKF: wrong filter_type/model_type pair")
  }
  
  if(model_type == "Lorenz05lin" & filter_type == "EKF"){
    print(model_type)
    print(filter_type)  
    stop("KF_EKF: wrong filter_type/model_type pair")
  }
  
  #---------------------------------------------------
  # Lorenz: From atmospheric time to Lorenz time
  
  if(model_type == "Lorenz05" | model_type == "Lorenz05lin") {
    dt_atm_h = dt /3600
    dt_Lorenz = dt_atm_h/6*0.05    # unitless, "Lorenz time" 
                                   # (6h atmospheric time ~ 0.05 Lorenz time units)
    Q_Lorenz=sd_noise^2 * diag(n)  # Lorenz's Q
  }
  
  Q_DSADM_diag = ssigma^2  /h *dt # model-error variances per model time step
  
  #---------------------------------------------------
  # main loop over MODEL time steps
  
  i_filter=0
  
  for(i in 1:ntime_model){    
    
    # (1) Fcst
    # (1.1) run deterministic fcst started from the previous anls
    
    if(model_type == "DSADM"){ 
      N_det=1
      XXf[,i] = dsadm_step(Xa, n, N_det, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = FALSE)

    }else if(model_type == "Lorenz05"){
      XXf[,i] = lorenz05_step(Xa, n, dt_Lorenz, F_Lorenz, J_Lorenz, rep(0,n))
      
    }else if(model_type == "Lorenz05lin"){
      XXf[,i] = lorenz05lin_step(Xa, X_ref, n, dt_Lorenz, F_Lorenz, J_Lorenz, rep(0,n))
    }
    
    # (1.2) fcst covs
    
    if(model_type == "DSADM"){
      # In the implicit DSADM time integration scheme, 
      # model error is added Before the mdl oprt is applied ==>
      
      CVM_forcing=diag( Q_DSADM_diag[,i] ) 
      AQ = A + CVM_forcing   # A is the previous-cycle anls-err CVM
                             # NB: CVM_forcing is not exactly Q 
      
      # B = F * AQ * F^T 
      # (1) PHI := F * AQ
      # (2) B = PHI * F^T = (F * PHI^T)^T = F * PHI^T
      
      if(filter_type == "KF"){
        
        #PHI = apply( AQ,    2, function(x) dsadm_step(x, n, 1, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], 
        #                                              Rem, forcing = FALSE) )
        #B   = apply( t(PHI),2, function(x) dsadm_step(x, n, 1, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], 
        #                                              Rem, forcing = FALSE) ) 
        
        PHI = dsadm_step(AQ,     n, n, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = FALSE)
        B   = dsadm_step(t(PHI), n, n, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = FALSE)
      
      }else if(filter_type == "EKF"){  # for testing only
        PHI = apply( AQ,     2, 
                     function(x) ApplyJacobian_fd( dsadm_step, Xa, XXf[,i], x, eps,
                                                   n, n, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], 
                                                   Rem, forcing = FALSE) )
        B   = apply( t(PHI), 2, 
                     function(x) ApplyJacobian_fd( dsadm_step, Xa, XXf[,i], x, eps,
                                                   n, n, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], 
                                                   Rem, forcing = FALSE) )
      }
      

            
    }else if(model_type == "Lorenz05"){
      
      # In the Lorenz model, system noise is added after the fcst ==>
      # # B = F * A * F^T  + Q
      # (1) PHI := F * AQ
      # (2) B = PHI * F^T = (F * PHI^T)^T = F * PHI^T
      
      PHI = apply( A,  2, 
                   function(x) ApplyJacobian_fd(lorenz05_step, Xa, XXf[,i], x, eps,
                                                n, dt_Lorenz, F_Lorenz, J_Lorenz, rep(0,n)) )
      P   = apply( t(PHI), 2, 
                   function(x) ApplyJacobian_fd(lorenz05_step, Xa, XXf[,i], x, eps,
                                                n, dt_Lorenz, F_Lorenz, J_Lorenz, rep(0,n)) )
      P = (P + t(P)) /2 # eliminate computational non-symmetry
      B = P + Q_Lorenz
      
      
      
    }else if(model_type == "Lorenz05lin"){
      
      if(filter_type == "KF"){
        
        PHI = apply(AQ,    2, function(x) lorenz05lin_step(x, Xref, n, 
                                          dt_Lorenz, F_Lorenz, J_Lorenz, rep(0,n)) )
        P   = apply(t(PHI),2, function(x) lorenz05lin_step(x, Xref, n,  
                                          dt_Lorenz, F_Lorenz, J_Lorenz, rep(0,n)) )
        B = P + Q_Lorenz
      }
    }
    
    
    # (2) Anls
    
    # Separate model time steps when the anls is to be or not to be performed.
    # ANLS are to be done at    t=stride*k +1,    where k=1,2,3,...
    # Therefore at the anls times, t-1 should divisible by stride:
    
    if(((i-1) %% stride) != 0){  # no anls, continue fcst
      Xa = XXf[,i]
      A = B

    }else{                       # perform anls
      BHT  = B[             , ind_obs_space]  # B*H^T
      HBHT = B[ind_obs_space, ind_obs_space]  # H*B*H^T
      HBHTpR = HBHT + R
      K = BHT %*% solve(HBHTpR)
      
      Xa = XXf[,i] + K %*% (OBS[,i] - XXf[ind_obs_space,i])
    
      A = B - K%*%B[ind_obs_space,] # (I-KH)B
      
      # use  the Joseph form: 
      ## A=(I-KH)B(I-KH)^T + KRK^T 
      ## --> Yields the same results.
      
      #ImKH=diag(n) - K %*% H
      #A=ImKH %*% B %*% t(ImKH) + K %*% R %*% t(K)
    
      # store BB
      
      if(predict_BB_KF & (i-1) %% stride == 0){
      i_filter=i_filter +1
      
      BB[,,i_filter] = B 
      } 
      
      # Averaging of B, A
      
      B_mean = B_mean + B
      A_mean = A_mean + A
      
    }
    
    XXa[,i] = Xa
    #AA[,,i] = A
    
  }  # end time loop
  
  B_mean = B_mean / ntime_filter
  A_mean = A_mean / ntime_filter
  
  return(list(XXa   = XXa  [, ind_time_anls], 
              XXf   = XXf  [, ind_time_anls], 
              B_mean = B_mean,
              A_mean = A_mean,
              BB = BB 
              #AA = AA[,,ind_time_anls]
              ))
}
