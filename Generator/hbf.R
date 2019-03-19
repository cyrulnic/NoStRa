
HHBEF <- function(ntime_filter, n, dt, stride, ind_obs_space, ind_time_anls, Rem, 
                UU, rrho, nnu, ssigma,  
                F_Lorenz, J_Lorenz, sd_noise,
                R_diag, m, OBS, 
                X_flt_start, Ba_start, Xae_start, B_clim,
                N, w_cvr, w_evp10, 
                inflation, spa_shift_max_Bf, spa_shift_max_S, C_lclz,
                model_type){
  #-----------------------------------------------------------------------------------
  # HHBEF: a Generalized HBEF, which includes EnKF, Var, EnVar, HBEF as special cases
  # and blends sample cvm with any or all of
  # (i) B_clim
  # (ii) time-smoothed covs
  # (iii) space-smoothed covs
  # As compared to the HBEF, in the HHBEF:
  # 
  # 1) The scnd filter treats B rather than P and Q separately,
  # 2) The fcst step of the scnd filter involves a "regression to the mean"
  #    ("climatological") background-error covariance matrix:
  #-------------------------------------------------------------------    
  #  B_f(t) = w_cvr*B_clim + (1-w_cvr)*B_a(t-1)               (1)
  #-------------------------------------------------------------------  
  #  whr B_clim is the static time-mean climatological CVM,
  #      B_a(t-1) is the final estimate of B at the previous time step, and
  #      w_cvr weighs the static CVM B_clim vs the recent-past CVM B_a(t-1).
  #  
  #  The anls step of the scnd flt is 
  #-------------------------------------------------------------------  
  # B_a(t) = (1-w_evp)*B_f(t) + w_evp*S(t)                   (2)
  #-------------------------------------------------------------------  
  # whr 
  # w_evp = N /(theta + N), find from 
  # w_evp10 = 10 / (theta +10)  ==>
  # theta = 10/w_evp10 - 10  ==>
  #---------------------------------- 
  # w_evp = N /(10/w_evp10 - 10 + N)
  #----------------------------------
  # w_pve=1-w_evp 
  #  whr N is the ensm size.
  # 
  # Note that the "anls-fcst" cycle starts here from an timestep-0 field X_flt_start
  # and one cycle is {(i) fcst, (ii) anls}.
  # 
  # NB: We assume that the model-error and obs-err statistics are perfect,
  #     therefore, the deterministic (control) forecast is preferred here 
  #     over the ensemble mean -- both in the fcst and anls.
  # 
  #--------------------------------------------------
  # Args:
  # 
  # ntime_filter - nu of ANLS (filter) time steps 
  #                (onefilter time step is
  #               a multiple (stride) of the model time step)
  # dt - MODEL time step (atmospheric time), sec.
  # stride - nu of model time steps between consecutive analyses
  # ind_obs_space - vector of indices of the state vector, where (in space) OBS are present
  # ind_time_anls - model time steps at which the anls is to be performed 
  #                (= 1, 1 + stride, 1 + 2*stride, ...)
  # 
  # UU, rrho, nnu,  ssigma - scnd flds for DSADM
  # F_Lorenz, J_Lorenz, sd_noise - Lorenz-2005 params
  # 
  # R_diag - diagonal (a vector) of the obs-err CVM
  # m - distance between adjacent obs in space (in grid meshes, integer)
  # OBS - obs at ALL MODEL time steps at the obs locations defined by ind_obs_space
  # X_flt_start - at time step 1, the fcst is to be started from X_flt_start
  # Ba_start - "anls" B at  the virtual time step 0
  # Xae_start - the "anls" ensemble at the virtual time step 0 
  #           (from which the 1st ensm fcst starts)
  # B_clim - static B
  # N - ensm size
  # w_cvr - in computing the prior Bf(t): the relative weight of the static "Climatological" CVM 
  #         vs the evolved Recent past CVM Ba(t-1), see Eq(1) above
  # w_evp10 - relative weight of the (localized) Ensemble sample CVM S vs the Prior Bf ("ensm-vs-prior"),
  #           see Eq(2) above 
  # inflation - covariance inflation coeficient 
  #       (defined as the multiplier of the fcts-ensm perturbations, i.e.
  #       the covariances are effectively multiplied by inflation^2)
  # spa_shift_max_Bf, spa_shift_max_S - number of spatial shifts to be made in one of the 
  #          two directionsin order to spatially smooth the cvm (Bf and S, resp.)
  #          with an even triangular weighting function
  # C_lclz - localization matrix 
  # model_type = "DSADM" or "Lorenz05" or "Lorenz05lin"
  # 
  # return: arrays at the anls times only.
  # 
  # M Tsyrulnikov (current code owner),
  # A Rakitko
  # 
  # Mar 2019
  #-----------------------------------------------------------------------------------
  
  ntime_model = ntime_filter * stride  # nu of filter (anls) time steps
  
  # Checks (only DSADM model permitted at the moment)
  
  if(model_type != "DSADM"){ # & model_type != "Lorenz05" & model_type != "Lorenz05lin"){
    print(model_type)  
    stop("HHBEF: wrong model_type")
  }
  
  #---------------------------------------------------
  # Preliminaries
  
  h = 2*pi*Rem/n                    # spatial mesh size
  
  XXf = matrix(NA, nrow = n, ncol = ntime_model)
  XXa = matrix(NA, nrow = n, ncol = ntime_model)
  B_mean = matrix(0,  nrow = n, ncol = n)
  #BBf = array(NA,dim=c(n, n, ntime_model))
  #BBa = array(NA,dim=c(n, n, ntime_model))
  #SS  = array(NA,dim=c(n, n, ntime_model))
  
  # Obs related variables
  
  n_obs=length(ind_obs_space)        # number of obs
  H = matrix(0, nrow=n_obs, ncol=n)  # obs oprt
  for (i in 1:n_obs){
    H[i, ind_obs_space[i]] = 1
  }
  
  R=diag(R_diag)                     # obs-err CVM
  sqrt_R_diag=sqrt(R_diag)           # obs-err st.devs
  
  # Starting (step-0) conditions
  
  Xa  = X_flt_start # the 1st deterministic fcst starts at time step 0 from this field
  Xae = Xae_start   # step-0 anls ensm
  Ba  = Ba_start    # step-0 anls estimate of B(t=0)
  
  w_evp = N / (10/w_evp10 - 10 + N) # ensm-vs-prior weight in computing Ba
  
  #---------------------------------------------------
  # The main loop over MODEL time steps
  
  ntm10 =floor(ntime_model /10)
  ntm100=floor(ntime_model /100)
  
  for(i in (1:ntime_model)){
    
    if(i %% ntm10 == 0){
      message(i / ntm100)
    }
    
    #-------------------------
    # (1) Fcst
    # (1.1) run deterministic fcst started from the previous anls
    
    N_det=1
    Xf = dsadm_step(Xa, n, N_det, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = FALSE)
     
    # (1.2) Fcst ensm
    
    Xfe = dsadm_step(Xae, n, N, dt, UU[,i], rrho[,i], nnu[,i], ssigma[,i], Rem, forcing = TRUE)
  
    
    # Separate model time steps i when the anls is to be or not to be performed:
    # ANLS are to be done at    i=stride*k +1,    where i=0,1,2,...
    # Therefore at the anls times, i-1 should divisible by stride:
    
    if(((i-1) %% stride) != 0){  # no anls, continue fcst
      
      Xa  = Xf
      Xae = Xfe
      
      
    }else{                       # time step when anls is performed
      
      #-------------------------
      # (1.3) Scnd flt B: fcst, compute the Prior B, i.e. Bf
      
      Bf = w_cvr*B_clim + (1-w_cvr)*Ba
      
      # Spa smoo Bf
      
      if(spa_shift_max_Bf > 0){
        C=Bf
        spa_shift_max = spa_shift_max_Bf
        Bf = SpaSmooCVM(n, C, spa_shift_max)
      }
      
      #-------------------------
      # (2) Anls
      
      # (2.0.0) Inflation
      
      dXfe=(Xfe - rep(Xf, N)) * inflation  # inflated fcst-ensm perturbations
      Xfe_inflated = rep(Xf, N) + dXfe
      
      # (2.0.1) Sample CVM
      
      S=dXfe %*% t(dXfe) /N  # non-shifted sample CVM
      
      # (2.0.2) Spa smoo S
      
      if(spa_shift_max_S > 0){
        C=S
        spa_shift_max = spa_shift_max_S
        S = SpaSmooCVM(n, C, spa_shift_max)
      }
      
      # (2.0.3) Covariance Localization

      S_lclz = S * C_lclz             # localization
      
      # (2.1) Scnd flt: B: anls, compute the Posterior B, i.e. Ba
      
      Ba = (1-w_evp)*Bf + w_evp*S_lclz
      
      # (2.3) KF anls with Ba
      
      # Kalman gain
      BHT  = Ba[             , ind_obs_space]  # B*H^T
      HBHT = Ba[ind_obs_space, ind_obs_space]  # H*B*H^T
      HBHTpR = HBHT + R
      K = BHT %*% solve(HBHTpR)
      
      # Deterministic anls
      Xa = Xf + K %*% (OBS[,i] - Xf[ind_obs_space])

      # (2.4) Ensm anls
      # (2.4.1) Simulate obs errs
      
      obsN01 = rnorm(n_obs*N, mean=0, sd=1) # N(0,1) noise
      simOBS_err = matrix(obsN01*sqrt_R_diag, nrow=n_obs, ncol=N)
      
      # (2.4.2) Generate perturbed obs
      
      simOBS = OBS[,i] + simOBS_err
      
       
      # (2.4.3) Anls ensm
      
      Xae = Xfe_inflated + K %*% (simOBS - Xfe_inflated[ind_obs_space,])
      
      # Averaging of B
      B_mean = B_mean + Ba
      
    }                # end anls time step
  
    #-------------------------
    # Store arrays
    
    XXf[,i] = Xf
    XXa[,i] = Xa
    #BBf[,,i] = Bf
    #BBa[,,i] = Ba
    #SS[,,i] = S
    
  }  # end time loop
  
  B_mean = B_mean / ntime_filter
  
  return(list(XXf = XXf[, ind_time_anls], 
              XXa = XXa[, ind_time_anls],
              B_mean = B_mean
              #SS  = SS [,,ind_time_anls], 
              #BBa = BBa[,,ind_time_anls],
              #BBf = BBf[,,ind_time_anls]
  ))

}
