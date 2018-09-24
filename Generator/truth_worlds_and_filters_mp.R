# Run DSADM, Lorenz-2005, Lorenz-2005-TLM, and 
# filters: KF, EKF, EnKF, Var, EnVar, HBEF, and  
# HBF, which includes all the above filters as special cases.
# Write output files to be examined by other programs.
# 
# A Rakitko
# M Tsyrulnikov (current code owner)
# 17 Jul 2018

library(plot3D)
library(limSolve)

source("ext2mod_params4dsadm_statio.R")
source("ext2mod_params4transform.R")
source('dsadm_start_field_cvm.R')
source('dsadm_step.R')
source('dsadm_generator.R')
source('predictSpatialFldCVM.R')
source('kf_ekf.R')
source('hbf.R')
source('lclz_matrix.R')
source('symm_pd_mx_sqrt.R')
source("gfunction.R")
source("transform_function.R")
source("rmse.R")
source('lorenz05.R')
source('lorenz05_step.R')
source('lorenz05lin.R')
source('lorenz05lin_step.R')
source('ApplyJacobian_fd.R')


path = './Out'

# Check that the two dirs do exist: ./Out and ./Out/Data.
# Otherwise, uncomment the following two lines:
#dir.create(path)
#dir.create(paste0(path,'/DATA'))
# -----------

config <- read.table('./config.txt', sep = ';')

# NB: time_filter is the number of assimilation steps (analyses), whereas
#     time_model is the number of model time steps
#
#     dt_h is the model time step (in hours), so that
#     the time interval between the consecutive analyses is dt_h*stride (h)

mode            = config[config$V1 == "mode", 2]
n               = config[config$V1 == "n", 2]
stride          = config[config$V1 == "stride", 2]   # model time steps in one asml cycle
time_filter     = config[config$V1 == "time_filter", 2] # nu of filter time steps (analyses)
dt_h            = config[config$V1 == "dt_h", 2]
U_mean          = config[config$V1 == "U_mean", 2]
V_char          = config[config$V1 == "V_char", 2]
L_mult          = config[config$V1 == "L_mult", 2]
L_perturb_mult  = config[config$V1 == "L_perturb_mult", 2]
sd_x            = config[config$V1 == "sd_x", 2]
sd_U            = config[config$V1 == "sd_U", 2]
kappa_rho       = config[config$V1 == "kappa_rho", 2]
pi_rho          = config[config$V1 == "pi_rho", 2]
kappa_nu        = config[config$V1 == "kappa_nu", 2]
pi_nu           = config[config$V1 == "pi_nu", 2]
kappa_sigma     = config[config$V1 == "kappa_sigma", 2]
M               = config[config$V1 == "M", 2] # number of replicates (worlds)
N               = config[config$V1 == "N", 2]
m               = config[config$V1 == "m", 2] #OBS GRID MESH SIZE
sqrt_R          = config[config$V1 == "sqrt_R", 2]
lclz_mult       = config[config$V1 == "lclz_mult", 2]
inflation       = config[config$V1 == "inflation", 2]
seed            = config[config$V1 == "seed", 2]
perform_kf_ekf  = config[config$V1 == "perform_kf_ekf", 2]
KF_SaveClim     = config[config$V1 == "KF_SaveClim", 2]
perform_HBF     = config[config$V1 == "perform_HBF", 2]
HBF_SelfClim    = config[config$V1 == "HBF_SelfClim", 2]
w_cvr           = config[config$V1 == "w_cvr", 2]
w_evp10         = config[config$V1 == "w_evp10", 2]
F_Lorenz        = config[config$V1 == "F_Lorenz", 2]
J_Lorenz        = config[config$V1 == "J_Lorenz", 2]        
sd_noise_Lorenz = config[config$V1 == "sd_noise_Lorenz", 2]
model_type      = config[config$V1 == "model_type", 2]

# seeds

seed_for_secondary_fields = seed
seed_for_filters          = seed + 12345

# model selector 

if(model_type == 1) {
  model_type="DSADM"
  
} else if(model_type == 2) {
  model_type="Lorenz05"
  
} else if(model_type == 3) {
  model_type="Lorenz05lin"
  
} else{
  print(model_type)  
  stop("KF_EKF: wrong model_type")
}

# KF filter type selector

if(perform_kf_ekf == 1) {
  filter_type="KF"
  
} else if(perform_kf_ekf == 2) {
  filter_type="EKF"
  
}else{
  filter_type="None"
}

#------------------------------------------------------
# mode-dependent parameters
# mode=0: single mdl&flt runs (produce Fields only), 
# mode=1: predicted BBx (fieds' Covariances), 
# mode=2: predicted BB_KF (KF's background-error Covariances),
# mode=3: worlds-averaged Fields' Covariances 
# mode=4: worlds-averaged HBF's   Covariances 

# By default and with mode=0: 
# switch off time-specific Bx & B_flt computations (both recurrent and by worlds-averaging)

if(mode == 0) {
  predict_BBx   = FALSE
  predict_BB_KF = FALSE
  worlds           = FALSE
  worldsAve_BBx    = FALSE
  worldsAve_BB_HBF = FALSE
}

# mode>0: ensure that ntime is small

if(mode > 0){           # single runs of the truth and the filters
  #if(time_model > 400) time_model = 400
} 

# mode=1: predict BBx

if(mode == 1){          # estm CVM of x by Predicting them
  M=1
  predict_BBx   = TRUE
  predict_BB_KF = FALSE
  worlds           = FALSE
  worldsAve_BBx    = FALSE
  worldsAve_BB_HBF = FALSE
  
  perform_kf_ekf=-1
  perform_HBF=-1
}

# mode=2: predict BB_KF

if(mode == 2){          # estm CVM of KF' background errors by Predicting them
  M=1
  predict_BBx   = FALSE
  predict_BB_KF = TRUE
  worlds           = FALSE
  worldsAve_BBx    = FALSE
  worldsAve_BB_HBF = FALSE
  
  perform_kf_ekf=1
  perform_HBF=-1
  
}

# mode=3: estm BBx by worlds-averaging

if(mode == 3 & M > 1){ # estm BBx_worldsAve (by averaging over worlds)
  predict_BBx   = FALSE
  predict_BB_KF = FALSE
  worlds           = TRUE
  worldsAve_BBx    = TRUE
  worldsAve_BB_HBF = FALSE
  
  perform_kf_ekf=-1
  perform_HBF=-1
}

# mode=4: estm BB_HBF by worlds-averaging

if(mode == 4 & M > 1){ # estm BB_HBF_worldsAve (by averaging over worlds)
  predict_BBx   = FALSE
  predict_BB_KF = FALSE
  worlds           = TRUE
  worldsAve_BBx    = FALSE
  worldsAve_BB_HBF = TRUE
  
  perform_kf_ekf=-1
  perform_HBF=1
}

#------------------------------------------------------
# time_model & synonyms

time_model = time_filter * stride  # nu of model time steps

time=time_model
ntime_model=time_model
ntime_cycles=time_filter
ntime_filter=time_filter
ntime_model_time_steps=time_model
dt=dt_h *3600 # s
dt_model = dt
dt_filter = dt * stride


#------------------------------------------------------
# calculated basic parameters

Rem       = 6.37e6           #     Earth radius, m
a         = 1 / sqrt(2*pi*Rem)
delta_s   = 2*pi*Rem/n       # spatial mesh, m
L         = L_mult*delta_s   # mean spatial len scale of x, m
L_perturb = L*L_perturb_mult # spatial len scale of the scnd flds, m
T         = L / V_char       # mean time scale for x, s

ExtClim   = -HBF_SelfClim    # +-1; B_clim is taken from a KF run: short (-1) or long (+1)
# NB: SelfClim for time_filter=10,000 seems to be almost as good as a 100,000 B_clim.

#------------------------------------------------------
# Derived-type variables containing external parameters
# 
# NB: U_i      is the same for all tertiary==pre-secondary fields and the unperturbed model.
#     V_char_i is the same for all tertiary==pre-secondary fields and the unperturbed model.
#     
# 1) Tertiary fields.
# Each tertiary==pre-secondary field theta_i (i=1,2,3,4) is characterized by the four External Parameters:
# (U_i, L_i, V_char_i, SD_i),
# which are to be converted to the four Model Parameters:
#  (U_i, rho_i, nu_i, sigma_i)
#  
# (i=1) theta_1* = U*

tert1_extpar = list(U=U_mean, L=L_perturb, V_char=V_char, SD=sd_U)

# (i=2) theta_2* =  rho*

sd_tert2=log(kappa_rho)
tert2_extpar = list(U=U_mean, L=L_perturb, V_char=V_char, SD=sd_tert2)


# (i=3) theta_3* =  nu*

sd_tert3=log(kappa_nu)
tert3_extpar = list(U=U_mean, L=L_perturb, V_char=V_char, SD=sd_tert3)


# (i=4) theta_4* = sigma*

sd_tert4=log(kappa_sigma)
tert4_extpar = list(U=U_mean, L=L_perturb, V_char=V_char, SD=sd_tert4)

# 2) unperturbed model for x

x_unpert_extpar = list(U=U_mean, L=L, V_char=V_char, SD=sd_x)

#---------
# From ext params (U, L, V_char, SD) to model params (U, rho, nu, sigma)

tert1_modpar = ext2mod_params4dsadm_statio(tert1_extpar, Rem)
tert2_modpar = ext2mod_params4dsadm_statio(tert2_extpar, Rem)
tert3_modpar = ext2mod_params4dsadm_statio(tert3_extpar, Rem)
tert4_modpar = ext2mod_params4dsadm_statio(tert4_extpar, Rem)
x_unpert_modpar = ext2mod_params4dsadm_statio(x_unpert_extpar, Rem)

#---------
# Transform params:

epsilon_rho = ext2mod_params4transform(pi_rho, kappa_rho)
epsilon_nu  = ext2mod_params4transform(pi_nu,  kappa_nu )

#------------------------------------------------------
# Obs-err CVM (diagonal) & lclz radius

R_diag = rep(sqrt_R^2, n%/%m+min(1,n%%m)) # vector on the diagonal of the R mx
L_loc  = lclz_mult*L                      # lclz radius, m

#------------------------------------------------------
# initialize the output variable "filter"

filter = list() # the output variable

# create list for parameters (for output file)

parameters = list()

parameters$mode                 <- mode
parameters$n                    <- n
parameters$stride               <- stride 
parameters$time_filter          <- time_filter
parameters$dt_filter            <- dt_filter  # filter time step, s
parameters$U_mean               <- U_mean
parameters$V_char               <- V_char
parameters$sd_x                 <- sd_x
parameters$kappa_rho            <- kappa_rho
parameters$pi_rho               <- pi_rho
parameters$kappa_nu             <- kappa_nu
parameters$pi_nu                <- pi_nu
parameters$sd_U                 <- sd_U
parameters$sqrt_R               <- sqrt_R
parameters$L_loc                <- L_loc
parameters$kappa_sigma          <- kappa_sigma
parameters$M                    <- M
parameters$mesh_obs             <- m
parameters$L_mean               <- L
parameters$L_perturb            <- L_perturb
parameters$N                    <- N
parameters$inflation            <- inflation 
parameters$seed                 <- seed
parameters$perform_kf_ekf       <- perform_kf_ekf
parameters$KF_SaveClim          <- KF_SaveClim
parameters$perform_HBF          <- perform_HBF
parameters$HBF_SelfClim         <- HBF_SelfClim
parameters$w_cvr                <- w_cvr  
parameters$w_evp10              <- w_evp10
parameters$F_Lorenz             <- F_Lorenz
parameters$J_Lorenz             <- J_Lorenz
parameters$model_type           <- model_type

# Store the Parameters in the output "filter" variable

filter$parameters = parameters


set.seed(seed_for_secondary_fields)

#------------------------------------------------------
# DSADM: generate SECONDARY fields.
# Technique: generate TERTIARY fields and (point-wise) transform them to the SECONDARY fields.
# NB: Only ONE realization of each scnd field is generated: N_scnd=1

if(model_type == "DSADM"){
  
  message("Generate SECONDARY fields")
  
  N_scnd=1
  
  # Generate tert1(t,s)
  
  if(tert1_extpar$SD != 0){
    start_field = dsadm_start_field_cvm(tert1_modpar, n, Rem)$x_start
    
    tert_field  = dsadm_generator(as.matrix(start_field, nrow=n, ncol=N_scnd), 
                      n, N_scnd, ntime_model, dt, 
                      matrix(tert1_modpar$U,     nrow = n, ncol = time_model),
                      matrix(tert1_modpar$rho,   nrow = n, ncol = time_model),
                      matrix(tert1_modpar$nu,    nrow = n, ncol = time_model),
                      matrix(tert1_modpar$sigma, nrow = n, ncol = time_model), Rem)
 
  # Transform tert1(t,s) --> UU(t,s)
  
    UU = x_unpert_modpar$U + tert_field
  }else{
    UU = matrix(tert1_modpar$U, nrow = n, ncol = time_model)
  }
  
  #  Generate tert2(t,s)
  
  if(tert2_extpar$SD != 0){
    start_field = dsadm_start_field_cvm(tert2_modpar, n, Rem)$x_start
    
    tert_field  = dsadm_generator(start_field, n, N_scnd, ntime_model, dt, 
                      matrix(tert2_modpar$U,     nrow = n, ncol = time_model),
                      matrix(tert2_modpar$rho,   nrow = n, ncol = time_model),
                      matrix(tert2_modpar$nu,    nrow = n, ncol = time_model),
                      matrix(tert2_modpar$sigma, nrow = n, ncol = time_model), Rem)
    
  # Transform tert2(t,s) --> rho(t,s)  
    
    rrho = x_unpert_modpar$rho * transform_function(tert_field, epsilon_rho)
  }else{
    rrho = matrix(x_unpert_modpar$rho, nrow = n, ncol = time_model)
  }

  #  Generate tert3(t,s)
  
  if(tert3_extpar$SD != 0){
    start_field = dsadm_start_field_cvm(tert3_modpar, n, Rem)$x_start
    
    tert_field  = dsadm_generator(start_field, n, N_scnd, ntime_model, dt, 
                    matrix(tert3_modpar$U,     nrow = n, ncol = time_model),
                    matrix(tert3_modpar$rho,   nrow = n, ncol = time_model),
                    matrix(tert3_modpar$nu,    nrow = n, ncol = time_model),
                    matrix(tert3_modpar$sigma, nrow = n, ncol = time_model), Rem)
  # Transform tert3(t,s) --> nu(t,s)
    
    nnu = x_unpert_modpar$nu * transform_function(tert_field, epsilon_nu)
  }else{
    nnu = matrix(x_unpert_modpar$nu, nrow = n, ncol = time_model)
  }
  
  # Generate tert4(t,s) 
  
  if(tert4_extpar$SD != 0){
    start_field = dsadm_start_field_cvm(tert4_modpar, n, Rem)$x_start
    
    tert_field  = dsadm_generator(start_field, n, N_scnd, ntime_model, dt, 
                    matrix(tert4_modpar$U,     nrow = n, ncol = time_model),
                    matrix(tert4_modpar$rho,   nrow = n, ncol = time_model),
                    matrix(tert4_modpar$nu,    nrow = n, ncol = time_model),
                    matrix(tert4_modpar$sigma, nrow = n, ncol = time_model), Rem)
  # Transform tert4(t,s) --> sigma(t,s)
    
    ssigma = x_unpert_modpar$sigma * transform_function(tert_field)
  }else{
    ssigma = matrix(x_unpert_modpar$sigma, nrow = n, ncol = time_model)
  }
}

#------------------------------------------------------
# Setup params, ini field, and noise for Lorenz05

if(model_type == "Lorenz05" | model_type == "Lorenz05lin"){
  
  dt_atm_h = min(3, dt_h)  #h
  if(F_Lorenz > 32 & n > 60) dt_atm_h = dt_atm_h /2 # for stability, may appear to be reduced
  
  dt_Lorenz = dt_atm_h/6*0.05 # unitless, "Lorenz time" (6h atmospheric time ~ 0.05 Lorenz time units)

  seed_ini_cond_Lorenz=seed
  seed_noise_Lorenz   =seed *12.3456
  
  # Specify ini cond that *minimize* the initial transient
  
  assumed_mean = 1.2 * F_Lorenz^(1/3) # from Lorenz-05, p.1577, bef Eq(5)
  assumed_ms = assumed_mean * F_Lorenz # from Lorenz-05, p.1577, bef Eq(5)
  assumed_rms = sqrt(assumed_ms)
  assumed_var=assumed_ms - assumed_mean^2
  assumed_sd = sqrt(assumed_var)
  
  set.seed(seed_ini_cond_Lorenz)
  X1=rnorm(n, mean=assumed_mean, sd=assumed_sd) # ini condition
  
  # smooth ini cond to reduce the Lorenz-05's initial transient period
  
  nsweep=ceiling( 3*(60/n) * J_Lorenz)
  xx2=X1
  
  if(nsweep >0){
    
    for(sweep in 1:nsweep){
      xx1=xx2
      
      for (i in 1:n){
        im1=i-1
        if(im1 < 1) im1=im1+n
        
        ip1=i+1
        if(ip1 > n) ip1=ip1-n
        
        xx2[i]=(xx1[im1] + 2* xx1[i] + xx1[ip1]) /4
      }
    }
    
    X1=xx2 * max(nsweep/3, 1) # *nsweep/3 -- because smoothing reduces the magnitude
  }
  
  if(model_type == "Lorenz05lin"){
    X1_lin=X1 /5  # ini field to start the fcst t=1:2 in the filter;  5 is smth >1
  }
  
  # Specify system noise (model error) space-time field
  #  for Lorenz 05 & Lorenz05lin
  # Note that sd_noise_Lorenz is specified per 6h atm time, 2pi/60 mesh size
  # while the noise is white both in space and in time.
  # The whiteness implies that 
  # 
  # Var x_dscr(t,s) = 1/(delta_t * delta_s) ~ n/delta_t
  # 
  # For the specified  delta_t  and  n  grid points, we obtain
  # 
  # Var x_dscr = Var x_dscr_6h_60points * (n/60) / (delta_t / delta_t_6h)  ==>
  # sd_noise = sd_noise_Lorenz * sqrt( (n/60) / (delta_t / delta_t_6h) )

  dt_h_6h = 6 # h
  sd_noise = sd_noise_Lorenz * sqrt( (n/60) / (dt_h / dt_h_6h) )
  
  set.seed(seed_noise_Lorenz)
  noise_Lorenz=rnorm(time_model*n, mean=0, sd=sd_noise)
  
  noise_Lorenz=matrix(noise_Lorenz, nrow=n)  # model error (system noise)
  noise_Lorenz_zero=noise_Lorenz
  noise_Lorenz_zero[,]=0   # for the reference trajectory
}

#------------------------------------------------------
# Specify OBS network in space & time

ind_obs_space = seq(1, n, m)             # spatial obs network
ind_time_anls = seq(1, time_model, stride) # model time steps when anls is to be done

#------------------------------------------------------
# generate one-world TRUTH

N_truth=1  # number of realizations when generating the truth

# 1) Ini conditions
#    X_mdl_start  is the state vector at t=1

if(model_type == "DSADM"){
  
  Start = dsadm_start_field_cvm(x_unpert_modpar, n, Rem)
  X_mdl_start = Start$x_start
  Bx_start = Start$CVM
  
}else if(model_type == "Lorenz05" | model_type == "Lorenz05lin"){
  X_mdl_start=X1
}

# 2) run the model

if(model_type == "DSADM"){
  message("Main model run")
  X_true_mdl_steps = dsadm_generator(X_mdl_start, n, N_truth, ntime_model, dt, 
                                     UU, rrho, nnu, ssigma, Rem)
  
  if(predict_BBx){
    message("predict BBx")
    BBx = predictSpatialFldCVM(ntime_filter, n, dt, stride, Rem,
                               UU, rrho, nnu, ssigma,
                               Bx_start)
    filter$BBx   =  BBx  # anls steps
  }

  
}else if(model_type == "Lorenz05"){ # NLIN mdl
  X_true_mdl_steps = lorenz05(X_mdl_start, n, time_model, dt_Lorenz, 
                              F_Lorenz, J_Lorenz, noise_Lorenz)  # perturbed

}else if(model_type == "Lorenz05lin"){ # LINEARIZED mdl

  # (1) Reference trajectory for lorenz05lin  (non-perturbed !) and
  #     lin-lorenz05 TRUE PERTURBATION trajectory (UU)
  # Let the initial prtbn be =0
  
  n_pert=1 # one truth
  U_start=matrix(0, nrow=n, ncol=n_pert) # ncol=n_pert is one perturbation
  noise_spatim=array(noise_Lorenz, dim=c(n,n_pert, time_model))
  
  LIN = lorenz05lin(X_mdl_start, n, U_start, n_pert, time_model, dt_Lorenz, F_Lorenz, J_Lorenz, noise_spatim)
  UU              = LIN$UUU[,n_pert,]
  X_ref_mdl_steps = LIN$XX_ref
  
  # (3) The lin-lorenz05 TRUTH (X_ref_mdl_steps + UU)
  
  X_true_mdl_steps = X_ref_mdl_steps + UU
}

# TRUTH evaluated at the anls times

X_true = X_true_mdl_steps[,ind_time_anls] # truth at the anls times
TRUTH_RMSE = rmse(X_true[,], 0)

ind_time_flt_plot = seq(from=1, to=min(400, time_model), by=stride)

image2D(X_true_mdl_steps[,ind_time_flt_plot], main="X_true", xlab="Space", ylab="Time")
image2D(rrho  [,ind_time_flt_plot], main="rho  ", xlab="Space", ylab="Time")
image2D(nnu   [,ind_time_flt_plot], main="nu   ", xlab="Space", ylab="Time")
image2D(ssigma[,ind_time_flt_plot], main="sigma", xlab="Space", ylab="Time")
image2D(X_true_mdl_steps[,ind_time_flt_plot], main="X_true", xlab="Space", ylab="Time")

if(time_model > 10000){
  ind_time_flt_plot = seq(from=1, to=time_model, by=stride*100)
  
  image2D(X_true_mdl_steps[,ind_time_flt_plot], main="X_true long", xlab="Space", ylab="Time")
  image2D(rrho  [,ind_time_flt_plot], main="rho long ", xlab="Space", ylab="Time")
  image2D(nnu   [,ind_time_flt_plot], main="nu long ", xlab="Space", ylab="Time")
  image2D(ssigma[,ind_time_flt_plot], main="sigma long", xlab="Space", ylab="Time")
  image2D(X_true_mdl_steps[,ind_time_flt_plot], main="X_true", xlab="Space", ylab="Time")
  
}

max(abs(X_true_mdl_steps))
mean(abs(X_true_mdl_steps))

#-----------------------------------------------
# Store the Fields in the output "filter" variable
# All fields are written at the anls times only.

if(model_type == "DSADM"){
  
  filter$rrho   = rrho  [,ind_time_anls]
  filter$nnu    = nnu   [,ind_time_anls]
  filter$UU     = UU    [,ind_time_anls]
  filter$ssigma = ssigma[,ind_time_anls]
}

filter$X_true_anlstimes = X_true

#------------------------------------------------------
# Generate OBS

n_obs=n%/%m+min(1,n%%m)
OBS_NOISE = matrix(rnorm(n_obs*time_model, mean=0, sd=sqrt_R),
                   nrow=n_obs,
                   ncol=time_model)  # at the moment, OBS are generated at all mdl time steps (to be changed)
OBS = X_true_mdl_steps[ind_obs_space,] + OBS_NOISE

#------------------------------------------------------
# One-world KF/EKF



if(perform_kf_ekf > 0){
  
  message("Run KF")
  
  X_flt_start = X_mdl_start         # the 1st FILTER fcst starts from the truth, X_flt_start
  A_start = diag(n);   A_start[,]=0 # correspondingly, the time=0 anls is error-free
  
  KF_res = KF_EKF(ntime_filter, n, dt, stride, ind_obs_space, ind_time_anls, Rem,
                      UU, rrho, nnu, ssigma,  
                      F_Lorenz, J_Lorenz, sd_noise,
                      R_diag, m, OBS, 
                      X_flt_start, A_start, 
                      model_type, filter_type,
                      predict_BB_KF)
  
  image2D(KF_res$XXf[,1:min(200, time_filter)], main="XXf")
  image2D(KF_res$XXa[,1:min(200, time_filter)], main="XXa")
  image2D(X_true        [,1:min(200, time_filter)], main="X_true")
  
  TRUTH_RMSE = rmse(X_true[,], 0)
  KF_fRMSE  = rmse(KF_res$XXf[,], X_true[,])
  KF_fRMSE
  KF_aRMSE  = rmse(KF_res$XXa[,], X_true[,]) 
  KF_aRMSE
  
  B_var_mean = sum(diag(KF_res$B_mean)) /n
  sqrt_B_var_mean=sqrt(B_var_mean)
  sqrt_B_var_mean
  
  A_var_mean = sum(diag(KF_res$A_mean)) /n
  sqrt_A_var_mean=sqrt(A_var_mean)
  sqrt_A_var_mean

  # Test (permanent) KF_EKF:
  
  print("KF_fRMSE: real fcst-err RMSE")
  KF_fRMSE
  print("sqrt_B_var_mean: fcst-err RMSE assumed by the filter")
  sqrt_B_var_mean
  print("-- sh.be close to each other")
  
  print("KF_aRMSE: real anls-err RMSE")
  KF_aRMSE
  print("sqrt__var_mean: anls-err RMSE assumed by the filter")
  sqrt_A_var_mean
  print("-- sh.be close to each other")
  
  # Save B_clim
  
  if(KF_SaveClim > 0){
    B_clim_long = KF_res$B_mean
    B_clim_and_params = list(B_clim_long=B_clim_long,
                             parameters=parameters)
    filename=paste0("B_clim_", seed, "_tenKappa", kappa_rho*10, 
                    "_tenLpert_mult", L_perturb_mult*10, ".RData")
    save(B_clim_and_params, file=filename)
  }
  
  filter$KF = KF_res
  
}

#-----------------------------------------------
# One-world HBF (Generalized HBEF, which includes EnKF, HBEF, Var, and EnVar as special cases)


if(perform_HBF == 1){

  message("Run HBF")
  
  set.seed(seed_for_filters)
  
  X_flt_start = X_mdl_start # the 1st flt fcst starts from the truth
  
  # select B_clim
  
  if(HBF_SelfClim > 0){
    
    B_clim = KF_res$B_mean  # KF'c CVM from the current run
    
  }else{
    
    ClimFile=paste0("B_clim_", seed, "_tenKappa", kappa_rho*10, 
                    "_tenLpert_mult", L_perturb_mult*10, ".RData")
    load(ClimFile, verbose= TRUE)
    
    B_clim = B_clim_and_params$B_clim_long
  }
  
  # Lclz
  
  L_loc  = lclz_mult*L                # lclz radius, m
  C_lclz = lclz_matrix(n, L_loc, Rem) # lclz mx
  
  # Initial conditions
  
  Ba_start = matrix(0, nrow=n, ncol=n) # start flt at t=0 from the truth, as in the KF
  A_start = Ba_start # assume no obs at t=0
  
  # Compute the step-0 anls ensm members as sqrt(A)*g, whr g~N(0,I)
  
  sqrt_A_start = symm_pd_mx_sqrt(A_start)$sq
  N01 = matrix(rnorm(n^2, mean=0, sd=1), nrow=n, ncol=N)
  dXae_start = sqrt_A_start %*% N01    # perturbations yet
  Xae_start = X_flt_start + dXae_start # ensm members
  
  HBF_res = HBF(ntime_filter, n, dt, stride, ind_obs_space, ind_time_anls, Rem, 
                UU, rrho, nnu, ssigma,  
                F_Lorenz, J_Lorenz, sd_noise,
                R_diag, m, OBS, 
                X_flt_start, Ba_start, Xae_start, B_clim,
                N, w_cvr, w_evp10, 
                C_lclz, inflation,
                model_type)

  HBF_one_world = HBF_res
  
  filter$HBF_one_world = HBF_one_world
  
  HBF_fRMSE  = rmse(HBF_res$XXf[,], X_true[,])
  HBF_fRMSE
  
  HBF_RMSE_assumed = sqrt( mean( diag( HBF_res$B_mean ) ) )
  HBF_RMSE_assumed
  
  HBF_aRMSE  = rmse(HBF_res$XXa[,], X_true[,]) 
  HBF_aRMSE

}

#------------------------------------------------------
# Worlds: generate TRUTH & perform filtering

if(worlds){
  
  message("Run worlds")
  
  X     = list()
  X_HBF = list()
  
  set.seed(seed_for_filters)
  
  for(iter in 1:M){  
    cat("\r",paste0(round(iter/M*100,0),'%'))
    
    # Gen truth (one world -- one truth)
    
    X_start    = dsadm_start_field_cvm(x_unpert_modpar, n, Rem)$x_start
    
    X_true_mdl_steps  = dsadm_generator(X_start, n, N_truth, ntime_model, dt, 
                                 UU, rrho, nnu, ssigma, Rem)  # model steps
    X_true = X_true_mdl_steps[,ind_time_anls] # truth at the anls times
    X[[iter]] = X_true  # only anls steps

    #save(X_true, file = paste0(path,'/DATA/truth_',iter, '.Rdata'))
    
    #-----------
    # HBF filtering
    
    if(perform_HBF == 1){
      
      # gen OBS
      OBS_NOISE = matrix(rnorm((n%/%m+min(1,n%%m))*time_model, mean=0, sd=sqrt_R),
                         nrow=n%/%m+min(1,n%%m), ncol=time_model)
      OBS       = X_true_mdl_steps[ind_obs_space,] + OBS_NOISE
      
      # Initial conditions
      
      Ba_start = matrix(0, nrow=n, ncol=n) # start flt at t=0 from the truth, as in the KF
      A_start = Ba_start # assume no obs at t=0
      
      # Compute the step-0 anls ensm members as sqrt(A)*N(0,I)
      
      sqrt_A_start = symm_pd_mx_sqrt(A_start)$sq
      N01 = matrix(rnorm(n^2, mean=0, sd=1), nrow=n, ncol=N)
      dXae_start = sqrt_A_start %*% N01    # perturbations yet
      
      X_flt_start = X[[iter]][,1]
      Xae_start   = X_flt_start + dXae_start # ensm members
      
      
      HBF_res = HBF(ntime_filter, n, dt, stride, ind_obs_space, ind_time_anls, Rem, 
                    UU, rrho, nnu, ssigma,  
                    F_Lorenz, J_Lorenz, sd_noise,
                    R_diag, m, OBS, 
                    X_flt_start, Ba_start, Xae_start, B_clim,
                    N, w_cvr, w_evp10, 
                    C_lclz, inflation,
                    model_type)
      
      #save(HBF_res, file = paste0(path,'/DATA/HBF_',iter, '.Rdata'))
      
      X_HBF[[iter]] = HBF_res$XXf  # only anls steps
      
    } # end if_HBF
  } # end loop_iter
  
} # end if_worlds

#-----------------------------------------------
# Estm B_x_true by averaging over the worlds

if(worldsAve_BBx){
  
  message('Covariance matrix of x: averaging over the worlds')
  
  X_arr <- array(NA, dim = c(M, n, time_filter))
  
  for(iter in 1:M){
    cat("\r",paste0(round(iter/M*100,0),'%'))
    X_arr[iter,,] <- X[[iter]]
  }
  
  BBx_worlds_ave = array(NA, dim = c(n, n, time_filter))

  for(t in 1:time_filter){
    BBx_worlds_ave[,,t] <- cov(X_arr[,,t])
  }
  
  image2D(BBx_worlds_ave[,,time_filter/2], main="CVM of x at t=time_filter/2")
  image2D(BBx_worlds_ave[,,time_filter],   main="CVM of x at t=time_filter")

  CVM_x_mean = apply(BBx_worlds_ave, c(1,2), mean)
  image2D(CVM_x_mean, main="Time mean spatial CVM of x")
  
  var_x_mean = mean(diag(CVM_x_mean))
  sqrt(var_x_mean)
  sd_x_true = sd(as.vector(X_true))  # one world only..
  sd_x_true  # sh.be close to sqrt(var_x_mean). OK.
  
  filter$BBx_worlds_ave = BBx_worlds_ave
}

#-----------------------------------------------
# Estm BB_HBF by averaging over the worlds

if(worldsAve_BB_HBF){
  
  message('HBF: background-error (sample covariance) matrix S and one-world flt stats')

  BB_HBF_worldsAve = array(0, dim = c(n, n, time_filter))
  
  for(iter in 1:M){
    #load(paste0(path,'/DATA/truth_',iter,'.Rdata'))  
    #load(paste0(path,'/DATA/HBF_',iter,'.Rdata'))  
    
    X_true = X[[iter]] 
    Xf_HBF = X_HBF[[iter]]
    
    for(i in 1:time_filter){
      BB_HBF_worldsAve[,,i] = BB_HBF_worldsAve[,,i] + 
        ((Xf_HBF[,i] - X_true[,i])) %*% t(Xf_HBF[,i] - X_true[,i])
    }
    cat("\r",paste0(round(iter/M*100,0),'%'))
  }
  
  BB_HBF_worldsAve = BB_HBF_worldsAve / M
  
  filter$BB_HBF_worldsAve = BB_HBF_worldsAve
} 

#-----------------------------------------------
# write output 

save=TRUE # FALSE TRUE

if(save){
  if(mode == 0) {                            # single-run model/filters output
    save(filter, file = paste0(path,"/fields_", seed, "_mode",mode, ".RData"))
  }
  
  if(mode > 0) {                             # Fields and Covs output
    save(filter, file = paste0(path,"/fields_covs_", seed,   "_mode",mode, ".RData"))
  }
}

#--------------------------------------------
# Write output txt file

outfile=paste0("./RMSE.txt")

if(KF_SaveClim > 0) outfile=paste0("./RMSE_long_", seed, ".txt")

unlink(outfile)
sink(outfile, append=TRUE)

cat("\n")

if(perform_kf_ekf > 0){
  
  cat("TRUTH_RMSE=")
  print(TRUTH_RMSE)
  
  cat("\n")
  
  cat("KF_fRMSE: real fcst-err RMSE")
  print(KF_fRMSE)
  cat("sqrt_B_var_mean: fcst-err RMSE assumed by the filter")
  print(sqrt_B_var_mean)
  cat("-- sh.be close to each other")
  
  cat("\n")
  
  cat("KF_aRMSE: real anls-err RMSE")
  print(KF_aRMSE)
  cat("sqrt__var_mean: anls-err RMSE assumed by the filter")
  print(sqrt_A_var_mean)
  cat("-- sh.be close to each other")
  
  cat("\n")
  
  cat("FG VARIANCE reduction in the anls")
  print(100*(KF_fRMSE^2 - KF_aRMSE^2) / KF_fRMSE^2)

}


cat("\n")

if(perform_HBF > 0){
  
  cat("TRUTH_RMSE=")
  print(TRUTH_RMSE)
  
  cat("\n")
  
  cat("HBF_fRMSE=")
  print(HBF_fRMSE)
  
  cat("\n")
  
  cat("HBF_RMSE_assumed=")
  print(HBF_RMSE_assumed)
  
  cat("\n")
  
  cat("HBF_aRMSE=")
  print(HBF_aRMSE)

  cat("\n")
  
  cat("(HBF_fRMSE - KF_fRMSE) / KF_fRMSE, %:")
  print(100*(HBF_fRMSE - KF_fRMSE) / KF_fRMSE)
  
  
}

sink()
#------------------------------------------------
