# Examine the output of the DSM on S1
# making use of the TRUE covariance matrix computed by:
# either (1) predicting the field and KF covariances
# or     (2) averaging over a large number of "worlds".
# 
# Which field (X_true or KF's background error field or HBF's background error field)
# is examined is decided according the the filter$parameters$mode switch:
# 
# mode=0: single mdl&flt runs (produce Fields only), 
# mode=1: predicted BBx (fieds' Covariances), 
# mode=2: predicted BB_KF (KF's background-error Covariances),
# mode=3: worlds-averaged Fields' Covariances 
# mode=4: worlds-averaged HBF's background-error Covariances 
#
# (there is no config file for this script)
#
# (1) Plot:
#    xi,
#    the secondary fields U, rho, nu, sigma, 
#    V, 
#    Lambda
#    lambda.
# 2) Calc: 
#    Variability (max/min and the ratio 90% / 10% quantile) in
#      V, Lambda, microscale
#    crl(theta, V), where theta=U,rho,nu, sigma
#    crl(V, Lambda)
#    ...
#
# M Tsyrulnikov
# 17 Jul 2018


library(MASS)
library(stats)
library(plot3D)

source('est_crf_timser.R')
source('symm_cvm_row.R')


universe= "All" # Rho  Nu  U  Sigma All

datafile=paste0("./fields_covs.RData")


load(datafile) # contains:
#           parameters:
#             dim=n, delta_t=delta_t_model, 
#             stride defines the asml time step: Ta=delta_t * stride
#             ...
#           all scnd fields: 
#             Rho, Nu, U, Sigma
#           Cov_mat=CVM_truth[1:n,1:n, 1:ntime] 


n=filter$parameters$n
ntime_filter=filter$parameters$time_filter
stride=filter$parameters$stride
ntime=ntime_filter   # only anls-time instants are analyzed

M=filter$parameters$M
seed=filter$parameters$seed
mode=filter$parameters$mode

message("mode=", mode)

# mode:
# mode=0: single mdl&flt runs (Fields only available), 
# mode=1: predicted BBx (fieds' Covariances), 
# mode=2: predicted BB_KF (KF's background-error Covariances),
# mode=3: worlds-averaged Fields' Covariances 
# mode=4: worlds-averaged HBF's   Covariances 

# Primary (xi) and secondary (U, rho, nu, sigma) fields -- at ALL mdl time steps

# Field and COVs selector:
# mode=0: no covs to analyse; stop execution.
# mode=1: xi=the field x_true itself, B=B_x_true  Predicted.
# mode=3: xi=the field x_true itself, B=B_KF      Predicted.
# mode=2: xi=KF's  background error field, B=B_x_true worlds-averaged
# mode=4: xi=HBF's background error field, B=B_HBF    worlds-averaged


if(mode == 0) {
  
  field="x_true"
  stop("mode=0. No covs to examine. Stop.")
  
}else if(mode == 1){
 
  field="x_true"
  xi=filter$X_true_anlstimes # the true field
  B=filter$BBx               # predicted x_true covs
 
}else if(mode == 2){
  
  field="KF"
  xi=filter$KF$XXf - filter$X_true_anlstimes # the FK's background-error field
  B=filter$KF$BB             # predicted KF background-error covs
  
}else if(mode == 3){
  
  field="x_true"
  xi=filter$X_true_anlstimes # the true field
  B=filter$BBx_worlds_ave    # estimated by world-averaging x_true covs
  
}else if(mode == 4){
  
  field="HBF"
  xi=filter$HBF_one_world$XXf - filter$X_true_anlstimes # the FK's background-error field
  B=filter$BB_HBF_worldsAve
}

# scnd fields

rho=filter$rrho
nu=filter$nnu
sigma=filter$ssigma
U=filter$UU


# Find out whether x_true is stationary

stationarity=FALSE
if(sd(as.vector(rho))<1e-16 & sd(as.vector(nu))<1e-16 & sd(as.vector(U))<1e-16 & sd(as.vector(sigma))<1e-16){
  stationarity=TRUE  
}


if(universe == "Rho"){
  theta_name="rho" 
  theta=rho
  
} else if (universe == "Nu"){
  theta_name="nu"
  theta=nu
  
} else if (universe == "U"){
  theta_name="U"
  theta=U
  
} else if (universe == "Sigma"){
  theta_name="sigma"
  theta=sigma
  
} else if (universe == "All"){
  theta_name="all"
  theta=rho # !!!!!!! just to fill in smth
} 


Rem=6.37*10^6 # m
rekm=Rem/10^3 # km

mesh    = 2*pi/n  # rad
mesh_km = mesh * rekm # grid spacing, km
grid_km = c(0:(n-1))*mesh_km
step_h=filter$parameters$dt / 3600
tgrid_h=c(1:ntime)*step_h

V = matrix(0, nrow=n, ncol=ntime)

for (i in (1:n)){
  V[i,]=B[i,i,]
}
std_st = sqrt(V)

macroscale_st_km  <- matrix(0, nrow=n, ncol=ntime)
microscale_st_km  <- matrix(0, nrow=n, ncol=ntime)

CRM=B  # initialize

for (t in (1:ntime)){
  CRM[,,t]=diag(1/std_st[,t]) %*% B[,,t] %*% diag(1/std_st[,t])  # Crl Mx
  for (ii in (1:n)){
    macroscale_st_km[ii, t]  = sum(CRM[ii,,t])*mesh_km/2
    iim1=ii-1
    if(ii==1) iim1=n
    iip1=ii+1
    if(ii==n) iip1=1
    microscale_st_km[ii, t] = 1/sqrt( (-CRM[ii,iim1,t] + 2*CRM[ii,ii,t] - CRM[ii,iip1,t])/mesh^2 ) * rekm
  }
}  

Lambda=macroscale_st_km
lambda=microscale_st_km

# Spatially averaged fields

V_spaAve = apply(V, 2, mean)
Lambda_spaAve = apply(Lambda, 2, mean)

#------------------------------------------------
# Plots: theta, xi, V, Lambda, lambda

ntime_plot=400
nt=floor(ntime/ntime_plot)
nt=min(3,nt)

for(it in 1:nt){
  
  t1=(it-1)*ntime_plot +1
  t2=t1 + ntime_plot -1
  tt=t1:t2
  
  if(theta_name == "all"){
    
    namefile=paste0(field, "_rho_",t1,".pdf")
    #png(namefile, width=5.1, height=5.1, units = "in", res=300)
    pdf(namefile, width=7, height=7)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(rho[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
            main=bquote("Secondary field"~rho), 
            cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
    dev.off()
    
    namefile=paste0(field, "_nu_",t1,".pdf")
    #png(namefile, width=5.1, height=5.1, units = "in", res=300)
    pdf(namefile, width=7, height=7)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(nu[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
            main=bquote("Secondary field"~nu),
            cex.main=1.7, cex.axis=1.5, cex.lab=1.7,) 
    dev.off()
    
    namefile=paste0(field, "_sigma_",t1,".pdf")
    #png(namefile, width=5.1, height=5.1, units = "in", res=300)
    pdf(namefile, width=7, height=7)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(sigma[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
            main=bquote("Secondary field"~sigma),
            cex.main=1.7, cex.axis=1.5, cex.lab=1.7,) 
    dev.off()
    
    namefile=paste0(field, "_U_",t1,".pdf")
    #png(namefile, width=5.1, height=5.1, units = "in", res=300)
    pdf(namefile, width=7, height=7)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(U[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
            main=bquote("Secondary field"~U),
            cex.main=1.7, cex.axis=1.5, cex.lab=1.7,) 
    dev.off()
    
  }else{
    
    #namefile=paste0(field, theta_name,"_",t1,".png")
    #png(namefile, width=5.1, height=5.1, units = "in", res=300)
    #par(mai=c(1.2,1.2,0.7,0.7))
    #image2D(theta[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
    #        main=paste(theta_name), 
    #        cex.main=1.7, cex.axis=1.5, cex.lab=1.7,) 
    #dev.off()
    
  }
  
  if(stationarity){
    fieldtitle="Stationary field "
  }else{
    fieldtitle="Non-stationary field "
  }
  
  namefile=paste0(field, "_xi_",t1,".pdf")
  #png(namefile, width=5.1, height=5.1, units = "in", res=300)
  pdf(namefile, width=7, height=7)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(xi[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
          main=fieldtitle,
          #main="Non-stationary field",
          cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
  dev.off()
  
  
  namefile=paste0(field, "_V_",t1,".pdf")
  #png(namefile, width=5.1, height=5.1, units = "in", res=300)
  pdf(namefile, width=7, height=7)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(log10(V[,tt]), x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
          main=bquote(log[10]~"(Var"~xi~")"), 
          cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
  dev.off()
  
  
  namefile=paste0(field, "_Lambda_",t1,".pdf")
  #png(namefile, width=5.1, height=5.1, units = "in", res=300)
  pdf(namefile, width=7, height=7)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(Lambda[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
          main=as.expression( bquote(Lambda) ), 
          cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
  dev.off()
  
  # spa Ave plots
  
  namefile=paste0(field, "_V_SpaAve",t1,".pdf")
  #png(namefile, width=5.1, height=5.1, units = "in", res=300)
  pdf(namefile, width=7, height=7)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(log10(V_spaAve[tt]),  xlab="Time, h",
          main=bquote(log[10]~"(Var_spaAve"~xi~")"), 
          cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
  dev.off()
  
  
  namefile=paste0(field, "_Lambda_SpaAve",t1,".pdf")
  #png(namefile, width=5.1, height=5.1, units = "in", res=300)
  pdf(namefile, width=7, height=7)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(Lambda_spaAve[tt], xlab="Time, h",
          main=bquote("(Lambda_spaAve)"), 
          cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
  dev.off()
  
  
  
  
  #namefile=paste0("Lambda_mic_",t1,".png")
  #png(namefile, width=5.1, height=5.1, units = "in", res=300)
  #par(mai=c(1.2,1.2,0.7,0.7))
  #image2D(lambda[,tt], x=grid_km, y=tgrid_h[tt], xlab="Space, km", ylab="Time, h",
  #        main=paste("lambda"), 
  #        cex.main=1.7, cex.axis=1.5, cex.lab=1.7,)
  #dev.off()
  
}

#------------------------------------------------
# Selected i,t spatial crl stats

SpaCrl_plots=TRUE # TRUE FALSE
if(SpaCrl_plots){
  numplot=5
  for (iplot in c(1:numplot)) {
    
    ncurve=30
    
    tt <- sample(c(1:ntime), ncurve, replace = FALSE, prob = NULL) # select random time instant
    ii <- sample(c(1:n),     ncurve, replace = FALSE, prob = NULL) # select random spatial point
    
    j_center=floor(n/2) # for plotting
    
    mn=0
    mx=0
    for (icurve in 1:ncurve){
      mn <- min( mn, min( CRM[ii[icurve], , tt[icurve]] ) )
      mx <- max( mx, max( CRM[ii[icurve], , tt[icurve]] ) )
    }
    
    t=tt[1] # 1st curve
    i=ii[1] # 1st curve
    
    row=CRM[i,,t]
    srow=symm_cvm_row(row, n, i)
    
    namefile=paste0(field, "_SpaCrl_t",iplot,".pdf")
    #png(namefile, width=7.1, height=5.1, units = "in", res=300)
    pdf(namefile, width=7, height=5)
    par(mai=c(1.2,1.2,0.7,0.7))
    plot(grid_km[]-grid_km[j_center], srow, 
         xlab="Distance, km", ylab="Correlation", 
         type="l", ylim=c(mn,mx), 
         cex.main=1.7, cex.axis=1.5, cex.lab=1.7,,
         main=paste("Spatial forecast-error correlation") )
    
    for(icurve in 2:ncurve){
      i=ii[icurve]
      t=tt[icurve]
      row=CRM[i,,t]
      srow=symm_cvm_row(row, n, i)
      lines(grid_km[]-grid_km[j_center], srow, type="l")
    }
    
    abline(h=0)
    dev.off()
    
  } # End for(iplot in c(1:numplot))  
}   # End if(SpaCrl_plots)


mx_mac <- max(macroscale_st_km[,t])
mx_std=max(std_st[,t])

plot (macroscale_st_km[,t]/mx_mac, ylim=c(0,1.1), type="l", col="green",
      main=("GREEN: Rel. Len.scale | RED: Rel. std"))
lines(std_st[,t]/mx_std, type="l", col="red")


plot (microscale_st_km[,t]/macroscale_st_km[,t], type="l", col="green",
      main="Rel. curv. radius at the origin")

#------------------------------------------------

ind_max=which(Lambda == max(Lambda), arr.ind=TRUE)
plot(B[,ind_max[1],ind_max[2]])

ind_min=which(Lambda == min(Lambda), arr.ind=TRUE)
plot(B[,ind_min[1],ind_min[2]])
                            
#------------------------------------------------
# Cross-field crls

if(!stationarity){
  crl_Lambda_V=cor(as.vector(Lambda), as.vector(V), method="spearman")
  
  crl_sigma_V=cor(as.vector(sigma), as.vector(V), method="spearman")
  crl_rho_V=cor(as.vector(rho), as.vector(V), method="spearman")
  crl_nu_V=cor(as.vector(nu), as.vector(V), method="spearman")
}

#-------------
# Assess L:
# LL=sqrt(nu/rho).
# To account for negative rho, let 
# LL=sqrt(nu/(rho - rho_min + abs(rho_min)/2),
# the same with nu.

rho_min=min(rho)
nu_min =min(nu)
rho_regul=rho - rho_min + abs(rho_min)/2
nu_regul =nu  - nu_min  + abs(nu_min) /2
LL=sqrt(nu_regul/rho_regul)

#-------------
# Assess V:
# We have the eqn for the statio model
# V = a^2 sigma^2/(2*rho) * sum_{m=-n/2}^{n/2} 1/[1 + (Lm/R)^2].
# We simplify it by replacing the sum with the integral:
# V ~ a^2 sigma^2/(2*rho) * int_{-Ln/(2R)}^{Ln/(2R)} dx /[1 + x^2] =
#   a^2 sigma^2/(2*rho) * 2*arctg(Ln/(2R)).
# With large enough n, we obtain
#
# V ~ a^2 sigma^2/(2*rho) * pi
#
# a=1/sqrt(2*pi*R)  ==>
# a^2*pi = 1/(2*R)
#.............................
# V ~ 1/(2*R) sigma^2/(2*rho) 
#.............................

VV=1/(2*Rem) * sigma^2 / (2 * rho_regul)

if(!stationarity){
  crl_LL_Lambda=cor(as.vector(LL), as.vector(Lambda), method="spearman")
  crl_ll_lambda=cor(as.vector(LL), as.vector(lambda), method="spearman")
  crl_VV_V=cor(as.vector(VV), as.vector(V), method="spearman")
  
  crl_Lambda_V=cor(as.vector(Lambda), as.vector(V), method="spearman")
  crl_lambda_V=cor(as.vector(lambda), as.vector(V), method="spearman")
}

micro_d_macro_scale_mean=mean(lambda / Lambda)



#------------------------------------------------
# time shifted cross-field crls
# V should be delayed wrt VV
# L should be delayed wrt LL

# (V)

SHIFT=FALSE
if(SHIFT){
  nsh=5
  
  ccor_VV_V = c(1:nsh)
  ccor_VV_V[]=0
  
  for (sh in (1:nsh)) {
    ntime_minus_sh=ntime - sh
    ccor_VV_V[sh] = cor(as.vector(V[,(sh+1):ntime]), as.vector(VV[,1:ntime_minus_sh]),
                        method="spearman") # spearman better than pearson,
  }  
  
  ccor_VV_V
  crl_VV_V_abs_max=max(abs(ccor_VV_V))
  shift_VV_V_opt=which(abs(ccor_VV_V) == abs(crl_VV_V_abs_max))
  crl_VV_V_max=ccor_VV_V[shift_VV_V_opt]
  shift_VV_V_opt
  crl_VV_V_max 
  
  # (L)
  
  ccor_LL_Lambda = c(1:nsh)
  ccor_LL_Lambda[]=0
  
  for (sh in (1:nsh)) {
    ntime_minus_sh=ntime - sh
    ccor_LL_Lambda[sh] = cor(as.vector(Lambda[,(sh+1):ntime]), as.vector(LL[,1:ntime_minus_sh]),
                             method="spearman") # spearman better than pearson,
  }  
  
  ccor_LL_Lambda
  crl_LL_Lambda_abs_max=max(abs(ccor_LL_Lambda))
  shift_LL_Lambda_opt=which(abs(ccor_LL_Lambda) == abs(crl_LL_Lambda_abs_max))
  crl_LL_Lambda_max=ccor_LL_Lambda[shift_LL_Lambda_opt]
  shift_LL_Lambda_opt
  crl_LL_Lambda_max
}

# --> time delay not significant.
#------------------------------------------------
# In the U-universe, check if V & Lambda_xi are correlated with 
# dU/dx (x is  an integer in the centered finite difference)
# dU/dt (t is  an integer in the centered finite difference)

dUdx = matrix(0, nrow=n, ncol=ntime)
dUdt = matrix(0, nrow=n, ncol=ntime)
             
if(universe == "U") {

  for (t in c(2:(ntime-1))) {
    for (i in c(1:n)) {
      im1=i-1
      ip1=i+1
      if(i == 1) im1=n
      if(i == n) ip1=1
      dUdx[i,t]=(U[ip1,t] - U[im1,t])/2
      dUdt[i,t]=(U[i,t+1] - U[i,t-1])/2
    }
  }
  
# smooth dUdt
  
  nsmoo=64 # >0 (16 is good, 64 is even better - for Lambda and V)
  dUdx_=dUdx # initialize

  for (smoo in c(1:nsmoo)){
    dUdx=dUdx_
  
    for (t in c(2:(ntime-1))) {
      for (i in c(1:n)) {
        im1=i-1
        ip1=i+1
        if(i == 1) im1=n
        if(i == n) ip1=1
        dUdx_[i,t]=( dUdx[i,t]+ 
                  0.5*(dUdx[ip1,t]+dUdx[im1,t]+dUdx[i,t-1]+dUdx[i,t+1]) )/3
      }
    }
  }
 
  image2D(dUdx_)
  crl_Lambda_dUdx_smoothed =
    cor( as.vector(Lambda          [,2:(ntime-1)]) , 
         as.vector(dUdx_[,2:(ntime-1)]), method="pearson" ) # spearman  pearson
  crl_lambda_dUdx_smoothed =
    cor( as.vector(microscale_st_km[,2:(ntime-1)]) , 
         as.vector(dUdx_[,2:(ntime-1)]), method="pearson" )
  crl_V_dUdx_smoothed =
    cor( as.vector(V               [,2:(ntime-1)]) , 
         as.vector(dUdx_[,2:(ntime-1)]), method="pearson" )
  
}

# ==> Both V & Lambda positively correlate with a smoothed dU/dx,
#     with nsmoo=64, the crls are both =0.6.
#------------------------------------------------
#------------------------------------------------
# Write the output file


outfile=paste0("out_",field,".txt")
unlink(outfile)
sink(outfile)


cat("ntime=")
print(ntime)

cat("M=")
print(M)

cat("max(V) / min(V)=")
print(max(V)/min(V))

cat("max(Lambda) / min(Lambda)=")
print(max(Lambda)/min(Lambda))

cat("max(lambda) / min(lambda)=")
print(max(lambda)/min(lambda))

V_decile_ratio=quantile(V[, 1:ntime], probs=0.9) / quantile(V[, 1:ntime], probs=0.1)
cat("V_decile_ratio=")
print(V_decile_ratio)

Lambda_decile_ratio=quantile(Lambda[, 1:ntime], probs=0.9) / quantile(Lambda[, 1:ntime], probs=0.1)
cat("Lambda_decile_ratio=")
print(Lambda_decile_ratio)

cat("\n")

if(!stationarity){
  cat("crl_Lambda_V=")
  print(crl_Lambda_V)
  
  cat("crl_rho_V=")
  print(crl_rho_V)
  
  cat("crl_nu_V=")
  print(crl_nu_V)
  
  cat("crl_sigma_V=")
  print(crl_sigma_V)
  
  cat("crl_LL_Lambda=")
  print(crl_LL_Lambda)
  
  cat("crl_ll_lambda=")
  print(crl_ll_lambda)
  
  cat("crl_VV_V=")
  print(crl_VV_V)
  
}

if(SHIFT){
  cat("Time shifted: crl_VV_V_max=")
  print(crl_VV_V_max)
  
  cat("Time shifted: crl_LL_Lambda_max=")
  print(crl_LL_Lambda_max)
}

cat("micro_d_macro_scale_mean=")
print(micro_d_macro_scale_mean)

cat("\n")


if(universe == "U") {
  cat("nsmoo=")
  print(nsmoo)
  
  cat("crl_Lambda_dUdx_smoothed=")
  print(crl_Lambda_dUdx_smoothed)
  
  cat("crl_lambda_dUdx_smoothed=")
  print(crl_lambda_dUdx_smoothed)
  
  cat("crl_V_dUdx_smoothed=")
  print(crl_V_dUdx_smoothed)
}



sink()
print("Normal finish")
#------------------------------------------------
