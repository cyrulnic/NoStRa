# Examine the output from one world of the DSADM on S1:
# 1) Calc the time-mean temporal and the spatial crfs of xi
# 3) Gau q-q plot for xi (unconditionally on the scnd flds)
# 
# All fields in the input "filter" variable are written at the anls times only.
#
# M Tsyrulnikov
# 8 Jul 2018

library(MASS)
library(stats)
library(plot3D)
library(nortest)
library(moments)

source('est_crf_timser.R')
source('est_crf_timser_cyc.R')
source('fft_Rspe_rearrange.R')

datafile=paste0("./fields.RData")


Rem=6.37*10^6 # m
rekm=Rem/10^3 # km

# Specify the field to be examined

exam = "x_true" # x_true  or FG errors: xi_KF or xi_HBF

#------------------------------------------------------------------------
# Read fields to be examined

load(datafile, verbose=TRUE) # contains:
#           parameters:
#             dim=n, dt_model, 
#             stride defines the asml time step: Ta=dt * stride
#             ...
#           all scnd fields: 
#             rrho, nnu, UU, ssigma
#           x_true
#           filter$HBF_one_world$XXf
#           filter$KF$XXf

n=filter$parameters$n
ntime=filter$parameters$time_filter

#ntime=min(ntime, 5000)  ##   TMP !!!

step_h=filter$parameters$dt_filter / 3600
#M=filter$parameters$M
seed=filter$parameters$seed

# Final piece of setup

mesh <- 2*pi/n  # rad
mesh_km = mesh * rekm # grid spacing, km
grid_km = c(0:(n-1))*mesh_km

tgrid_h=c(1:ntime)*step_h
ntime_plot=min(5000, ntime)
ntime4plots_steps=150 # floor(4000/step_h)
ntime4plots_steps=min(ntime,ntime4plots_steps)

#---------------------------------------------
# fill in scnd flds

rho    = filter$rrho
nu     = filter$nnu
sigma  = filter$ssigma
U      = filter$UU  

#---------------------------------------------
# fill in the field xi to be examined

if(exam == "x_true"){         # x_true 
  
  xi=                           filter$x_anls_times
    
}else if(exam == "xi_KF"){    # KF's  FG error field 
  
  xi=filter$KF$XXf -            filter$x_anls_times
 
}else if(exam == "xi_HBF"){   # HBF's FG error field 
  
  xi=filter$HBF_one_world$XXf - filter$x_anls_times
} 

#------------------------------------------------------------------------
# dummies

universe= "All" # Rho  Nu  U  Sigma All

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

#------------------------------------------------------------------------
# Diagnose the most general field stats

xi_spamean=apply(xi,2, mean)
xi_spasd  =apply(xi,2, sd)

plot(xi_spamean[1:min(ntime,2000)], main="spamean[t]")
plot(xi_spasd  [1:min(ntime,2000)], main="spasd[t]")

mean(xi[1:n,1:ntime])
sd(xi[1:n,1:ntime])
median(abs(xi[1:n,1:ntime]))

#------------------------------------------------------------------------
# Plots: xi and the scnd flds

plots_fld = TRUE #  TRUE FALSE

if(plots_fld){
  nt1=1
  if(ntime < nt1*ntime4plots_steps) nt1=1
  dt1=floor(ntime/nt1)
  
  for (it1 in 1:nt1){
    t1=(it1 -1)*dt1 +1 
    t2=t1 + ntime4plots_steps -1
    
    pngname=paste0(exam, "_",t1, ".png")
    png(pngname, width=5.1, height=5.1, units = "in", res=300)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(xi[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
            xlab="Space, km", ylab="Time, h",
            #main=as.expression( bquote(xi) ), 
            main=exam,
            cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
    dev.off()
    
    
    pngname=paste0("rho_",t1, ".png")
    png(pngname, width=5.1, height=5.1, units = "in", res=300)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(rho[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
            xlab="Space, km", ylab="Time, h",
            main=as.expression( bquote(rho) ), 
            cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
    dev.off()
    
    
    pngname=paste0("nu_",t1, ".png")
    png(pngname, width=5.1, height=5.1, units = "in", res=300)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(nu[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
            xlab="Space, km", ylab="Time, h",
            main=as.expression( bquote(nu) ), 
            cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
    dev.off()
    
    pngname=paste0("sigma_",t1, ".png")
    png(pngname, width=5.1, height=5.1, units = "in", res=300)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(sigma[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
            xlab="Space, km", ylab="Time, h",
            main=as.expression( bquote(sigma) ), 
            cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
    dev.off()
    
    pngname=paste0("U_",t1, ".png")
    png(pngname, width=5.1, height=5.1, units = "in", res=300)
    par(mai=c(1.2,1.2,0.7,0.7))
    image2D(U[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
            xlab="Space, km", ylab="Time, h",
            main=paste('U'), 
            cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
    dev.off()
  
  }
  
  # long plot
  
  lent=min(ntime, 1000)
  indt=seq(from=1, to=ntime, length.out =lent)
  
  pngname=paste0(exam, "_long.png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(xi[,indt], x=grid_km, y=tgrid_h[indt], 
          xlab="Space, km", ylab="Time, h",
          #main=as.expression( bquote(xi) ), 
          main=exam,
          cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  dev.off()
  
} # end if(plots_fld)

#------------------------------------------------------------------------
# NGaussianity of  xi

plots_nGau = TRUE # TRUE FALSE
if(plots_nGau){
  
  xx=as.vector(xi)
  nn=length(xx)
  rarefy=16 # >1 to save CPU time while uniformly sampling the whole spacetimea 64 is fast & OK.
  ind=seq(from=1, to=nn, by=rarefy)
  yy=xx[ind]
  
  pngname=paste0("QQplot_xi_perturb_",theta_name,".png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  qqnorm(yy, main = paste("Normal q-q plot: xi"),
         xlab="Theoretical quantile", ylab="Sample quantile", 
         cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  qqline(yy)
  dev.off()
  
  kurtosis_xi=kurtosis(yy)
  skewness_xi=skewness(xx)

} # end if(plots_nGau)

#------------------------------------------------
# Spatial crf estm (direct)
# (tested OK)
#  Turn off because the fft-based crf (see below) is cheaper to compute.

estm_spacrf=TRUE # FALSE
if(estm_spacrf){
  nt=ntime /1
  sscrf=matrix(0, nrow=nt, ncol=n)
  
  mxi=mean(xi)
  stdxi=sd(as.vector(xi))
  
  for (t in 1:nt){
    sscrf[t,]=est_crf_timser_cyc(xi[,t], mxi, stdxi)$crf
  }
  scrf=colMeans(sscrf)
  
  dsmax=n/2 # min(20, n/2)
  mn=min(scrf)
  mx=max(scrf)
  
  pngname=paste0("SpaCrl_mean_",exam,".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(grid_km[1:(dsmax+1)] - grid_km[1], scrf[1:(dsmax+1)], ylim=c(mn,mx), type="o",
       xlab="Distance, km", ylab="Correlation",
       main=paste("Mean spa crl (direct estm) ", exam))
  abline(h=0)
  dev.off()
}

#------------------------------------------------------------------------
# Spatial spe estm

nt=ntime /1
sspe=matrix(0, nrow=n, ncol=nt)
mxi=mean(xi)

for (t in 1:nt){
  sspe[,t] = abs( fft(xi[,t] - mxi, inverse=FALSE) ) ^2
}
spaspe_full=Re( rowMeans(sspe) )
#spa_spe_symm=fft_Rspe_rearrange(spaspe)
spaspe=c(1:((n/2)+1)) # init

spaspe[1]=spaspe_full[1]
for (m in (1:(n/2))){ # m is the wvn
  ifull1=m+1
  ifull2=n - m +1
  spaspe[m+1]=( spaspe_full[ifull1] + spaspe_full[ifull2] ) /2
}

spaspe_grid=(1:((n/2)+1)) -1

mn=min(spaspe)
mx=max(spaspe)

pngname=paste0("SpaSpe_mean_",exam,".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(spaspe_grid, spaspe, ylim=c(mn,mx), type="o",
     xlab="Global wavenumber", ylab="Spectrum",
     main=paste("Mean spa spectrum. ", exam))
dev.off()

mn=min(log(spaspe))
mx=max(log(spaspe))

pngname=paste0("LogSpaSpe_mean_",exam,".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(log(spaspe_grid +1), log(spaspe), ylim=c(mn,mx), type="o",
     xlab="Log(wavenumber +1)", ylab="Log(Spectrum)",
     main=paste("Log Mean spa spectrum. ", exam))
dev.off()


#------------------------------------------------------------------------
# Spatial crf from the estimated spatial spectrum 
# (tested OK)

spacrf=Re( fft(spaspe_full, inverse=TRUE) / n)
spacrf=spacrf/spacrf[1]

dsmax=n/2 # min(20, n/2)
mn=min(spacrf)
mx=max(spacrf)

pngname=paste0("SpaCrl_mean2_",exam,".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(grid_km[1:(dsmax+1)] - grid_km[1], spacrf[1:(dsmax+1)], ylim=c(mn,mx), type="o",
     xlab="Distance, km", ylab="Correlation",
     main=paste("Mean spa crl (from spe) ", exam))
abline(h=0, lwd=1)
dev.off()

#------------------------------------------------------------------------
# Time crf estm

TIM_CRF=TRUE # FALSE TRUE
if(TIM_CRF){
  dtmax=min(ntime/2,200)
  ttcrf=matrix(0, nrow=n, ncol=dtmax+1)
  
  half_t_Sample=FALSE
  
  if(half_t_Sample){
    nt1=2
    dt1=ntime/nt1
    
    for (it1 in 1:nt1){
      t1=(it1 -1)*dt1 +1 
      t2=t1 + ntime/nt1 -1
      for (i in 1:n){
        ttcrf[i,]=est_crf_timser(xi[i,t1:t2], dtmax)$crf
      }
      tcrf=colMeans(ttcrf)
      
      mn=min(tcrf)
      mx=max(tcrf)
      
      pngname=paste0("TimCrl_mean_",t1,"_",exam,".png")
      png(pngname, width=7.1, height=5.1, units = "in", res=300)
      par(mai=c(1.2,1.2,0.7,0.7))
      plot(tgrid_h[1:(dtmax+1)] - grid_km[1], tcrf, ylim=c(mn,mx),
           xlab="Time, h", ylab="Correlation",
           main=paste("Half-ntime mean temporal correlation. ", exam))
      dev.off()
    }
  }
  
  t1=2
  t2=ntime
  
  for (i in 1:n){
    ttcrf[i,]=est_crf_timser(xi[i,t1:t2], dtmax)$crf
  }
  tcrf=colMeans(ttcrf)
  
  mn=min(tcrf)
  mx=max(tcrf)
  
  pngname=paste0("TimCrl_mean_",exam,".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(tgrid_h[1:(dtmax+1)] - grid_km[1], tcrf, ylim=c(mn,mx), type="o",
       xlab="Time, h", ylab="Correlation",
       main=paste("Mean temporal correlation. ", exam))
  dev.off()
  
}

#------------------------------------------------
#------------------------------------------------
# Write output file

outfile=paste0("out_xi_",exam,".txt")
unlink(outfile)
sink(outfile, append=TRUE)

cat("\n")

if(plots_nGau){
  cat("kurtosis_xi=")
  print(kurtosis_xi)
  
  cat("skewness_xi=")
  print(skewness_xi)
  
}

if(TIM_CRF){
  cat("\n")
  
  cat("tcrf[1:25]=")
  print(tcrf[1:25])
  
}


sink()
#------------------------------------------------
