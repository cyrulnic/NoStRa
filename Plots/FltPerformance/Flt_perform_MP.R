# Filters' performance wrt KF:
# r(f) := [RMSE(f)- RMSE(KF)] / RMSE(KF)

# 1) As functions of the nstatio-ty strength

xx <- c(0,1,2,3) # 4 categories: 0=statio, 1=weakly nstatio, 2=medium nstatio (BASIC), 3=strong nststio
nx=length(xx)

# Plot FG RMSEs wrt KF
#  so that the rel err of KF is =0

Var=  0.01*c(0.0002, 1.35, 2.51, 5.99)
enkf= 0.01*c(2,      2.38, 2.92, 4.51)
ev=   0.01*c(0.0002, 0.98, 1.54, 3.46)
hbef= 0.01*c(0.4,    1.06, 1.27, 1.78)
hhbef=0.01*c(0.002,  0.82, 1.16, 1.78)

mn=min(min(Var), min(enkf), min(ev), min(hbef), min(hhbef))
mn=0
mx=max(max(Var), max(enkf), max(ev), max(hbef), max(hhbef))

namefile=paste0("RMSE_flt_NstatioStrength.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, Var, type='l' ,col='black', lwd=3, lty=3,
     xlab="Strength of non-stationarity", 
     xaxp=c(0,3,3),
     ylab="Relative deterioration wrt KF", ylim=c(mn,mx),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7,
     main=paste("Filters' FG RMSEs relative to KF")
     )
lines(xx, enkf, col="orange", lwd=2, lty=2)
lines(xx, ev, col="green", lwd=1)
lines(xx, hbef, col="purple", lwd=2)
lines(xx, hhbef, col="blue", lwd=3)

leg.txt<-c('(Var-KF)/KF','(EnKF-KF)/KF', '(EnVar-KF)/KF', '(HBEF-KF)/KF', '(HHBEF-KF)/KF')
leg.col<-c("black",  "orange", "green", "purple", "blue")
legend("topleft", inset=0, leg.txt, col=leg.col, lwd=c(3,2,1,2,3), lty=c(3,2,1,1,1),  
       pt.lwd=3, cex=1.3, pt.cex=3, bg="white")

dev.off()


#----------------------------------------------
# 2) As functions of the nstatio-ty time (and space) scale

xx <- c(1,2,3) # 3 categories: L_perturb=L*1, L_perturb=L*2 (BASIC), L_perturb=L*3
nx=length(xx)

# Plot FG RMSEs wrt KF
#  so that the rel err of KF is =0

Var  =0.01* c(2.04, 2.51, 3.51)
enkf =0.01* c(2.76, 2.92, 3.70)
ev   =0.01* c(1.34, 1.54, 2.45)
hbef =0.01* c(1.35, 1.27, 1.33)
hhbef=0.01* c(1.12, 1.16, 1.33)



mn=min(min(Var), min(enkf), min(ev), min(hbef), min(hhbef))
mn=0
mx=max(max(Var), max(enkf), max(ev), max(hbef), max(hhbef))
mx=0.045


namefile=paste0("RMSE_flt_Lperturb.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, Var, type='l' ,col='black', lwd=3, lty=3,
     xlab="Time scale of nonstationarity", 
     xaxp=c(0,3,3),
     ylab="Relative deterioration wrt KF", ylim=c(mn,mx*1.00),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7,
     main=paste("Filters' FG RMSEs relative to KF")
)
lines(xx, enkf, col="orange", lwd=2, lty=2)
lines(xx, ev, col="green", lwd=1)
lines(xx, hbef, col="purple", lwd=2)
lines(xx, hhbef, col="blue", lwd=3)

leg.txt<-c('(Var-KF)/KF','(EnKF-KF)/KF', '(EnVar-KF)/KF', '(HBEF-KF)/KF', '(HHBEF-KF)/KF')
leg.col<-c("black",  "orange", "green", "purple", "blue")
legend("topleft", inset=0, leg.txt, col=leg.col, lwd=c(3,2,1,2,3), lty=c(3,2,1,1,1),  
       pt.lwd=3, cex=1.3, pt.cex=3, bg="white")
dev.off()




#----------------------------------------------
# 3) As functions of the ensm size N

xx <- c(10,20,30) # 3 categories: N=10 (BASIC) , N=20, and N=30
nx=length(xx)

# Plot FG RMSEs wrt KF
#  so that the rel err of KF is =0

Var=  0.01* c(2.51, 2.51, 2.51)
enkf= 0.01* c(2.92, 1.93, 1.06)
ev=   0.01* c(1.54, 1.30, 0.83)
hbef= 0.01* c(1.27, 1.08, 0.55)
hhbef=0.01* c(1.16, 0.92, 0.49)


mn=min(min(Var), min(enkf), min(ev), min(hbef), min(hhbef))
mn=0
mx=max(max(Var), max(enkf), max(ev), max(hbef), max(hhbef))

namefile=paste0("RMSE_flt_N.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, Var, type='l' ,col='black', lwd=3, lty=3,
     xlab="Ensemble size", 
     xaxp=c(10,30,2),
     ylab="Relative deterioration wrt KF", ylim=c(mn,mx*1.00),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7,
     main=paste("Filters' FG RMSEs relative to KF")
)
lines(xx, enkf, col="orange", lwd=2, lty=2)
lines(xx, ev, col="green", lwd=1)
lines(xx, hbef, col="purple", lwd=2)
lines(xx, hhbef, col="blue", lwd=3)

leg.txt<-c('(Var-KF)/KF','(EnKF-KF)/KF', '(EnVar-KF)/KF', '(HBEF-KF)/KF', '(HHBEF-KF)/KF')
leg.col<-c("black",  "orange", "green", "purple", "blue")
legend("bottomleft", inset=0, leg.txt, col=leg.col, lwd=c(3,2,1,2,3), lty=c(3,2,1,1,1),  
       pt.lwd=3, cex=1.3, pt.cex=3, bg="white")
dev.off()





#----------------------------------------------
# 4) As functions of stride

xx <- c(6,12,24) # 3 categories: stride=1, stride=2 (BASIC) , stride=4
nx=length(xx)

# Plot FG RMSEs wrt KF
#  so that the rel err of KF is =0

Var  =0.01* c(2.69, 2.51, 4.14)
enkf =0.01* c(4.12, 2.92, 2.16)
ev   =0.01* c(2.15, 1.54, 1.65)
hbef =0.01* c(1.55, 1.27, 0.28)
hhbef=0.01* c(1.55, 1.16, 0.27)


mn=min(min(Var), min(enkf), min(ev), min(hbef), min(hhbef))
mn=0
mx=max(max(Var), max(enkf), max(ev), max(hbef), max(hhbef))

namefile=paste0("RMSE_flt_stride.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, Var, type='l' ,col='black', lwd=3, lty=3,
     xlab="Interval between observations, h", 
     xaxp=c(6,24,3),
     ylab="Relative deterioration wrt KF", ylim=c(mn,mx*1.00),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7,
     main=paste("Filters' FG RMSEs relative to KF")
)
lines(xx, enkf, col="orange", lwd=2, lty=2)
lines(xx, ev, col="green", lwd=1)
lines(xx, hbef, col="purple", lwd=2)
lines(xx, hhbef, col="blue", lwd=3)

leg.txt<-c('(Var-KF)/KF','(EnKF-KF)/KF', '(EnVar-KF)/KF', '(HBEF-KF)/KF', '(HHBEF-KF)/KF')
leg.col<-c("black",  "orange", "green", "purple", "blue")
legend("bottomleft", inset=0, leg.txt, col=leg.col, lwd=c(3,2,1,2,3), lty=c(3,2,1,1,1),  
       pt.lwd=3, cex=1.3, pt.cex=3, bg="white")
dev.off()

#----------------------------------------------
#----------------------------------------------
# Compare 
# a) the gain provided by B_clim   on top of the current ensemble covs,
# b) the gain provided by B_recent on top of the current ensemble covs
# -- both as functions of stride.

xx <- c(1,2,3) # 3 categories: stride=1, stride=2 (BASIC), stride=4
nx=length(xx)

# Plot r(ev)/r(enkf) and r(hbef)/r(enkf) (the lower the better)
# Recall that
# r(f) := [RMSE(f)- RMSE(KF)] / RMSE(KF)


gclim  =c(0.48, 0.47, 0.24)
grecent=c(0.62, 0.57, 0.46)


mn=0
mx=max(max(gclim), max(grecent))

namefile=paste0("gain_clim_recent_stride.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, gclim, type='l' ,col='green', lwd=1.5, lty=1,
     xlab="Observations stride", 
     xaxp=c(1,3,2),
     ylab="Roles (the higher the better)", ylim=c(mn,mx*1.00),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7,
     main=expression( paste( 'Roles of B'[clim], ' (thin green) and B'[recent], ' (thick red)') )
)
lines(xx, grecent, col="red", lwd=3)

leg.txt<-(c('[r(EnKF) - r(EnVar)] / r(EnKF)','[r(EnKF) - r(HBEF)] / r(EnKF)'))
leg.col<-c("green",  "red")
legend("bottomleft", inset=0, leg.txt, col=leg.col, lwd=c(1.5,3), lty=c(1,1),  
       pt.lwd=3, cex=1.3, pt.cex=3, bg="white")

dev.off()
