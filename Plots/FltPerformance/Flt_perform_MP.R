# Filters' performance wrt KF:
# r(f) := [RMSE(f)- RMSE(KF)] / RMSE(KF)

# 1) As functions of the nstatio-ty strength

#------------
# 1a) EnKF + x

xx <- c(0,1,2,3) # 4 categories: 0=statio, 1=weakly nstatio, 2=medium nstatio (BASIC), 3=strong nststio
nx=length(xx)

# Plot FG RMSEs wrt KF
#  so that the rel err of KF is =0

enkf =0.01*c(2.04,   2.38, 2.92, 4.51)
ec   =0.01*c(0.0,    0.98, 1.54, 3.46)
es   =0.01*c(0.41,   1.26, 1.69, 3.36)
et   =0.01*c(0.28,   1.06, 1.27, 1.78)
hhbef=0.01*c(0.0,    0.70, 0.88, 1.78)


mn=min(min(enkf), min(ec), min(es), min(et), min(hhbef))
mn=0
mx=max(max(enkf), max(ec), max(es), max(et), max(hhbef))

namefile=paste0("RMSE_flt_NstatioStrength_Eplus.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, enkf, type='b',col='black', lwd=3, lty=1, pch=4, cex=1.8,
     xlab="Strength of non-stationarity", 
     xaxp=c(0,3,3),
     ylab="Relative performance score", ylim=c(mn,mx),
     #main=paste("Filters' FG RMSEs relative to KF"),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7
     )
lines(xx, ec, col="orange", type='b', lwd=3, pch=3, cex=1.8)
lines(xx, es, col="green", type='b', lwd=3, pch=0, cex=1.8)
lines(xx, et, col="purple", type='b', lwd=3, pch=2, cex=1.8)
lines(xx, hhbef, col="blue", type='b', lwd=3, pch=8, cex=1.8)

leg.txt<-c('EnKF','EnKF +C (=EnVar)', 'EnKF +S', 'EnKF +T (=HBEF)', 'HHBEF')
leg.col<-c("black",  "orange", "green", "purple", "blue")
legend("topleft", inset=0, leg.txt, col=leg.col, lwd=c(3,3,3,3,3), lty=c(1,1,1,1,1),  
       pch=c(4,3,0,2,8), 
       pt.lwd=3, cex=1.3, pt.cex=1.8, bg="white")

dev.off()


#------------
# 1b) HHBEF - x

h  =0.01*c(0.0,    0.70, 0.88, 1.78)
h_c=0.01*c(0.28,   0.95, 1.04, 1.78)
h_s=0.01*c(0.0,    0.82, 1.16, 1.78)
h_t=0.01*c(0.0,    0.94, 1.28, 3.2)


mn=min(min(h), min(h_c), min(h_s), min(h_t))
mn=0
mx=max(max(h),max(h_c), max(h_s), max(h_t))

namefile=paste0("RMSE_flt_NstatioStrength_Hminus.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, h, type='b',col='blue', lwd=3, lty=1, pch=8, cex=1.8,
     xlab="Strength of non-stationarity", 
     xaxp=c(0,3,3),
     ylab="Relative performance score", ylim=c(mn,mx),
     #main=paste("Filters' FG RMSEs relative to KF"),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7
)
lines(xx, h_c, col="gold2",        type='b', lwd=3, pch=3, cex=1.8)
lines(xx, h_s, col="springgreen3", type='b', lwd=3, pch=0, cex=1.8)
lines(xx, h_t, col="violetred",    type='b', lwd=3, pch=2, cex=1.8)

leg.txt<-c('HHBEF','HHBEF -C', 'HHBEF -S', 'HHBEF -T')
leg.col<-c("blue",  "gold2", "springgreen3",  "violetred")
legend("topleft", inset=0, leg.txt, col=leg.col, lwd=c(3,3,3,3), lty=c(1,1,1,1),  
       pch=c(8,3,0,2), 
       pt.lwd=3, cex=1.3, pt.cex=1.8, bg="white")

dev.off()


#----------------------------------------------
# 2) As functions of the nstatio-ty time (and space) scale

xx <- c(1,2,3) # 3 categories: L_perturb=L*1, L_perturb=L*2 (BASIC), L_perturb=L*3
nx=length(xx)

# Plot FG RMSEs wrt KF
#  so that the rel err of KF is =0

enkf =0.01* c(2.76, 2.92, 3.7)
ec   =0.01* c(1.34, 1.54, 2.45)
es   =0.01* c(1.80, 1.69, 2.11)
et   =0.01* c(1.3,  1.27, 1.27)
hhbef=0.01* c(0.99, 0.88, 1.11)



mn=0
mx=max(max(enkf), max(ec), max(es), max(et), max(hhbef))
mx=0.045



namefile=paste0("RMSE_flt_Lperturb_Eplus.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
plot(xx, enkf, type='b',col='black', lwd=3, lty=1, pch=4, cex=1.8,
     xlab="Time scale of nonstationarity", 
     xaxp=c(0,3,3),
     ylab="Relative performance score", ylim=c(mn,mx),
     #main=paste("Filters' FG RMSEs relative to KF"),
     cex.main=1.7, cex.axis=1.5, cex.lab=1.7
)
lines(xx, ec, col="orange", type='b', lwd=3, pch=3, cex=1.8)
lines(xx, es, col="green", type='b', lwd=3, pch=0, cex=1.8)
lines(xx, et, col="purple", type='b', lwd=3, pch=2, cex=1.8)
lines(xx, hhbef, col="blue", type='b', lwd=3, pch=8, cex=1.8)

leg.txt<-c('EnKF','EnKF +C (=EnVar)', 'EnKF +S', 'EnKF +T (=HBEF)', 'HHBEF')
leg.col<-c("black",  "orange", "green", "purple", "blue")
legend("topleft", inset=0, leg.txt, col=leg.col, lwd=c(3,3,3,3,3), lty=c(1,1,1,1,1),  
       pch=c(4,3,0,2,8), 
       pt.lwd=3, cex=1.3, pt.cex=1.8, bg="white")

dev.off()



