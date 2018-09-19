est_crf_timser = function(xx, dtmax){
  # xx[1:ntime]
  
  ntime=length(xx)
  nt=ntime-dtmax
  ind=c(1:nt)
  
  cvf=c(0:dtmax)
  crf=cvf

  for(dt in (0:dtmax)){
    x1=xx[ind]
    x2=xx[ind+dt]
    m1=mean(x1)
    m2=mean(x2)
    x1=x1-m1
    x2=x2-m2
    
    s1=sd(x1)
    s2=sd(x2)
    cvf[dt+1]=sum(x1*x2) /(nt-1)
    crf[dt+1]=cvf[dt+1] / (s1*s2)
  }
  
  return(list("cvf"=cvf, "crf"=crf))
}


est_crosscrf_timser = function(xx, yy, dtmax){
  # xx[1:ntime]
  # yy[1:ntime]
  # NB: Positive shifts dt!:
  #   crl(xx(t), yy(t+dt))
  #   (i.y. yy is later than xx)
  
  nt=ntime-dtmax
  ind=c(1:nt)
  
  cvf=c(0:dtmax)
  crf=cvf
  
  for(dt in (0:dtmax)){
    x1=xx[ind]
    x2=yy[ind+dt]
    m1=mean(x1)
    m2=mean(x2)
    x1=x1-m1
    x2=x2-m2
    
    s1=sd(x1)
    s2=sd(x2)
    cvf[dt+1]=sum(x1*x2) /(nt-1)
    crf[dt+1]=cvf[dt+1] / (s1*s2)
  }
  
  return(list("cvf"=cvf, "crf"=crf))
}