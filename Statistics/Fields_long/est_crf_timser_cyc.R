est_crf_timser_cyc = function(x, m=NULL, stdev=NULL){

  # Estm crf from one realization of the random process on the circle, x[1:n],
  # With small n=length(x), m and stdev should be supplied to get correct results,
  #  otherwise these are evaluated here.
  #  
  # M Tsy 2018
  
  n=length(x)  # n is the number of grid points on the circle
  ind1=c(1:n)
  
  cvf=c(1:n) # init
  crf=cvf
  
  if(is.null(m)) m=mean(x)
  
  if(is.null(stdev)){
    stdev=sd(x)
  }
  
  v = stdev^2
  dx = x-m

  for(ds in (0:(n-1))){
    x1=dx[ind1]
    ind2=(ds:(ds+n-1)) %%n +1  # ind1 shifted cyclically by ds
    x2=dx[ind2]

    cvf[ds+1]=sum(x1*x2) /(n-1)
    crf[ds+1]=cvf[ds+1] / v
  }
  
  return(list("cvf"=cvf, "crf"=crf))
}


