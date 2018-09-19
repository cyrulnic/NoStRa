dsadm_start_field_cvm = function(modpar, n, Rem){
  
  #--------------------------------------------------------------------------
  # Generate the initial field z(s) and its CVM 
  # for the stationary advection-diffusion-decay model 
  # with constant coefficients: modpar$rho, modpar$nu, modpar$sigma
  # (note that modpar$U does NOT influence the stationary spectrum and therefore is
  # not used here)
  #
  # Args:
  # 
  # modpar - model params list (U, rho, nu, sigma)
  # n - number if spatial-grid points
  # Rem - Earth radius, m
  # 
  # Return:
  # x_start - initial field on the n-grid
  # CVM - initial field's cov mx
  # 
  # Tested OK.
  # M Tsy 2018 Jul
  #--------------------------------------------------------------------------
  
  a = 1 / sqrt(2*pi*Rem)
  
  rho=modpar$rho
  nu=modpar$nu
  sigma=modpar$sigma
  
  # Statio spectrum:
  # bm=a^2 *sigma^2 / 2 / {rho + nu*m^2/Rem^2 }
  
  mm=seq(from=-(n/2)+1, to= n/2, by=1)
  nmR2pr=nu*(mm/Rem)^2 +rho
  
  bm <- (a*sigma)^2 /2  / nmR2pr  # spectrum from m=-n/2 to m=n/2
  
  # convert the ordinary spectrum to that assumed by fft:
  # shift argument by pi/2: from wvns in (-pi/2 to pi/2) to (0,pi)

  bm_fft=c(1:n)
  bm_fft[1 : ((n/2)+1)] = bm[(n/2) : n]
  bm_fft[((n/2)+1) : n] = rev(bm[((n/2)+1) : n])
  
  # Simulate field:
  # 
  # z_0 \sim {\cal N}(0,1) \cdot \sqrt{b_0}
  # z_{n/2} \sim {\cal N}(0,1) \cdot \sqrt{b_{n/2}}
  # 
  # For m=1,\dots, (n/2)-1},
  # v_m \sim {\cal N}(0,1) \cdot \sqrt{b_m/2}
  # w_m \sim {\cal N}(0,1) \cdot \sqrt{b_m/2}
  # z_m=v_m + \i w_m
  # 
  # Negative wvns:
  # z_{-m} = Conj(z_m)
  # But z_{-m} == z_{n-m} ==>
  # Compute
  # z_{n-m} = Conj(z_m)
  
  zz = c(1:n)   # In the zz vector, its index=m+1
  zz[]=0+0i     # init cplx
  
  # fill in zz (spec coeffs for z) for m \in [0, n/2]

  for (index in (1:((n/2)+1))){
    VAR=bm_fft[index]
    SD=sqrt(VAR)

    if(index == 1 | index == (n/2)+1){   # zz[m] is real for m=0 and m=n/2
      SD_real = SD
      SD_imag = 0
    }else{
      SD_real = SD/sqrt(2)
      SD_imag = SD/sqrt(2)
    }
    
    v = rnorm(1, mean=0, sd=SD_real)
    w = rnorm(1, mean=0, sd=SD_imag)
    zz[index] = v + 1i*w
  }
  
  # Fill in the rest of zz: z_{n-m} = Conj(z_m)
  
  zz[((n/2)+2) : n] = rev(  Conj( zz[2:(n/2)]  )  )

  x_start = fft(zz, inverse = TRUE)
  
  # tested: Im(x_start)=0
  
  x_start = Re(x_start)
  
  #------------------------------------------
  # CVM
  
  cvf = Re(fft(bm_fft, inverse = TRUE))
  
  CVM = matrix(NA, nrow=n, ncol=n)
  
  for (i in 1:n){
    for (j in 1:n){
      d=abs(i-j)
      if(d > n/2) d = n - d
      CVM[i,j] = cvf[d+1]
    }
  }
  
  #------------------------------------------
  # tested: mean((start_field)^2) --> sum(bm[]) for a large number of trials.

  return(list(x_start=x_start,
              CVM=CVM
              ))
}


