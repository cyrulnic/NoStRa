SpaSmooCVM = function(n, C, spa_shift_max){
  
  #--------------------------------------------------
  # Perform Spatial Smoothing of a CVM for a 1D random field on a regular grid on S1.
  # C_smoo = sum_{s=-smax}^{smax} F^s * C * F^{-s}
  # The smoothing weighting funuction is triangular, even, having itx max at the shift zero,
  # and vanishing at |shift| = spa_shift_max+1
  #
  # Args:
  # n - grid size
  # C - the CVN to be smoothed
  # spa_shift_max=smax - number of spatial shifts to be made in one of the 
  #          two directionsin order to spatially smooth the cvm:
  #   for shift in (-spa_shift_max, spa_shift_max) 
  #     for m in (1,N)  # ensm member
  #       dXfe_shifted[m] = FF^{shift} dXfe[m]  # FF is the forward-shift oprt
  #       CC_shifted[shift] = dXfe_shifted[m] %*% t(dXfe_shifted[m]) /N
  #     end for
  #   end for
  #   C_smoo=sum weight_shift[]*C[], 
  #  weight_shift[1:nshift] - the weights that sum up to 1 and have a triangular shape
  #                  with a max at shift=0 and become zero at shift = +-(spa_shift_max+1) 
  # nshift = 2*spa_shift_max +1
  # 
  # return C_smoo
  # 
  # M Tsy 2019 Mar
  #--------------------------------------------------

  if(spa_shift_max == 0){
    
    C_smoo = C
  
  }else{
    
    nshift = 2*spa_shift_max +1
    
    # specify weight_shift: triangular shape & sum up to 1
    
    weight_shift = c(1:nshift)
    weight_shift[]=0
    weight_shift[1:(spa_shift_max+1)] = c(1:(spa_shift_max+1)) / (spa_shift_max +1)
    weight_shift[(spa_shift_max+1):nshift] = rev(weight_shift[1:(spa_shift_max+1)])
    tmp=sum(weight_shift)
    weight_shift = weight_shift /tmp
    
    # Calc shifted covs
    
    ind_spa = c(1:n)
    
    shift_F0 = spa_shift_max + 1   # no shift
    
    CC_shifted = array(0, dim = c(n,n,nshift))
    CC_shifted[,,shift_F0] = C
    
    for (s in 1:spa_shift_max) {
      ind_Fs = (ind_spa - (s+1))%%n +1  # forw s-shift (indexes)
      ind_Bs = (ind_spa + (s-1))%%n +1  # bckw s-shift (indexes)
      
      C_Fs = C[ind_Fs, ind_Fs]
      C_Bs = C[ind_Bs, ind_Bs]
      
      shift_Fs = spa_shift_max + 1 + s
      CC_shifted[,,shift_Fs] = C_Fs
      
      shift_Bs = spa_shift_max + 1 - s
      CC_shifted[,,shift_Bs] = C_Bs
    }
    
    C_smoo = apply(CC_shifted, c(1,2), function(x) sum(x*weight_shift))
    
  }
  
  return(C_smoo)
}
