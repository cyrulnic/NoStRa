symm_pd_mx_sqrt = function(A){
  
  #------------------------------------------------------------
  # Compute symmetric pos.definite SQRT of the symm pos-def mx A
  # if there are negative ei-values, these are coerced to be 0.
  # 
  # Return:
  # 
  # 1) $sq      = sqrt(A)
  # 2) $evalmin = min(eval) before setting it to 0.
  # 
  # M Tsy 2017
  #------------------------------------------------------------
    
  ei=eigen(A)
  eval=ei$values
  evalmin=min(eval)
  
  if(is.complex(eval) == TRUE){
    message("symm_pd_mx_sqrt")
    print(eval)
    stop("Non-symmetric imput mx")
  }
  
  ind=which(eval < 0, arr.ind = TRUE)
  eval[ind]=0
  
  evec=ei$vectors
  
  sq = evec %*% diag(sqrt(eval)) %*% t(evec)
  
  return(list("sq"=sq, "evalmin"=evalmin))
}

  

