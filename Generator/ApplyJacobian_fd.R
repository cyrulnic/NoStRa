ApplyJacobian_fd = function(F, x0, x1, u, eps, ...){
  
  # Approx computation of the Jacobian J of function F(x) applied to u --
  # using finite differencing.
  # Applies J to u.
  # Linearization around the trajectory 
  # x1=F(x0)
  # 
  # Method:
  # 
  # J=dF/dx|_x0
  # F(x0 + eps*u) - F(x0) = eps* J*u + o(eps)  ==>
  # 
  # v=J*u \approx (F(x0 + eps*u) - F(x0)) /eps == (F(x0 + eps*u) - x1) /eps   (*)
  # 
  # Args:
  # 
  # F - function name. 
  #   NB: The 1st argument of F should be x (x0 in this case).
  # x0 - point around which the linearization is done
  # x1 - should be =F(x0)
  # u - "perturbation to which J is to be applied
  # eps - small real number, multiplies u in calculating the finite diference
  # ... - arguments of F after x
  # 
  # M Tsyrulnikov
  # Jan 2018
  #browser()
  
  v = (F(x0 + eps*u, ...) - x1) /eps
  return(v)
}

