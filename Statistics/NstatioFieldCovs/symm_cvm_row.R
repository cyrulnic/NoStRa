symm_cvm_row = function(row, n, i) {
  # From the i-th row of the cvm (of size n),
  # Compute the vector 
  # srow[1:n] such that :
  # 1) it is the cyclic rotation of the original mx row
  # 2) the i-th position of the original row (ie at the diagonal of the cvm)
  #    corresponds to the n/2 (center) position in srow
  # In other words, we just rotate the original row by the angle 2pi/n * (n/2 -i)
  
  srow=row # init
  shift = floor(n/2) -i
  
  for (js in 1:n){ # js is the position in the shifted ("symmetric") row
    
    j=(js - shift -1) %% n  +1  # j is the position in the original mx row
    
    srow[js] = row[j]
  }
  return(srow)
}