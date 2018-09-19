ext2mod_params4dsadm_statio = function(extpar, Rem) {
  
  #--------------------------------------------------------------------------
  # Converts the list of the Stationary-DSADM's External Parameters
  # 
  # extpar$U=U, m/s
  # extpar$L=L, m
  # extpar$V_char=V_char, m/s
  # extpar$SD=SD
  # 
  # to the list of the Stationary-DSADM's Model Parameters (which is returned)
  # 
  # modpar$U=U
  # modpar$rho=rho
  # modpar$nu=nu
  # modpar$sigma=sigma
  # 
  # Rem - Earth radius, m
  # 
  # Return:
  # 
  # modpar
  # 
  # M Tsy 2018 June
  #--------------------------------------------------------------------------
  
  modpar=list() # init 
  
  L=extpar$L
  V_char=extpar$V_char
  
  T=L/V_char
  
  # Find rho from T:
  # rho = 1/T * S2 / S1
  
  mm=seq(from=-(n/2)+1, to= n/2, by=1)
  LmdR2p1=(L*mm/Rem)^2 +1
  S2=sum( LmdR2p1^(-2) )
  S1=sum( LmdR2p1^(-1) )
  rho=S2/S1/T
  
  # Find nu=rho*L^2
  
  nu=rho*L^2
  
  # Find sigma:
  # SD^2 = a^2 * sigma^2 / (2 rho) * S1 , for S1 see above  ==>
  # sigma = SD/a * sqrt(2*rho/S1)
  
  a = 1 / sqrt(2*pi*Rem)
  sigma = extpar$SD/a * sqrt(2*rho/S1)
  
  # output variable:
  
  modpar$U=extpar$U
  modpar$rho=rho
  modpar$nu=nu
  modpar$sigma=sigma
  
  return(modpar) # contains U, rho, nu, sigma
}