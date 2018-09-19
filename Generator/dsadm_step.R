dsadm_step = function(x_start, n, N, dt, U, rho, nu, sigma, Rem=6.37e6, forcing = TRUE){
  
  # Integrate DSADM FOR ONE TIME STEP using an implicit scheme
  # (1-step prediction): generate an ensemble of fcsts from the enssemble of the initial conditons.
  # 
  # Can be used as the FCST model (ie unforced) or as 
  # the generating model of TRUTH (forced), or as 
  # the generating model for an ensemble member (possible truth, forced).
  # In the two latter (forced) cases, the forcing is generated within this routine
  #  as the white noise multiplied by sigma(t,s).
  # By default, forcing is ON.
  #-------------------------------------------
  # Numerical scheme (see technical doc for details): 
  # Each time step, solve the system of lin algebr eqs for x=x_fcst
  # 
  # G * x = x_start  + forc                            (1)
  # 
  # NB: noise[1:n] is added to x_start, ie before (not after) the time step.
  # 
  # forc(s) = dt*sigma(s)*alpha(s),
  # alpha(s_j) is the discretized spatio-temporal white noise (j=1,...,n):
  # 
  # alpha(s_j) = N(0, VAR=1/(ds*dt)),
  # 
  # whr ds=h=2*pi*Rem /n, 
  #     dt is the time step, so
  # SD(alpha_j) = sqrt(1/(ds*dt))  ==>
  # 
  # forc[j] = dt*sigma[j] * sqrt(1/(ds*dt)) * N(0,1)   (2)
  # 
  # G[,] is specified below (the M mx).
  #
  #------
  #  To reduce the number of multiplications by dt, Eq(1) is rewritten as
  #
  #  Mx=y     (2)
  #
  #  whr M = G/dt and 
  #  y = (x_start  + forc) /dt == x_start/dt + noise
  #  noise[j] = sigma[j] * sqrt(1/(ds*dt)) * N(0,1)
  #
  #  M is the cyclic tridiag mx:
  #  The i-th row of M contains the following non-zero entries: 
  #  
  #  a_i=M[i,im1] (the lower sub-diagonal)
  #  b_i=M[i,i]   (the main      diagonal)
  #  c_i=M[i,ip1] (the upper sub-diagonal)
  #  
  #  whr 
  #  
  #  im1=i-1 unless i=1, in which case im1=n
  #  ip1=i+1 unless i=n, in which case ip1=1
  #  
  #  a,b,c are evaluated as follows (an upwind scheme):
  #  
  #  If U_i >0, then
  #  
  #  a_i = -U_i/dx - nu_i /dx^2
  #  b_i = 1/dt + U_i/dx +rho_i + 2* nu_i /dx^2
  #  c_i = -nu_i /dx^2
  #
  #  If U_i <0, then
  #  
  #  a_i = -nu_i /dx^2
  #  b_i = 1/dt - U_i/dx +rho_i + 2* nu_i /dx^2
  #  c_i = U_i/dx - nu_i /dx^2
  #  
  #-------------------------------------------
  # Arguments:
  #
  # x_start[1:n, 1:N] - ensemble of initial fields
  # n - dim-ty of the state vector x
  # N - ensm size
  # dt - time step, s
  # U[1:n], rho[1:n], nu[1:n], sigma[1:n] - coefficient fields of the DSADM
  # Rem - Earth radius, m
  # forcing - logical switch: if TRUE, then forcing is computed here and added to the solution,
  #                           if FALSE, the pure FCST is computed.
  # 
  # Return: the forecast x_fcst[1:n, 1:N].
  #
  # M Tsyrulnikov
  # June 2018
  #***************************************************************
  
  ## test part1 --------------------
  ## Specify the following rho,no,sigma such that SD(x)=10
  ## U is arbitrary
  #library(limSolve)
  #Rem = 6.37e6 
  #n=10
  #N=2
  #x_start=matrix(0, nrow=n, ncol=N)
  #dt=21600
  #forcing=TRUE
  #U=rep(10, n) ; U[n]=-U[n] ; U[n/2]=-U[n/2]
  #rho=rep(1.56241e-06, n)
  #nu= rep(17380852, n)
  #sigma=rep(46.6087, n)
  ## end test part1 -----------------
  
  #------------------------------
  # Checks
  
  if(nrow(as.matrix(x_start)) != n |
     length(U)       != n |
     length(rho)     != n |
     length(nu)      != n |
     length(sigma)   != n ) {
    stop("dsadm_step: Wrong length of an input vector")
  }
  
  #------------------------------
  
  h=2*pi*Rem /n
  
  #  Mx=y
  #------------------------------
  # (1) The RHS
  #  y = (x_start  + noise) /dt == x_start/dt + sigma(s)* sqrt(1/(h*dt)) *N(0,1)
  
  if(forcing){    # forcing is ON: add dscr white noise scaled by sigma
    
    noise=matrix( sigma[] * sqrt(1/(h*dt)) * rnorm(n*N,0,1) , nrow=n, ncol=N)
    
  } else{         # forcing is OFF
    
    noise=matrix(0, nrow=n, ncol=N)
  }

  y=x_start/dt + noise  # mx of right-hand sides of the system to be solved
  
  #------------------------------
  # (2) M: the matrix of the system to be solved (its 3 diagonals aa, bb, cc)
  # 
  #..........................................
  #    |b_1 c_1 0     ...  0     0     a_1  |
  #    |a_2 b_2 c_2 0 ...  0     0     0    |
  # M= |...                                 |
  #    |0   0   0 ...      a_n-1 b_n-1 c_n-1|
  #    |c_n 0   0     ...  0     a_n   b_n  |
  #..........................................
  
  aa = c(1:n) # init
  bb=aa
  cc=aa
  
  Udh = U[]/h
  nudh2 = nu[] /h^2
  hinv=1/dt
  
  # Upwind fin differences. Fill in aa,bb,cc.

  ind_U_posit = U > 0
  
  aa[ind_U_posit]  = -Udh[ind_U_posit] - nudh2[ind_U_posit]
  aa[!ind_U_posit] =                   - nudh2[!ind_U_posit]
  
  bb = hinv + abs(Udh) + rho + 2*nudh2
  
  cc[ind_U_posit] =                    - nudh2[ind_U_posit]
  cc[!ind_U_posit] = Udh[!ind_U_posit] - nudh2[!ind_U_posit]
  
  ## test  part2 --------------------
  ## Find x by using the standard non-sparse solve() (and forming M explicitly)
  #M=matrix(0, nrow=n, ncol=n)
  #ii=c(1:n)
  #iip1 = (1:n)%%n+1
  #iim1 = (1:n-2)%%n+1
  #for (i in 1:n){
  #  im1=iim1[i]
  #  ip1=iip1[i]
  #  M[i, im1] =aa[i]
  #  M[i, i  ] =bb[i]
  #  M[i, ip1] =cc[i]
  #}
  #x_=solve(M,y)
  #y_=M %*% x_
  #norm(y-y_, "f")  # tested OK (discrep 8e-20)
  ## end test  part2 -----------------
  
  #------------------------------
  # Find the solution using the superposition principle:
  # 
  # With MM := M[-1,-1], yy ":= y[-1], xx := x[-1], chi := x[1]:
  # #..........................................
  #    |b_1 | c_1 0     ...  0     0     a_1  |
  #    |______________________________________|
  #    |a_2 | b_2 c_2 0 ...  0     0     0    |
  # M= |...                                   |
  #    |0   | 0   0 ...      a_n-1 b_n-1 c_n-1|
  #    |c_n | 0   0     ...  0     a_n   b_n  |
  #..........................................
  # (1) solve
  #          MM*uu = yy 
  # and
  #          MM*vv = phi,
  # whr 
  #          phi = -(aa[2], rep(0,n-3), cc[n])
  # (2) find chi:
  #          chi = (y[1] - cc[1]*uu[1] - aa[1]*uu[n-1]) / (bb[1] + cc[1]*vv[1]  + aa[1]*vv[n-1])
  # (3) compute 
  #          xx = uu[] + chi*vv[]
  # (4)
  #          x=c(chi,xx[])
  
  # 3 diagonals of MM
  
  diam1 = aa[3:n]
  dia   = bb[2:n]
  diap1 = cc[2:(n-1)]

  yy=y[-1,]
  phi = -c(aa[2], rep(0,n-3), cc[n])

  uu = Solve.tridiag (diam1, dia, diap1, yy)  # yy and uu are (n-1)*N matrices
  
  ## test  part3 --------------------
  ## The true solution to MM*u=yy is
  #MM = M[-1,-1]
  #uu_ = solve(MM,yy)
  #norm(as.matrix(uu_), "f")
  #norm(as.matrix(uu), "f")
  #norm(as.matrix(uu-uu_), "f")  # tested OK (discrep =0)
  ## end test  part3 ----------------
  
  vv = Solve.tridiag (diam1, dia, diap1, phi) # phi and vv are (n-1)-vectors
  
  chi = (y[1,] - cc[1]*uu[1,] - aa[1]*uu[n-1,]) / (bb[1] + cc[1]*vv[1]  + aa[1]*vv[n-1]) # N-vector
  
  xx = uu + vv %*% t(chi)
  
  x=rbind(chi,xx)
  
  ## test  part4 --------------------
  #norm(as.matrix(x-x_), "f")  # tested OK (discrep =1e-15)
  ## end test  part4 ----------------
  ## 
  #------------------------------
  
  x_fcst=x
  return(x_fcst)
}
