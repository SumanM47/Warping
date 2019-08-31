## Code for full method

Forecast_func_f <- function(Y_in, X_in, s, sp, cell, Xf, Yf = NULL,
                            iters = 1000000, burn = 600000, thin = 40, step_int = 30, ridge = FALSE,
                            L = 15, 
                            knots_param,
                            rho_a = 0.95, rho_x = 0.95,
                            a_med = 0.1,
                            sig_0_a = 0.01, sig_0_b = 0.01,
                            rho_a0 = 0.95,
                            a_y = 0.01, b_y = 0.01){
  
  source(paste(getwd(),"/auxiliary/Rfunctions.R",sep=""))
  require(Rcpp)
  require(RcppArmadillo)
  sourceCpp(paste(getwd(),"/auxiliary/cppfuncs.cpp",sep=""))
  require(spdep)
  
  Hour_Y <- Y_in
  Hour_X <- X_in
  
  r_add <- 0
  if(ridge){r_add <- 1e-5}
  
  min_lon <- min(c(s[,1],sp[,1]))
  max_lon <- max(c(s[,1],sp[,1]))
  min_lat <- min(c(s[,2],sp[,2]))
  max_lat <- max(c(s[,2],sp[,2]))
  
  nx <- length(unique(sp[,1]))
  ny <- length(unique(sp[,2]))
  
  d1 <- (max_lon - min_lon)/nx
  d2 <- (max_lat - min_lat)/ny
  
  ksx <- knots_param[1]
  ksy <- knots_param[2]
  
  
  
  ## Initial Values & Bookkeeping
  
  n_y <- nrow(Y_in)
  n_x <- nrow(X_in)
  
  KS <- ksx*ksy
  non_mis_y_ind <- which(!is.na(Hour_Y),arr.ind = TRUE)
  mis_y_ind <- which(is.na(Hour_Y),arr.ind = T)
  my_ind <- which(is.na(Hour_Y))
  mis_B_ind <- NULL
  if(length(mis_y_ind)>0){
  mis_B_ind <- cbind(mis_y_ind,1)
  if(L>1){
  for(ll in 2:L){
    mis_B_ind <- rbind(mis_B_ind,cbind(mis_y_ind,ll))
  }}
  }
  
  nhnmnum <- apply(Hour_Y,1,function(x){sum(!is.na(x))})
  s_y <- s
  s_x <- sp
  
  
  x.proj <- unique(sp[,1])
  y.proj <- unique(sp[,2])
  m1 <- length(x.proj)
  m2 <- length(y.proj)
  mm <- ncol(Y_in)
  np <- nrow(Xf)*ncol(Xf)
  
  Hour_smoothX <- Hour_X[cell,]
  
  
  w1 <- 2*pi*(seq(0,m1-1,1))/m1        
  w2 <- 2*pi*(seq(0,m2-1,1))/m2
  w1 <- matrix(w1,m1,m2,byrow=FALSE)
  w2 <- matrix(w2,m1,m2,byrow=TRUE)
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  
  nhx <- ncol(Hour_X)
  nhxf <- ncol(Xf)
  nsx <- nrow(Y_in)
  nsxf <- nrow(Xf)
  
  
  nbx <- nhx*nsx
  nbxf <- nhxf*nsxf
  
  Bptilde <- array(0,c(m1*m2,ncol(Hour_X),L))
  Bpftilde <- array(0,c(m1*m2,ncol(Xf),L))
  Xmat <- matrix(1,n_y*ncol(Hour_X),L+1)
  Xfmat <- matrix(1,n_x*ncol(Xf),L+1)
  
  for(nh in 1:ncol(Hour_X)){
    z <- fft(matrix(Hour_X[,nh],nx,ny),inverse=TRUE) 
    for(l in 1:L){
      W           <- dbinom(l-1,L-1,delta/(2*pi))
      x           <- fft(z*W)
      tempbas <- Re(x)/(m1*m2)
      Bptilde[,nh,l] <- as.vector(tempbas)
      Xmat[(((nh-1)*n_y + 1): (nh*n_y)),(l+1)] <- as.vector(tempbas)[cell]
    }
  }
  for(nh in 1:ncol(Xf)){
    zf <- fft(matrix(Xf[,nh],nx,ny),inverse=TRUE)
    for(l in 1:L){
      W           <- dbinom(l-1,L-1,delta/(2*pi))
      x           <- fft(zf*W)
      tempbas <- Re(x)/(m1*m2)
      Bpftilde[,nh,l] <- as.vector(tempbas)
      Xfmat[(((nh-1)*n_x + 1): (nh*n_x)),(l+1)] <- as.vector(tempbas)
    }
  }
  
  Btilde <- subarr(Bptilde,cell)
  Btilde[mis_B_ind] <- 0
  
  
  A_mat <- matrix(0,KS,2)
  
  BB <- basis(s_x,ksx,ksy)
  AS_x <- BB$B
  AS_y <- submat(AS_x,cell)
  ws_x <- warp(sp,AS_x,A_mat)
  ws_y <- warp(s,AS_y,A_mat)
  
  knotsx <- BB$knotsx
  knotsy <- BB$knotsy
  knotsS <- as.matrix(expand.grid(knotsx,knotsy))
  
  
  temp_ob <- dnearneigh(knotsS,0,1)
  temp1 <- nb2listw(temp_ob)
  temp2 <- listw2sn(temp1)
  el <- cbind(temp2$from,temp2$to)
  C <- matrix(0,KS,KS)
  C[el] <- 1
  DC <- diag(rowSums(C))
  mat_as_inv <- DC - rho_a*C
  C_mat_as_inv <- chol(mat_as_inv)
  logdet_mat_as_inv <- 2*sum(log(diag(C_mat_as_inv)))
  
  mat_a0_inv <- DC - rho_a0*C
  C_mat_a0_inv <- chol(mat_a0_inv)
  hlogdet_mat_a0_inv <- sum(log(diag(C_mat_a0_inv)))
  
  
  D_x_inv <- diag(L+1);
  if(L > 1){
    mat_x_inv <- 2*diag(L)
    mat_x_inv[1,1] <- mat_x_inv[L,L] <- 1
    diag(mat_x_inv[-1,]) <- -rho_x
    diag(mat_x_inv[,-1]) <- -rho_x
    D_x_inv[-1,-1] <- mat_x_inv
  }
  
  Sig_2 <- diag(2)
  logdet_Sig_2 <- 0
  Sig_2_inv <- solve(Sig_2)
  
  
  
  Sig_a_inv <- kronecker(Sig_2_inv,mat_as_inv)
  
  
  match_ind_x <- roundnmatch(ws_x,d1=d1,d2=d2,nx=nx,ny=ny,min_lon=min_lon,min_lat=min_lat)
  match_ind_y <- roundnmatch(ws_y,d1=d1,d2=d2,nx=nx,ny=ny,min_lon=min_lon,min_lat=min_lat)
  
  warp_Xf <- Xfmat
  warp_Xmat <- Xmat
  warp_Xfs <- matrix(1,nrow(Y_in)*ncol(Xf),L+1)
  Btilde <- subarr(Bptilde,match_ind_y)
  Btilde[mis_B_ind] <- 0
  warp_Xmat[,-1] <- matrix(Btilde,ncol = L)
  warp_Xmat[my_ind,] <-0
  
  
  scale_a <- a_med
  
  sigma_a <- scale_a
  
  sigma_beta_y <- 10
  
  tempfit <- lm(as.vector(Hour_Y)~warp_Xmat+0)
  b0 <- rep(0,KS)
  sigma_y <- summary(tempfit)$sigma^2
  
  tau.a <- rep(0.04,KS)
  tau.va <- tau.x <- 0.5
  tau.sig_a <- 0.5
  tau.a0 <- 0.5
  
  acc.a <- rep(0,length(tau.a))
  acc.va <- acc.x <- acc.sig_a <- 0
  acc.a0 <- 0
  
  keep.a <- array(0,dim=c((iters-burn)/thin,KS,2))
  keep.sigma_a <- rep(0,(iters-burn)/thin)
  keep.sigma_y <- rep(0,(iters-burn)/thin)
  keep.rho_a <- rep(0,(iters - burn)/thin)
  keep.rho_x <- rep(0,(iters - burn)/thin)
  keep.rho_a0 <- rep(0,(iters - burn)/thin)
  keep.sigma_0 <- rep(0,(iters - burn)/thin)
  keep.b0 <- matrix(0,(iters - burn)/thin,KS)
  forecast <- matrix(0,n_x,ncol(Xf))
  warphat <- matrix(0,nrow(sp),2)
  keep.warp_sig <- array(0,c((iters - burn)/thin,nrow(sp),2))
  if(!is.null(Yf)){forekeepmean <- list(); forekeepsd <- list()}
  
  bb.mn <- l1 <- l2 <- linf <- tempforecastvec <- MSE <- MAD <- 0
  
  Hour_Y0 <- Hour_Y
  Hour_Y0[which(is.na(Hour_Y),arr.ind = T)] <- 0
  Hour_Z <- Hour_Y - as.vector(AS_y%*%b0)
  Hour_Z0 <- Hour_Z
  Z_mis_ind <- which(is.na(Hour_Z),arr.ind = T)
  Hour_Z0[Z_mis_ind] <- 0
  
  ## GO!
  
  for(i in 1:iters){
    
    
    
    ## sigma_0
    
    b0_quad <- crossprod(C_mat_a0_inv%*%b0)
    sigma_0 <- 1/rgamma(1,sig_0_a + 0.5*n_y, sig_0_b + 0.5*b0_quad)
    
    ## rho_a0
    
    l_rho_a0 <- 0.5*b0_quad/sigma_0 + hlogdet_mat_a0_inv + 9*log(rho_a0)
    
    va0 <- log(rho_a0/(1-rho_a0))
    can_va0 <- va0 + tau.a0*rnorm(1)
    can_rho_a0 <- 1/(1+exp(-can_va0))
    
    can_mat_a0_inv <- DC - can_rho_a0*C
    can_C_mat_a0_inv <- chol(can_mat_a0_inv)
    can_hlogdet_mat_a0_inv <- sum(log(diag(can_C_mat_a0_inv)))
    
    can_l_rho_a0 <- 0.5*crossprod(can_C_mat_a0_inv%*%b0)/sigma_0 + can_hlogdet_mat_a0_inv + 9*log(can_rho_a0)
    
    a_l_rho_a0 <- can_l_rho_a0 - l_rho_a0 - can_va0 - 2*log(1+exp(-can_va0)) + va0 + 2*log(1+exp(-va0))
    if(log(runif(1)) < a_l_rho_a0){
      rho_a0 <- can_rho_a0
      va0 <- can_va0
      mat_a0_inv <- can_mat_a0_inv
      C_mat_a0_inv <- can_C_mat_a0_inv
      hlogdet_mat_a0_inv <- can_hlogdet_mat_a0_inv
      acc.a0 <- acc.a0 + 1/step_int
    }
    
    ## b0
    
    XtX <- crossprod(warp_Xmat)
    BtKtX <- t(AS_y)%*%(cbind(nhnmnum,apply(Btilde,c(1,3),"sum")))
    Sig_temp <- XtX + (1/sigma_beta_y)*D_x_inv
    C_Sig_temp <- chol(Sig_temp)
    lterm <- crossprod(forwardsolve(t(C_Sig_temp),t(BtKtX)))
    
    Sigma_new_inv <- (mat_a0_inv/sigma_0) + (crossprod(sqrt(nhnmnum)*AS_y)/sigma_y) - (lterm/sigma_y)
    
    XtY_0 <- crossprod(warp_Xmat,as.vector(Hour_Y0))
    BtKtY <- as.vector(t(AS_y)%*%rowSums(Hour_Y,na.rm = T))
    
    Sigma_new <- solve(Sigma_new_inv)
    
    MM1 <- (BtKtY - as.vector(BtKtX%*%(chol2inv(C_Sig_temp)%*%XtY_0)))/sigma_y
    MM2 <- t(chol(Sigma_new))%*%rnorm(KS)
    
    b0 <- as.vector(Sigma_new%*%MM1 + MM2)
    
    
    int_0 <- as.vector(AS_y%*%b0)
    Hour_Z <- Hour_Y - int_0
    Hour_Z0 <- Hour_Z
    Hour_Z0[Z_mis_ind] <- 0
    
    ## rho_a
    
    l_rho_a <- -0.5*(crossprod(C_mat_as_inv%*%A_mat[,1]) + crossprod(C_mat_as_inv%*%A_mat[,2]))/sigma_a + logdet_mat_as_inv + 9*log(rho_a)
    
    va <- log(rho_a/(1-rho_a))
    can_va <- va + tau.va*rnorm(1)
    can_rho_a <- 1/(1+exp(-can_va))
    
    can_mat_as_inv <- DC - can_rho_a*C
    can_C_mat_as_inv <- chol(can_mat_as_inv)
    can_logdet_mat_as_inv <- 2*sum(log(diag(can_C_mat_as_inv)))
    
    can_l_rho_a <- -0.5*(crossprod(can_C_mat_as_inv%*%A_mat[,1]) + crossprod(can_C_mat_as_inv%*%A_mat[,2]))/sigma_a + can_logdet_mat_as_inv + 9*log(can_rho_a)
    
    a_l_rho_a <- can_l_rho_a - l_rho_a - can_va - 2*log(1+exp(-can_va)) + va + 2*log(1+exp(-va))
    if(log(runif(1)) < a_l_rho_a){
      rho_a <- can_rho_a
      va <- can_va
      mat_as_inv <- can_mat_as_inv
      C_mat_as_inv <- can_C_mat_as_inv
      logdet_mat_as_inv <- can_logdet_mat_as_inv
      acc.va <- acc.va + 1/step_int
    }
    
    
    ## Sig_2
    
    temp_as <- 0.5*(crossprod(C_mat_as_inv%*%A_mat[,1]) + crossprod(C_mat_as_inv%*%A_mat[,2]))
    
    l_sig_a <- -KS*log(sigma_a) -temp_as/sigma_a - (2.58^2)*0.5*(sigma_a^2)
    
    can_sigma_a <- exp(log(sigma_a) + tau.sig_a*rnorm(1))
    can_l_sig_a <- -KS*log(can_sigma_a) -temp_as/can_sigma_a - (2.58^2)*0.5*(can_sigma_a^2)
    
    a_sig_a <- can_l_sig_a - l_sig_a - log(sigma_a) + log(can_sigma_a)
    if(log(runif(1)) < a_sig_a){
      sigma_a <- can_sigma_a
      acc.sig_a <- acc.sig_a + 1/step_int
    }
    
    
    
    logdet_Sig_y_Inv <- -as.numeric(determinant((1/sigma_beta_y)*D_x_inv + XtX,log = T)$modulus)
    quad <- sum(Hour_Z0^2) - crossprod(forwardsolve(t(C_Sig_temp),crossprod(warp_Xmat,as.vector(Hour_Z0))))
    
    temp <- 0.5*(logdet_Sig_y_Inv - quad/sigma_y)
    
    
    ## rho_x
    
    l_rho_x <- temp + 9*log(rho_x)
    
    vx <- log(rho_x/(1-rho_x))
    can_vx <- vx + tau.x*rnorm(1)
    can_rho_x <- 1/(1+exp(-can_vx))
    
    can_mat_x_inv <- mat_x_inv
    diag(can_mat_x_inv[-1,]) <- -can_rho_x
    diag(can_mat_x_inv[,-1]) <- -can_rho_x
    can_D_x_inv <- D_x_inv
    can_D_x_inv[-1,-1] <- can_mat_x_inv
    
    
    can_logdet_Sig_y_Inv <- -as.numeric(determinant((1/sigma_beta_y)*can_D_x_inv + XtX,log = T)$modulus)
    can_quad <- sum(Hour_Z0^2) - crossprod(forwardsolve(t(chol(XtX + (1/sigma_beta_y)*can_D_x_inv)),crossprod(warp_Xmat,as.vector(Hour_Z0))))
    
    can_temp <- 0.5*(can_logdet_Sig_y_Inv - can_quad/sigma_y)
    
    can_l_rho_x <- can_temp + 9*log(can_rho_x)
    
    a_l_rho_x <- can_l_rho_x - l_rho_x - can_vx - 2*log(1+exp(-can_vx)) + vx + 2*log(1+exp(-vx))
    if(log(runif(1)) < a_l_rho_x){
      rho_x <- can_rho_x
      vx <- can_vx
      temp <- can_temp
      mat_x_inv <- can_mat_x_inv
      D_x_inv <- can_D_x_inv
      logdet_Sig_y_Inv <- can_logdet_Sig_y_Inv
      quad <- can_quad
      acc.x <- acc.x + 1/step_int
    }
    
    ## a
    
    for(k in 1:KS){
      
      temp_a <- 0.5*(crossprod(C_mat_as_inv%*%A_mat[,1]) + crossprod(C_mat_as_inv%*%A_mat[,2]))
      
      l_a <- temp - temp_a/sigma_a 
      
      can_a <- A_mat
      
      can_a[k,] <- A_mat[k,] + tau.a[k]*rnorm(2)
      
      
      ww_y <- s_y+AS_y%*%can_a
      can_match_ind_y <- roundnmatch(ww_y,d1,d2,nx,ny,min_lon=min_lon,min_lat=min_lat)
      can_Btilde <- subarr(Bptilde,can_match_ind_y)
      can_Btilde[mis_B_ind] <- 0
      can_warp_Xmat <- warp_Xmat
      can_warp_Xmat[,-1] <- matrix(can_Btilde,ncol = L)
      can_warp_Xmat[my_ind,1] <- 0
      
      can_XtX <- crossprod(can_warp_Xmat)
      can_logdet_Sig_y_Inv <- -as.numeric(determinant((1/sigma_beta_y)*D_x_inv + can_XtX,log = T)$modulus)
      can_quad <- sum(Hour_Z0^2) - crossprod(forwardsolve(t(chol(can_XtX + (1/sigma_beta_y)*D_x_inv)),crossprod(can_warp_Xmat,as.vector(Hour_Z0))))
      
      can_temp <- 0.5*(can_logdet_Sig_y_Inv - can_quad/sigma_y)
      can_temp_a <- 0.5*(crossprod(C_mat_as_inv%*%can_a[,1]) + crossprod(C_mat_as_inv%*%can_a[,2]))
      
      can_l_a <- can_temp - can_temp_a/sigma_a
      
      a_l_a <- can_l_a - l_a
      if(log(runif(1)) < a_l_a){
        A_mat <- can_a
        ws_y <- ww_y
        match_ind_y <- can_match_ind_y
        acc.a[k] <- acc.a[k] + 1/step_int
        temp <- can_temp
      }
    }
    
    
    
    Btilde <- subarr(Bptilde,match_ind_y)
    Btilde[mis_B_ind] <- 0
    warp_Xmat[,-1] <- matrix(Btilde,ncol = L)
    warp_Xmat[my_ind,1] <- 0
    XtX <- crossprod(warp_Xmat)
    logdet_Sig_y_Inv <- -as.numeric(determinant((1/sigma_beta_y)*D_x_inv + XtX,log = T)$modulus)
    C_mid <- chol(XtX + (1/sigma_beta_y)*D_x_inv)
    XtY <- crossprod(warp_Xmat,as.vector(Hour_Z0))
    quad <- sum(Hour_Z0^2) - crossprod(forwardsolve(t(C_mid),XtY))
    
    
    
    ## sigma_y
    
    sigma_y <- 1/rgamma(1,a_y+0.5*nrow(non_mis_y_ind), b_y + 0.5*quad)
    
    
    ## Step_Size
    
    if(i <= burn/4 & i%%step_int == 0 & i > 20){
      
      tau.a <- ifelse(acc.a < 0.3, tau.a*0.8, ifelse(acc.a > 0.5, tau.a*1.2,tau.a))
      tau.va <- ifelse(acc.va < 0.3, tau.va*0.8, ifelse(acc.va > 0.5, tau.va*1.2, tau.va))
      tau.x <- ifelse(acc.x < 0.3, tau.x*0.8, ifelse(acc.x > 0.5, tau.x*1.2, tau.x))
      tau.sig_a <- ifelse(acc.sig_a < 0.3, tau.sig_a*0.8, ifelse(acc.sig_a > 0.5, tau.sig_a*1.2, tau.sig_a))
      tau.a0 <- ifelse(acc.a0 < 0.3, tau.a0*0.8, ifelse(acc.a0 > 0.5, tau.a0*1.2, tau.a0))
      
      acc.a <- rep(0,length(tau.a))
      acc.va <- 0
      acc.x <- 0
      acc.sig_a <- 0
      acc.a0 <- 0
      
    }
    
    ## Store
    if(i%%thin==0 & i > burn){
      index <- (i - burn)/thin
      keep.sigma_a[index] <- sigma_a
      keep.a[index,,] <- A_mat
      keep.sigma_y[index] <- sigma_y
      keep.rho_a[index] <- rho_a
      keep.rho_x[index] <- rho_x
      keep.rho_a0[index] <- rho_a0
      keep.sigma_0[index] <- sigma_0
      keep.b0[index,] <- b0
    }
    
    
    if(i > burn & i%%thin == 0){
      ws_x <- s_x + AS_x%*%A_mat
      keep.warp_sig[(i - burn)/thin,,] <- ws_x - s_x
      match_ind_x <- roundnmatch(ws_x,d1,d2,nx,ny,min_lon=min_lon,min_lat=min_lat)
      warphat <- warphat + ws_x/((iters - burn)/thin)
      
      
      Bpftildenew <- subarr(Bpftilde,match_ind_x)
      Bpp <- subarr(Bpftilde,match_ind_y)
      warp_Xf[,-1] <- matrix(Bpftildenew,ncol = L)
      warp_Xfs[,-1] <- matrix(Bpp,ncol = L)
      rm("Bpftildenew","Bpp")
      Cm <- chol2inv(C_mid)
      temppred <- matrix(warp_Xf%*%(Cm%*%XtY),n_x,ncol(Xf))
      if(!is.null(Yf)){d0 <- rowSums((AS_y%*%Sigma_new)*AS_y,na.rm = T);AAA <- .subset(temppred,cell,1:ncol(Xf)); forekeepmean[[(i-burn)/thin]] <- as.vector(AAA + int_0); Res <- (AAA + int_0 - Yf); MSE <- MSE + mean(Res^2,na.rm = T); MAD <- MAD + mean(abs(Res),na.rm=T); forekeepsd[[(i-burn)/thin]] <- sqrt(rep(d0,nhxf) + sigma_y*(1+as.vector(get_sig2(warp_Xfs,Cm))))}
      
      
      ## Get beta_0f
      
      b0f <- as.vector(AS_x%*%b0)
      #beta_0f <- rep(b0f,nhxf)
      
      
      tempforecastvec <- tempforecastvec + temppred + b0f
      rm("temppred","AAA")
      
      
    }
    
  }
  
  Forecast <- tempforecastvec/((iters-burn)/thin)
  qarr <- apply(keep.warp_sig,2:3,"quantile",prob  = c(0.0125,0.9875))
  non_sig_ind <- which(qarr[1,,1] <= 0 & 0 <= qarr[2,,1] & qarr[1,,2] < 0 & 0 < qarr[2,,2])
  sig_ind <- setdiff(1:nrow(sp),non_sig_ind)
  keep.warp_y <- .subset(keep.warp_sig,1:((iters-burn)/thin),cell,1:2)
  
 
  
  if(!is.null(Yf)){Yfvec <- as.vector(Yf); MSE <- MSE/((iters - burn)/thin);  MAD <- MAD/((iters - burn)/thin); Cov <- mean(sapply(1:length(Yfvec),function(j){if(!is.na(Yfvec[j])){ttm <- sapply(forekeepmean,.subset2,j); tts <- sapply(forekeepsd,.subset2,j); q <- quantile(ttm + tts*rnorm((iters-burn)/thin),c(0.025,0.975)); return(as.numeric(q[1] <= Yfvec[j] & Yfvec[j] <= q[2]));}else{return(NA)}}),na.rm = T); CRPS <- mean(sapply(1:length(Yfvec),function(j){if(!is.na(Yfvec[j])){forec <- sapply(forekeepmean,.subset2,j); mean(abs(forec - Yfvec[j])) - 0.5*mean(sapply(1:((iters-burn)/thin),function(kk){mean(abs(forec[-kk] - forec[kk]))}))}else{return(NA)}}),na.rm = T)}else{MSE <- Cov <- MAD <- CRPS <- NA}      
  ll <- list("Sigma_a" = keep.sigma_a,"a" = keep.a,"sigma_y" = keep.sigma_y, "forecast" = Forecast, "coverage" = Cov, "MSE" = MSE, "warphat" = warphat, "rho_a" = keep.rho_a, "rho_x" = keep.rho_x, "rho_a0" = keep.rho_a0, "sigma_0" = keep.sigma_0, "b0" = keep.b0, "sig_ind" = sig_ind, "non_sig_ind" = non_sig_ind, "MAD" = MAD, "CRPS" = CRPS, "forekeepmean" = forekeepmean, "warp_y" = keep.warp_y)
  return(ll)
  
}