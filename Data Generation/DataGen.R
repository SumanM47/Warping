## Function to compute the indices of the nearest grid points 
## for each location for a given set of coordinates
roundnmatch <- function(s,d1,d2,nx,ny,min_lon,min_lat){
  s1 <- (s[,1] - min_lon)/d1
  s2 <- (s[,2] - min_lat)/d2
  s1 <- ceiling(s1)
  s2 <- ceiling(s2)
  s1[s1 < 1] <- 1
  s1[s1 > nx] <- nx
  s2[s2 < 1] <- 1
  s2[s2 > ny] <- ny
  match_ind <- (s2-1)*nx + s1
}



## Function to generate warped and smoothed image of a given image X with coordinates S

## Also generates a noisy image of the warped and smoothed image 
## on (possibly) select subset of points on the grid

## Default set to translation warp and smoothing with L=10
## Default set to produce the noisy image of the whole picture

gen_data <- function(nrep=1,
                     X,S,
                     cells=1:nrow(S),
                     L=10, 
                     ws = S+0.1,
                     pmiss=0,
                     beta_0_y = 1.5, beta_1_y = 0.25, sig_y = 1){
  
  ## Regulating proprtion of missingness
  if(pmiss<0 | pmiss>=1){stop("pmiss must be between 0 and 1")}
  
  
  ## Bookkeeping
  nsy <- length(cells)
  nt <- ncol(X)
  nsx <- nrow(X)
  ntoty <- nsy*nt
  
  
  ## Components for the rounding and smoothing functions
  min_lon <- min(S[,1])
  max_lon <- max(S[,1])
  min_lat <- min(S[,2])
  max_lat <- max(S[,2])
  
  m1 <- length(unique(S[,1]))
  m2 <- length(unique(S[,2]))
  
  d1 <- (max_lon - min_lon)/m1
  d2 <- (max_lat - min_lat)/m2
  
  ## Smoothing coefficients
  b_true <- sort(1+rnorm(L),decreasing = T)
  
  
  ##smoothing
  
  w1 <- 2*pi*(seq(0,m1-1,1))/m1        
  w2 <- 2*pi*(seq(0,m2-1,1))/m2
  w1 <- matrix(w1,m1,m2,byrow=FALSE)
  w2 <- matrix(w2,m1,m2,byrow=TRUE)
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  
  Bptilde <- array(0,dim = c(m1*m2,ncol(X),L))
  
  ## Creating sub-processes
  for(nh in 1:nt){
    z <- fft(matrix(X[,nh],m1,m2),inverse=TRUE) 
    for(l in 1:L){
      W           <- dbinom(l-1,L-1,delta/(2*pi))
      x           <- fft(z*W)
      tempbas <- Re(x)/(m1*m2)
      Bptilde[,nh,l] <- as.vector(tempbas)
    }
  }
  
  smoothX <- apply(Bptilde, 2, function(z){z%*%b_true})
  
  ## Rounding the warped locations to grid points
  match_ind <- roundnmatch(ws,d1,d2,m1,m2,min_lon,min_lat)
  X_true <- .subset(smoothX,match_ind,1:nt)
  
  ## Subsetting the X image for creating noisy Y image
  X_true_sub <- .subset(X_true,cells,1:nt)
  Y_true <- array(0,c(nsy,nt,nrep))
  
  
  ## Creating scaled, translated and noisy images
  ## adding missingness
  for(i in 1:nrep){
    Y_t <- beta_0_y + beta_1_y*X_true_sub + matrix(rnorm(ntoty),nsy,nt)
    if(pmiss>0){
      Y_tv <- as.vector(Y_t)
      Y_tv[sample(1:ntoty,floor(pmiss*ntoty),replace = F)] <- NA
      Y_t <- matrix(Y_t,nsy,nt)
    }
    Y_true[,,i] <- Y_t
  }
  
  return(list("Y_arr" = Y_true, "X_true" = X_true))
}



## Function to generate X on a regular grid and S given grid sizes

## Default grid limits set to be [0,1]^2
## Defaults to exponential covariance with no nugget
## Default mean 0

## Uses grf from geoR

genXS <- function(num=1,
		 numx, numy, 
                 xlim=c(0,1),ylim=c(0,1),
                 mn=0,
                 sdev=1,
                 range=0.1,nu=0.5){
  
  ## Handling inputs
  nin=numx*numy
  lmn <- length(mn)
  if(lmn!=1 & lmn!=nin){stop("mn must be scalar or be a vector of length numx*numy")}
  
  require(fields)
  x <- seq(xlim[1],xlim[2],length.out=numx)
  y <- seq(ylim[1],ylim[2],length.out=numy)

  S <- as.matrix(expand.grid(x,y))
  X <- matrix(0,numx*numy,num)

  grid <- list(x=x,y=y)
  obj <- matern.image.cov(grid=grid,theta=range,smoothness=nu,setup=TRUE)
  for(ind in 1:num){
     X[,ind] <- mn + sdev*as.vector(sim.rf(obj))
  }
  
  return(list("X"=X,"S"=S))
}