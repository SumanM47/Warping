warp_bs <- function(s,ksx,ksy,A_vec){
  if(length(A_vec)!=2*ksx*ksy){stop("A_vec must have length twice ksx*ksy")}
  require(splines)
  b1 <- splines::bs(unique(s[,1]),ksx)
  b2 <- splines::bs(unique(s[,2]),ksy)
  out <- outer(b1,b2)
  Bmat <- matrix(NA,nrow(s),ksx*ksy)
  for(i in 1:ksy){
    for(j in 1:ksx){
      Bmat[,(i-1)*ksx + j] <- out[,j,,i]
    }
  }
  return(s + Bmat%*%matrix(A_vec,ncol=2))
}

warp_diffeomorph <- function(s,theta){
  x <- s[,1]
  y <- s[,2]
  s1 <- (x - min(x))/(max(x)-min(x))
  s2 <- (y - min(y))/(max(y)-min(y))
  theta1 <- theta[1]
  theta2 <- theta[2]
  w1 <- s1 - 2*theta1*s2*sin(s1)*sin(s2)*(cos(pi*s1) + 1)*(cos(pi*s2) + 1)
  w2 <- s2 - 2*theta2*s1*sin(s1)*sin(s2)*cos(s1*pi/2)*cos(pi*s2*3/2)
  return(cbind((max(x)-min(x))*w1 + min(x),(max(y) - min(y))*w2 + min(y)))
}

warp_trans <- function(s,tra){
  ltra <- length(tra)
  if(ltra<1 | ltra>2){stop("tra must be a scalar or vector of length 2")}
  if(ltra==1){tra = rep(tra,2)}
  return(t(t(s)+tra))
}