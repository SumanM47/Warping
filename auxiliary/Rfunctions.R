basis <- function(s,ksx,ksy){
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
  return(list("B" = Bmat,"knotsx" = 1:ksx,"knotsy" = 1:ksy))
}

warp <- function(s,AS,A_mat){
  return(s + AS%*%A_mat)
}

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