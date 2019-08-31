rm(list=ls())

## Reqiured packages for the exercise
library(fields)
library(Rcpp)
library(RcppArmadillo)
library(spdep)
library(splines)
library(ggplot2)
library(viridis)

## Data Generation----

## Importing neccesary functions for data generation

source("Data Generation/DataGen.R")

## Importing pre-defined warping functions to use to generate data

source("Data Generation/warp_fn.R")


## Generating X and S

##Choose gridsize
xsize <- 75
ysize <- 75
nhrs <- 24

## Default 
XSobj <- genXS(num=nhrs,numx=xsize,numy=ysize) ## generates from an exponential covariance on a 50x50 grid on [0,1]^2 for 24 hours

## Other versions, uncomment and change inputs as you like
# XSobj <- genXS(num=nhrs,numx=xsize,numy=ysize,
#                xlim=c(-1,1),ylim=c(-1,1)) ## changing the grid boundaries

# XSobj <- genXS(num=nhrs,numx=xsize,numy=ysize,
#                mn=1,
#                sdev=1.5,
#                range=0.05,nu=1.5) ## changing mean and covariance of the process


X <- XSobj$X
S <- XSobj$S

## Setting the number of monitor stations to be 1/25-th of total points
## Getting the locations randomly
cells <- sample(1:(xsize*ysize),floor(xsize/5)*floor(ysize/5),replace=F)

## Coordinates for the monitor locations
s <- S[cells,]



## Using pre-defined warp functions
## Returns warped location matrix

## Using Identity warp
# ws <- S

## Using Translation warp
ws <- warp_trans(S,0.12)

## Another use of translation warp, uncomment to use
# ws <- warp_trans(S,c(0.1,-0.12))

## Diffeomorphism warp, uncomment to use
# ws <- warp_diffeomorph(S,c(0.5,0.1))

## B-spline based warp, uncomment to use
# A_vec <- rnorm(2*floor(xsize/5)*floor(ysize/5))
# ws <- warp_bs(S,floor(xsize/5),floor(ysize/5),A_vec) ## Keeping the number of basis functions same as number of monitors


## Generating Y and warped and smoothed X

## Default definition, uncomment to change inputs and use
# data_obj <- gen_data(nrep=1,
#                      X,S,
#                      cells=1:nrow(S),
#                      L=10,
#                      ws = S+0.1,
#                      pmiss=0,
#                      beta_0_y = 1.5, beta_1_y = 0.25, sig_y = 1)


data_obj <- gen_data(X=X,S=S,
                     cells=cells,
                     ws=ws)

## Adding 5% missingness, uncomment to use
# data_obj <- gen_data(X=X,S=S,
#                      cells=cells,
#                      ws=ws,
#                      pmiss=0.05)

Y <- data_obj$Y_arr[,,1]
X_true <- data_obj$X_true

rm("XSobj")

## Data Handling----

## Divinding in Training and testing parts
train_ind <- 1:ceiling(0.75*nhrs)
test_ind <- (ceiling(0.75*nhrs)+1):nhrs

Y_in <- Y[,train_ind]
X_in <- X[,train_ind]

Yf <- Y[,test_ind]
Xf <- X[,test_ind]


## Application----

## Loading the function for applying full model

source("Methods/Full_Model.R")


## Default definition, uncomment to change and use

# out <- Forecast_func_f(Y_in, X_in, s, sp, cell, Xf, Yf = NULL,
#                        iters = 1000000, burn = 600000, thin = 40, step_int = 30, ridge = FALSE,
#                        L = 15, 
#                        knots_param,
#                        rho_a = 0.95, rho_x = 0.95,
#                        a_med = 0.1,
#                        sig_0_a = 0.01, sig_0_b = 0.01,
#                        rho_a0 = 0.95,
#                        a_y = 0.01, b_y = 0.01)


out <- Forecast_func_f(Y_in=Y_in, X_in=X_in, s=s, sp=S, cell=cells, Xf=Xf, Yf=Yf,
                       iters=50000, burn=20000, thin=3, step_int=30,
                       L=15,
                       knots_param=c(floor(xsize/5),floor(ysize/5)))


## Performace metrics----
print(perform_metrics <- c("MSE"=out$MSE, "MAD"=out$MAD, "Coverage"=out$coverage, "CRPS"=out$CRPS))



## Warp Estimates----
## Red means significant
## Green means non-significant

## Thinning out to make the picture presentable
ind <- which(1:(xsize*ysize)%%43 == 0)
sig_ind <- intersect(ind,out$sig_ind)
non_sig_ind <- intersect(ind,out$non_sig_ind)

plot(S[ind,1],S[ind,2],col = "white", xlab = "longitude", ylab = "lattitude")
arrows(S[sig_ind,1], S[sig_ind,2], out$warphat[sig_ind,1], out$warphat[sig_ind,2], col = "red", angle = 30, length = 0.05)
arrows(S[non_sig_ind,1], S[non_sig_ind,2], out$warphat[non_sig_ind,1], out$warphat[non_sig_ind,2], col = "forestgreen", angle = 30, length = 0.05)


## Forecasts----

## New forecasts with monitor station data overlayed
fff <- out$forecast
ii <- 1 ## index for forecast hour
  ggplot()+
    geom_point(aes(x = S[,1], y= S[,2], col = fff[,ii]),data = as.data.frame(fff[,ii]))+
    geom_point(aes(x = s[,1], y = s[,2], col = Yf[,ii]), data = as.data.frame(Yf[,ii]), alpha = 1, size = 4)+
    scale_color_gradientn(colors = viridis(25))

  
## Diagnostics----
  
  ## All commented out, uncomment to use
  
  ## coefficients for the warping function
  # pdf("a.pdf")
  # par(mfrow=c(3,2))
  # for(kk in 1:dim(out$a)[2]){
  #   for(kkk in 1:2){
  #     plot(out$a[,kk,kkk],type = "l", main = paste("coeff (",kk,",",kkk,")",sep = ""))
  #     acf(out$a[,kk,kkk],lag.max = 40, main = paste("coeff (",kk,",",kkk,")",sep = ""))
  #   }
  # }
  # dev.off()
   
  ## coefficients for spatially varying intercept
  # pdf("b0.pdf")
  # par(mfrow=c(3,2))
  # for(inn in 1:ncol(out$b0)){
  #   plot(out$b0[,inn], type = "l", main = paste("b0[",inn,"]",sep=""))
  #   acf(out$b0[,inn], lag.max = 40, main = paste("b0[",inn,"]",sep=""))
  # }
  # dev.off()
   
  ## w(s) - s
  # pdf("warpfn.pdf")
  # par(mfrow=c(2,1))
  # for(inn in 1:dim(out$warp_y)[2]){
  #   plot(out$warp_y[,inn,1], type = "l", main = expression(paste("w"[1],"(s) - s"[1])), xlab = "iterations", ylab = "")
  #   plot(out$warp_y[,inn,2], type = "l", main = expression(paste("w"[2],"(s) - s"[2])), xlab = "iterations", ylab = "")
  # }
  # dev.off()
   
  ## others
  # pdf("others.pdf")
  # par(mfrow = c(2,2))
  # plot(out$Sigma_a,type = "l", main = "sigma_a")
  # acf(out$Sigma_a,lag.max = 40, main = "sigma_a")
  # plot(out$sigma_y,type = "l", main = "sigma_y")
  # acf(out$sigma_y,lag.max = 40, main = "sigma_y")
  # plot(out$rho_a,type = "l", main = "rho_a")
  # acf(out$rho_a,lag.max = 40, main = "rho_a")
  # plot(out$rho_x,type = "l", main = "rho_x")
  # acf(out$rho_x,lag.max = 40, main = "rho_x")
  # plot(out$rho_a0,type = "l", main = "rho_a0")
  # acf(out$rho_a0,lag.max = 40, main = "rho_a0")
  # plot(out$sigma_0,type = "l", main = "sigma_0")
  # acf(out$sigma_0,lag.max = 40, main = "sigma_0")
  # dev.off()
  