library(stinepack)


plotIvols_svibss<- function (ivolData, paths=1e4, n=100 ,
                             bssMatrix=NULL, bssplotscale=0.5,  kappa=1,
                             sviMatrix = NULL, slices = NULL) 
{
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  colnum <- sqrt(nSlices * 2)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  atmVol <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  
  i<-1
  
  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(bidVol)
    kmin <- min(k[include])
    kmax <- max(k[include])
    ybottom <- 0.8 * min(bidVol[include])
    ytop <- 1.2 * max(askVol[include], na.rm = T)
    xrange <- c(kmin, kmax)
    yrange <- c(0.1, ytop)
    plot(k, bidVol, col = "red", pch = 18, cex = 0.5, xlim = xrange, 
         ylim = yrange, main = paste("T =", format(t, digits = 2, 
                                                   nsmall = 2)), xlab = "Log-Strike", ylab = "Implied Vol.")
    par(new = T)
    plot(k, askVol, col = "blue", pch = 18, cex = 0.5, xlim = xrange, 
         ylim = yrange, main = NA, xlab = NA, ylab = NA)
    
    
    
    ### plot svi fit
    if ((!is.null(sviMatrix))) {
      vol <- function(k) {
        sqrt(svi(sviMatrix[slice, ], k)/t)
      }
      par(new = T)
      curve(vol(x), from = kmin, to = kmax, col = "orange", 
            lwd = 2, add = T)
    }
    
    
    ### plot bss fit
    if (! is.null(bssMatrix)){
      N<-n
      if(t<0.1) {N<- max(N, 400)}
      
      finalPrices<- hybridScheme(bssMatrix[i,])(paths,N,kappa, t)
      # print(bssMatrix[slice,])
      bsskmin<- kmin
      bsskmax<- kmax
      if (!is.na(bssplotscale)) {
        bsskmin<-bsskmin*bssplotscale
        bsskmax<-bsskmax*bssplotscale
      }
      curve(ImpliedVol(x, finalPrices, t), from = bsskmin, to = bsskmax, col='dark green',lwd=3,  add=TRUE)
    }
    
    
    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVol[i] <- volInterp(0)
    atmSkew[i] <- (volInterp(0.01) - volInterp(-0.01))/0.02
    
    i<- i+1
  }
  par(mfrow = c(1, 1), mex = 1)
  par(new = F)
  return(list(expiries = expDates[slices], atmVol = atmVol, atmSkew = atmSkew))
}



plotIvols_BSS<- function (ivolData, paths, n ,bssMatrix, kappa,  sviMatrix = NULL, slices = NULL) 
{
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  colnum <- sqrt(nSlices * 2)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  atmVol <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(bidVol)
    kmin <- min(k[include])
    kmax <- max(k[include])
    ybottom <- 0.8 * min(bidVol[include])
    ytop <- 1.2 * max(askVol[include], na.rm = T)
    xrange <- c(kmin, kmax)
    yrange <- c(0.1, ytop)
    plot(k, bidVol, col = "red", pch = 18, cex = 0.5, xlim = xrange, 
         ylim = yrange, main = paste("T =", format(t, digits = 2, 
                                                   nsmall = 2)), xlab = "Log-Strike", ylab = "Implied Vol.")
    par(new = T)
    plot(k, askVol, col = "blue", pch = 18, cex = 0.5, xlim = xrange, 
         ylim = yrange, main = NA, xlab = NA, ylab = NA)
    
    finalPrices<- hybridScheme(bssMatrix[slice,])(paths,n,kappa, t)
    curve(ImpliedVol(x, finalPrices), from = kmin, to = kmax, col='dark green', add=TRUE)
    
    
    if ((!is.null(sviMatrix))) {
      vol <- function(k) {
        sqrt(svi(sviMatrix[slice, ], k)/t)
      }
      par(new = T)
      curve(vol(x), from = kmin, to = kmax, col = "orange", 
            lwd = 2, add = T)
    }
    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVol[slice] <- volInterp(0)
    atmSkew[slice] <- (volInterp(0.01) - volInterp(-0.01))/0.02
  }
  par(mfrow = c(1, 1), mex = 1)
  par(new = F)
  return(list(expiries = expDates, atmVol = atmVol, atmSkew = atmSkew))
}


