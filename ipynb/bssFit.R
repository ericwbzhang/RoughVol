# Fit BSS independently to each slice

bssSqFit<- function (ivolData, paths=5e3, n=400 , kappa=1, kscale= .4,
                     guess=list(S0=1, eta=2 , xi=0.24^2, alpha=-0.44, rho=-.9 ) ,
                     slices= NULL){
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
  
  S0<- rep(1, nSlices)
  alpha<- numeric(nSlices)
  eta<-   numeric(nSlices)
  rho<-   numeric(nSlices)
  xi<-    numeric(nSlices)
  
  TEXP<- expDates[slices]
  alpha1<- numeric(nSlices)
  eta1<-   numeric(nSlices)
  rho1<-   numeric(nSlices)
  xi1<-    numeric(nSlices)

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

    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVol[i] <- volInterp(0)
    N<-n;
    if (t<0.1) { N<- max(N,400)}
    
    xi_obj<- function (param_xi, param, k_in, midVol_in){
      param$xi<-param_xi
      bssVol<- bssiv(bssparams = param, paths = paths, n = N, kappa = kappa, k = k_in, t = t )
      res<- sum((bssVol- midVol_in)^2) *1e4
      
      return (res)
    }
    
    eta_obj<- function (param_eta, param, k_in, midVol_in){
      param$eta<- param_eta
      bssVol<- bssiv(bssparams = param, paths = paths, n = N, kappa = kappa, k = k_in, t = t )
      res<- sum((bssVol- midVol_in)^2) *1e4
      
      return (res)
    }
    
    rho_obj<-function(param_rho, param, k_in, midVol_in){
      param$rho<- param_rho
      bssVol<- bssiv(bssparams = param, paths = paths, n = N, kappa = kappa, k = k_in, t = t )
      res<- sum((bssVol- midVol_in)^2) *1e4
      
      return (res)
    }
    
    bssparam<- guess
    
    k_low<- kmin*kscale
    k_high<- kmax*kscale
    idx<- which(kIn> k_low & kIn<k_high)
    midVol_in<- volIn[idx]
    k_in<- kIn[idx]
    
    ## First, optimize xi
    set.seed(9081)
    lower<-(atmVol[i]-0.03)^2
    upper<- (atmVol[i]+0.03)^2
    fit_xi<- optimize(f = xi_obj, param= bssparam, k_in= k_in, midVol_in= midVol_in, 
                      lower = lower, upper = upper, maximum = F)
    print(fit_xi)
    if(fit_xi$minimum<= lower + 0.0001 || fit_xi$minimum >= upper-0.0001){
      message("Warning! The boundary of Xi is hit.")
    }
    bssparam$xi<- fit_xi$minimum
    fit_xi$param<- bssparam
    
    ## Second, optimize eta, with optimized xi
    set.seed(9081)
    lower<- 1.8
    upper<- 2.5
    fit_eta<- optimize(f = eta_obj, param= bssparam,k_in=k_in, midVol_in=midVol_in, 
                       lower = lower, upper = upper, maximum = F)
    print(fit_eta)
    if(fit_eta$minimum<=lower+0.0001 || fit_eta$minimum >= upper-0.0001){
      message('Warning! The boundary of Eta is hit.')
    }
    bssparam$eta<- fit_eta$minimum
    fit_eta$param<- bssparam
    
    ## Last, optimize rho, with optimized xi and eta
    set.seed(9081)
    lower<- -0.95
    upper<--.6
    
    fit_rho<- optimize(f= rho_obj, param= bssparam, k_in= k_in, midVol_in=midVol_in, 
                       lower = lower, upper= upper, maximum = F)
    print(fit_rho)
    if(fit_rho$minimum< lower+0.0001 || fit_rho$minimum > upper-0.0001){
      message('Warning! The boundary of Rho is hit.')
    }
    bssparam$rho<- fit_rho$minimum
    fit_rho$param<- bssparam
    
    if (!(fit_xi$objective>= fit_eta$objective/3 && fit_eta$objective>= fit_rho$objective/3)){
      message ('The optimization kernal shows suspicious contradicts with t= ', t)
    }
    
    alpha[i]<- bssparam$alpha
    rho[i]  <- bssparam$rho
    eta[i]  <- bssparam$eta
    xi[i]   <- bssparam$xi
    
    bestfit<- fit_xi
    if (fit_eta$object< bestfit$object) {bestfit<- fit_eta}
    if (fit_rho$object< bestfit$object) {bestfit<- fit_rho}
    
    alpha1[i]<- bestfit$param$alpha
    rho1[i]<-   bestfit$param$rho
    eta1[i]<-   bestfit$param$eta
    xi1[i]<-    bestfit$param$xi
    
    i<- i+1
  }

  res<- list(normal= data.frame(texp= TEXP, S0=S0, xi= xi, alpha= alpha, eta= eta, rho=rho),
             best= data.frame(texp= TEXP, S0= S0, xi= xi1, alpha= alpha1, eta= eta1, rho=rho1))
  return (res)
  
}
