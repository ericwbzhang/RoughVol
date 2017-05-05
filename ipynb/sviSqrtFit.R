# Fit the SVI-Sqrt parametrization
# This version gets rid of alpha and fits only slices beyond a certain cutoff expiration.
# Other slices are obtained by projection

sviSqrt <- function(sviSqrtParams,k,w0){
rho <- sviSqrtParams$rho;
eta <- sviSqrtParams$eta;
w <- w0/2*(1+rho*eta/sqrt(w0)*k+sqrt((eta/sqrt(w0)*k+rho)^2+1-rho^2));
return(w);
}

sviSqrtParams <- list(rho=-0.7,eta=1);

# Compute w0 from implied vol smiles
computeW0 <- function(ivolData){

expDates <- unique(ivolData$Texp);
nSlices <- length(expDates);
w0 <- numeric(nSlices);

for (slice in 1:nSlices){
    t <- expDates[slice];
    texp <- ivolData$Texp;
    bidVol <- bidVols[texp==t];
    askVol <- askVols[texp==t];
    pick <- (bidVol>0)&(askVol>0)&(bidVol<askVol);
    pick <- !is.na(pick);
    midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
    f <- (ivolData$Fwd[texp==t])[1];
    k <- log(ivolData$Strike[texp==t]/f)[pick]; 
    w0[slice] <- t*stinterp(k,midVar,0)$y;
}
return(w0);
};

# Now fit SVI-Sqrt
sviSqrtFit <- function(ivolData){

bidVols <- as.numeric(ivolData$Bid);
askVols <- as.numeric(ivolData$Ask);
expDates <- unique(ivolData$Texp);
nSlices <- length(expDates);
slices <- 1:nSlices;
sviMatrix <- array(dim=c(nSlices,5));
#w0 <- computeW0(ivolData);
nrows <- dim(ivolData)[1];
midV <- numeric(nrows)*NA;
kk <- numeric(nrows)*NA;
ww0 <- numeric(nrows)*NA;

# Compute w0, keeping k and midVar as we go
for (slice in 1:nSlices){
    t <- expDates[slice];
    texp <- ivolData$Texp;
    bidVol <- bidVols[texp==t];
    askVol <- askVols[texp==t];
    pick <- (bidVol>0)&(askVol>0);
    pick <- !is.na(pick)&(bidVol<askVol);
    midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
    f <- (ivolData$Fwd[texp==t])[1];
    k <- log(ivolData$Strike[texp==t]/f)[pick]; 
    w0 <- t*stinterp(k,midVar,0)$y;
    # Now put in correct place in columns
    ww0[texp==t][pick]<- w0;
    midV[texp==t][pick]<- midVar;
    kk[texp==t][pick]<- k;    
}

tcutoff <- min(0.1,max(expDates));

# Define objective function
obj <- function(sviSqrtParams){
    tmp1 <- as.list(sviSqrtParams);
    names(tmp1) <- c("rho","eta");
    sviSqrtVar <- sviSqrt(tmp1,kk,ww0)/texp;
    tmp <- sum(((midV-sviSqrtVar)^2)[texp >= tcutoff],na.rm=T);
    return(tmp);
};

sviSqrtGuess <- list(rho = - 0.7, eta = 1.0); 

fit <- optim(sviSqrtGuess,obj,method="L-BFGS-B",lower=c(-.999,-Inf),upper=c(+.999,+Inf));
res <-as.list(fit$par);

# Now convert the result to an SVI matrix
sel <- !is.na(ww0);
w0r <- unique(ww0[sel]);
rho <- rep(res$rho,nSlices);
a <- w0r/2*(1-rho^2);
gg <- res$eta/sqrt(w0r);
b <- w0r/2*gg;
m <- -rho/gg;
sig <- sqrt(1-rho^2)/gg;

tmp <- data.frame(a,b,sig,rho,m);

return(tmp);

}
