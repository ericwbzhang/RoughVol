# Fit SVI independently to each slice

sviFit <- function(ivolData){

bidVols <- as.numeric(ivolData$Bid);
askVols <- as.numeric(ivolData$Ask);
expDates <- unique(ivolData$Texp);
nSlices <- length(expDates);
slices <- 1:nSlices;
sviMatrix <- array(dim=c(nSlices,5));

######################################
for (slice in slices){
t <- expDates[slice];
texp <- ivolData$Texp;
bidVol <- bidVols[texp==t];
askVol <- askVols[texp==t];
pick <- (bidVol>0)&(askVol>0);
pick <- !is.na(pick);
midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
f <- (ivolData$Fwd[texp==t])[1];
k <- log(ivolData$Strike[texp==t]/f)[pick]; # Plot vs log-strike

sviGuess <- list(a = mean(midVar,na.rm=T), b = 0.1, sig = 0.1, rho = - 0.7, m = 0); 

# Define objective function
obj <- function(sviparams){
    tmp1 <- as.list(sviparams);
    names(tmp1) <- c("a","b","sig","rho","m")
    sviVar <- svi(tmp1,k);
    minVar <- tmp1$a+tmp1$b*tmp1$sig*sqrt(abs(1-tmp1$rho^2));
    negVarPenalty <- min(100,exp(-1/minVar));
    tmp <- sum((midVar-sviVar)^2,na.rm=T)+negVarPenalty;
    return(tmp*10000);
};

fit <- optim(sviGuess,obj);
if(abs(fit$par["rho"])>.999){
fit <- optim(fit$par,obj,method="L-BFGS-B",lower=c(-10,0,0.00000001,-.999,-10),upper=c(+10,100,100,+.999,+10));}
sviMatrix[slice,] <- fit$par*c(t,t,1,1,1);

                    }# end of for loop
######################################
sviMatrix <- as.data.frame(sviMatrix);
colnames(sviMatrix) <- c("a","b","sig","rho","m");
return(sviMatrix);
}# end of function
