# Function to check sviMatrix for arbitrage
# We do this by checking for real roots, slice by slice pairwise



sviCalendarArbitrageCheck <- function(sviMatrix){

n <- dim(sviMatrix)[1];
arbtot <- 0;

arbcheck <- function(i){sviRoots(sviMatrix[i:(i+1),])$crossedness}; # Recall that sviRoots returns discriminants when there is a true real root.
tmp <- sum(sapply(1:(n-1),arbcheck));

return(tmp);
}

sviSliceArbitrageCheck <- function(sviMatrix){

sviMatrix <- as.data.frame(sviMatrix);
n <- dim(sviMatrix)[1];

g <- function(sviparams,k){

    a <- sviparams$a;
    b <- sviparams$b;
    sig <- sviparams$sig;
    rho <- sviparams$rho;
    m <- sviparams$m;
    
    discr <- sqrt((k-m)*(k-m) + sig*sig);
    w <- a + b *(rho*(k-m)+ discr);
    dw <- b*rho + b *(k-m)/discr;
    d2w <- b*sig^2/(discr*discr*discr);
    
    return(1 - k*dw/w + dw*dw/4*(-1/w+k*k/(w*w)-4) +d2w/2);
}

arbitrageableSlices <- numeric(0);

for (slice in 1:n){
    res <- optimize(function(x){g(sviMatrix[slice,],x)},interval=c(-2,2));
    if(res$objective < 0) {arbitrageableSlices <- c(arbitrageableSlices,slice)};
}

return(arbitrageableSlices);

}


sviSliceArbitrageCompute <- function(sviMatrix){

sviMatrix <- as.data.frame(sviMatrix);
n <- dim(sviMatrix)[1];

g <- function(sviparams,k){

    a <- sviparams$a;
    b <- sviparams$b;
    sig <- sviparams$sig;
    rho <- sviparams$rho;
    m <- sviparams$m;
    
    discr <- sqrt((k-m)*(k-m) + sig*sig);
    w <- a + b *(rho*(k-m)+ discr);
    dw <- b*rho + b *(k-m)/discr;
    d2w <- b*sig^2/(discr*discr*discr);
    
    return(1 - k*dw/w + dw*dw/4*(-1/w+k*k/(w*w)-4) +d2w/2);
}

sliceArb <- NULL;

for (slice in 1:n){
    res <- optimize(function(x){g(sviMatrix[slice,],x)},interval=c(0,2));
    kArb <- res$minimum
    gArb <- min(res$objective,0)
    sliceArb <- rbind(sliceArb,c(slice,kArb,gArb))
}

sliceArb <- as.data.frame(sliceArb)
names(sliceArb) <- c("slice","kMin","gArb")


return(as.data.frame(sliceArb));

}