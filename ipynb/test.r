rm(list=ls())
library(MASS)
library(parallel)

bstar <- function(k, alpha){
  result <- ((k^(alpha+1) - (k-1)^(alpha+1)) / (alpha+1)) ^ (1/alpha)
  return(result)
}

BSImpliedVolCall <- function(S0, K, T, r, C)
{
  nK <- length(K);
  sigmaL <- rep(1e-10,nK);
  CL <- BSFormula(S0, K, T, r, sigmaL);
  sigmaH <- rep(10,nK);
  CH <- BSFormula(S0, K, T, r, sigmaH);
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    CM <- BSFormula(S0, K, T, r, sigma);
    CL <- CL + (CM < C)*(CM-CL);
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL);
    CH <- CH + (CM >= C)*(CM-CH);
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH);
  }
  return(sigma);
  
}

MyImpliedVolCall <- function(S0, K, T, r, C){
  sigmaL <- 1e-3
  sigmaH <- 10
  while(sigmaH-sigmaL > 1e-10){
    sigma <- (sigmaL + sigmaH)/2
    if(BSFormula(S0, K, T, r, sigma)==C) {break}
    else if(BSFormula(S0, K, T, r, sigma)>C) {sigmaH <- sigma}
    else{
      sigmaL <- sigma
    }
  }
  return(sigma)
}

BSFormula <- function(S0, K, T, r, sigma)
{
  x <- log(S0/K)+r*T;
  sig <- sigma*sqrt(T);
  d1 <- x/sig+sig/2;
  d2 <- d1 - sig;
  pv <- exp(-r*T);
  return( S0*pnorm(d1) - pv*K*pnorm(d2));
}

hybridScheme <- function(params){
  
  S0 <- params$S0
  xi <- params$xi
  eta <- params$eta
  alpha <- params$alpha
  rho <- params$rho
  
  #n = 100
  #kappa = 1
  covMatrix <- function(n, kappa){
    ### generate the covariance matrix
    ### @param n: the distance between two points is 1/n
    ### @param kappa: how many power low terms to include around zero
    
    sigma <- matrix(0, nrow=kappa+1, ncol=kappa+1)
    sigma[1,1] <- 1/n
    if(kappa==0){
      return (sigma)
    }
    else{
      for(j in 2:(kappa+1)){
        sigma[1, j] <- ((j-1)^(alpha+1)-(j-2)^(alpha+1)) / ((alpha+1)*n^(alpha+1))
        sigma[j, 1] <- ((j-1)^(alpha+1)-(j-2)^(alpha+1)) / ((alpha+1)*n^(alpha+1))
      }
      for(j in 2:(kappa+1)){
        for(k in 2:(kappa+1)){
          if(j==k){
            sigma[j,k] <- ((j-1)^(2*alpha+1)-(j-2)^(2*alpha+1)) / ((2*alpha+1)*n^(2*alpha+1))
          }
          
        }
      }
      return(sigma)
    }
  }
  
  Simulation <- function(n, kappa, T, W, Z, Gamma, tseq){
    #print(W)
    if(kappa==0){Y2 <- convolve(Gamma, rev(W), type="open")[1:floor(n*T)]}
    else{Y2 <- convolve(Gamma, rev(W[,1]), type="open")[1:floor(n*T)]}
    
    Y1 <- rep(0, floor(n*T))
    
    for(i in 1:floor(n*T)){
      Y1[i] <-0
      if (kappa!=0 ){
        for (k in 1:min(i,kappa)){
          Y1[i] = Y1[i] + W[i+1-k,k+1]
        }
      }
    }
    Y <- Y1+Y2 ## The simulated series of main interst
    v <- xi*exp(eta*sqrt(2*alpha+1)*Y - tseq)
    v <- c(xi, v[1:length(v)-1])
    S <- S0 * exp(sum(v^0.5*Z) - 1/2*sum(v)/n)
    return(S)
  }
  
  MC <- function(N, n, kappa, T){
    Gamma <- sapply(seq(1:floor(n*T)), function(x){(bstar(x, alpha)/n)^alpha}) 
    if(kappa!=0) {Gamma[1:kappa] <- 0}
    tseq <- eta*eta/2*sapply(seq(1:floor(n*T)),function(x){(x/n)^(2*alpha+1)})
    steps <- floor(n*T)
    W <- mvrnorm(steps*N, mu=rep(0, kappa+1), Sigma=covMatrix(n, kappa))
    Wperp <- rnorm(steps*N, sd=sqrt(1/n))
    Z <- rho * W[,1] + sqrt(1-rho*rho)*Wperp
    return(sapply(seq(1:N), function(loopNum){Simulation(
      n, kappa, T, W[(1+(loopNum-1)*steps):(loopNum*steps),],Z[(1+(loopNum-1)*steps):(loopNum*steps)], Gamma, tseq)}))
    #return(unlist(mclapply(seq(1:N), function(loopNum){Simulation(
    #  n, kappa, T, W[(1+(loopNum-1)*steps):(loopNum*steps),],Z[(1+(loopNum-1)*steps):(loopNum*steps)], Gamma, tseq)}, mc.cores=4)))
  }
  return(MC)
}

impvol <- function(k, st, T){
  payoff <- (st > exp(k)) * (st - exp(k))
  return(BSImpliedVolCall(1, exp(k), T, 0, mean(payoff)))
}




params <- list(S0=1, xi=0.235^2, eta=1.9, alpha=-0.43, rho=-0.9)
T <- 1
n <- 10 # step length is 1/n
paths <- 1e5
kappa <-1
set.seed(9081)

Rprof()
finalPrices <- hybridScheme(params)(paths, n, kappa, T)
Rprof(NULL)
summaryRprof()

summary(finalPrices)



vol <- function(k, finalP){sapply(k, function(x){impvol(x, finalP, T)})}
payoff <- function(k, finalP){sapply(k, function(x){mean((finalP > exp(x)) * (finalP - exp(x)))})}

k<-5
clr<- rainbow(k)
finalP<- list()

curve(vol(x, finalPrices), from=-0.35, to=0.35, col='black', xlab="Log strike", ylab="Implied Vol", ylim= c(0, 0.4) ,lwd=2)
