# setwd("/Users/aliasgarhaji/ali/study/semester7/MTH9903/code/mth9903/a")

source("BlackScholes.R");

# This module takes in option metrics data from WRDS and returns the Ivols in the standard format.

getDays <- function(s,e) {
  start <- as.Date(as.character(s),format="%Y%m%d");
  expiry <- as.Date(as.character(e),format="%Y%m%d");
  return( difftime(expiry,start,units="days") );
}

generateOptionMetricsIvols <- function(spxData){
  
  spxData <- data.frame(days=getDays(spxData$date, spxData$exdate) , spxData );
  out2 <- NULL;
  
  for (numdays in sort(unique(spxData$days))) { #One expiration at a time
    
    vdc <- spxData[(spxData$days==numdays)&(spxData$cp_flag=="C"),];# Select calls for just this expiration
    vdp <- spxData[(spxData$days == numdays)&(spxData$cp_flag=="P"),];# Select puts for just this expiration
    expiration <- unique(spxData[(spxData$days == numdays),]$exdate);
    
    ##get the put , call strikes and then the unique of the two.
    callStrikes <- sort(unique(vdc$strike_price));
    putStrikes <- sort(unique(vdp$strike_price));
    strikes <- callStrikes[callStrikes%in%putStrikes];
    
    ##setup for calibration.
    nK <- length(strikes);
    vols <- numeric(nK);
    imid <- numeric(nK);
    ca <- numeric(nK);
    cb <- numeric(nK);
    pa <- numeric(nK);
    pb <- numeric(nK);
    
    if((nK>=6)&(numdays>0)){
      cbb <- vdc$best_bid;
      pbb <- vdp$best_bid;
      cba <- vdc$best_offer;
      pba <- vdp$best_offer;
      
      for (i in 1:nK){
        
        k <- strikes[i];
        cb[i] <- mean(cbb[vdc$strike_price ==k]);
        pb[i] <- mean(pbb[vdp$strike_price ==k]);
        ca[i] <- mean(cba[vdc$strike_price ==k]);
        pa[i] <- mean(pba[vdp$strike_price ==k]);
        
        ## this is so that we can run put-call parity.
        ## C - P = PV * ( F - K );
        ibid <- cb[i]-pb[i];
        iask <- ca[i]-pa[i];
        imid[i] <- (ibid+iask)/2;
        
      }
      
      pvGuess <- 1;
      fGuess <- mean(imid+strikes);
      nearTheMoneyStrikes <- strikes[order(abs(imid))][1:6];
      
      include <- (strikes%in%nearTheMoneyStrikes); #Optimize only on near-the-money strikes
      obj <- function(params){
        f <- params[1]; pv <- params[2];
        ifit <- pv*(f-strikes); # This is a vector length nk
        errmid <- (ifit-imid)*include;
        return(sum(errmid^2));
      }
      fit <- optim(c(fGuess,pvGuess),obj,method="L-BFGS-B", lower=c(min(strikes),0.5),upper=c(max(strikes),2));
      ffit <- fit$par[1];
      pvfit <- fit$par[2];
      # Get implieds for OTM options
      texp <- numdays/365.25;
      ivolcbid <- BSImpliedVolCall(ffit, strikes, texp, 0,cb/pvfit);
      ivolcask <- BSImpliedVolCall(ffit, strikes, texp, 0,ca/pvfit);
      ivolpbid <- BSImpliedVolPut(ffit, strikes, texp, 0,pb/pvfit);
      ivolpask <- BSImpliedVolPut(ffit, strikes, texp, 0,pa/pvfit);
      
      ivolbid <- ivolcbid*(strikes>ffit)+ivolpbid*(strikes<=ffit); # This version outputs OTM vols
      ivolask <- ivolcask*(strikes>ffit)+ivolpask*(strikes<=ffit); # This version outputs OTM vols
      
      callBid <- BSFormula(ffit, strikes, texp, 0, ivolbid);
      callAsk <- BSFormula(ffit, strikes, texp, 0, ivolask);
      exclude <- (cb==0)|(pb==0)
      callMid <- (callBid+callAsk)/2;
      
      out <- data.frame(expiration,texp,strikes,ivolbid,ivolask,ffit,callMid);
      out$callMid[exclude]<-NA;
      out$ivolbid[exclude]<-NA;
      colnames(out) <- c("Expiry","Texp","Strike","Bid","Ask","Fwd","CallMid");
      out2 <- rbind(out2,out);
    } # end of if(nk>6)

  } # end biggest for.
  out2$Bid[(out2$Bid<10^-8)]<-NA;
  out3 <- out2[order(out2$Texp,out2$Strike),];# Sort columns for output
  return(out3);
  
} # End of function
