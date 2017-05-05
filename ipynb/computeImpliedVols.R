# This module takes in Rich Holowczak's option quote data and returns implied volatilities in a standard format.

mthCodes <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
"T", "U", "V", "W", "X");

findMonth <- function(mthCode){(which(mthCodes==mthCode)-1)%%12+1;}
callPut <- function(mthCode){((which(mthCodes==mthCode))<=12);}

generateIvols <- function(spxData){

  startDate <- as.Date(unique(spxData$Date),format="%m/%d/%y");
  year <- spxData$Expiration_Year+2000;
	mthCode <-spxData$Expiration_Month_Code;
	month <- sapply(mthCode,findMonth);
	day <- spxData$Expiration_Day_of_Month-1;	
	expiry <- as.Date(paste(year,month,day,sep="/"));
	days <- difftime(expiry,startDate);
	t <- as.numeric(days)/365.25; # Better code would compute business time
	cp <- sapply(mthCode,callPut); # Call or put?
	out2 <- NULL;

for (numdays in sort(unique(days))){ #One expiration at a time
		
		vdc <- spxData[(days==numdays)&cp,];# Select just this expiration
		vdp <- spxData[(days == numdays)&(!cp),];# Select just this expiration
		expiration <- unique(expiry[(days == numdays)]);
		callStrikes <- sort(unique(vdc$Strike_Price));
		putStrikes <- sort(unique(vdp$Strike_Price));
		strikes <- callStrikes[callStrikes%in%putStrikes];
		nK <- length(strikes);
		vols <- numeric(nK);
		imid <- numeric(nK);
		ca <- numeric(nK);
		cb <- numeric(nK);
		pa <- numeric(nK);
		pb <- numeric(nK);
		
		if((nK>=6)&(numdays>0)){
    			cbb <- vdc$Option_Bid_Price;
    			pbb <- vdp$Option_Bid_Price;
    			cba <- vdc$Option_Offer_Price;
    			pba <- vdp$Option_Offer_Price;
    			
    			for (i in 1:nK){
    				
    				k <- strikes[i];
    				cb[i] <- mean(cbb[vdc$Strike_Price ==k]);
    				pb[i] <- mean(pbb[vdp$Strike_Price ==k]);
    				ca[i] <- mean(cba[vdc$Strike_Price ==k]);
    				pa[i] <- mean(pba[vdp$Strike_Price ==k]);
    				
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
		} # end of if()
	} # End of for{}
out2$Bid[(out2$Bid<10^-8)]<-NA;
out3 <- out2[order(out2$Texp,out2$Strike),];# Sort columns for output
return(out3);
} # End of function

