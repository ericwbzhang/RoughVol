########################################################################################################
# The SVI Jump Wings parametrization

# Convert SVI matrix to SVI-JW
sviToJw <- function(sviMatrix,texp){

a <- sviMatrix$a;
b <- sviMatrix$b;
sig <- sviMatrix$sig;
rho <- sviMatrix$rho;
m <- sviMatrix$m;

vt <- (a+b*(-rho*m+sqrt(m^2+sig^2)))/texp; # Recall that svi parameters a and b are in w terms
bhat <- sqrt(1/(vt*texp))*b;
psit <- bhat/2*(-m/sqrt(m^2+sig^2)+rho);
pt <- bhat*(1-rho);
ct <- bhat*(1+rho);
varmint <- (a+b*abs(sig)*sqrt(1-rho^2))/texp;

tmp <- data.frame(vt,psit,pt,ct,varmint,texp);
return(tmp);
}

# Convert SVI matrix to SVI-JW
jwToSvi <- function(jwMatrix){

vt <- jwMatrix$vt;
psit <- jwMatrix$psit;
pt <- jwMatrix$pt;
ct <- jwMatrix$ct;
varmint <- jwMatrix$varmint;
texp <- jwMatrix$texp;

sqrtw <- sqrt(vt*texp);
bhat <- (pt+ct)/2;
b <- bhat*sqrtw;
rho <- 1- pt/bhat; 
bet <- (rho-2*psit/bhat);
alpha <- sign(bet)*sqrt(1/bet^2-1);

m <- (vt-varmint)*texp/b/(-rho+sign(alpha)*sqrt(1+alpha^2)-alpha*sqrt(1-rho^2));
sig <- alpha*m;
a <- varmint*texp - b*sig*sqrt(1-rho^2);

tmp <- data.frame(a,b,sig,rho,m);
return(tmp);
}
