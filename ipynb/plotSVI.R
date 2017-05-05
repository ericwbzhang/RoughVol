# Function to plot SVI surface

plotSVI <- function(sviMatrix,xrange=NULL,yrange=NULL){

m <- dim(sviMatrix)[1];
clr <- rainbow(m);

if (is.null(xrange)) xrange <- c(-1,1);
x1 <- xrange[1]; x2 <- xrange[2];
topleft <- svi(sviMatrix[m,],x1);
topright <- svi(sviMatrix[m,],x2);
if (is.null(yrange)) yrange <- c(0,max(topleft,topright));

for (i in 1:m){
sviparams <- sviMatrix[i,];
curve(svi(sviparams,x),from=x1,to=x2,col=clr[i],xlab="Log-strike k",ylab=expression(w=sigma^2*T),ylim=yrange);
par(new=T)
}
par(new=F)
}
