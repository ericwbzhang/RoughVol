# SVI quartic root finder

# The SVI formula
svi <- function(sviparams,k){
  a <- sviparams$a 
  b <- sviparams$b 
  sig <- sviparams$sig 
  rho <- sviparams$rho 
  m <- sviparams$m 
  return(a + b *(rho*(k-m)+ sqrt((k-m)*(k-m) + sig*sig))) 
}

# Function to find roots and compute crossedness of two slices
# This function takes two SVI total variance slices, the earlier slice in row 1, the later
# slice in row 2.  The columns are {a,b,sigma,rho,m}.
# The function returns the positions of the points where the lines on the total variance
# cross and the crossedness of the two slices.

sviRoots <- function(sviData){

a1 <- sviData$a[1] 
b1 <- sviData$b[1] 
s1 <- sviData$sig[1] 
r1 <- sviData$rho[1] 
m1 <- sviData$m[1] 

a2 <- sviData$a[2] 
b2 <- sviData$b[2] 
r2 <- sviData$rho[2] 
m2 <- sviData$m[2] 
s2 <- sviData$sig[2] 

# The standard form of the quartic is q4 x^4 + q3 x^3 +q2 x^2 +q1 x + q0 == 0
# Note that multiplying all the coefficients qi by a factor should have no effect on the roots.

q2 <- 1000000 * -2 * (-3 * b1 ^ 4 * m1 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 + 4 * b1 ^ 2 * b2 ^ 2 * m1 * m2 + 
            b1 ^ 2 * b2 ^ 2 * m2 ^ 2 - 3 * b2 ^ 4 * m2 ^ 2 + 6 * b1 ^ 4 * m1 ^ 2 * r1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r1 ^ 2 + 4 * b1 ^ 2 * b2 ^ 2 * m1 * m2 * r1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r1 ^ 2 - 3 * b1 ^ 4 * m1 ^ 2 * r1 ^ 4 - 6 * b1 ^ 3 * b2 * m1 ^ 2 * r1 * r2 - 
            6 * b1 ^ 3 * b2 * m1 * m2 * r1 * r2 - 6 * b1 * b2 ^ 3 * m1 * m2 * r1 * r2 - 
            6 * b1 * b2 ^ 3 * m2 ^ 2 * r1 * r2 + 6 * b1 ^ 3 * b2 * m1 ^ 2 * r1 ^ 3 * r2 + 
            6 * b1 ^ 3 * b2 * m1 * m2 * r1 ^ 3 * r2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r2 ^ 2 + 
            4 * b1 ^ 2 * b2 ^ 2 * m1 * m2 * r2 ^ 2 + b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r2 ^ 2 + 6 * b2 ^ 4 * m2 ^ 2 * r2 ^ 2 - 
            3 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r1 ^ 2 * r2 ^ 2 - 12 * b1 ^ 2 * b2 ^ 2 * m1 * m2 * r1 ^ 2 * r2 ^ 2 - 
            3 * b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r1 ^ 2 * r2 ^ 2 + 6 * b1 * b2 ^ 3 * m1 * m2 * r1 * r2 ^ 3 + 
            6 * b1 * b2 ^ 3 * m2 ^ 2 * r1 * r2 ^ 3 - 3 * b2 ^ 4 * m2 ^ 2 * r2 ^ 4 - 
            a1 ^ 2 * (b1 ^ 2 * (-1 + 3 * r1 ^ 2) - 6 * b1 * b2 * r1 * r2 + b2 ^ 2 * (-1 + 3 * r2 ^ 2)) - 
            a2 ^ 2 * (b1 ^ 2 * (-1 + 3 * r1 ^ 2) - 6 * b1 * b2 * r1 * r2 + b2 ^ 2 * (-1 + 3 * r2 ^ 2)) - 
            2 * a2 * (3 * b1 ^ 3 * m1 * r1 * (-1 + r1 ^ 2) - b1 ^ 2 * b2 * (2 * m1 + m2) * (-1 + 
                3 * r1 ^ 2) * r2 - 3 * b2 ^ 3 * m2 * r2 * (-1 + r2 ^ 2) + b1 * b2 ^ 2 * (m1 + 2 * m2) * 
               r1 * (-1 + 3 * r2 ^ 2)) + 2 * a1 * (3 * b1 ^ 3 * m1 * r1 * (-1 + r1 ^ 2) - 
              b1 ^ 2 * b2 * (2 * m1 + m2) * (-1 + 3 * r1 ^ 2) * r2 - 3 * b2 ^ 3 * m2 * r2 * (-1 + 
                r2 ^ 2) + b1 * b2 ^ 2 * (m1 + 2 * m2) * r1 * (-1 + 3 * r2 ^ 2) + 
              a2 * (b1 ^ 2 * (-1 + 3 * r1 ^ 2) - 6 * b1 * b2 * r1 * r2 + b2 ^ 2 * (-1 + 3 * r2 ^ 2))) - 
            b1 ^ 4 * s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * s1 ^ 2 + b1 ^ 4 * r1 ^ 2 * s1 ^ 2 - 2 * b1 ^ 3 * b2 * r1 * r2 * 
             s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * r2 ^ 2 * s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * s2 ^ 2 - b2 ^ 4 * s2 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * r1 ^ 2 * s2 ^ 2 - 2 * b1 * b2 ^ 3 * r1 * r2 * s2 ^ 2 + b2 ^ 4 * r2 ^ 2 * s2 ^ 2)

q4 <- 1000000 * (b1 ^ 4 * (-1 + r1 ^ 2) ^ 2 - 4 * b1 ^ 3 * b2 * r1 * (-1 + r1 ^ 2) * r2 - 
            4 * b1 * b2 ^ 3 * r1 * r2 * (-1 + r2 ^ 2) + b2 ^ 4 * (-1 + r2 ^ 2) ^ 2 + 
            2 * b1 ^ 2 * b2 ^ 2 * (-1 - r2 ^ 2 + r1 ^ 2 * (-1 + 3 * r2 ^ 2)))

q0 <- 1000000 * (a1 ^ 4 + a2 ^ 4 + b1 ^ 4 * m1 ^ 4 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 + b2 ^ 4 * m2 ^ 4 - 
            2 * b1 ^ 4 * m1 ^ 4 * r1 ^ 2 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 * r1 ^ 2 + b1 ^ 4 * m1 ^ 4 * r1 ^ 4 + 
            4 * b1 ^ 3 * b2 * m1 ^ 3 * m2 * r1 * r2 + 4 * b1 * b2 ^ 3 * m1 * m2 ^ 3 * r1 * r2 - 
            4 * b1 ^ 3 * b2 * m1 ^ 3 * m2 * r1 ^ 3 * r2 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 * r2 ^ 2 - 
            2 * b2 ^ 4 * m2 ^ 4 * r2 ^ 2 + 6 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 * r1 ^ 2 * r2 ^ 2 - 
            4 * b1 * b2 ^ 3 * m1 * m2 ^ 3 * r1 * r2 ^ 3 + b2 ^ 4 * m2 ^ 4 * r2 ^ 4 + 
            4 * a2 ^ 3 * (b1 * m1 * r1 - b2 * m2 * r2) - 4 * a1 ^ 3 * (a2 + b1 * m1 * r1 - 
              b2 * m2 * r2) + 2 * b1 ^ 4 * m1 ^ 2 * s1 ^ 2 - 2 * b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * s1 ^ 2 - 
            2 * b1 ^ 4 * m1 ^ 2 * r1 ^ 2 * s1 ^ 2 + 4 * b1 ^ 3 * b2 * m1 * m2 * r1 * r2 * s1 ^ 2 - 
            2 * b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r2 ^ 2 * s1 ^ 2 + b1 ^ 4 * s1 ^ 4 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * s2 ^ 2 + 
            2 * b2 ^ 4 * m2 ^ 2 * s2 ^ 2 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r1 ^ 2 * s2 ^ 2 + 
            4 * b1 * b2 ^ 3 * m1 * m2 * r1 * r2 * s2 ^ 2 - 2 * b2 ^ 4 * m2 ^ 2 * r2 ^ 2 * s2 ^ 2 - 
            2 * b1 ^ 2 * b2 ^ 2 * s1 ^ 2 * s2 ^ 2 + b2 ^ 4 * s2 ^ 4 + 4 * a2 * (b1 * m1 * r1 - b2 * m2 * r2) * 
             (-2 * b1 * b2 * m1 * m2 * r1 * r2 + b1 ^ 2 * (m1 ^ 2 * (-1 + r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + r2 ^ 2) - s2 ^ 2)) - 4 * a1 * (a2 + b1 * m1 * r1 - 
              b2 * m2 * r2) * (a2 ^ 2 - 2 * b1 * b2 * m1 * m2 * r1 * r2 + 2 * a2 * (b1 * m1 * r1 - 
                b2 * m2 * r2) + b1 ^ 2 * (m1 ^ 2 * (-1 + r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + r2 ^ 2) - s2 ^ 2)) + 2 * a2 ^ 2 * 
             (-6 * b1 * b2 * m1 * m2 * r1 * r2 + b1 ^ 2 * (m1 ^ 2 * (-1 + 3 * r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + 3 * r2 ^ 2) - s2 ^ 2)) + 
            2 * a1 ^ 2 * (3 * a2 ^ 2 - 6 * b1 * b2 * m1 * m2 * r1 * r2 + 6 * a2 * (b1 * m1 * r1 - 
                b2 * m2 * r2) + b1 ^ 2 * (m1 ^ 2 * (-1 + 3 * r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + 3 * r2 ^ 2) - s2 ^ 2)))


q3 <- 1000000 * -4 * (b1 ^ 4 * m1 * (-1 + r1 ^ 2) ^ 2 - b1 ^ 3 * r1 * (-1 + r1 ^ 2) * 
             (a1 - a2 + b2 * (3 * m1 + m2) * r2) + b2 ^ 3 * (-1 + r2 ^ 2) * 
             ((a1 - a2) * r2 + b2 * m2 * (-1 + r2 ^ 2)) + b1 * b2 ^ 2 * r1 * 
             (a1 - 3 * a1 * r2 ^ 2 - b2 * (m1 + 3 * m2) * r2 * (-1 + r2 ^ 2) + 
              a2 * (-1 + 3 * r2 ^ 2)) + b1 ^ 2 * b2 * ((a1 - a2) * (-1 + 3 * r1 ^ 2) * r2 + 
              b2 * (m1 + m2) * (-1 - r2 ^ 2 + r1 ^ 2 * (-1 + 3 * r2 ^ 2))))

q1 <- 1000000 * 4 * (-(b1 ^ 4 * m1 ^ 3) + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 + b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 - 
            b2 ^ 4 * m2 ^ 3 + 2 * b1 ^ 4 * m1 ^ 3 * r1 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 * r1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 * r1 ^ 2 - b1 ^ 4 * m1 ^ 3 * r1 ^ 4 - b1 ^ 3 * b2 * m1 ^ 3 * r1 * r2 - 
            3 * b1 ^ 3 * b2 * m1 ^ 2 * m2 * r1 * r2 - 3 * b1 * b2 ^ 3 * m1 * m2 ^ 2 * r1 * r2 - b1 * b2 ^ 3 * m2 ^ 3 * r1 * r2 + b1 ^ 3 * b2 * m1 ^ 3 * r1 ^ 3 * r2 + 3 * b1 ^ 3 * b2 * m1 ^ 2 * m2 * 
             r1 ^ 3 * r2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 * r2 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 * r2 ^ 2 + 
            2 * b2 ^ 4 * m2 ^ 3 * r2 ^ 2 - 3 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 * r1 ^ 2 * r2 ^ 2 - 
            3 * b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 * r1 ^ 2 * r2 ^ 2 + 3 * b1 * b2 ^ 3 * m1 * m2 ^ 2 * r1 * r2 ^ 3 + 
            b1 * b2 ^ 3 * m2 ^ 3 * r1 * r2 ^ 3 - b2 ^ 4 * m2 ^ 3 * r2 ^ 4 + a1 ^ 3 * (b1 * r1 - b2 * r2) + 
            a2 ^ 3 * (-(b1 * r1) + b2 * r2) + a2 ^ 2 * (b1 ^ 2 * (m1 - 3 * m1 * r1 ^ 2) + 3 * b1 * b2 * (m1 + m2) * r1 * r2 + b2 ^ 2 * m2 * (1 - 3 * r2 ^ 2)) + 
            a1 ^ 2 * (b1 ^ 2 * (m1 - 3 * m1 * r1 ^ 2) + 3 * b1 * r1 * (-a2 + b2 * (m1 + m2) * r2) + 
              b2 * (3 * a2 * r2 + b2 * (m2 - 3 * m2 * r2 ^ 2))) - b1 ^ 4 * m1 * s1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m2 * s1 ^ 2 + b1 ^ 4 * m1 * r1 ^ 2 * s1 ^ 2 - b1 ^ 3 * b2 * m1 * r1 * r2 * s1 ^ 2 - 
            b1 ^ 3 * b2 * m2 * r1 * r2 * s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * m2 * r2 ^ 2 * s1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m1 * s2 ^ 2 - b2 ^ 4 * m2 * s2 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 * r1 ^ 2 * s2 ^ 2 - 
            b1 * b2 ^ 3 * m1 * r1 * r2 * s2 ^ 2 - b1 * b2 ^ 3 * m2 * r1 * r2 * s2 ^ 2 + 
            b2 ^ 4 * m2 * r2 ^ 2 * s2 ^ 2 + a2 * (b1 ^ 2 * b2 * r2 * (m1 ^ 2 * (-1 + 3 * r1 ^ 2) + 
                2 * m1 * m2 * (-1 + 3 * r1 ^ 2) - s1 ^ 2) + b1 ^ 3 * r1 * (-3 * m1 ^ 2 * 
                 (-1 + r1 ^ 2) + s1 ^ 2) + b2 ^ 3 * r2 * (3 * m2 ^ 2 * (-1 + r2 ^ 2) - s2 ^ 2) + 
              b1 * b2 ^ 2 * r1 * (m1 * m2 * (2 - 6 * r2 ^ 2) + m2 ^ 2 * (1 - 3 * r2 ^ 2) + s2 ^ 2)) + 
            a1 * (3 * a2 ^ 2 * (b1 * r1 - b2 * r2) + a2 * (2 * b1 ^ 2 * m1 * (-1 + 3 * r1 ^ 2) - 
                6 * b1 * b2 * (m1 + m2) * r1 * r2 + 2 * b2 ^ 2 * m2 * (-1 + 3 * r2 ^ 2)) + 
              b1 ^ 3 * r1 * (3 * m1 ^ 2 * (-1 + r1 ^ 2) - s1 ^ 2) + b1 ^ 2 * b2 * r2 * ( 
                m1 * m2 * (2 - 6 * r1 ^ 2) + m1 ^ 2 * (1 - 3 * r1 ^ 2) + s1 ^ 2) + 
              b1 * b2 ^ 2 * r1 * (2 * m1 * m2 * (-1 + 3 * r2 ^ 2) + m2 ^ 2 * (-1 + 3 * r2 ^ 2) - 
                s2 ^ 2) + b2 ^ 3 * r2 * (-3 * m2 ^ 2 * (-1 + r2 ^ 2) + s2 ^ 2)))

term16 <- (2 * q2 ^ 3 + 27 * q3 ^ 2 * q0 - 72 * q4 * q2 * q0 - 9 * q3 * q2 * q1 + 27 * q4 * q1 ^ 2)

term21 <- (q2 ^ 2 / 4 + 3 * q4 * q0 - 3 * q3 * q1 / 4)

term1sq <- -256 * term21 ^ 3 + term16 ^ 2

term1 <- sqrt(term1sq+0*1i)  # Note use of complex arithmetic in R

term23 <- (term16 + term1) ^ (1/3)

term22 <- 3*q4*term23 

temp1 <- (4 * 2 ^ (1 / 3) * term21)
temp2 <- (3 * 2 ^ (1 / 3) * q4)
temp3 <- q3 ^ 2 / (4 * q4 ^ 2) - (2 * q2) / (3 * q4)
temp4 <- temp1/ term22 + term23/temp2 

rr <- sqrt(temp3+temp4) 

temp5 <- q3 ^ 2 / (2 * q4 ^ 2) - (4 * q2) / (3 * q4)
temp6 <- (-q3 ^ 3 / 4 + q4 * q3 * q2 - 2 * q4 ^ 2 * q1) / (q4 ^ 3)

ee <- q3 ^ 2 / (2 * q4 ^ 2) - (4 * q2) / (3 * q4) - (4 * 2 ^ (1 / 3) * term21) / term22 - term23 / (3 * 2 ^ (1 / 3) * q4) -
(-q3 ^ 3 / 4 + q4 * q3 * q2 - 2 * q4 ^ 2 * q1) / (q4 ^ 3 * rr) 

dd <- q3 ^ 2 / (2 * q4 ^ 2) - (4 * q2) / (3 * q4) - (4 * 2 ^ (1 / 3) * term21) / term22 - term23 / (3 * 2 ^ (1 / 3) * q4) + 
(-q3 ^ 3 / 4 + q4 * q3 * q2 - 2 * q4 ^ 2 * q1) / (q4 ^ 3 * rr) 

temp7 <- -q3 / (4 * q4) 

# Potential roots are given by
roots <- c(
-q3/(4*q4)+rr/2+sqrt(dd)/2,
-q3/(4*q4)+rr/2-sqrt(dd)/2,
-q3/(4*q4)-rr/2+sqrt(ee)/2,
-q3/(4*q4)-rr/2-sqrt(ee)/2) 

# Need to check these are really roots
kr <- roots*(abs(Im(roots)) < 10^(-10)) 
test <- function(k){(a1 + b1 * (r1 * (k - m1) + sqrt((k - m1) ^ 2 + s1 ^ 2))) - 
            (a2 + b2 * (r2 * (k - m2) + sqrt((k - m2) ^ 2 + s2 ^ 2)))} 
            
num <- which(abs(test(kr))<10^(-10)) # Which potential root is actually a root?

roots <- sort(Re(kr[num]))# Roots in ascending order
nRoots <- length(roots)

crossedness <- 0

if(nRoots>1){ midPoints <- (roots[1:(nRoots-1)]+roots[2:nRoots])/2 } else {midPoints<-c()}  
if(nRoots>0){
  samplePoints <- c(roots[1]-1,midPoints,roots[nRoots]+1)  # Choose some sensible sampling points
	sviShort <- svi(sviData[1,],samplePoints) 
	sviLong <- svi(sviData[2,],samplePoints) 
	crossedness <- max(c(0,max(sviShort-sviLong)))  # Maximal amount of crossing bounded below by zero
	} 

return(list(roots=roots,crossedness=crossedness))   # Function returns the sum of the discriminants for true roots 
}

# Example (Remove comments to run)
# slice1 <- c(1.8,.8,0,-.5,0)
# slice2 <- c(1,1,1,-.5,0)
# testData <- as.data.frame(rbind(slice1,slice2))
# colnames(testData)<- c("a","b","sig","rho","m")
# sviRoots(testData)

