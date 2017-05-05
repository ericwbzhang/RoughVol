rm(list=ls())

library(stinepack)

# download.file(url="http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip",
#               destfile="oxfordRvData.zip")
# unzip(zipfile="oxfordRvData.zip")
# OxfordManRVRaw.data <- read.csv("OxfordManRealizedVolatilityIndices.csv", stringsAsFactors = F)
# save(OxfordManRVRaw.data, file='OxfordManRVRaw.rData')

#load("spxOptionMetrics.rData")
# load('OxfordManRVRaw.rData')

source("BlackScholes.R")
source("optionMetricsToIvols.R")
source("sviFit0.R")
source("svi.R")
source('sviarbitrage.R')
source('sviRoots.R')
source('sviSqrtFit.R')
source("computeImpliedVols.R")
source("sviFitQuarticRoots.R")
source("sviVolSurface.R")
source('plotIvols_svibss.R')
source('hybrid_bss.R')
source('bssFit.R')
source('bss.R')
source('bssiv.R')


# download.file(url="http://mfe.baruch.cuny.edu/wp-content/uploads/2014/09/spxOptionMetrics.rData_.zip", destfile="spxOptionMetrics.rData.zip")
# unzip(zipfile="spxOptionMetrics.rData.zip")


load("spxOptionMetrics.rData")


spxData110915$strike_price <- spxData110915$strike_price/1000
spxOptData <- generateOptionMetricsIvols(spxData110915)
# spxOptData <- spxOptData[spxOptData$Texp>0.04 ,  ]


paths<-1e3*5
n<-400


# alpha<- rep(-0.43, 14)
# S0<- rep(1,14)
# eta<- rep(1.9,14)
# rho<- rep(-.9, 14)
# xi<- rep(0.235^2,14)
# 
# 
# bssmatrix<- data.frame(S0=S0, alpha= alpha, eta=eta, xi=xi, rho=rho)

res<- plotIvols_svibss(ivolData = spxOptData)
plot(log(res$expiries), log(res$atmVol))
fit.lm<- lm(log(res$atmVol[-1])~ log(res$expiries[-1]))
summary(fit.lm)
abline(fit.lm)


bssmatrix<- bssSqFit(spxOptData, paths = 5e3, n = 400)
save(bssmatrix, file= 'spx20110915bss.rData')

# svifitSqrt <- sviSqrtFit(spxOptData)
# svifitQR<- sviFitQR(ivolData = spxOptData, sviGuess = svifitSqrt ) 

load('spx20110915bss.rData')

plotIvols_svibss(ivolData = spxOptData, paths= 1e5, n= 400, bssplotscale = .6, 
                 kappa = 1, bssMatrix = bssmatrix$best  )



