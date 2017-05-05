
bss<- function (bssparams, k, paths, n, kappa, t){
  finalPrices<- hybridScheme(bssparams)(paths, n, kappa, t)
  return (ImpliedVol(k, finalPrices)^2*t)
}
