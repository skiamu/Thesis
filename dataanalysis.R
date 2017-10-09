# test script to get used to the ghyp package
library(R.matlab)
library(ghyp)
load("mcshapiro.test.RData")
data(smi.stocks)
# save stock data in a .mat file 
filename <- paste("stocks", ".mat", sep="")
writeMat(filename, smi.stocks=smi.stocks)
# import data red in Matlab
Returns <- readMat("Yahoo.mat"); Returns <- Returns$Returns
# check gaussianty --> reject H0
mcshapiro.test(Returns)

# fitting NIG distribution
fitted.NIG <- fit.NIGmv(data = Returns,silent = T)
logLik(fitted.NIG)
# fitting normal distribution
fitted.gauss <- fit.gaussmv(data = Returns)
logLik(fitted.gauss)
# fitting ghyp distribution
fitted.ghyp <- fit.ghypmv(data = Returns,silent = T)
logLik(fitted.ghyp)
# fitted t distribution
fitted.t <- fit.tmv(data = Returns,silent = T)
logLik(fitted.t)


lik.ratio.test(fitted.ghyp,fitted.NIG, conf.level = 0.95)
lik.ratio.test(fitted.ghyp,fitted.gauss, conf.level = 0.95)
lik.ratio.test(fitted.ghyp,fitted.t,conf.level = 0.95)

# filename <- paste("dist", ".mat", sep="")
# writeMat(filename, param = coef(fitted.NIG))
pairs(fitted.ghyp, cex = 0.5, legend.cex = 0.5, nbins = 50)



