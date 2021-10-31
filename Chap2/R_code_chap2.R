
######################################
#    R code for chapter 2
######################################
# Simulation of an ARMA(1,1) and ARCH(1)
library(fGarch) 
set.seed(12345)
arma11=arima.sim(model=list(ar=0.4, 
              ma=0.6), n=300) # Xt=0.4 X(t-1)+e(t)+0.6 e(t-1)
arch1=garchSim(garchSpec(model=list(omega=0.4,alpha=0.7, beta=0)),
               n=300, extended = T)
layout(matrix(c(1,2,1,3),nc=2))
plot(arma11, col=2,xlab="", ylab="", main="Simulation ARMA(1,1)")
plot(arch1$garch, type="l", col=2, main="Simulation ARCH(1)", 
       xlab="", ylab="")
plot(arch1$sigma^2, type="l", col=2, 
     main="Variance conditionnelle", xlab="", ylab="")

##  pacf of BTC returns

library(quantmod); library(zoo)
btc=getSymbols("BTC-USD", src="yahoo", from="2014-09-17",
               to="2021-10-20",auto.assign = FALSE)
closedAdj=zoo(na.omit(btc[,6])); rt=diff(log(closedAdj))
pacf(rt^2,na.action = na.pass, col=4, xlab="")

# Information criteria

library(rugarch)   # charger le package
spec1=ugarchspec(variance.model = list(garchOrder=c(1,0)),
                 mean.model = list(armaOrder=c(0,0)))
spec4=ugarchspec(variance.model = list(garchOrder=c(4,0)),
                 mean.model = list(armaOrder=c(0,0)))
fit1=ugarchfit(spec1,rt)  # Estimation
fit4=ugarchfit(spec4,rt)
t(infocriteria(fit1))  # critères d'information (normalisés)
t(infocriteria(fit4))

# Results of BTC returns estimation (ARCH(1)) with normal
fit1@fit$matcoef
show(fit1)

# histogram

hist(rt,col=4, prob=T, breaks = 50) # histogramme or rt
lines(density(rt), col=2,lwd=3) # density estimation (kernel)
curve(dnorm(x, mean(rt),sd(rt)), lwd=3, col="darkorchid3", add=T) # normal density
legend("topleft",bty="n",lwd=3, col=c(2,"darkorchid3"),
       legend=c("Kernel","Normal"))

# Results of BTC returns estimation (ARCH(1)) with student

specSt=ugarchspec(list(garchOrder=c(1,0)), list(armaOrder=c(0,0)),distribution.model = "std")
fitst=ugarchfit(specSt,rt)
show(fitst)




