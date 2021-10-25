##########################################################
#         R code for chapter 1
##########################################################

# Exercce

psi_i=ARMAtoMA(ar=0.4,lag.max = 10)  # représentation de Wald
psi_i

ARMAacv(ar=0.4,n=10)  

######################

#  Graph of Bitcoin

library(quantmod); library(zoo)
btc=getSymbols("BTC-USD", src="yahoo", from="2014-09-17",
               to="2021-10-20",auto.assign = FALSE)
closedAdj=zoo(btc[,6]); rt=diff(log(closedAdj))
par(mfrow=c(1,2), mar=c(2.5,3.8,1,1))
plot(closedAdj, col=4, xlab=""); plot(rt,col=4, xlab="")


# ACF and PACF

par(mfrow=c(2,2), mar=c(2.5,3.8,1,1))
acf(closedAdj, 40,na.action = na.pass, col=4, xlab="") 
acf(rt, 40,na.action = na.pass, col=4, xlab="")
pacf(closedAdj,40,na.action = na.pass,col=4, xlab="")
pacf(rt,40,na.action = na.pass,col=4, xlab="")

# ARCH LM Test of Bitcoin returns

et=rt-mean(rt,na.rm=T)
archLM(et, 12)

# reset test

reset.test(na.omit(rt), p=2,1)

# keenan test

keenan.test(na.omit(rt),2)

# McLeod test

arma_pq=arima(rt,c(1,0,1))
resid2=arma_pq$residuals^2
m=1L:7L; Q=NULL ; pval=NULL
for (i in m) { Q[i]=Box.test(resid2,i,"Ljung-Box")$statistic
pval[i]= Box.test(resid2,i,"Ljung-Box")$p.value
}
tab=rbind(m,Q,p.value=pval)
print(tab, digits=4)

# BDS test

mod=arima(rt,c(1,0,1))
residd=c(mod$residuals) ; 
fNonlinear::bdsTest(na.omit(residd))@test

# PR Test

NTS::PRnd(na.omit(residd))


















