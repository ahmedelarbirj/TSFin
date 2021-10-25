#####################################################
#        R functions
# Author: Hamrita Mohamed Essaied (mhamrita@gmail.com)
#
#####################################################

##########################################
#  Theoritical Autocovariance function
##########################################

ARMAacv=function(ar=0,ma=0,n=10){
  p=c(1,ARMAtoMA(ar,ma,n+10000))
  gg=NULL
  for(ii in 1:n) gg[ii]=sum(p[1:10000]*p[(ii+1):(10000+ii)])
  c(sum(p^2),gg)
}

###################################################
#   Tests for (non) linearity
###################################################

##  ARCH LM Test

archLM=function(x,p=1, disp=T){
  # x: time series;  p: lag order; disp: display or not the result
  mat=embed(x^2,(p+1))
  reg=summary(lm(mat[,1]~mat[,-1]))
  r2=reg$r.squared
  statistic=r2*length(resid(reg))
  cv=qchisq(0.95,p)
  pvalue=1-pchisq(statistic,p)
  if(disp){
    cat("================================================\n")
    cat("=======           ARCH LM test         =========\n")
    cat("================================================\n")
    cat("====   H0: no ARCH effects             =========\n\n")
    
    cat(" Chi-squared: ", statistic, "  df: ", p,  "  p.value: ", pvalue, 
        "  critical value 5%: ", cv, " \n")
    cat(" F-statistic", reg$fstatistic[1]," df:  ", reg$fstatistic[-1],"p.value: ",
        1-pf(reg$fstatistic[1], reg$fstatistic[2],reg$fstatistic[3]), 
        " critical value 5%: ", qf(0.95,reg$fstatistic[2], reg$fstatistic[3]), "\n")
  }
  return(invisible(list("Chi-squared"=statistic, "p-value"=pvalue,
                        "df"=p)))
}

##  RESET Test

reset.test=function(x,p=1,s=1){
  # Estimate of AR(p)
  xx=embed(x,(p+1))
  mod1=lm(xx[,1]~xx[,-1])
  scr0=sum(mod1$residuals^2)
  xx.hat=mod1$fitted.values
  et=mod1$residuals;  TT=length(xx.hat)
  # Estimate of model 2
  xxx=matrix(0,nr=TT, nc=s)
  for(ii in 1:s)  xxx[,ii]=xx.hat^(s+1)
  X=cbind(xx[,-1],xxx)
  mod2=lm(et~X);  scr1=sum(mod2$residuals^2)
  df2=TT-p-s-1 ; df1=s+p+1
  F_stat=df2/df1*(scr0-scr1)/scr1
  pval=pf(F_stat,df1,df2,lower.tail = F)
  cval=qf(0.95,df1,df2)
  dec=ifelse(pval < 0.05, "H0 is rejected", "H0 is accepted")
  cat("=======================================\n")
  cat("        RESET test                     \n")
  cat("=======================================\n\n")
  cat("      H0: the model is linear          \n\n")
  cat("F-statistic: ", F_stat, "df: ", c(df1,df2), " p-value: ", pval, "CV 5%: ", cval,"\n")
  cat("Decision: ", dec , "\n")
  
  return(invisible(list(Fstatistic=F_stat,df=c(df1,df2),p.value=pval)))
}

## Keenan Test

keenan.test=function(x,p=1){
  # Estimate of AR(p)
  xx=embed(x,(p+1))
  mod1=lm(xx[,1]~xx[,-1])
  xx.hat=mod1$fitted.values
  et=mod1$residuals
  # Estimate of model 2
  # (i)
  mod21=lm(xx.hat^2~xx[,-1])
  vt=mod21$residuals
  mod22=lm(et ~ -1+vt)
  F_stat=summary(mod22)$fstatistic[1]
  df=summary(mod22)$fstatistic[-1]
  pval=pf(F_stat,df[1],df[2],lower.tail = F)
  cval=qf(0.95,df[1],df[2])
  dec=ifelse(pval < 0.05, "H0 is rejected", "H0 is accepted")
  cat("=======================================\n")
  cat("        Keenan test                    \n")
  cat("=======================================\n")
  cat("      H0: the model is linear          \n\n")
  print(summary(mod22)$coefficients)
  cat("\n\nF-statistic: ", F_stat, "df: ", c(df), " p-value: ", pval, "CV 5%: ", cval,"\n")
  cat("Decision: ", dec , "\n")
  
  return(invisible(list(Fstatistic=F_stat,df=c(df),p.value=pval)))
}

