rm(list=ls())

library("MASS")
library("gmm")

ff25_exret<- as.matrix(data.frame(read.table("LL_portfolios.txt",header=F,sep="",na.strings="na")))
returns<-100*ff25_exret[2:141,]
fact<-data.frame(read.table("LL_factors.txt",header=T,sep="",na.strings="na"))
cay<-scale(fact$cay[1:140], scale=FALSE) #demeaned cay
cons<-fact$cons[2:141]
mkt<-fact$mkt[2:141]

#fact.list<-100*cbind(cay, cons, 10*cay*cons) #scaling as was done in LL2001
fact.list<-100*cbind(cons) #scaling as was done in LL2001
n.obs<-dim(returns)[1]
n.port<-dim(returns)[2]
factors<-array(fact.list, c(n.obs, length(fact.list)/n.obs))
n.fact<-dim(factors)[2]

################## first stage  %%%%%%%%%%%
b.factors<-cbind(array(1, c(n.obs,1)), factors)
b.betas<-t(ginv(t(b.factors)%*%b.factors) %*% t(b.factors) %*% returns)
ts.res<-returns-b.factors %*% t(b.betas)
Sigma.r<-var(ts.res)
Sigma.f<-var(factors)
b.Sigma.f<-rbind(array(0, c(1, n.fact+1)), cbind(array(0,c(n.fact,1)), Sigma.f))

tstats<-c()
for (i in 1:25){
  ols<-lm(returns[,i] ~ cons)
  print(summary(ols))
}
tstats

################## second stage  %%%%%%%%%%%
b.betas[,1]<-1
av.ret<-array(colMeans(returns), c(dim(returns)[2],1))
FM.est<-ginv((t(b.betas)%*% b.betas)) %*% t(b.betas) %*% av.ret
FM.fit.ret<-b.betas %*% FM.est
FM.R2<-1-var(av.ret-FM.fit.ret)/var(av.ret)

################## standard error for FM: OLS  %%%%%%%%%%%
P<-ginv((t(b.betas) %*% b.betas)) %*% t(b.betas)
FM.vcov.ols<-P %*% Sigma.r %*% t(P)/n.obs
FM.std.ols<-sqrt(diag(FM.vcov.ols)) 
FM.t.ols<-FM.est/FM.std.ols
FM.pval.ols<-array(0, c((n.fact+1),1))
for (i in 1:(n.fact+1))  {FM.pval.ols[i,1]<-2*min(pnorm(FM.t.ols[i,1]), 1-pnorm(FM.t.ols[i,1]))}

################# standard error for FM: Shanken   %%%%%%%%%%%%%%%%
mult<-1+t(FM.est[2:(n.fact+1),1]) %*% ginv(Sigma.f) %*% FM.est[2:(n.fact+1),1]
FM.vcov.sh<-as.numeric(mult) *P %*% Sigma.r %*% t(P)/n.obs +b.Sigma.f/n.obs
FM.std.sh<-sqrt(diag(FM.vcov.sh)) 
FM.t.sh<-FM.est/FM.std.sh
FM.pval.sh<-array(0, c((n.fact+1),1))
for (i in 1:(n.fact+1))  {FM.pval.sh[i,1]<-2*min(pnorm(FM.t.sh[i,1]), 1-pnorm(FM.t.sh[i,1]))}

######################## printing out FM estimation ########################
print("Simple Fama-MacBeth estimation:")
output<-cbind(FM.est, FM.t.ols,   FM.t.sh)
colnames(output) <- c("lambda", "t_ols", "t_shanken" )
print(output)
cat("R-squared is ", round(FM.R2,2),"\n")

######################  plotting cross-sectional fit ########################
ylim.xr=c(min(av.ret,FM.fit.ret)-0.1, max(av.ret,FM.fit.ret) +0.1)
xlim.xr=c(min(av.ret,FM.fit.ret)-0.1, max(av.ret,FM.fit.ret) +0.1)

labels_ff25<-c('11','12','13','14','15', '21','22','23','24','25', 
               '31','32','33','34','35', '41','42','43','44','45', '51','52','53','54','55')

plot(av.ret,FM.fit.ret, ylim=ylim.xr, xlim=xlim.xr, col="white",
     main="Cross-sectional model fit, Fama-MacBeth", xlab="Average returns", ylab="Fitted returns")
text(av.ret,FM.fit.ret, labels=labels_ff25) 
abline(a=0,b=1)

ylim.xr=c(min(av.ret)-0.1, max(av.ret) +0.1)
xlim.xr=c(min(b.betas[,2])-0.1, max(b.betas[,2]) +0.1)

labels_ff25<-c('11','12','13','14','15', '21','22','23','24','25', 
               '31','32','33','34','35', '41','42','43','44','45', '51','52','53','54','55')

plot(b.betas[,2], av.ret, ylim=ylim.xr, xlim=c(2.5,5.5), col="white",
     main="Cross-sectional model fit, Fama-MacBeth", xlab="Betas", ylab="Expected returns")
text(av.ret,FM.fit.ret, labels=labels_ff25) 
#abline(a=0,b=1)

########## a moment condition that includes factor means too ################
g<-function(theta, x){
  r<-x[,1:25]
  f<-x[,26:(dim(x)[2])]
  if (length(f)==140){f<-array(f, c(140,1))}
  n.f<-dim(f)[2]
  n.t<-dim(f)[1]
  n.p<-dim(r)[2]
  lambda.0<-theta[1]
  lambda<-array(theta[2:(n.f+1)], c(n.f,1))
  f.mean<-array(theta[(n.f+2):(2*n.f+1)], c(1, n.f))
  f.demeaned<-f-matrix(rep(f.mean,each=n.t),nrow=n.t)
  rf<-array(0, c(n.t, n.p))
  for (i in 1:n.f){
    rf<-rf + r * f.demeaned[,i] *lambda[i] }
  e.r<-r-lambda.0 - rf
  e<-cbind(e.r, f.demeaned)
  return(e)
}

################# Optimal GMM #######################
x<-as.matrix(cbind(returns, factors))
res<- gmm(g, x = x, t0 = c(FM.est, colMeans(factors)), wmatrix=c("optimal"), type=c("twoStep"), vcov=c("HAC"))
gmm.op.lambda<-res$coefficients
gmm.op.t<-gmm.op.lambda/sqrt(diag(res$vcov))
f.mean<-array(gmm.op.lambda[(n.fact+2):(2*n.fact+1)], c(1, n.fact))
factors.demeaned<- factors - matrix(rep(f.mean,each=n.obs),nrow=n.obs)
rf<-array(0, c(n.obs, n.port))
for (i in 1:n.fact){
  rf<-rf + returns * factors.demeaned[,i] *gmm.op.lambda[i+1] }
gmm.op.fit.ret<- colMeans(gmm.op.lambda[1] + rf)
gmm.op.R2<-1-var(av.ret-gmm.op.fit.ret)/var(av.ret)

######################## printing out GMM, W=optimal,  estimation ########################
print("GMM with optimal weight matrix:")
output<-cbind(gmm.op.lambda, gmm.op.t)
colnames(output) <- c("lambda", "t_gmm")
print(output)
cat("R-squared is ", round(gmm.op.R2,2),"\n")

######################  plotting cross-sectional fit ########################
ylim.xr=c(min(av.ret,gmm.op.fit.ret)-0.1, max(av.ret,gmm.op.fit.ret) +0.1)
xlim.xr=c(min(av.ret,gmm.op.fit.ret)-0.1, max(av.ret,gmm.op.fit.ret) +0.1)

plot(av.ret,gmm.op.fit.ret, ylim=ylim.xr, xlim=xlim.xr, col="white",
     main="Cross-sectional model fit, GMM, Optimal W", xlab="Average returns", ylab="Fitted returns")
text(av.ret,gmm.op.fit.ret, labels=labels_ff25)
abline(a=0,b=1)




################## fixed weight marix GMM #################
x<-as.matrix(cbind(returns, factors))
res<- gmm(g, x = x, t0 = c(FM.est, colMeans(factors)), wmatrix=c("ident"), type=c("twoStep"), vcov=c("HAC"))
gmm.lambda<-res$coefficients
gmm.t<-gmm.lambda/sqrt(diag(res$vcov))
f.mean<-array(gmm.lambda[(n.fact+2):(2*n.fact+1)], c(1, n.fact))
factors.demeaned<- factors - matrix(rep(f.mean,each=n.obs),nrow=n.obs)
rf<-array(0, c(n.obs, n.port))
for (i in 1:n.fact){
  rf<-rf + returns * factors.demeaned[,i] *gmm.lambda[i+1] }
gmm.fit.ret<- colMeans(gmm.lambda[1] + rf)
gmm.R2<-1-var(av.ret-gmm.fit.ret)/var(av.ret)

######################## printing out GMM, W=identity,  estimation ########################
print("GMM with identity weight matrix:")
output<-cbind(gmm.lambda, gmm.t)
colnames(output) <- c("lambda", "t_gmm")
print(output)
cat("R-squared is ", round(gmm.R2,2),"\n")

######################  plotting cross-sectional fit ########################
ylim.xr=c(min(av.ret,gmm.fit.ret)-0.1, max(av.ret,gmm.fit.ret) +0.1)
xlim.xr=c(min(av.ret,gmm.fit.ret)-0.1, max(av.ret,gmm.fit.ret) +0.1)

plot(av.ret, gmm.fit.ret, ylim=ylim.xr, xlim=xlim.xr, col="white",
     main="Cross-sectional model fit, GMM, W=I", xlab="Average returns", ylab="Fitted returns")
text(av.ret, gmm.fit.ret, labels=labels_ff25)
abline(a=0,b=1)