### R code from vignette source 'german.Rnw'

###################################################
### code chunk number 1: german.Rnw:24-25
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: german.Rnw:30-34
###################################################
library("mpath")
library("zic")
library("pscl")
data(docvisits)


###################################################
### code chunk number 3: german.Rnw:36-37
###################################################
barplot(with(docvisits, table(docvisits)), ylab="Frequency",xlab="Doctor office visits")


###################################################
### code chunk number 4: german.Rnw:41-44
###################################################
dt <- docvisits[,-(2:3)]
tmp <- model.matrix(~age30*health+age35*health+age40*health+age45*health+age50*health+age55*health+age60*health, data=dt)[,-(1:9)]
dat <- cbind(dt, tmp)


###################################################
### code chunk number 5: german.Rnw:47-52 (eval = FALSE)
###################################################
## m1 <- zeroinfl(docvisits~.|., data=dat, dist="negbin")
## summary(m1)
## cat("loglik of zero-inflated model", logLik(m1))
## cat("BIC of zero-inflated model", AIC(m1, k=log(dim(dat)[1])))
## cat("AIC of zero-inflated model", AIC(m1))


###################################################
### code chunk number 6: german.Rnw:55-59 (eval = FALSE)
###################################################
## fitbe <- be.zeroinfl(m1, data=dat, dist="negbin", alpha=0.01,  trace=FALSE)
## summary(fitbe)
## cat("loglik of zero-inflated model with backward selection", logLik(fitbe))
## cat("BIC of zero-inflated model with backward selection", AIC(fitbe, k=log(dim(dat)[1])))


###################################################
### code chunk number 7: german.Rnw:62-63 (eval = FALSE)
###################################################
##     fit.lasso <- zipath(docvisits~.|.,data = dat, family = "negbin", nlambda=100, lambda.zero.min.ratio=0.001, maxit.em=300, maxit.theta=25, theta.fixed=FALSE, trace=FALSE, penalty="enet", rescale=FALSE)


###################################################
### code chunk number 8: german.Rnw:66-69 (eval = FALSE)
###################################################
## minBic <- which.min(BIC(fit.lasso))
## coef(fit.lasso, minBic)
## cat("theta estimate", fit.lasso$theta[minBic])


###################################################
### code chunk number 9: german.Rnw:72-73 (eval = FALSE)
###################################################
## se(fit.lasso, minBic, log=FALSE)


###################################################
### code chunk number 10: german.Rnw:76-79 (eval = FALSE)
###################################################
## AIC(fit.lasso)[minBic]
## BIC(fit.lasso)[minBic]
## logLik(fit.lasso)[minBic]


###################################################
### code chunk number 11: german.Rnw:82-89 (eval = FALSE)
###################################################
## n <- dim(dat)[1]
## K <- 10
## foldid <- split(sample(1:n), rep(1:K, length = n))
## fitcv <- cv.zipath(docvisits ~ . | ., data = dat, family = "negbin", nlambda=100, 
## lambda.count=fit.lasso$lambda.count[1:30], lambda.zero= fit.lasso$lambda.zero[1:30],
## maxit.em=300, maxit.theta=1, theta.fixed=FALSE, trace=FALSE, penalty="enet", rescale=FALSE, foldid=foldid)
## cat("max of cv loglik", max(fitcv$cv))


###################################################
### code chunk number 12: german.Rnw:94-106 (eval = FALSE)
###################################################
##     fit.mcp <- zipath(docvisits~.|.,data = dat, family = "negbin", gamma.count=2.7, gamma.zero=2.7, lambda.count=fit.lasso$lambda.count[1:30], lambda.zero= fit.lasso$lambda.zero[1:30], maxit.em=300, maxit.theta=10, theta.fixed=FALSE, penalty="mnet")
## minBic <- which.min(BIC(fit.mcp))
## coef(fit.mcp, minBic)
## cat("theta estimate", fit.mcp$theta[minBic])
## se(fit.mcp, minBic, log=FALSE)
## AIC(fit.mcp)[minBic]
## BIC(fit.mcp)[minBic]
## logLik(fit.mcp)[minBic]
## fitcv <- cv.zipath(docvisits ~ . | ., data = dat, family = "negbin", gamma.count=2.7, gamma.zero=2.7, 
## lambda.count=fit.lasso$lambda.count[1:30], lambda.zero= fit.lasso$lambda.zero[1:30],
## maxit.em=300, maxit.theta=1, theta.fixed=FALSE, trace=FALSE, penalty="mnet", rescale=FALSE, foldid=foldid)
## cat("max of cv loglik", max(fitcv$cv))


###################################################
### code chunk number 13: german.Rnw:110-122 (eval = FALSE)
###################################################
##     fit.scad <- zipath(docvisits~.|.,data = dat, family = "negbin", gamma.count=2.5, gamma.zero=2.5, lambda.count=fit.lasso$lambda.count[1:30], lambda.zero= fit.lasso$lambda.zero[1:30], maxit.em=300, maxit.theta=10, theta.fixed=FALSE, penalty="snet")
## minBic <- which.min(BIC(fit.scad))
## coef(fit.scad, minBic)
## cat("theta estimate", fit.scad$theta[minBic])
## se(fit.scad, minBic, log=FALSE)
## AIC(fit.scad)[minBic]
## BIC(fit.scad)[minBic]
## logLik(fit.scad)[minBic]
## fitcv <- cv.zipath(docvisits ~ . | ., data = dat, family = "negbin", gamma.count=2.5, gamma.zero=2.5, 
## lambda.count=fit.lasso$lambda.count[1:30], lambda.zero= fit.lasso$lambda.zero[1:30],
## maxit.em=300, maxit.theta=1, theta.fixed=FALSE, trace=FALSE, penalty="snet", rescale=FALSE, foldid=foldid)
## cat("max of cv loglik", max(fitcv$cv))


