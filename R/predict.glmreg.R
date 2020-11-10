predict.glmreg <- function(object,newx,newoffset, which=1:length(object$lambda), type=c("link","response","class","coefficients","nonzero"), na.action=na.pass, ...){
 type=match.arg(type)
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
     }
  else{
if(!is.null(object$terms)){
 mf <- model.frame(delete.response(object$terms), newx, na.action = na.action, xlev = object$xlevels)
 newx <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
 newx <- newx[, -1, drop = FALSE] ### remove the intercept
}
}
  b0=as.matrix(object$b0)
  rownames(b0)="(Intercept)"
  nbeta=rbind2(b0,as.matrix(object$beta))[,which]
    vnames=dimnames(nbeta)[[1]]
    lambda=object$lambda[which]
  nlambda <- length(lambda)
  if(type=="coefficients")return(nbeta)
  if(type=="nonzero"){
  #if(nlambda>1) return(nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE))
  if(nlambda>1) return(eval(parse(text="glmnet:::nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE)")))
  else return(which(abs(nbeta[-1]) > 0))
}
  famtype <- switch(object$family,
					"gaussian"=1,
					"binomial"=2,
					"poisson"=3,
					"negbin"=4)
  n <- dim(newx)[1]
  m <- dim(newx)[2]
  if(is.list(newx))
      newx <- matrix(unlist(newx), ncol=m)
      #newx <- as.matrix(as.data.frame(lapply(newx, as.numeric)))
  if(dim(object$beta)[1] != m) stop("number of columns in newx should be the same as number of rows in beta\n")
  if(object$is.offset)
    if(is.null(newoffset))
      stop("offset is used in the estimation but not provided in prediction\n")
    else offset <- newoffset
  else offset <- rep(0, n)
  res <- .Fortran("pred",
				  n=as.integer(n),
				  m=as.integer(m),
				  nlambda=as.integer(nlambda),
				  x=as.double(newx),
				  b=as.double(object$beta[,which]),
				  a0=as.double(object$b0[which]),
				  offset=as.double(offset),
                  family=as.integer(famtype),
				  eta = as.double(matrix(0,n,nlambda)),
				  mu = as.double(matrix(0,n,nlambda)),
				  PACKAGE="mpath")
  eta <- matrix(res$eta, ncol=nlambda)
  mu <- matrix(res$mu, ncol=nlambda)
#  to be considered below: 02/28/2020
#  if(object$is.offset){
#	  if(is.null(newoffset)) stop("newoffset is required\n")
#          else 
#     {
#      eta <- eta + newoffset
#      mu <- eta
#     }
#  }
  colnames(eta) <- colnames(mu) <- colnames(object$beta[,which])
  #pihat <- exp(eta)/(1+exp(eta))
    if (object$family=="gaussian" | type=="link") return(eta)
    if (type=="response") return(mu)
	if (object$family=="binomial" & type=="class") return(eta>0)
  }
coef.glmreg <- function(object,which=1:length(object$lambda),...)
  {
  b0=matrix(object$b0, nrow=1)
  rownames(b0)="(Intercept)"
  nbeta=rbind2(b0,as.matrix(object$beta))
    return(nbeta[,which])
  }
