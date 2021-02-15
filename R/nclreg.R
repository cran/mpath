### penalized nonconvex loss based linear models (NCL) for regression and classification
nclreg <- function(x, ...) UseMethod("nclreg")

nclreg.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(nclreg.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

nclreg.formula <- function(formula, data, weights, offset=NULL, contrasts=NULL, ...){
    ## extract x, y, etc from the model formula and frame
    if(!attr(terms(formula, data=data), "intercept"))
        stop("non-intercept model is not implemented")
    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms") # allow model.frame to have updated it

    Y <- model.response(mf, "any") # e.g. factors are allowed
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    ## null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    weights <- as.vector(model.weights(mf))
    if(!length(weights)) weights <- rep(1, nrow(mf))
    else if(any(weights < 0)) stop("negative weights not allowed")
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if(length(weights) != length(Y))
        stop("'weights' must be the same length as response variable")

    offset <- as.vector(model.offset(mf))
    if(!is.null(offset)) {
        if(length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
    }

    RET <- nclreg_fit(X[,-1], Y, weights=weights, offset=offset, ...)
    RET$call <- match.call()
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    class(RET) <- "nclreg"
    RET
}
nclreg.matrix <- function(x, y, weights, offset=NULL, ...){
    RET <- nclreg_fit(x, y, weights, offset=offset, ...)
    RET$call <- match.call()
    return(RET)
}

### compute second derivative of robust loss, cf: Appendix in Wang (2016), Quadratic majorization for nonconvex loss with applications to the boosting algorithm
compute.2d <- function(y, f, s, family=c("clossR", "closs", "gloss", "qloss")){
    family <- match.arg(family)
    if(family=="clossR") u <- y-f
    else u <- y*f
    if(family %in% c("clossR", "closs")){
        b <- -(u^2-2*u-s^2+1)/(s^4*exp((1-u)^2/(2*s^2)))
        cval <- 1
        if(family=="closs")
            cval <- 1/(1 - exp(-1/(2*s^2)))
        cval*b
    }
    else if(family=="gloss"){
        2^s*s*exp(u)*(1+exp(u))^(-s-2)*(s*exp(u)-1)
    }
    else if(family=="qloss")
        sqrt(2/pi)*u/s^3*exp(-u^2/(2*s^2))
}

nclreg_fit <- function(x,y, weights, offset, rfamily=c("clossR", "closs", "gloss", "qloss"), s=NULL, fk=NULL, iter=10, reltol=1e-5, penalty=c("enet","mnet","snet"), nlambda=100, lambda=NULL, type.path=c("active", "nonactive", "onestep"), decreasing=FALSE, lambda.min.ratio=ifelse(nobs<nvars,.05, .001),alpha=1, gamma=3, standardize=TRUE, intercept=TRUE, penalty.factor = NULL, maxit=1000, type.init=c("bst", "ncl", "heu"), mstop.init=10, nu.init=0.1, eps=.Machine$double.eps, epscycle=10, thresh=1e-6, trace=FALSE){
### compute h value
    compute.h <- function(rfamily, y, fk_old, s, B){
        if(rfamily=="clossR")
            h <- gradient(family=rfamily, u=y-fk_old, s=s)/B+fk_old
        else if(rfamily %in% c("closs", "gloss", "qloss")) 
            h <- -y*gradient(family=rfamily, u=y*fk_old, s=s)/B+fk_old
        h
    }
    call <- match.call()
    rfamily <- match.arg(rfamily)
    penalty <- match.arg(penalty)
    type.path <- match.arg(type.path)
    type.init <- match.arg(type.init)
    if(type.path=="active")
	    active <- 1
    else active <- 0
    if (!is.null(lambda) && length(lambda) > 1 && all(lambda == cummin(lambda))){
	    decreasing <- TRUE
#	    if(type.path=="active"){
#	    decreasing <- FALSE
#	      warnings("choose type.path='nonactive' with increasing sequence of lambda, or let lambda=rev(lambda) as computed below\n")
#	    lambda <- rev(lambda) 
#       }
    }
    else if(!is.null(lambda) && length(lambda) > 1 && all(lambda == cummax(lambda)))
	    decreasing <- FALSE
#    else if(is.null(lambda) && decreasing && type.path=="active")
#	    stop("set decreasing=FALSE or type.path='nonactive'")
    if(rfamily %in% c("closs", "gloss", "qloss"))
        if(!all(names(table(y)) %in% c(1, -1)))
            stop("response variable must be 1/-1 for family ", rfamily, "\n")
    nm <- dim(x)
    nobs <- n <- nm[1]
    if(length(y)!=n) stop("length of y is different from row of x\n")
    nvars <- m <- nm[2]
    B <- bfunc(family=rfamily, s=s)
    if(missing(weights)) weights=rep(1,nobs)
    weights <- as.vector(weights)
    w <- weights/sum(weights)
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    ## check weights and offset
    if( !is.null(weights) && any(weights < 0) ){
        stop("negative weights not allowed")
    }
    if (is.null(offset)){
         is.offset <- FALSE
         offset <- rep.int(0, nobs)
     }else is.offset <- TRUE
    pentype <- switch(penalty,
                      "enet"=1,
                      "mnet"=2,
                      "snet"=3)

    if(is.null(penalty.factor))
        penalty.factor <- rep(1, nvars)
    if(all(penalty.factor==0)){	
        lambda <- rep(0, nvars) 
        penalty.factor <- rep(1, nvars)
    }
    xold <- x
    if(standardize){
         tmp <- stan(x, weights)
         x <- tmp$x
         meanx <- tmp$meanx
         normx <- tmp$normx
    }
    if(is.null(s)){
        warning("missing nonconvex loss tuning parameter s value, and a default value is provided\n")
        s <- switch(rfamily,       
                    "closs"=1,
                    "gloss"=1,
                    "qloss"=2,
                    "clossR"=1
                    )
    }
    penfac <- penalty.factor/sum(penalty.factor) * nvars
    zscore <- function(rfamily, RET, s){
        if(rfamily=="clossR"){
            z <- gradient(family=rfamily, u=y-RET$fitted.values, s=s)
            scores <- abs(crossprod(x, w*z))/(penfac*alpha)
        } 
        else if(rfamily %in% c("closs", "gloss", "qloss")){ 
            z <- gradient(family=rfamily, u=y*RET$fitted.values, s=s)
            scores <- abs(crossprod(x, w*(y*z)))/(penfac*alpha) ### note y=1/-1, thus can be further simplified without y
        }
    }
### initiate predictive value
    start <- NULL
    if(is.null(fk) || is.null(lambda)){
        if(type.init %in% c("ncl", "heu")){ ### use ncl function to generate intercept-only model
            RET <- ncl_fit(x=matrix(1, ncol=1, nrow=length(y)), y=y, iter=10000, reltol=1e-20, weights=weights, s=s, rfamily=rfamily, trace=FALSE)
            if(type.init=="ncl") start <- c(coef(RET), rep(0, nvars))
### it is similar to the following optimization results
                                        #fn <- function(b) sum(loss(y, f=b, cost, family = rfamily, s=s))
                                        #res <- optim(0, fn, method="SANN")
                                        #RET$fitted.values <- res$par
                                        #RET$h <- compute.h(rfamily, y, RET$fitted.values, s, B) 
### end of intercept computing
### check KKT condition for intercept only model, the following value should be close to 0
### if(rfamily=="gloss")
### mean(w*s*2^s*y*exp(y*coef(RET))*(exp(y*coef(RET))+1)^(-s-1)
            else if(type.init=="heu"){ ### heuristic for choosing starting values
                v <- zscore(rfamily, RET, s)
                ix <- which(v >= quantile(v, 0.9))
                b0.1 <- coef(RET)
                beta.1 <- rep(0, nvars)
                beta.1[ix] <- 1
                start <- c(b0.1, beta.1)
                RET$fitted.values <- x %*% beta.1 + b0.1
                RET$h <- compute.h(rfamily, y, RET$fitted.values, s, B) 
            }
        }
        else if(type.init=="bst") ### use bst function to generate results
        {
	    RET <- bst(x, y, family=rfamily, ctrl = bst_control(mstop=mstop.init, nu=nu.init, s=s, intercept=TRUE))
            RET$fitted.values <- RET$yhat
            RET$h <- compute.h(rfamily, y, RET$fitted.values, s, B) 
	    start <- c(attributes(coef(RET))$intercept, coef(RET))
        }
    }
    else {
        RET <- NULL
        RET$fitted.values <- fk
    }
    cost <- 0.5 # the optimization problem doesn't use cost argument, but loss function is from bst package
    los <- 2*loss(y, f=RET$fitted.values, cost, family = rfamily, s=s, fk=RET$fitted.values)
    if(is.null(lambda)){
        h <- RET$h                               
### If standardize=TRUE, x has already been standardized. Therefore, in the sequel, we don't standardize x again
### method A, to obtain lambda values from fitting the penalized regression. 
        lambda <- glmreg_fit(x=x*sqrt(B), y=h*sqrt(B), weights=weights, offset=offset, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, intercept=intercept, penalty.factor = penalty.factor, maxit=1, eps=eps, family="gaussian", penalty=penalty)$lambda
### with type.init, add two different lambda sequences and make lambda values flexible to have different solution paths
        #if(type.init %in% c("bst", "heu")){
            if(!decreasing)#{### solution path backward direction
                lambda <- rev(lambda)
        #}
                                        #cat("lambda in method A", lambda[1], "\n")
### method B is the same as method A to obtain lmax 
                                        #w <- weights/sum(weights)
                                        #        h1 <- h*sqrt(B)
                                        #        x1 <- x*sqrt(B)
                                        #	lmax <- max(abs(crossprod(x1,w*(h1-mean(h1))))/(penfac*alpha))
                                        #        cat("lambda in method B", lmax, "\n")
### method C considers the KKT conditions for the original nonconvex loss
                                        #assuming    RET <- ncl(y~1, data=data.frame(y, 1), iter=10000, del=1e-20, weights=weights, s=s, rfamily=rfamily, trace=FALSE)
					#if(rfamily=="clossR"){
                                        #z <- gradient(family=rfamily, u=y-RET$fitted.values, s=s)
                                        #lmax <- max(abs(crossprod(x, w*z))/(penfac*alpha))
                                        # } 
                                        #else if(rfamily %in% c("closs", "gloss", "qloss")){ 
                                        #    z <- gradient(family=rfamily, u=y*RET$fitted.values, s=s)
                                        #    lmax <- max(abs(crossprod(x, w*(y*z)))/(penfac*alpha)) ### note y=1/-1, thus can be further simplified without y
                                        #}
					#cat("lambda in method C", lmax, "\n")
### method D is the same as method C, but with the expressed gradients formula
                                        #assuming    RET <- ncl(y~1, data=data.frame(y, 1), iter=10000, del=1e-20, weights=weights, s=s, rfamily=rfamily, trace=FALSE)
					#b0 <- coef(RET)
					#if(rfamily=="clossR")
					#lmax <- max(abs(crossprod(x, w/s^2*(y-b0)*exp(-(y-b0)^2/(2*s^2)))/(penfac*alpha)))
                                        # if(rfamily=="closs"){
                                        #    cval <- 1/(1 - exp(-1/(2*s^2)))
                                        #    lmax <- max(abs(crossprod(x, w*cval/s^2*y*(y*coef(RET)-1)*exp(-(1-y*coef(RET))^2/(2*s^2)))/(penfac*alpha)))
                                        #}
                                        # else if(rfamily=="gloss")
                                        #    lmax <- max(abs(crossprod(x, w*s*2^s*y*exp(y*coef(RET))*(exp(y*coef(RET))+1)^(-s-1))/(penfac*alpha)))
                                        # else if(rfamily=="qloss")
                                        #       lmax <- max(abs(crossprod(x, w*sqrt(2/pi)/s*y*exp(-coef(RET)^2/(2*s^2)))/(penfac*alpha)))
					#cat("lambda in method D", lmax, "\n")
                                        #        lpath <- seq(log(lmax), log(lambda.min.ratio * lmax), length.out=nlambda)
                                        #        lambda <- exp(lpath)
    }
    else
        nlambda <- length(lambda)
    beta <- matrix(0, ncol=nlambda, nrow=m)
    b0 <- rep(0, nlambda)
    ### initialize start values -begin
    fk_old <- RET$fitted.values
    h <- compute.h(rfamily, y, fk_old, s, B)
    if(any(is.nan(h))){
      cat('h has NaN, fk_old used to compute h at beginning', fk_old)
      h <- y
    }
    tmp <- init(w, h, offset, family="gaussian")
    mustart <- tmp$mu
    etastart <- tmp$eta
    ### initialize the start values -end
    stopit <- FALSE
### for each element of the lambda sequence, iterate until convergency
    typeA <- function(beta, b0){
        i <- 1
        los <- pll <- matrix(NA, nrow=iter, ncol=nlambda)
        while(i <= nlambda){
            if(trace) message("loop in lambda:", i, "\n")
            if(trace) {
                cat("Quadratic majorization iterations ...\n")
            }
            k <- 1
            d1 <- 10
            while(d1 > reltol && k <= iter){
                fk_old <- RET$fitted.values
                h <- compute.h(rfamily, y, fk_old, s, B)
                if(any(is.nan(h))){ # exit loop
                    stopit <- TRUE
                    break
                }
		RET <- glmreg_fit(x=x*sqrt(B), y=h*sqrt(B), weights=weights, offset=offset, lambda=lambda[i],alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, intercept=intercept, penalty.factor = penalty.factor, maxit=maxit, eps=eps, family="gaussian", penalty=penalty, start=start)
		RET$b0 <- RET$b0/sqrt(B)
### for LASSO, the above two lines are equivalent to the next line
                                        #RET <- glmreg_fit(x=x, y=h, weights=weights, lambda=lambda[i]/B,alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, penalty.factor = penalty.factor, maxit=maxit, eps=eps, family="gaussian", penalty=penalty)
                fk <- RET$fitted.values <- predict(RET, newx=x)
                start <- coef(RET)
                los[k, i] <- 2*mean(loss(y, f=fk, cost, family = rfamily, s=s, fk=NULL))
###penalized loss value for beta
                penval <- .Fortran("penGLM",
                                   start=as.double(RET$beta),
                                   m=as.integer(m),
                                   lambda=as.double(lambda[i]*penalty.factor),
                                   alpha=as.double(alpha),
                                   gam=as.double(gamma),
                                   penalty=as.integer(pentype),
                                   pen=as.double(0.0),
                                   PACKAGE="mpath")$pen
                if(standardize)  ### lambda value, hence penval, depends on whether standardize is TRUE/FALSE
                    pll[k, i] <- los[k, i] + n*penval
                else pll[k, i] <- los[k, i] + penval
                d1 <- sum((fk_old - fk)^2)
                                        #d1 <- sum((fk_old - fk)^2)/sum(fk_old^2) ### this can cause a problem if fk_old is zero
                if(trace) cat("\n  iteration", k, ": relative change of fk", d1, ", robust loss value", los[k, i], ", penalized loss value", pll[k, i], "\n")
                if(trace) cat("  d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
                k <- k + 1
            }
            if(!stopit){
                beta[,i] <- as.vector(RET$beta)
                b0[i] <- RET$b0
                i <- i + 1
            }
            else i <- nlambda + 1
        }
        list(beta=beta, b0=b0, RET=RET, risk=los, pll=pll)
    }
### only for direction="bwd", cycle through only on active set. for each element of the lambda sequence, iterate until convergency
    typeB <- function(beta, b0){
	i <- 1
	x.act <- x ### current active set
	start.act <- start
	m.act <- dim(x.act)[2]
	penalty.factor.act <- penalty.factor
        los <- pll <- matrix(NA, nrow=iter, ncol=nlambda)
        while(i <= nlambda){
            if(trace) message("loop in lambda:", i, ", lambda=", lambda[i], "\n")
            if(trace) {
                cat("Quadratic majorization iterations ...\n")
            }
            k <- 1
            d1 <- 10
            while(d1 > reltol && k <= iter){
                fk_old <- RET$fitted.values
                h <- compute.h(rfamily, y, fk_old, s, B)
	    	if(any(is.nan(h))){ # exit loop 
                   cat('h has NaN, in typeB fk_old used to compute h', fk_old)
	    		stopit <- TRUE
                    break
                }
               RET <- .Fortran("glmreg_fit_fortran",
	                        x=as.double(x.act), 
				y=as.double(h), 
				weights=as.double(weights),
			        n=as.integer(n),
			        m=as.integer(m.act),	
				start=as.double(start.act), 
				etastart=as.double(etastart),
	                        mustart=as.double(mustart),
			       	offset=as.double(offset),
			       	nlambda=as.integer(1),
			       	lambda=as.double(lambda[i]/B), 
				alpha=as.double(alpha),        
			        gam=as.double(gamma), 
				rescale=as.integer(0),
			       	standardize=as.integer(0),
			       	intercept=as.integer(intercept),
				penaltyfactor=as.double(penalty.factor.act),
				thresh=as.double(thresh), 
				epsbino=as.double(0), 
				maxit=as.integer(maxit),
			        eps=as.double(eps), 
				theta=as.double(0),
			       	family=as.integer(1), 
				penalty=as.integer(pentype), 
				trace=as.integer(trace),
				beta=as.double(matrix(0, ncol=1, nrow=m)),
				b0=as.double(0),
				yhat=as.double(rep(0, n)),
                satu=as.integer(0),
				PACKAGE="mpath")
     # epsbino, theta are not used in the above Fortran call with family="gaussian"
#		RET <- glmreg_fit(x=x.act*sqrt(B), y=h*sqrt(B), weights=weights, offset=offset, lambda=lambda[i],alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, penalty.factor = penalty.factor.act, maxit=maxit, eps=eps, family="gaussian", penalty=penalty, start=start.act)
#		RET$b0 <- RET$b0/sqrt(B)  ### this line should be removed after calling glmreg_fit_fortran
### for LASSO, the above two lines are equivalent to the next line if and only if standardize=FALSE
                                        #RET <- glmreg_fit(x=x, y=h, weights=weights, lambda=lambda[i]/B,alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, penalty.factor = penalty.factor, maxit=maxit, eps=eps, family="gaussian", penalty=penalty)
#                fk <- RET$fitted.values <- predETt(RET, newx=x.act)
#                start.act <- coef(RET)
                fk <- RET$fitted.values <- RET$yhat
                start.act <- c(RET$b0, RET$beta)
                etastart <- RET$yhat
		mustart <- RET$mustart
	  	los[k, i] <- 2*mean(loss(y, f=fk, cost, family = rfamily, s=s, fk=NULL))
###penalized loss value for beta
if(any(is.na(RET$beta))){
                cat('rfamily=', rfamily, 's=', s, 'B=', B, "\n")
		cat('fk_old', fk_old, "\n")
	        cat("h=", h, "\n")
	        cat("start.act=", start.act, "\n")
	        cat("mustart=", mustart, "\n")
	        cat("etastart=", etastart, "\n")
}
                penval <- .Fortran("penGLM",
                                   start=as.double(RET$beta),
                                   m=as.integer(m.act),
                                   lambda=as.double(lambda[i]/B*penalty.factor.act),
                                   alpha=as.double(alpha),
                                   gam=as.double(gamma),
                                   penalty=as.integer(pentype),
                                   pen=as.double(0.0),
                                   PACKAGE="mpath")$pen
     		if(standardize)  ### lambda value, hence penval, depends on whether standardize is TRUE/FALSE
                    pll[k, i] <- los[k, i] + n*penval
                else pll[k, i] <- los[k, i] + penval
                d1 <- sum((fk_old - fk)^2)
                                        #d1 <- sum((fk_old - fk)^2)/sum(fk_old^2) ### this can cause a problem if fk_old is zero
                if(trace) cat("\n  iteration", k, ": relative change of fk", d1, ", robust loss value", los[k, i], ", penalized loss value", pll[k, i], "\n")
                if(trace) cat("  d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
                k <- k + 1
            }
            if(!stopit){
                tmp <- as.vector(RET$beta)
                activeset <- which(abs(tmp) > 0)
                if(length(activeset)==0)    break
                beta[activeset,i] <- tmp[activeset]
		start.act <- start.act[which(abs(start.act) > 0)]
		x.act <- x[, activeset] #update active set
		if(length(activeset)==1)
                    x.act <- matrix(x.act, ncol=1)
		m.act <- length(activeset) #update active set
	        penalty.factor.act <- penalty.factor[activeset]
                                        # beta[,i] <- as.vector(RET$beta)
                b0[i] <- RET$b0
      		i <- i + 1
            }
            else i <- nlambda + 1
        }
        list(beta=beta, b0=b0, RET=RET, risk=los, pll=pll)
    }
#typeBB was converted from typeB, works for both "bwd" and "fwd".
### these number conversion should be consistent to those in Fortran subroutines src/compute_h and loss
    rfamilytype <- switch(rfamily,       
                    "clossR"=11,
                    "closs"=12,
                    "gloss"=13,
                    "qloss"=14
                    )
    typeBB <- function(beta, b0){
        if(type.path=="active" && decreasing)
	    RET <- .Fortran("nclreg_ad",
			    x=as.double(x), 
			    y=as.double(y),
			    weights=as.double(weights),
			    n=as.integer(n),
			    m=as.integer(m),
			    start=as.double(start), 
			    etastart=as.double(etastart),
			    mustart=as.double(mustart),
			    offset=as.double(offset),
			    iter=as.integer(iter),
			    nlambda=as.integer(nlambda), 
			    lambda=as.double(lambda), 
			    alpha=as.double(alpha),
		        gam=as.double(gamma), 
			    standardize=as.integer(0), 
			    intercept=as.integer(intercept), 
			    penaltyfactor=as.double(penalty.factor),
			    maxit=as.integer(maxit), 
			    eps=as.double(eps), 
			    epscycle=as.double(epscycle), 
			    penalty=as.integer(pentype), 
			    trace=as.integer(trace), 
			    del=as.double(reltol), 
			    rfamily=as.integer(rfamilytype),
			    B=as.double(B), 
			    s=as.double(s),
			    thresh=as.double(thresh),
			    #cost=as.double(cost),
			    decreasing=as.integer(decreasing),
			    beta=as.double(beta),
		        b0=as.double(b0), 
			    yhat=as.double(RET$fitted.values), 
			    los=as.double(rep(0, nlambda)),
			    pll=as.double(rep(0, nlambda)),
			    nlambdacal=as.integer(0),
			    PACKAGE="mpath")
    else
	    RET <- .Fortran("nclreg_fortran",
			    x=as.double(x), 
			    y=as.double(y),
			    weights=as.double(weights),
			    n=as.integer(n),
			    m=as.integer(m),
			    start=as.double(start), 
			    etastart=as.double(etastart),
			    mustart=as.double(mustart),
			    offset=as.double(offset),
			    iter=as.integer(iter),
			    nlambda=as.integer(nlambda), 
			    lambda=as.double(lambda), 
			    alpha=as.double(alpha),
		        gam=as.double(gamma), 
			    standardize=as.integer(0), 
			    intercept=as.integer(intercept), 
			    penaltyfactor=as.double(penalty.factor),
			    maxit=as.integer(maxit), 
			    eps=as.double(eps), 
			    epscycle=as.double(epscycle), 
			    penalty=as.integer(pentype), 
			    trace=as.integer(trace), 
			    del=as.double(reltol), 
			    rfamily=as.integer(rfamilytype),
			    B=as.double(B), 
			    s=as.double(s),
			    thresh=as.double(thresh),
			    #cost=as.double(cost),
			    decreasing=as.integer(decreasing),
			    active=as.integer(active),
			    beta=as.double(beta),
		        b0=as.double(b0), 
			    yhat=as.double(RET$fitted.values), 
			    los=as.double(rep(0, nlambda)),
			    pll=as.double(rep(0, nlambda)),
			    nlambdacal=as.integer(0),
			    PACKAGE="mpath")
       	    list(beta=matrix(RET$beta, ncol=nlambda), b0=RET$b0, RET=RET, risk=RET$los, pll=RET$pll, nlambdacal=RET$nlambdacal)
    }
### update for one element of lambda depending on increasing sequence of lambda (last element of lambda) or decreasing sequence of lambda (then first element of lambda) in each MM iteration, and iterate until convergency of prediction. Then fit a solution path based on the sequence of lambda. Experiment code. Argument direction has been removed.
    typeC <- function(beta, b0){
        if(trace) {
            cat("Quadratic majorization iterations ...\n")
        }
        los <- pll <- rep(NA, nlambda)
        k <- 1
        fk <- RET$fitted.values
        d1 <- 10
        lam <- lambda[ifelse(decreasing, nlambda, 1)]
	while(d1 > reltol && k <= iter){
	    fk_old <- fk
            h <- compute.h(rfamily, y, fk_old, s, B)
	    RET <- glmreg_fit(x=x*sqrt(B), y=h*sqrt(B), weights=weights, offset=offset, lambda=lam, alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, intercept=intercept, penalty.factor = penalty.factor, maxit=maxit, eps=eps, family="gaussian", penalty=penalty, start=start)
            RET$b0 <- RET$b0/sqrt(B)
	    fk <- predict(RET, newx=x)
            start <- coef(RET)
            d1 <- sum((fk_old - fk)^2)
            if(trace) cat("\n  iteration", k, ": relative change of fk", d1, "\n")
            if(trace) cat("  d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
            k <- k + 1
        }
### fit a solution path    
        RET <- glmreg_fit(x=x*sqrt(B), y=h*sqrt(B), weights=weights, offset=offset, lambda=lambda, alpha=alpha,gamma=gamma, rescale=FALSE, standardize=FALSE, intercept=intercept, penalty.factor = penalty.factor, maxit=maxit, eps=eps, family="gaussian", penalty=penalty, start=start)
        RET$b0 <- RET$b0/sqrt(B)
        for(i in 1:nlambda){
            los[i] <- 2*weighted.mean(loss(y, f=RET$fitted.values[,i], cost, family = rfamily, s=s, fk=NULL), w=w)
###penalized loss value for beta
            penval <- .Fortran("penGLM",
                               start=as.double(RET$beta),
                               m=as.integer(m),
                               lambda=as.double(lambda[i]*penalty.factor),
                               alpha=as.double(alpha),
                               gam=as.double(gamma),
                               penalty=as.integer(pentype),
                               pen=as.double(0.0),
                               PACKAGE="mpath")$pen
            if(standardize)  ### lambda value, hence penval, depends on whether standardize is TRUE/FALSE
                pll[i] <- los[i] + n*penval
            else pll[i] <- los[i] + penval
	}
        beta <- RET$beta
        b0 <- RET$b0
        list(beta=beta, b0=b0, RET=RET, risk=los, pll=pll)
    }
#    if(type.path=="active") tmp <- typeB(beta, b0)
# typeA and typeB are combined into typeBB
    if(type.path=="active" || type.path=="nonactive") tmp <- typeBB(beta, b0)
#    if(type.path=="active") tmp <- typeBB(beta, b0)
#    else if(type.path=="nonactive") tmp <- typeA(beta, b0)
    else if(type.path=="onestep") tmp <- typeC(beta, b0)
    beta <- tmp$beta
    b0 <- tmp$b0
    RET <- tmp$RET
    RET$risk <- tmp$risk
    if(standardize) {
        beta <- beta/normx
        b0 <- b0 - crossprod(meanx, beta)
    }
    if(is.null(colnames(x))) varnames <- paste("V", 1:ncol(x), sep="")
    else varnames <- colnames(x)
    dimnames(beta) <- list(varnames, round(lambda, digits=4))
    RET$beta <- beta
    RET$b0 <- matrix(b0, nrow=1)
    RET$x <- xold
    RET$y <- y
    RET$call <- call
    RET$lambda <- lambda
    RET$nlambda <- nlambda
    RET$penalty <- penalty
    RET$s <- s
    RET$risk <- tmp$risk
    RET$pll <- tmp$pll
    RET$nlambdacal <- tmp$nlambdacal
    RET$rfamily <- RET$family <- rfamily
    RET$type.init <- type.init
    RET$mstop.init <- mstop.init
    RET$nu.init <- nu.init
    RET$decreasing <- decreasing
    RET$type.path <- type.path
    RET$is.offset <- is.offset
    class(RET) <- "nclreg"
    RET
}

predict.nclreg <- function(object, newdata=NULL, newy=NULL, newoffset=NULL, which=1:object$nlambda, type=c("response","class","loss", "error", "coefficients", "nonzero"), na.action=na.pass, ...){
    type=match.arg(type)
    if(is.null(newdata)){
        if(!match(type,c("coefficients", "nonzero"),FALSE))stop("You need to supply a value for 'newdata'")
        ynow <- object$y
    }
    else{
        if(!is.null(object$terms)){
            mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$xlevels)
            ynow <- model.frame(object$terms, newdata)[,1] ### extract response variable
            newdata <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
            if(!is.null(ynow) && !is.null(newy))
                warnings("response y is used from newdata, but newy is also provided. Check if newdata contains the same y as newy\n")
        }
        else ynow <- newy
    }
    if(type=="coefficients") return(coef.glmreg(object)[,which])
    if(type=="nonzero"){
        nbeta <- object$beta[,which]
        if(length(which)>1) return(eval(parse(text="glmnet:::nonzeroCoef(nbeta[,,drop=FALSE],bystep=TRUE)")))
        else return(which(abs(nbeta) > 0))
    }
    if(is.null(newdata))
        newx <- as.data.frame(object$x)
    else newx <- as.data.frame(newdata)
    if(dim(newx)[2]==dim(object$beta)[1]) ### add intercept
        newx <- cbind(1, newx)
    if(object$is.offset)
     if(is.null(newoffset))
       stop("offset is used in the estimation but not provided in prediction\n")
     else offset <- newoffset
   else offset <- rep(0, length(ynow))
    res <- offset + as.matrix(newx) %*% rbind(object$b0, object$beta)
    if(type=="response") return(res[, which])
    if(type %in% c("loss", "error") && is.null(ynow)) stop("response variable y missing\n")
    if(type=="loss"){
        tmp <- rep(NA, length(which))
        for(i in 1:length(which))
            tmp[i] <- mean(loss(ynow, res[,which[i]], family = object$family, s=object$s))
        return(tmp)
    }
    if(type=="error"){
        tmp <- rep(NA, length(which))
        for(i in 1:length(which))
            tmp[i] <- evalerr(object$family, ynow, res[,which[i]])
        return(tmp)
    }
}


coef.nclreg <- function(object, ...)
    coef.glmreg(object)

### standardize predictor variables to have mean 0 and mean value of sum of squares = 1
stan <- function(x, weights){
    w <- weights/sum(weights)
    xx <- x
    meanx <- drop(w %*% x)
    x <- scale(x, meanx, FALSE) # centers x
    x <- sqrt(w) * x
    one <- rep(1, dim(x)[1])
    normx <- sqrt(drop(one %*% (x^2)))
    x <- scale(xx, meanx, normx)
    list(x=x, meanx=meanx, normx=normx) 
}
