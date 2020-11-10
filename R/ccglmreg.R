### penalized optimization for regression and classification
ccglmreg <- function(x, ...) UseMethod("ccglmreg")

ccglmreg.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(ccglmreg.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

ccglmreg.formula <- function(formula, data, weights, offset=NULL, contrasts=NULL, ...){
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

    RET <- ccglmreg_fit(X[,-1], Y, weights=weights, offset=offset, ...)
    RET$call <- match.call()
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    class(RET) <- "ccglmreg"
    RET
}
ccglmreg.matrix <- function(x, y, weights, offset=NULL, ...){
    RET <- ccglmreg_fit(x, y, weights, offset=offset, ...)
    RET$call <- match.call()
    return(RET)
}

ccglmreg_fit <- function(x,y, weights, offset, cfun="ccave", dfun="gaussian", 
                      s=NULL, delta=0.1, fk=NULL, iter=10, reltol=1e-5, penalty=c("enet","mnet","snet"), nlambda=100, lambda=NULL, type.path=c("active", "nonactive"), decreasing=TRUE, lambda.min.ratio=ifelse(nobs<nvars,.05, .001),alpha=1, gamma=3, rescale=TRUE, standardize=TRUE, intercept=TRUE, penalty.factor = NULL, maxit=1000, type.init=c("bst", "co", "heu"), init.family=NULL, mstop.init=10, nu.init=0.1, eps=.Machine$double.eps, epscycle=10, thresh=1e-6, parallel=FALSE, n.cores=2, theta, trace=FALSE, tracelevel=1){
### compute h value
    compute.h <- function(rfamily, y, fk_old, s, B){
        if(rfamily=="clossR")
            h <- gradient(family=rfamily, u=y-fk_old, s=s)/B+fk_old
        else if(rfamily %in% c("closs", "gloss", "qloss"))
            h <- -y*gradient(family=rfamily, u=y*fk_old, s=s)/B+fk_old
        h
    }
    dfunold2 <- dfun
    dfunold <- dfun[[1]]
    ### there will be no gaussianC to distinguish
    if(dfunold %in%c("gaussian","gaussianC") || dfunold %in%c(1,4)) rescale <- FALSE
    cfunold <- cfun
    cfun <- cfun2num(cfun)
    dfun <- dfun2num(dfunold)
    if(!(cfun %in% 1:8)) stop("cfun is not implemented\n")
    if(!(dfun %in% c(1, 4, 5, 8, 9))) stop("dfun is not implemented\n")
    if(dfun==1 | dfun==4) dfunnew <- "gaussian"
    else dfunnew <- dfunold
     if(dfunold=="poisson"){
         theta <- 1 ### not used anyway
         family <- 3
     }else if(dfunold=="negbin"){
         if(missing(theta)) stop("theta has to be provided for family=negbin()")
         family <- 4
     }else theta <- 0 # not used
    #else if(dfun==5) dfunnew <- "binomial"
    #else stop("not implemented yet for dfun=", dfun, "\n")
    call <- match.call()
    penalty <- match.arg(penalty)
    type.path <- match.arg(type.path)
    type.init <- match.arg(type.init)
    if(type.path=="active")
        active <- 1
    else active <- 0
    if (!is.null(lambda) && length(lambda) > 1 && all(lambda == cummin(lambda))){
        decreasing <- TRUE
    }
    else if(!is.null(lambda) && length(lambda) > 1 && all(lambda == cummax(lambda)))
        decreasing <- FALSE
    if(!is.null(init.family)) rfamily <- init.family else{
                                                         if(cfun==4 && dfun==1) rfamily <- "clossR"
                                                         else if(cfun==4 &&     dfun==4) rfamily <- "closs"
                                                         else rfamily <- "gaussian"
                                                     }
    if(dfun %in% 4:6)
        y <- y2num(y)
    if(!is.matrix(x)) x <- matrix(x)
    nm <- dim(x)
    nobs <- n <- nm[1]
    if(length(y)!=n) stop("length of y is different from row of x\n")
    nvars <- m <- nm[2]
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
    if(is.null(s) || is.na(s)) s <- assign_s(cfun, y) else check_s(cfun, s)
    if(cfun==6)
        if(s > 1) delta <- (s-1)/2 else
                                       if(s==1) delta <- 0 else{
                                                               if(is.null(delta)) stop("delta must be provided")
                                                               if(delta <= 0) stop("delta must be positive")
                                                           }
    penfac <- penalty.factor/sum(penalty.factor) * nvars
    zscore <- function(rfamily, RET, s){
        if(rfamily=="clossR"){
            z <- gradient(family=rfamily, u=y-RET$fitted.values, s=s)
            scores <- abs(crossprod(x, w*z))/(penfac*alpha)
        } 
        else if(rfamily=="closs"){ 
            z <- gradient(family=rfamily, u=y*RET$fitted.values, s=s)
            scores <- abs(crossprod(x, w*(y*z)))/(penfac*alpha) ### note y=1/-1, thus can be further simplified without y
        }
    }
### initiate predictive value
    start <- NULL
    if(is.null(fk) || is.null(lambda)){
        if(type.init %in% c("co", "heu")){ ### use co function to generate intercept-only model
            RET <- ccglm(y~1, data=data.frame(cbind(y, rep(1, n))), iter=10000, reltol=1e-20, weights=weights, s=s, cfun=cfun, dfun=dfunold2, init.family=init.family, trace=FALSE)
            if(type.init=="co") start <- c(coef(RET), rep(0, nvars))
### end of intercept computing
            else if(type.init=="heu"){ ### heuristic for choosing starting values
                v <- zscore(rfamily, RET, s)
                ix <- which(v >= quantile(v, 0.9))
                b0.1 <- coef(RET)
                beta.1 <- rep(0, nvars)
                beta.1[ix] <- 1
                start <- c(b0.1, beta.1)
                RET$fitted.values <- x %*% beta.1 + b0.1
            }
        }
        else if(type.init=="bst") ### use bst function to generate results
        {
            RET <- bst(x, y, family=rfamily, ctrl = bst_control(mstop=mstop.init, nu=nu.init, s=s, intercept=TRUE))
            RET$fitted.values <- RET$yhat
            RET$weights_update <- weights_update <- weights
            #.Fortran("update_wt",
            #                                                 n=as.integer(n),
            #                                                 weights=as.double(weights),
            #                                                 y=as.double(y),
            #                                                 f=as.double(RET$yhat),
            #                                                 cfun=as.integer(cfun),
            #                                                 dfun=as.integer(dfun),
            #                                                 s=as.double(s),
            #                                                 delta=as.double(delta),
            #                                                 weights_update=as.double(rep(0, n)),
            #                                                 PACKAGE="mpath")$weights_update
            start <- c(attributes(coef(RET))$intercept, coef(RET))
        }
    }
    else {
        RET <- NULL
        RET$fitted.values <- fk
    }
    #los <- .Fortran("loss2",
    #                n=as.integer(n),
    #                y=as.double(y),
    #                f=as.double(RET$fitted.values),
    #                weights=as.double(weights/n),
    #                cfun=as.integer(cfun),
    #                dfun=as.integer(dfun),
    #                s=as.double(s),
    #                delta=as.double(delta),
    #                los=as.double(0.0),
    #                PACKAGE="mpath")$los
    if(dfun==5) ytmp <- (y+1)/2 else ytmp <- y
    if(is.null(lambda)){
### method A, to obtain lambda values from fitting the penalized regression. 
        lambda <- try(glmreg_fit(x, ytmp, weights=RET$weights_update, offset=offset, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, alpha=alpha,gamma=gamma, rescale=FALSE, standardize=standardize, intercept=intercept, penalty.factor = penalty.factor, maxit=1, eps=eps, family=dfunnew, penalty=penalty)$lambda)  ### changed 8/17/2020
        if(inherits(lambda, "try-error"))
          stop("Initial value can't compute penalty lambda values. Possible reason: initial weights are all zero. Try to enlarge s value, or change type.init/init.family\n")
### with type.init, add two different lambda sequences and make lambda values flexible to have different solution paths
        if(!decreasing)#{### solution path backward direction
            lambda <- rev(lambda)
    }
    nlambda <- length(lambda)
    beta <- matrix(0, ncol=nlambda, nrow=m)
    fitted <- matrix(NA, ncol=nlambda, nrow=n)
    b0 <- rep(0, nlambda)
    weights_cc <- matrix(0, nrow=n, ncol=nlambda)
### initialize start values -begin
### Todo: initial values?
    mustart <- rep(0, n)
    etastart <- rep(0, n)
    fk_old <- RET$fitted.values
### check next line if update is needed
    if(type.init %in% c("bst", "co") && dfun==5)
        tmp <- init(RET$weights_update/sum(RET$weights_update), ytmp, offset, family=dfunnew)
    else if(type.init%in%c("bst","co")) ##it seems we don't need tmp next line
        tmp <- init(RET$weights_update, ytmp, offset, family=dfunnew)
    mustart <- tmp$mu
    etastart <- tmp$eta     
### initialize the start values -end
    stopit <- FALSE
    ### create a function for parallel computing can lead to object-not-found problem
    #typeP <- function(x, y, weights, start, etastart, mustart, offset,lambda,alpha,gamma,intercept,penalty.factor,maxit,eps, pentype, trace, iter, reltol, cfun, dfun, s, thresh, delta, m){
    if(isTRUE(parallel)){
    i <- 1
    cl <- parallel::makeCluster(n.cores, outfile="")
         registerDoParallel(cl)
         fitall <- foreach(i=1:nlambda, .packages=c("mpath")) %dopar%{
            RET <- .Fortran("ccglmreg_onelambda",
                            x_act=as.double(x), 
                            y=as.double(y), 
                            weights=as.double(weights),
                            n=as.integer(n),
                            m_act=as.integer(m),	
                            start_act=as.double(start), 
                            etastart=as.double(etastart),
                            mustart=as.double(mustart),
                            yhat=as.double(rep(0, n)),
                            offset=as.double(offset),
                            lambda_i=as.double(lambda[i]), 
                            alpha=as.double(alpha),        
                            gam=as.double(gamma), 
			                rescale=as.integer(rescale),
                            #standardize=as.integer(0),
                            standardize=as.integer(standardize),
                            intercept=as.integer(intercept),
                            penaltyfactor_act=as.double(penalty.factor),
                            maxit=as.integer(maxit),
                            eps=as.double(eps), 
                            theta=as.double(theta),
                            penalty=as.integer(pentype), 
                            trace=as.integer(trace),
                            iter=as.integer(iter),
                            del=as.double(reltol),
                            cfun=as.integer(cfun),
                            dfun=as.integer(dfun),
                            s=as.double(s),
                            thresh=as.double(thresh), 
                            beta_1=as.double(rep(0, m)),
                            b0_1=as.double(0),
                            fk=as.double(rep(0, n)),
                            delta=as.double(delta),
                            weights_update=as.double(rep(0, n)),
                            PACKAGE="mpath")
            list(beta=RET$beta, b0=RET$b0, yhat=RET$yhat, weights_update=RET$weights_update)   
        }
    parallel::stopCluster(cl)
    RET <- fitall[[nlambda]]    
    for(k in 1:nlambda){
             beta[,k] <- fitall[[k]]$beta
             b0[k] <- fitall[[k]]$b0
             fitted[,k] <- fitall[[k]]$yhat
             weights_cc[,k] <- fitall[[k]]$weights_update
         }
        tmp <- list(beta=beta, b0=b0, RET=RET, fitted=fitted, weights_cc=weights_cc)
    }
    typeA <- function(beta, b0){
        if(dfun %in% c(1, 4)) dfuntmp <- 1 else if(dfun==5) dfuntmp <- 2
        else if(dfun==8) dfuntmp <- 3 else if(dfun==9) dfuntmp <- 4
        else if(dfun==6)
                                        if(all(x[,1]==1)) xtmp <- x[,-1] else xtmp <- x
        i <- 1
        los <- pll <- matrix(NA, nrow=iter, ncol=nlambda)
        weights_cc <- matrix(NA, nrow=n, ncol=nlambda)
        if(trace && tracelevel==2) tracel <- 1 else tracel <- 0
        while(i <= nlambda){
            if(trace) message("\nloop in lambda:", i, ", lambda=", lambda[i], "\n")
             if(trace) {
                 cat(" COCO iterations ...\n")
             }
            k <- 1
            d1 <- 10
            weights_update <- weights
            satu <- 0
            while(d1 > reltol && k <= iter && satu==0){
                fitted.values <- fk
                RET <- .Fortran("glmreg_fit_fortran",
	                        x=as.double(x), 
				y=as.double(ytmp), 
				weights=as.double(weights_update),
			        n=as.integer(n),
			        m=as.integer(m),	
				start=as.double(start), 
				etastart=as.double(etastart),
	                        mustart=as.double(mustart),
			       	offset=as.double(offset),
			       	nlambda=as.integer(1),
			       	lambda=as.double(lambda[i]), 
				alpha=as.double(alpha),        
			        gam=as.double(gamma), 
				rescale=as.integer(rescale),
			    #standardize=as.integer(0), ### changed 08/17/2020
			   	standardize=as.integer(standardize),
			       	intercept=as.integer(intercept),
				penaltyfactor=as.double(penalty.factor),
				thresh=as.double(thresh), 
				epsbino=as.double(0), 
				maxit=as.integer(maxit),
			        eps=as.double(eps), 
				theta=as.double(theta),
			       	family=as.integer(dfuntmp),
				penalty=as.integer(pentype), 
				trace=as.integer(tracel),
				beta=as.double(matrix(0, ncol=1, nrow=m)),
				b0=as.double(0),
				yhat=as.double(rep(0, n)),
                                satu=as.integer(0),
				PACKAGE="mpath")
                satu <- RET$satu
                if(satu==1) warnings("saturated binomial model")
                fk <- RET$yhat
                etastart <- RET$etastart
                mustart <- RET$mustart
                if(dfun%in%c(1,4,5)){
                weights_update <- .Fortran("update_wt",
                                           n=as.integer(n),
                                           weights=as.double(weights),
                                           y=as.double(y),
                                           f=as.double(etastart),
                                           cfun=as.integer(cfun),
                                           dfun=as.integer(dfun),
                                           s=as.double(s),
                                           delta=as.double(delta),
                                           weights_update=as.double(rep(0, n)),
                                           PACKAGE="mpath")$weights_update
                los[k, i] <- .Fortran("loss2",
                                            n=as.integer(n),
                                            y=as.double(y),
                                            f=as.double(etastart),
                                            weights=as.double(weights),
                                            cfun=as.integer(cfun),
                                            dfun=as.integer(dfun),
                                            s=as.double(s),
                                            delta=as.double(delta),
                                            los=as.double(0.0),
                                            PACKAGE="mpath")$los
                }
                else{
                weights_update <- compute_wt3(y,mustart,theta,weights,cfun,family,s,delta) #for poisson/negbin
                los[k, i] <- sum(loss3(y,mustart,theta,weights,cfun,family,s,delta)$tmp)
                } 
                ### changed the following line 08/17/2020
                if(dfun!=5) start <- c(RET$b0, RET$beta)
                #if(!dfun %in% c(5, 8, 9)) start <- c(RET$b0, RET$beta)
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
                 if(standardize)  ### lambda value, hence penval, depends on whether standardize is TRUE/FALSE - to be confirmed
                     pll[k, i] <- los[k, i] + n*penval
                 else pll[k, i] <- los[k, i] + penval
                 if(k > 1)
                 d1 <- abs((pll[k, i] - pll[k-1, i])/pll[k-1, i])
                 if(trace) cat("\n  iteration", k, ": relative change of fk", d1, ", robust loss value", los[k, i], ", penalized loss value", pll[k, i], "\n")
                 if(trace) cat("  d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
                k <- k + 1
            }
                beta[,i] <- RET$beta
                b0[i] <- RET$b0
                weights_cc[,i] <- weights_update
                i <- i + 1
        }
        list(beta=beta, b0=b0, RET=RET, risk=los, pll=pll, weights_cc=weights_cc)
    }
    typeBB <- function(beta, b0){
        if(trace && tracelevel==2) tracel <- 1 else tracel <- 0
        if(type.path=="active" && decreasing)
            RET <- .Fortran("ccglmreg_ad",
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
			    rescale=as.integer(rescale),
			    standardize=as.integer(standardize), 
			    #standardize=as.integer(0), ### changed 8/17/2020
			    intercept=as.integer(intercept), 
			    penaltyfactor=as.double(penalty.factor),
			    maxit=as.integer(maxit), 
			    eps=as.double(eps), 
			    theta=as.double(theta), 
			    epscycle=as.double(epscycle), 
			    penalty=as.integer(pentype), 
			    trace=as.integer(tracel), 
			    del=as.double(reltol), 
			    cfun=as.integer(cfun),
			    dfun=as.integer(dfun),
			    s=as.double(s),
			    thresh=as.double(thresh),
			    decreasing=as.integer(decreasing),
			    beta=as.double(beta),
                            b0=as.double(b0), 
			    yhat=as.double(RET$fitted.values), 
			    los=as.double(rep(0, nlambda)),
			    pll=as.double(rep(0, nlambda)),
			    nlambdacal=as.integer(0),
			    delta=as.double(delta),
			    weights_cc=as.double(matrix(0, nrow=n, ncol=nlambda)),
			    PACKAGE="mpath")
        else
	    RET <- .Fortran("ccglmreg_fortran",
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
			    rescale=as.integer(rescale),
			    standardize=as.integer(standardize), 
			    #standardize=as.integer(0), #changed 8/17/2020
			    intercept=as.integer(intercept), 
			    penaltyfactor=as.double(penalty.factor),
			    maxit=as.integer(maxit), 
			    eps=as.double(eps), 
			    theta=as.double(theta), 
			    penalty=as.integer(pentype), 
			    trace=as.integer(tracel), 
			    del=as.double(reltol), 
			    cfun=as.integer(cfun),
			    dfun=as.integer(dfun),
			    s=as.double(s),
			    thresh=as.double(thresh),
			    decreasing=as.integer(decreasing),
			    active=as.integer(active),
			    beta=as.double(beta),
                            b0=as.double(b0), 
			    yhat=as.double(RET$fitted.values), 
			    los=as.double(rep(0, nlambda)),
			    pll=as.double(rep(0, nlambda)),
			    nlambdacal=as.integer(0),
			    delta=as.double(delta),
			    weights_cc=as.double(matrix(0, nrow=n, ncol=nlambda)),
			    PACKAGE="mpath")
        list(beta=matrix(RET$beta, ncol=nlambda), b0=RET$b0, RET=RET, risk=RET$los, pll=RET$pll, nlambdacal=RET$nlambdacal, weights_cc=RET$weights_cc)
    }
        #tmp <- typeP(x, y, weights, start, etastart, mustart, offset,lambda[i],alpha,gamma,intercept,penalty.factor,maxit,eps, pentype, trace, iter, reltol, cfun, dfun, s, thresh, delta, m) 
    if(is.null(parallel)) tmp <- typeA(beta, b0) else 
    if(!parallel) tmp <- typeBB(beta, b0)
    beta <- tmp$beta
    b0 <- tmp$b0
    RET <- tmp$RET
    RET$weights_update <- tmp$weights_cc
    if(!is.null(parallel) && standardize && dfun%in%c(5,8,9)){
        tmp1 <- stan(x, weights)
        meanx <- tmp1$meanx
        normx <- tmp1$normx
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
    #RET$risk <- tmp$risk
    #RET$pll <- tmp$pll
    RET$nlambdacal <- tmp$nlambdacal
    RET$cfun <- cfunold
    RET$dfun <- dfunold
    RET$type.init <- type.init
    RET$mstop.init <- mstop.init
    RET$nu.init <- nu.init
    RET$decreasing <- decreasing
    RET$type.path <- type.path
    RET$is.offset <- is.offset
    RET$fitted.values <- fitted

    class(RET) <- "ccglmreg"
    RET
}

### compute classification/prediction error based on predictive values
co_evalerr <- function(dfun, y, yhat){
    if(dfun %in% c(1:3,8,9))
        mean((y - yhat)^2)
    else if(dfun %in% 4:7)
        (mean(y != sign(yhat)))
}

predict.ccglmreg <- function(object, newdata=NULL, weights=NULL, newy=NULL, newoffset=NULL, which=1:object$nlambda, type=c("link", "response","class","loss", "error", "coefficients", "nonzero"), na.action=na.pass, ...){
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
        #if(length(which)>1) return(nonzeroCoef(nbeta[,,drop=FALSE],bystep=TRUE))
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
    if(type=="link") return(res[, which])
    if(type %in% c("link", "response")) return(predict(object, newx=newx, newoffset=newoffset, which=which, type=type))
    if(type %in% c("loss", "error") && is.null(ynow)) stop("response variable y missing\n")
    if(type=="loss"){
        object$cfun <- cfun2num(object$cfun)
        object$dfun <- dfun2num(object$dfun)
        if(object$dfun %in% 4:5) ynow <- y2num(ynow)
        tmp <- rep(NA, length(which))
        n <- length(ynow)
        if(missing(weights)) weights <- rep(1/n, n)
        for(i in 1:length(which)){
            if(object$dfun %in% c(1,4,5)) 
                tmp[i] <- .Fortran("loss2",
                               n=as.integer(n),
                               y=as.double(ynow),
                               f=as.double(res[,which[i]]),
                               weights=as.double(weights),
                               cfun=as.integer(object$cfun),
                               dfun=as.integer(object$dfun),
                               s=as.double(object$s),
                               delta=as.double(object$delta),
                               los=as.double(0.0),
                               PACKAGE="mpath")$los
        else{
          mu <- predict(object, newx=newx, newoffset=newoffset,which=which[i], type="response")
          tmp[i] <- sum(loss3(object$y,mu,object$theta,weights,object$cfun,family,object$s,object$delta)$tmp)
            }
        }
        return(tmp)
    }
    if(type=="error"){
        object$dfun <- dfun2num(object$dfun)
        if(object$dfun %in% 4:5) ynow <- y2num(ynow)
        tmp <- rep(NA, length(which))
        for(i in 1:length(which))
            tmp[i] <- co_evalerr(object$dfun, ynow, res[,which[i]])
        return(tmp)
    }
}


coef.ccglmreg <- function(object, ...)
    coef.glmreg(object)

