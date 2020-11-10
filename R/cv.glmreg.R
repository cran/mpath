"cv.folds" <-
    function(n, folds = 10)
{
    split(sample(1:n), rep(1:folds, length = n))
}

cv.glmreg <- function(x, ...) UseMethod("cv.glmreg", x)

cv.glmreg.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(cv.glmreg.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

cv.glmreg.formula <- function(formula, data, weights, offset=NULL, contrasts=NULL, ...){
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
### End of addition 08/07/2012 

    RET <- cv.glmreg_fit(X[,-1], Y, weights, offset=offset, ...)
    RET$call <- match.call()
    return(RET)
}
cv.glmreg.matrix <- function(x, y, weights, offset=NULL, ...){
    RET <- cv.glmreg_fit(x, y, weights, offset=offset, ...)
    RET$call <- match.call()
    return(RET)
}

cv.glmreg_fit <- function(x, y, weights, offset, lambda=NULL, balance=TRUE, 
                          family=c("gaussian", "binomial", "poisson", "negbin"), type=c("loss", "error"), 
                          nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, trace=FALSE, parallel=FALSE,  
                          ...){
    call <- match.call()
    if(missing(foldid) && nfolds < 3)
        stop("smallest nfolds should be 3\n")
    family <- match.arg(family)
    type <- match.arg(type)
    nm <- dim(x)
    nobs <- n <- nm[1]
    nvars <- m <- nm[2]
    if(missing(weights)) weights <- rep(1, nobs)
                                        #    if(is.null(offset)) offset <- rep(0, nobs)
    K <- nfolds
    glmreg.obj <- glmreg_fit(x, y, weights, offset=offset, lambda=lambda, family=family, ...)
    lambda <- glmreg.obj$lambda
    nlambda <- length(lambda)
    if(missing(foldid)){
        if(family=="binomial" && balance)  
            invisible(capture.output(all.folds <- eval(parse(text="pamr:::balanced.folds(y, K)"))))
        else all.folds <- cv.folds(length(y), K)
    }
    else all.folds <- foldid
    if(parallel){
        cl <- eval(parse(text="parallel:::makeCluster(n.cores)"))
        registerDoParallel(cl)
        i <- 1  ###needed to pass R CMD check with parallel code below
                                        #residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
                                        #residmat should be a list since it is possible to generate different length of loss values, thus, cbind may fail
        residmat <- foreach(i=seq(K), .packages=c("mpath")) %dopar% {
            omit <- all.folds[[i]]
            if(!is.null(offset)){
                offsetnow <- offset[- omit] 
                newoffset <- offset[omit] 
            }else offsetnow <- newoffset <- NULL 
            fitcv <- glmreg_fit(x[ - omit,,drop=FALSE ], y[ -omit], weights=weights[- omit], offset=offsetnow, lambda=lambda, family=family, ...)
            if(type=="loss") logLik(fitcv, newx=x[omit,, drop=FALSE], y[omit], weights=weights[omit], offset=newoffset)
            else if(family=="binomial" && type=="error"){
                tmp <- predict(fitcv, newx=x[omit,, drop=FALSE], weights=weights[omit], offset=newoffset, type="class")!=y[omit]
                apply(tmp, 2, mean)
            }
        }
        eval(parse(text="parallel:::stopCluster(cl)"))
        residmat <- sapply(residmat, '[', seq(max(sapply(residmat, length))))
        }
        else{
            residmat <- matrix(NA, nlambda, K)
            for(i in seq(K)) {
                if(trace)
                    cat("\n CV Fold", i, "\n\n")
                omit <- all.folds[[i]]
                if(!is.null(offset)){
                    offsetnow <- offset[- omit] 
                    newoffset <- offset[omit] 
                }else offsetnow <- newoffset <- NULL 
                fitcv <- glmreg_fit(x[ - omit,,drop=FALSE ], y[ -omit], weights=weights[- omit], offset=offsetnow, lambda=lambda, family=family, ...)
                if(type=="loss") residmat[1:fitcv$nlambda, i] <- logLik(fitcv, newx=x[omit,, drop=FALSE], y[omit], weights=weights[omit], offset=newoffset)
                else if(family=="binomial" && type=="error"){
                tmp <- predict(fitcv, newx=x[omit,, drop=FALSE], weights=weights[omit], offset=newoffset, type="class")!=y[omit]
                residmat[1:fitcv$nlambda, i] <- apply(tmp, 2, mean)
                }
            }
        }
        cv <- apply(residmat, 1, mean)
        cv.error <- sqrt(apply(residmat, 1, var)/K)
        if(type=="loss") lambda.which <- which.max(cv) else
            if(type=="error" && family=="binomial") lambda.which <- which.min(cv)
        if(family=="binomial")
            good <- which(!is.na(cv))
        else
            good <- 1:nlambda
        obj<-list(fit=glmreg.obj, residmat=residmat[good,], lambda = lambda[good], cv = cv[good], cv.error = cv.error[good], foldid=all.folds, lambda.which= lambda.which, lambda.optim = lambda[lambda.which], type=type)
        if(plot.it) plot.cv.glmreg(obj,se=se)
        class(obj) <- "cv.glmreg"
        obj
    }

    "plot.cv.glmreg" <-
        function(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...){
            lambda <- x$lambda
            cv <- x$cv
            cv.error <- x$cv.error
            if(is.null(ylab) && x$fit$family=="binomial" && x$type=="error")
                ylab <- "misclassification error"
            else if(is.null(ylab))
                ylab <- "log-likelihood"
            if(all(lambda > 0)){
            plot(log(lambda), cv, type = "b", xlab = expression(log(lambda)), ylab= ylab, ylim = range(cv, cv + cv.error, cv - cv.error), main=main)
            if(se)
                error.bars(log(lambda), cv + cv.error, cv - cv.error,
                           width = width, col=col)
            }else{
            plot(lambda, cv, type = "b", xlab = expression(lambda), ylab= ylab, ylim = range(cv, cv + cv.error, cv - cv.error), main=main)
            if(se)
                error.bars(lambda, cv + cv.error, cv - cv.error,
                           width = width, col=col)
            }

            invisible()
        }

    predict.cv.glmreg=function(object, newx, ...){
      	predict(object$fit,newx,which=object$lambda.which,...)
    }

    coef.cv.glmreg=function(object,which=object$lambda.which,...){
        coef(object$fit,which=which,...)
    }
