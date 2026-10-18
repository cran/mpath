cv.nclreg <- function(x, ...) UseMethod("cv.nclreg", x)

cv.nclreg.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(cv.nclreg.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

cv.nclreg.formula <- function(formula, data, weights, offset=NULL, ...){
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

    RET <- cv.nclreg_fit(X[,-1], Y, weights,...)
    RET$call <- match.call()
    return(RET)
}
cv.nclreg.matrix <- function(x, y, weights, offset=NULL, ...){
    RET <- cv.nclreg_fit(x, y, weights, offset=offset, ...)
    RET$call <- match.call()
    return(RET)
}

cv.nclreg_fit <- function(x, y, weights, offset, lambda=NULL, balance=TRUE, 
                          rfamily=c("clossR", "closs", "gloss", "qloss"), s=1.5,  
                          nfolds=10, foldid, type = c("loss", "error"), plot.it=TRUE, se=TRUE, n.cores=NULL, trace=FALSE, parallel=FALSE,  
                          ...){
    call <- match.call()
    type <- match.arg(type)
    if(missing(foldid) && nfolds < 3)
        stop("smallest nfolds should be 3\n")
    rfamily <- match.arg(rfamily)
    nm <- dim(x)
    nobs <- n <- nm[1]
    nvars <- m <- nm[2]
    if(missing(weights)) weights <- rep(1, nobs)
    K <- nfolds
    nclreg.obj <- nclreg_fit(x, y, weights, offset=offset, lambda=lambda, rfamily=rfamily, s=s, ...)
    lambda <- nclreg.obj$lambda
    nlambda <- length(lambda)
    if(missing(foldid)){
        if(rfamily %in% c("closs", "gloss", "qloss") && balance)  
            invisible(capture.output(all.folds <- eval(parse(text="balanced.folds(y, K)"))))
        else all.folds <- cv.folds(length(y), K)
    }
    else all.folds <- foldid
    foldenv <- list2env(list(x=x, y=y, weights=weights, offset=offset,
                             lambda=lambda, rfamily=rfamily, s=s,
                             type=type, all.folds=all.folds,
                             dots=list(...)),
                        parent=asNamespace("mpath"))
    foldwork <- function(i) {
        omit <- all.folds[[i]]
        if(!is.null(offset)){
            offsetnow <- offset[-omit]
            newoffset <- offset[omit]
        }else offsetnow <- newoffset <- NULL
        fitcv <- do.call(nclreg_fit, c(list(x=x[-omit,,drop=FALSE],
                                            y=y[-omit],
                                            weights=weights[-omit],
                                            offset=offsetnow, s=s,
                                            lambda=lambda,
                                            rfamily=rfamily), dots))
        predict(fitcv, newdata=x[omit,,drop=FALSE], newy=y[omit],
                weights=weights[omit], newoffset=newoffset, type=type)
    }
    environment(foldwork) <- foldenv

    if(parallel){
        dc <- parallel::detectCores(logical=TRUE)
        if(is.na(dc)) dc <- 1L
        if(is.null(n.cores)) n.cores <- dc
        n.cores <- max(1L, min(as.integer(n.cores), K, dc))
        if(.Platform$OS.type == "unix"){
            fold.results <- parallel::mclapply(seq_len(K), foldwork,
                                               mc.cores=n.cores)
            err <- sapply(fold.results, inherits, "try-error")
            if(any(err)) stop(fold.results[[which(err)[1]]])
        } else {
            cl <- parallel::makeCluster(n.cores)
            on.exit(parallel::stopCluster(cl), add=TRUE)
            fold.results <- parallel::parLapply(cl, seq_len(K), foldwork)
        }
    } else {
        fold.results <- vector("list", K)
        for(i in seq_len(K)) {
            if(trace) cat("\n CV Fold", i, "\n\n")
            fold.results[[i]] <- foldwork(i)
        }
    }
    residmat <- matrix(NA_real_, nlambda, K)
    for(i in seq_len(K)) {
        res <- fold.results[[i]]
        if(!is.null(res)) residmat[seq_along(res), i] <- res
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.min(cv)
    obj<-list(fit=nclreg.obj, residmat=residmat, lambda=lambda, cv = cv, cv.error = cv.error, foldid=all.folds, lambda.which= lambda.which, lambda.optim = lambda[lambda.which])
    if(plot.it) plot.cv.nclreg(obj,se=se, ylab=type)
    class(obj) <- "cv.nclreg"
    obj
}

"plot.cv.nclreg" <-
    function(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...){
        lambda <- x$lambda
	cv <- x$cv
        cv.error <- x$cv.error
	plot(log(lambda), cv, type = "b", xlab = expression(log(lambda)), ylab= ylab, ylim = range(cv, cv + cv.error, cv - cv.error), main=main)
        if(se)
            error.bars(log(lambda), cv + cv.error, cv - cv.error,
                       width = width, col=col)

        invisible()
    }

coef.cv.nclreg=function(object,which=object$lambda.which,...){
    coef(object$fit,which=which,...)
}
