### SVM based on composite operator
cv.irsvm <- function(x, ...) UseMethod("cv.irsvm")

cv.irsvm.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(cv.irsvm.matrix(x = x, ...))
    if (inherits(x, "data.frame"))
        return(cv.irsvm.matrix(x = x, ...))
    if (inherits(x, "numeric") && is.null(dim(x)))
        return(cv.irsvm.matrix(x = as.matrix(x), ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

cv.irsvm.formula <- function(formula, data, weights, contrasts=NULL, ...){
    ## extract x, y, etc from the model formula and frame
    if(!attr(terms(formula, data=data), "intercept"))
        stop("non-intercept model is not implemented with model formula, try model matrix instead")
    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
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
    RET <- cv.irsvm_fit(X[,-1], Y, weights,...)
    RET$call <- match.call()
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    RET
}
cv.irsvm.matrix <- function(x, y, weights, ...){
    RET <- cv.irsvm_fit(x, y, weights,...)
    RET$call <- match.call()
    return(RET)
}

cv.irsvm_fit <- function(x,y, weights, cfun="ccave", s=c(1, 5), type=NULL, 
                         kernel      = "radial",
                                        #degree      = 3,
                         gamma       = 2^(-4:10),
                                        #coef0       = 0,
                         cost = 2^(-4:4),
                                        #class.weights = NULL,
                                        #cachesize   = 100,
                                        #tolerance   = 0.001,
                         epsilon     = 0.1,
                                        #shrinking   = TRUE,
                                        #probability = FALSE,
                                        #fitted      = TRUE,
                                        #subset,
                                        #na.action,
                         balance=TRUE,
                         nfolds=10, foldid, trim_ratio=0.9, n.cores=NULL,
                         ...){
    call <- match.call()
    ## determine model type
          if (is.null(type)) type <-
            if (is.factor(y)) "C-classification"
          else "eps-regression"
    sval <- s
    costseq <- cost
    if(kernel=="linear") gammaseq <- 0 else gammaseq <- gamma
    if(is.null(dim(x))) x <- as.matrix(x) 
    nm <- dim(x)
    nobs <- n <- nm[1]
    nvars <- m <- nm[2]
    ## avoid any problems with 1D or nx1 arrays by as.vector.
     if(missing(weights)) weights=rep(1,n)
     weights <- as.vector(weights)
     if(!is.null(weights) && !is.numeric(weights))
         stop("'weights' must be a numeric vector")
     ## check weights 
     if( !is.null(weights) && any(weights < 0) ){
         stop("negative weights not allowed")
     }
     if(length(weights) != length(y))
         stop("'weights' must be the same length as response variable")
    K <- nfolds
    if(missing(foldid)){
        if(type %in% c("C-classification", "nu-classification") && balance)
         invisible(capture.output(all.folds <- eval(parse(text="balanced.folds(y, K)"))))
        else all.folds <- cv.folds(n, K)
    }
    else all.folds <- foldid
    param.grid <- expand.grid(cost=costseq, gamma=gammaseq, s=sval)
    param.grid <- param.grid[c("s", "gamma", "cost")]
    ngrid <- nrow(param.grid)
    if(!ngrid) stop("no parameters are valid for SVM")
    foldenv <- list2env(list(x=x, y=y, weights=weights, cfun=cfun,
                             type=type, trim_ratio=trim_ratio,
                             all.folds=all.folds, param.grid=param.grid,
                             ngrid=ngrid, dots=list(...)),
                        parent=asNamespace("mpath"))
    foldwork <- function(job) {
        i <- (job - 1L) %/% ngrid + 1L
        grid.i <- (job - 1L) %% ngrid + 1L
        omit <- all.folds[[i]]
        params <- param.grid[grid.i,]
        fitcv <- try(do.call(irsvm, c(list(
            x=x[-omit,,drop=FALSE], y=y[-omit],
            weights=weights[-omit], cfun=cfun,
            cost=params$cost, gamma=params$gamma, s=params$s,
            type=type), dots)), silent=TRUE)
        if(inherits(fitcv, "try-error")) return(rep(NA_real_, 4L))
        pred <- predict(fitcv, x[omit,,drop=FALSE], weights=weights[omit])
        if(type %in% c("C-classification", "nu-classification")) {
            error <- mean(pred != y[omit])
        } else {
            loss <- (y[omit] - pred)^2
            loss.order <- order(loss, decreasing=FALSE)
            trimmed <- loss.order[1:(length(omit) * trim_ratio)]
            error <- mean(loss[trimmed])
        }
        c(params$s, params$gamma, params$cost, error)
    }
    environment(foldwork) <- foldenv

    jobs <- seq_len(K * ngrid)
    dc <- parallel::detectCores(logical=TRUE)
    if(is.na(dc)) dc <- 1L
    if(is.null(n.cores)) n.cores <- dc
    n.cores <- max(1L, min(as.integer(n.cores), K, dc))
    if(.Platform$OS.type == "unix"){
        fold.results <- parallel::mclapply(jobs, foldwork,
                                           mc.cores=n.cores)
        err <- sapply(fold.results, inherits, "try-error")
        if(any(err)) stop(fold.results[[which(err)[1]]])
    } else {
        cl <- parallel::makeCluster(n.cores)
        on.exit(parallel::stopCluster(cl), add=TRUE)
        fold.results <- parallel::parLapply(cl, jobs, foldwork)
    }
    residmat <- do.call(cbind, lapply(seq_len(K), function(i) {
        first <- (i - 1L) * ngrid + 1L
        last <- i * ngrid
        do.call(rbind, fold.results[first:last])
    }))
    kseq <- 4*(1:K)
    valid <- rowSums(!is.na(residmat[, kseq, drop=FALSE])) == K
    if(!any(valid)) stop("no parameters are valid for SVM")
    residmat <- residmat[valid,,drop=FALSE]
    err <- apply(residmat[,kseq,drop=FALSE], 1, mean) #average error cross K folds
    opt <- residmat[which.min(err), 1:3]
    if(kernel=="linear"){
        residmat <- residmat[, -2*(1:K)]
        colnames(residmat) <- rep(c("s", "cost", "error"), K)
    }else                                                                       colnames(residmat) <- rep(c("s", "gamma", "cost", "error"), K)
    RET <- list(residmat=residmat, cost=opt[3], gamma=opt[2], s=opt[1])
    RET
}
