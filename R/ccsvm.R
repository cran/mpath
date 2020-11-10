### SVM based on composite operator
ccsvm <- function(x, ...) UseMethod("ccsvm")

ccsvm.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(ccsvm.matrix(x = x, ...))
    if (class(x)=="data.frame")
        return(ccsvm.matrix(x = x, ...))
    if (class(x)=="numeric" && is.null(dim(x)))
        return(ccsvm.matrix(x = as.matrix(x), ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

ccsvm.formula <- function(formula, data, weights, contrasts=NULL, ...){
    ## extract x, y, etc from the model formula and frame
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
    RET <- ccsvm_fit(X, Y, weights,...)
    RET$call <- match.call()
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    class(RET) <- "wsvm"
    RET
}
ccsvm.matrix <- function(x, y, weights, ...){
    RET <- ccsvm_fit(x, y, weights,...)
    RET$call <- match.call()
    return(RET)
}

ccsvm_fit <- function(x,y, weights, cfun="ccave", s=NULL, delta=0.0001, type=NULL, 
               kernel      = "radial",
               #degree      = 3,
               #gamma       = if (is.vector(x)) 1 else 1 / ncol(x),
               #coef0       = 0,
               cost        = 1,
               #class.weights = NULL,
               #cachesize   = 100,
               #tolerance   = 0.001,
               epsilon     = 0.1,
               #shrinking   = TRUE,
               #probability = FALSE,
               #fitted      = TRUE,
               #subset,
               #na.action,
               iter=10, reltol=1e-5, trace=FALSE,
               ...){
    call <- match.call()
     ## determine model type
         if (is.null(type)) type <-
           if (is.factor(y)) "C-classification"
         else "eps-regression"
    cfunold <- cfun
    cfun <- cfun2num(cfun)
    dfun <- switch(type,
                   "C-classification"=6,
                   "nu-classification"=6,
                   "eps-regression"=2,
                   "nu-regression"=2)
    if(dfun==6){
           if (is.factor(y)) {
                 y <- as.integer(y)
             } else {
                     if(any(as.integer(y) != y))
                         stop("y variable has to be of factor or         integer type for classification\n")
                     y <- as.factor(y)
                     y <- as.integer(y)
             }
    if(length(unique(y))!=2) stop("y must be a binary variable for classification\n")
    else if(!all(names(table(y)) %in% c(1, -1))){
        y[y==1] <- -1
        y[y==2] <- 1
        } 
    }
    if(is.null(s)) s <- assign_s(cfun, y) else check_s(cfun, s)
    if(cfun==6)
           if(s > 1) delta <- (s-1)/2 else
           if(s==1) delta <- 0 else{
             if(is.null(delta)) stop("delta must be provided")
             if(delta <= 0) stop("delta must be positive")
           }
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    if(missing(weights)) weights=rep(1,length(y))
    weights <- as.vector(weights)
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    ## check weights 
    if( !is.null(weights) && any(weights < 0) ){
        stop("negative weights not allowed")
    }
    if(length(weights) != length(y))
        stop("'weights' must be the same length as response variable")
    n <- length(y)
    if(all(x[,1]==1)) xtmp <- x[,-1] else xtmp <- x
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\niterations ...\n")
    }
    los <- rep(NA, iter)
    weights_update <- weights
    subset <- 1:n
    while(d1 > reltol && k <= iter){
        RET <- wsvm(xtmp, y, weights*weights_update, type=type, kernel=kernel, cost=cost, epsilon=epsilon, cross=0, subset=subset, ...)
        pred <- predict(RET, xtmp, decision.values=TRUE)
        if(dfun==6){
        fitted.values <- attr(pred, "decision.values")
         if(y[subset[1]]==-1)
            fitted.values <- -fitted.values
        }else fitted.values <- pred
        weights_update <- .Fortran("update_wt_ccsvm",
                                   n=as.integer(n),
                                   weights=as.double(weights),
                                   y=as.double(y),
                                   f=as.double(fitted.values),
                                   cfun=as.integer(cfun),
                                   dfun=as.integer(dfun),
                                   s=as.double(s),
                                   eps=as.double(epsilon),
                                   delta=as.double(delta),
                                   weights_update=as.double(rep(0, n)),
                                   PACKAGE="mpath")$weights_update
        subset <- which(weights_update > 0)
        if(trace && kernel=="linear"){
            beta <- t(RET$coefs) %*% as.matrix(xtmp[RET$index,]) #kernel="linear" only: convert dual coef to RET$coefs to primal coef
                                        #RET$b0 <- -RET$rho #intercept
###for some reason, depending on the sign of the first element of y, in wsvm, version 1.7-5 or svm version 1.7.3, beta and b0 are minus corresponding values coef_ and intercetp_ in Python scikit.learn.svm.SVC version xxx. Likewise, fitted values have inverse sign
        ### furthermore, it is questionable on the beta value from wsvm sometimes, but fitted values seem ok.
            if(dfun==6 && y[1]==-1)
                beta <- -beta
            los[k] <- .Fortran("loss2_ccsvm",
                               n=as.integer(n),
                               y=as.double(y),
                               f=as.double(fitted.values),
                               weights=as.double(weights),
                               cfun=as.integer(cfun),
                               dfun=as.integer(dfun),
                               s=as.double(s),
                               eps=as.double(epsilon),
                               delta=as.double(delta),
                               los=as.double(0.0),
                               PACKAGE="mpath")$los
            los[k] <- los[k]+0.5/cost*sum(beta^2) ### only for kernel="linear"
            cat("\niteration", k, ", robust loss value", los[k], "\n") 
            if(k > 1){
                d1 <- abs((los[k]-los[k-1])/los[k-1])
                if(los[k] > los[k-1]){
                    k <- iter
                }
            }
        }else{
            if(k > 1){
                d1 <- sum((fitold-fitted.values)^2)/sum(fitted.values^2)
            }
            fitold <- fitted.values
        }
        k <- k + 1
	if(trace) cat("d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
    }
    RET$cfun <- cfunold
    RET$s <- s
    RET$delta <- delta
    RET$weights_update <- weights_update
    RET
}

