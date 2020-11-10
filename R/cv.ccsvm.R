### SVM based on composite operator
cv.ccsvm <- function(x, ...) UseMethod("cv.ccsvm")

cv.ccsvm.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(cv.ccsvm.matrix(x = x, ...))
    if (class(x)=="data.frame")
        return(cv.ccsvm.matrix(x = x, ...))
    if (class(x)=="numeric" && is.null(dim(x)))
        return(cv.ccsvm.matrix(x = as.matrix(x), ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

cv.ccsvm.formula <- function(formula, data, weights, contrasts=NULL, ...){
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
    RET <- cv.ccsvm_fit(X[,-1], Y, weights,...)
    RET$call <- match.call()
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    RET
}
cv.ccsvm.matrix <- function(x, y, weights, ...){
    RET <- cv.ccsvm_fit(x, y, weights,...)
    RET$call <- match.call()
    return(RET)
}

cv.ccsvm_fit <- function(x,y, weights, cfun="ccave", s=c(1, 5), type=NULL, 
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
                         nfolds=10, foldid, trim_ratio=0.9, n.cores=2,
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
         invisible(capture.output(all.folds <- eval(parse(text="pamr:::balanced.folds(y, K)"))))
        else all.folds <- cv.folds(n, K)
    }
    else all.folds <- foldid
    cl <- eval(parse(text="parallel:::makeCluster(n.cores)")) 
    registerDoParallel(cl)
    i <- ii <- jj <- l <- 1
    residmat <- foreach(i=seq(K), .combine=cbind)%:% 
        foreach(l=sval, .combine='rbind')%:%
        foreach(ii=gammaseq, .combine='rbind')%:%
        foreach(jj=costseq, .combine='rbind', .packages="mpath")   %dopar%{
    #residmat <- NULL; for(i in seq(K)) for(l in sval) for(ii in gammaseq) for(jj in costseq){
            omit <- all.folds[[i]]
            m <- try(ccsvm(x[-omit,,drop=FALSE], y[-omit], weights=weights[-omit], cfun=cfun, cost=jj, gamma=ii, s=l, type=type, ...), silent=TRUE)
            if(!inherits(m, "try-error")){
                new <- predict(m, x[omit,,drop=FALSE], weights=weights[omit])
                if(type%in%c("C-classification", "nu-classification")) errorTU <- mean(new!=y[omit])
                else {
                    mloss <- (y[omit]-new)^2
                    mloss1 <- order(mloss, decreasing=FALSE)
                    mloss2 <- mloss1[1:(length(omit)*trim_ratio)]
                    errorTU <- mean(mloss[mloss2]) #trimmed MSE
                }
                c(l, ii, jj, errorTU)
                #residmat <- rbind(residmat, c(l, ii, jj, errorTU))
            }
        }
    eval(parse(text="parallel:::stopCluster(cl)")) 
    kseq <- 4*(1:K)
    err <- apply(residmat[,kseq], 1, mean) #average error cross K folds
    opt <- residmat[which.min(err), 1:3]
    if(kernel=="linear"){
        residmat <- residmat[, -2*(1:K)]
        colnames(residmat) <- rep(c("s", "cost", "error"), K)
    }else                                                                       colnames(residmat) <- rep(c("s", "gamma", "cost", "error"), K)
    RET <- list(residmat=residmat, cost=opt[3], gamma=opt[2], s=opt[1])
    RET
}
