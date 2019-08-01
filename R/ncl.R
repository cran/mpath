### nonconvex linear (NCL) models for regression and classification
ncl <- function(x, ...) UseMethod("ncl")

ncl.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(ncl.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

ncl.formula <- function(formula, data, weights, offset=NULL, contrasts=NULL, 
                        x.keep=FALSE, y.keep=TRUE, ...){
    ## extract x, y, etc from the model formula and frame
                                        #if(!attr(terms(formula, data=data), "intercept"))
                                        #    stop("non-intercept model is not implemented")
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
                                        #if(!attr(terms(formula, data=data), "intercept"))
                                        #RET <- ncl_fit(X[,-1], Y, weights,...)
                                        #else
    RET <- ncl_fit(X, Y, weights,...)
    RET$call <- match.call()
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    if(x.keep) RET$x <- data[,colnames(data)%in%attr(mt, "term.labels")]
    if(y.keep) RET$y <- Y
    class(RET) <- "ncl"
    RET
}
ncl.matrix <- function(x, y, weights, offset=NULL, ...){
    RET <- ncl_fit(x, y, weights, offset=offset, ...)
    RET$call <- match.call()
    return(RET)
}

ncl_fit <- function(x,y, weights, offset=NULL, cost=0.5, rfamily=c("clossR", "closs", "gloss", "qloss"), s=NULL, fk=NULL, iter=10, reltol=1e-5, trace=FALSE){
    call <- match.call()
    rfamily <- match.arg(rfamily)
    if(is.null(s)) stop("s must be a number\n")
     if(rfamily %in% c("closs", "qloss", "clossR") &&  s <= 0) stop("s must be positive\n")
         else if(rfamily == "gloss" && s <=1) stop("s must > 1 for gloss\n")
    if(rfamily %in% c("closs", "gloss", "qloss"))
        if(!all(names(table(y)) %in% c(1, -1)))
            stop("response variable must be 1/-1 for family ", rfamily, "\n")
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    if(missing(weights)) weights=rep(1,length(y))
    weights <- as.vector(weights)
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    ## check weights and offset
    if( !is.null(weights) && any(weights < 0) ){
        stop("negative weights not allowed")
    }
    if(length(weights) != length(y))
        stop("'weights' must be the same length as response variable")
    famtype <- switch(rfamily,
                      "closs"="clossMM",
                      "gloss"="glossMM",
                      "qloss"="qlossMM",
                      "clossR"="clossRMM"
                      )
    if(is.null(fk)){
        bsttype <- switch(rfamily,
                          "closs"="closs",
                          "gloss"="gloss",
                          "qloss"="qloss",
                          "clossR"="clossR"
                          )
### initiate predictive value
        RET <- bst(x, y, family=bsttype, ctrl = bst_control(mstop=1, s=s))
        RET$fitted.values <- RET$yhat
    }
    else {
        RET <- NULL
        RET$fitted.values <- fk
    }
    los <- loss(y, f=RET$fitted.values, cost, family = rfamily, s=s, fk=RET$fitted.values)
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\nQuadratic majorization iterations ...\n")
        cat("\ninitial loss", mean(los), "\n")
    }
    los <- rep(NA, iter)
    B <- bfunc(family=rfamily, s=s)
    while(d1 > reltol && k <= iter){
        fk_old <- RET$fitted.values
                                        #if(rfamily=="closs")
                                        #fk_old<- pmin(1, pmax(-1, fk_old)) ### truncated at -1 or 1
	if(rfamily=="clossR"){
            h <- gradient(family=rfamily, u=y-fk_old, s=s)/B+fk_old
        }
        else{ 
	    z <- y*fk_old
	    h <- -y*(gradient(family=rfamily, u=z, s=s)/B-z)
        }
    RET <- lm.wfit(x, h, weights, offset)
	fk <- RET$fitted.values
                                        #if(rfamily=="closs")
                                        #fk <- pmin(1, pmax(-1, fk)) ### truncated at -1 or 1
	los[k] <- mean(loss(y, f=fk, cost, family = rfamily, s=s, fk=NULL))
	if(trace){
        }
#	d1 <- sum((fk_old - fk)^2)/max(1, sum(fk_old^2))
	if(trace) cat("\niteration", k, ": relative change of fk", d1, ", robust loss value", los[k], "\n") 
	if(k > 1){
            d <- abs(los[k]-los[k-1]/los[k-1])
            if(los[k] > los[k-1]){
		      k <- iter
            }
        }
        k <- k + 1
	if(trace) cat("d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
    }
    RET$h <- h
    RET$x <- x
    RET$y <- y
    RET$s <- s
    RET$call <- call
    RET$cost <- cost
    RET$rfamily <- RET$family <- rfamily
    RET
}

predict.ncl=function(object, newdata=NULL, type=c("response","class","loss", "error", "coefficients"), na.action=na.pass, ...){
    type=match.arg(type)
    if(is.null(newdata)){
        if(!match(type,c("coefficients"),FALSE))stop("You need to supply a value for 'newdata'")
        ynow <- object$y
    }
    else{
        if(!is.null(object$terms)){
            mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$xlevels)
            ynow <- model.frame(object$terms, newdata)[,1] ### extract response variable
	    newdata <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
        }
    }
    if(type=="coefficients") return(object$coef)
    if(is.null(newdata))
        res <- object$fitted.values
    else
        res <- newdata %*% object$coef
    if(type=="response") return(res)
    if(type=="loss") return(mean(loss(ynow, res, family = object$family, s=object$s)))
    if(type=="error") return(evalerr(object$family, ynow, res))
}
