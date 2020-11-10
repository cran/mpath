### composite optimization based glm
ccglm <- function(x, ...) UseMethod("ccglm")

ccglm.formula <- function(formula, data, weights, offset=NULL, contrasts=NULL, 
cfun="ccave", dfun=gaussian(), s=NULL, delta=0.1, fk=NULL, init.family=NULL, iter=10, reltol=1e-5, theta, x.keep=FALSE, y.keep=TRUE, trace=FALSE, ...){
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
    x <- X
    y <- Y
    call <- match.call()
    cfunold <- cfun
    dfunold <- dfun[[1]]
    if(dfunold=="poisson"){
        theta <- 1 ### not used anyway
        family <- 3 
    }else if(dfunold=="negbin"){
        if(missing(theta)) stop("theta has to be provided for family=negbin()")
        family <- 4
    } 
    cfun <- cfun2num(cfun)
    if(dfunold%in%c("gaussian","gaussianC","binomial","poisson","negbin")) dfunold <- dfun2num(dfunold)
                                        #if(dfun==7 && !is.null(offset))
                                        #    warning("dfun=", dfun, ", offset not implmented\n")
    if(is.null(s) || is.na(s)) s <- assign_s(cfun, y) else check_s(cfun, s)
    if(cfun==6)
        if(s > 1) delta <- (s-1)/2 else
                                       if(s==1) delta <- 0 else{
                                                               if(is.null(delta)) stop("delta must be provided")
                                                               if(delta <= 0) stop("delta must be positive")
                                                           }
                                        #for dfunold=4, y must be 1 or -1
                                        #for dfunold=5, y must be 0 or 1
    if(dfunold==4)
        y <- y2num(y)
    if(dfunold==5)
        y <- y2num4glm(y)
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
    n <- length(y)
    if (is.null(offset)){
        is.offset <- FALSE
        offset <- rep.int(0, n)
    }else is.offset <- TRUE
     if(!is.null(init.family)) rfamily <- init.family else{
                                                          if(cfun==4 && dfunold==1)  rfamily <- "clossR"
                                                          else if(cfun==4 &&         dfunold==4) rfamily <- "closs"
                                                          else rfamily <- "gaussian"
                                                      }
     if(is.null(fk)){
 ### initiate predictive value
         RET <- bst(x, y, family=rfamily, ctrl = bst_control(mstop=1, s=s))
         RET$linear.predictors <- RET$yhat
     }
     else {
         RET <- NULL
         RET$linear.predictors <- fk
     }
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\niterations ...\n")
    }
    los <- rep(NA, iter)
    ytmp <- y
    if(dfunold==5)
        ytmp[ytmp==0] <- -1
    weights_update <- weights
###initial values
    if(!dfunold %in% c(1,4,5,8,9)){
        envir <- list2env(list(weights_update=weights_update, offset=offset), parent=environment(formula))
        environment(formula) <- envir
        RET <- suppressWarnings(glm(formula, data, family=dfun, weights=weights_update, offset=offset))
    }
    while(d1 > reltol && k <= iter){
        if(dfunold %in% c(1,4,5)){
            weights_update <- .Fortran("update_wt",
                                       n=as.integer(n),
                                       weights=as.double(weights),
                                       y=as.double(ytmp),
                                       f=as.double(RET$linear.predictors),
                                       cfun=as.integer(cfun),
                                       dfun=as.integer(dfunold),
                                       s=as.double(s),
                                       delta=as.double(delta),
                                       weights_update=as.double(rep(0, n)),
                                       PACKAGE="mpath")$weights_update
            los[k] <- .Fortran("loss2",
                               n=as.integer(n),
                               y=as.double(ytmp),
                               f=as.double(RET$linear.predictors),
                               weights=as.double(weights/n),
                               cfun=as.integer(cfun),
                               dfun=as.integer(dfunold),
                               s=as.double(s),
                               delta=as.double(delta),
                               los=as.double(0.0),
                               PACKAGE="mpath")$los
        }
        else{
            tmp1 <- loss3(y,mu=RET$fitted.values,theta,weights,cfun,family,s,delta)
            weights_update <- compute_wt(tmp1$z, weights, cfun, s, delta)
            los[k] <- sum(tmp1$tmp)
        }
        envir <- list2env(list(weights_update=weights_update, offset=offset), parent=environment(formula))
        environment(formula) <- envir
        RET <- suppressWarnings(glm(formula, data, family=dfun, weights=weights_update, offset=offset))
	if(trace) cat("\niteration", k, ", robust loss value", los[k], "\n") 
	if(k > 1){
            d1 <- abs((los[k]-los[k-1])/los[k-1])
            if(los[k] > los[k-1]){
                k <- iter
            }
        }
        k <- k + 1
	if(trace) cat("d1=", d1, ", k=", k, ", d1 > reltol && k <= iter: ", (d1 > reltol && k <= iter), "\n")
    }
    RET$s <- s
    RET$call <- call
    RET$cfun <- cfunold
    RET$dfun <- dfunold
    RET$weights <- weights
    RET$weights_update <- weights_update
    RET$is.offset <- is.offset
    RET <- c(RET, list(formula=formula, terms = mt, data=data,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    if(x.keep) RET$x <- data[,colnames(data)%in%attr(mt, "term.labels")]
    if(y.keep) RET$y <- Y
    class(RET) <- c("ccglm", "glm")
    RET
}

predict.ccglm=function(object, newdata=NULL, weights=NULL, newoffset=NULL, type=c("link","response","class","loss", "error", "coefficients"), na.action=na.pass, ...){
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
    if(object$is.offset)
        if(is.null(newoffset))
            stop("offset is used in the estimation but not provided in prediction\n")
        else offset <- newoffset
    else offset <- rep(0, length(ynow))
    if(type=="coefficients") return(object$coef)
    if(is.null(newdata))
        res <- object$linear.predictors
    else
        res <- offset + newdata %*% object$coef
    if(type=="link") return(res)
    if(type %in% c("link", "response")) return(predict(object, type=type)) 
    object$dfun <- dfun2num(object$dfun)
    if(object$dfun %in% 4) ynow <- y2num(ynow)
    if(type=="loss"){
        n <- length(ynow)
        if(missing(weights)) weights <- rep(1/n, n)
        object$cfun <- cfun2num(object$cfun)
        if(object$dfun%in%c(1,4)) return(.Fortran("loss2",
                                                  n=as.integer(length(ynow)),
                                                  y=as.double(ynow),
                                                  f=as.double(res),
                                                  weights=as.double(weights),
                                                  cfun=as.integer(object$cfun),
                                                  dfun=as.integer(object$dfun),
                                                  s=as.double(object$s),
                                                  delta=as.double(object$delta),
                                                  los=as.double(0.0),
                                                  PACKAGE="mpath")$los)
        else 
            sum(loss3(object$y,mu=predict(object, data.frame(newdata), offset=offset, type="response"),object$theta,weights,object$cfun,family,object$s,object$delta)$tmp)
    }
    if(type=="error" && object$dfun %in%c(1,4)) 
        return(co_evalerr(object$dfun, ynow, res))
    if(type=="error" && object$dfun!=5) 
        return(co_evalerr(object$dfun, ynow, predict(object, data.frame(newdata), offset=offset, type="response")))
    else if(type=="class" && object$dfun==5)
        return(predict(object, data.frame(newdata), offset=offset, type="response") > 0)
    else if(type=="error" && object$dfun==5)
        return(ynow!=(predict(object, data.frame(newdata), offset=offset, type="response") > 0))
}

check_s <- function(cfun, s){
    if(cfun %in% c(1:6,8) && s <= 0) stop("s must > 0 for cfun=", cfun, "\n")
    else if(cfun==7 && s < 0) stop("s must be nonnegative for cfun=", cfun, "\n")
    if(cfun %in% c(2, 3) && s < 1) warnings("The program can crash if s value is too close to 0 for cfun=", cfun, "\n")
}

                                        #for cfun=1, 2, 3, choose s value to achieve 95% efficiency for Huber's method, Andrew's since and Tukey's biweight loss, respectively. For cfun=4 and dfun=1 or 4, s < 1 is nonconvex; For cfun=4 and dfun=4, s=1 is an approximation to the hinge loss. cf: Singh et al 2014, Pattern Recognition47:441-453, page 445. For cfun=5, s=4 (or mu in the paper) was used in Krause and Singer (2004) ICML. For cfun=6, s=2 was used in Li and Bradic (2018) JASA. For cfun=7 and dfun=6, s=-1 works best in Wu and Liu 2007, JASA 102(479):974-983, page 981 
assign_s <- function(cfun, y){
    if(cfun==1) return(1.345)
    if(cfun==2) return(1.339)
    if(cfun==3) return(4.685)
    if(cfun==4) return(0.5)
    if(cfun==5) return(4)
    if(cfun==6) return(1)
    if(cfun==7) return(1)
    if(cfun==8) return(1)
}

cfun2num <- function(cfun){
    switch(cfun,
           "hcave"=1,
           "acave"=2,
           "bcave"=3,
           "ccave"=4,
           "dcave"=5,
           "gcave"=6,
           "tcave"=7,
           "ecave"=8)
}
dfun2num <- function(dfun){
    switch(dfun,
           "gaussian"=1,
           "eps-regression"=2,
           "huber"=3,
           "gaussianC"=4,
           "binomial"=5,
           "hinge"=6,
           "exponential"=7,
           "poisson"=8,
           "negbin"=9)
}
### output: +1/-1
y2num <- function(y){
    if (is.factor(y)) {
        y <- as.integer(y)
    } else {
        if(any(as.integer(y) != y))
            stop("y variable has to be of factor or integer type for classification\n")
        y <- as.factor(y)
        y <- as.integer(y)
    }
    if(length(unique(y))!=2) stop("y must be a binary variable for classification\n")
    else if(!all(names(table(y)) %in% c(1, -1))){
        y[y==1] <- -1
        y[y==2] <- 1
    }
    y
}
### output: 0/1
y2num4glm <- function(y){
    if (is.factor(y)) {
        y <- as.integer(y)
    } else {
        if(any(as.integer(y) != y))
            stop("y variable has to be of factor or integer type for classification\n")
        y <- as.factor(y)
        y <- as.integer(y)
    }
    if(length(unique(y))!=2) stop("y must be a binary variable for classification\n")
    else if(!all(names(table(y)) %in% c(0, 1))){
        y[y==1] <- 0
        y[y==2] <- 1
    }
    y
}

###compute composite loss values for dfun=poisson and negbin
### output: z is a vector of maximum log-likelihood value, i.e., s(u) values 
### output: tmp is a vector of g(s(u)) values
loss3 <- function(y,mu,theta,weights,cfun,family,s,delta){
    n <- length(y)
    z <- tmp <- rep(NA, n)
    for(i in 1:n){
        z[i] <- .Fortran("loglikFor",
                         n=as.integer(1),
                         y=as.double(y[i]),
                         mu=as.double(mu[i]),
                         theta=as.double(theta),
                         w=as.double(1),
                         family=as.integer(family),
                         ll=as.double(0),
                         PACKAGE="mpath")$ll
        z[i] <- -z[i] #convert MLE of log-likelihood to minimization of negative log-likelihood
        tmp[i] <- weights[i] * .Fortran("compute_g",
                                        cfun=as.integer(cfun),
                                        n=as.integer(1),
                                        z=as.double(z[i]),
                                        s=as.double(s),
                                        delta=as.double(delta),
                                        gval=as.double(0),
                                        PACKAGE="mpath")$gval
    }
    list(z=z, tmp=tmp)
}

###compute updated weights for composite loss function with dfun=poisson and negbin 
compute_wt3 <- function(y,mu,theta,weights,cfun,family,s,delta){
    n <- length(y)
    z <- tmp <- rep(NA, n)
    for(i in 1:n){
        z[i] <- .Fortran("loglikFor",
                         n=as.integer(1),
                         y=as.double(y[i]),
                         mu=as.double(mu[i]),
                         theta=as.double(theta),
                         w=as.double(1),
                         family=as.integer(family),
                         ll=as.double(0),
                         PACKAGE="mpath")$ll
    }
    z <- -z #convert MLE of log-likelihood to minimization of negative log-likelihood
    compute_wt(z, weights, cfun, s, delta)
}


