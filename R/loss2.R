loss2 <- function(y, f, weights, cfun, dfun, s, delta=0.0001){
    cfun <- cfun2num(cfun)
    dfun <- dfun2num(dfun)
    #if(!dfun %in% c(1, 4:7)) stop("dfun is not implemented")
    if(!cfun %in% 1:8) stop("cfun is not implemented")
    check_s(cfun, s)
    if(dfun %in% 4:7 && !setequal(unique(y), c(-1, 1))) stop("y must be -1/1 for binary classification")
    n <- length(y)
    if(n!=length(f)) stop("y and f should have comparable length\n")
    if(missing(weights)) weights <- rep(1, length(y))
    if(n!=length(weights)) stop("y and weights should have comparable length\n")
    if(cfun==6)
           if(s > 1) delta <- (s-1)/2 else
           if(s==1) delta <- 0 else{
             if(is.null(delta)) stop("delta must be provided")
             if(delta <= 0) stop("delta must be positive")
           }
    .Fortran("loss2",
                                    n=as.integer(n),
                                    y=as.double(y),
                                    f=as.double(f),
                                    weights=as.double(weights),
                                    cfun=as.integer(cfun),
                                    dfun=as.integer(dfun),
                                    s=as.double(s),
                                    delta=as.double(delta),
                                    los=as.double(0.0),
                                    PACKAGE="mpath")$los
}

loss2_ccsvm <- function(y, f, weights, cfun, dfun, s, eps, delta=0.0001){
    check_s(cfun, s)
    cfun <- cfun2num(cfun)
    if(!cfun %in% 1:8) stop("cfun is not implemented")
    dfun <- switch(dfun,
                    "eps-regression"=2,
                    "nu-regression"=2,
                    "C-classification"=6,
                    "nu-classification"=6)
    if(!dfun %in% c(2, 6)) stop("dfun is not intended for loss2_ccsvm")
    if(dfun==2 && eps < 0) stop("eps should be >= 0 for dfun=2 (eps-regression)")
    n <- length(y)
    if(n!=length(f)) stop("y and f should have comparable length\n")
    if(missing(weights)) weights <- rep(1, length(y))
    if(n!=length(weights)) stop("y and weights should have comparable length\n")
      if(cfun==6)
           if(s > 1) delta <- (s-1)/2 else
           if(s==1) delta <- 0 else{
             if(is.null(delta)) stop("delta must be provided")
             if(delta <= 0) stop("delta must be positive")
           }
    .Fortran("loss2_ccsvm",
                                    n=as.integer(n),
                                    y=as.double(y),
                                    f=as.double(f),
                                    weights=as.double(weights),
                                    cfun=as.integer(cfun),
                                    dfun=as.integer(dfun),
                                    s=as.double(s),
                                    eps=as.double(eps),
                                    delta=as.double(delta),
                                    los=as.double(0.0),
                                    PACKAGE="mpath")$los
}
