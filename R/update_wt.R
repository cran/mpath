### this is different from src/update_wt in that the inputs are different 
update_wt <- function(y, ypre, weights, cfun, s, dfun, delta=0.0001){
    if(dfun %in% 4:7 && !setequal(unique(y), c(-1, 1))) stop("y must be -1/1 for binary classification")
    n <- length(y)
    if(missing(weights)) weights <- rep(1, n)
    if(cfun==6)
        if(s > 1) delta <- (s-1)/2 else
                                       if(s==1) delta <- 0 else{
                                                               if(is.null(delta)) stop("delta must be provided")
                                                               if(delta <= 0) stop("delta must be positive")
                                                           }
    .Fortran("update_wt",
             n=as.integer(n),
             weights=as.double(weights),
             y=as.double(y),
             f=as.double(ypre),
             cfun=as.integer(cfun),
             dfun=as.integer(dfun),
             s=as.double(s),
             delta=as.double(delta),
             weights_update=as.double(rep(0, n)),
             PACKAGE="mpath")$weights_update
}
