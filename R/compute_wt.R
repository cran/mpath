### this is different from R/update_wt in that the inputs are different 
compute_wt <- function(z, weights, cfun, s, delta=0.0001){
    if(any(z < 0)) stop("z must be non-negative")
    if(is.null(s)) s <- assign_s(cfun, z) else check_s(cfun, s)
    if(missing(weights)) weights <- rep(1, length(z))
    if(cfun==6)
           if(s > 1) delta <- (s-1)/2 else
           if(s==1) delta <- 0 else{
             if(is.null(delta)) stop("delta must be provided")
             if(delta <= 0) stop("delta must be positive")
           }
    n <- length(z)
    v <- .Fortran("compute_v",
             cfun=as.integer(cfun),
             n=as.integer(n),
             z=as.double(z),
             s=as.double(s),
             delta=as.double(delta),
             v=as.double(rep(0, n)),
             PACKAGE="mpath")$v
    -v*weights
}
