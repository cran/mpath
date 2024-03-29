#compute composite loss function cfun(z), where z is the value computed from convex function dfun
compute_g <- function(z, cfun, s, delta=0.0001){
    if(is.null(s)) s <- assign_s(cfun, z) else check_s(cfun, dfun=1, s) #for now forget about dfun
    if(cfun==6)
           if(s > 1) delta <- (s-1)/2 else
           if(s==1) delta <- 0 else{
             if(is.null(delta)) stop("delta must be provided")
             if(delta <= 0) stop("delta must be positive")
           }
    n <- length(z)
    .Fortran("compute_g",
             cfun=as.integer(cfun),
             n=as.integer(n),
             z=as.double(z),
             s=as.double(s),
             delta=as.double(delta),
             gval=as.double(rep(0, n)),
             PACKAGE="mpath")$gval
}
