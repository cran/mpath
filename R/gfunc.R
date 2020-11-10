gfunc <- function(mu, family, epsbino){
    n <- length(mu)
     .Fortran("gfunc",
              mu=as.double(mu),
              n=as.integer(n),
              family=as.integer(family),
              epsbino=as.double(epsbino),
              eta=as.double(rep(0, n)),
              PACKAGE="mpath")$eta
}
