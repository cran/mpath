cv.glmregNB <- function(formula, data, weights, offset=NULL, 
lambda=NULL, nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=NULL, trace=FALSE, parallel=FALSE, 
...){
    call <- match.call()
 mf <- Call <- match.call()
                                     m <- match(c("formula", "data", "subset", "weights", "na.action",
                                                "etastart", "mustart", "offset"), names(mf), 0)
                                            mf <- mf[c(1, m)]
                                            mf$drop.unused.levels <- TRUE
                                            mf[[1L]] <- as.name("model.frame")
                                            mf <- eval.parent(mf)
                                            Terms <- attr(mf, "terms")
                                            Y <- model.response(mf, "numeric")
  ## null model support
                                            X <- if (!is.empty.model(Terms)) model.matrix(Terms, mf, contrasts) else matrix(,NROW(Y),0)
  nobs <- n <- length(Y)
  nvars <- m <- dim(X) - 1
 weights <- model.weights(mf)
 offset <- as.vector(model.offset(mf))
 if(!length(weights)) weights <- rep(1, nrow(mf))
  if(any(weights < 0)) stop("negative weights not allowed")

  if(missing(foldid) && nfolds < 3)
    stop("smallest nfolds should be 3\n")
  K <- nfolds
  glmregNB.obj <- do.call("glmregNB", list(formula, data, weights,
                                              offset=offset, lambda=lambda,
                                              parallel=FALSE, ...))
  lambda <- glmregNB.obj$lambda
  nlambda <- length(lambda)
  theta0 <- glmregNB.obj$theta
  if(is.null(offset)) offset <- rep(0, nobs)
  if(missing(foldid)) all.folds <- cv.folds(n, K)
  else all.folds <- foldid

  foldenv <- list2env(list(formula=formula, data=data, X=X, Y=Y,
                           weights=weights, offset=offset, lambda=lambda,
                           nlambda=nlambda, theta0=theta0,
                           all.folds=all.folds, dots=list(...), trace=trace),
                      parent=asNamespace("mpath"))
  foldwork <- function(i) {
    omit <- all.folds[[i]]
    if(!is.null(offset)){
      offsetnow <- offset[-omit]
      newoffset <- offset[omit]
    }else offsetnow <- newoffset <- NULL
    fitcv <- do.call("glmregNB", c(list(formula, data[-omit,], weights[-omit],
                                         offset=offsetnow, nlambda=nlambda,
                                         lambda=lambda, theta.est=FALSE,
                                         theta0=theta0, parallel=FALSE,
                                         trace=trace), dots))
    fitcv$terms <- NULL
    logLik(fitcv, newx=X[omit,-1, drop=FALSE], Y[omit],
           weights=weights[omit], offset=newoffset)
  }
  environment(foldwork) <- foldenv

  if(parallel){
    dc <- parallel::detectCores(logical=TRUE)
    if(is.na(dc)) dc <- 1L
    if(is.null(n.cores)) n.cores <- dc
    n.cores <- max(1L, min(as.integer(n.cores), K, dc))
    if(.Platform$OS.type == "unix"){
      fold.results <- parallel::mclapply(seq_len(K), foldwork,
                                         mc.cores=n.cores)
      err <- sapply(fold.results, inherits, "try-error")
      if(any(err)) stop(fold.results[[which(err)[1]]])
    } else {
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl), add=TRUE)
      fold.results <- parallel::parLapply(cl, seq_len(K), foldwork)
    }
  } else {
    fold.results <- vector("list", K)
    for(i in seq_len(K)) {
      if(trace) cat("\n CV Fold", i, "\n\n")
      fold.results[[i]] <- foldwork(i)
    }
  }
  nresid <- max(c(nlambda, lengths(fold.results)))
  residmat <- matrix(NA_real_, nresid, K)
  for(i in seq_len(K)) {
    res <- fold.results[[i]]
    if(!is.null(res)) residmat[seq_along(res), i] <- res
  }
   cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.max(cv)
    obj<-list(fit=glmregNB.obj, residmat=residmat, lambda = lambda, cv = cv, cv.error = cv.error, foldid=all.folds, lambda.which= lambda.which, lambda.optim = lambda[lambda.which])
    class(obj) <- "cv.glmreg"
    if(plot.it) plot(obj,se=se)
    obj
  }

