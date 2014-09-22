"vuong.test" <- function(m1,s1,m2,s2, type=c("None", "AIC", "BIC"), digits=getOption("digits")){
  type <- match.arg(type)
  ## get predicted probabilities for both models
  m1y <- m1$y
  m2y <- m2$y
  m1n <- length(m1y)
  m2n <- length(m2y)
  if(m1n==0 | m2n==0)
    stop("could not extract dependent variables from models")
  
  if(m1n != m2n)
    stop(paste("models appear to have different numbers of observations\n",
               "model 1 has ",m1n," observations\n",
               "model 2 has ",m2n," observations\n",
               sep="")
         )
  
  if(any(m1y != m2y))
    stop(paste("models appear to have different values on dependent variables\n"))
  if(class(m1) %in% c("glmreg", "zipath") && length(s1) > 1)
    stop("length of s1 should be 1\n")
  if(class(m2) %in% c("glmreg", "zipath") && length(s2) > 1)
    stop("length of s2 should be 1\n")
  at <- unique(m1$y)
  if(any(class(m1) %in% c("glmreg","glmregNB","zipath")))
  k1 <- sum(abs(unlist(coef(m1, s1))) > 0)
  if(any(class(m2) %in% c("glmreg","glmregNB","zipath")))
  k2 <- sum(abs(unlist(coef(m2, s2))) > 0)
  if(class(m1)[1]=="glmregNB")
  k1 <- k1+1
  if(class(m2)[1]=="glmregNB")
  k2 <- k2+1
  if(class(m1)[1]=="zipath"){
  if(m1$family=="negbin" && !m1$theta.fixed)
  k1 <- k1+1
}
  if(class(m2)[1]=="zipath"){
  if(m2$family=="negbin" && !m2$theta.fixed)
  k2 <- k2+1
}
  if(any(class(m1)%in% "glm")){
  k1 <- length(coef(m1))
  if(class(m1)[1]=="negbin")
  k1 <- k1+1
  m1 <- conv2glmreg(m1)
  s1 <- 1
  }
  if(any(class(m2)[1]=="glm")){
  k2 <- length(coef(m2))
  if(class(m2)[1]=="negbin")
  k2 <- k2+1
  m2 <- conv2glmreg(m2)
  s2 <- 1
  }
  if(class(m1)[1]=="zeroinfl"){
  k1 <- length(coef(m1))
  if(m1$dist=="negbin")
  k1 <- k1+1
  m1 <- conv2zipath(m1)
  s1 <- 1
  }
  if(class(m2)[1]=="zeroinfl"){
  k2 <- length(coef(m2))
  if(m2$dist=="negbin")
  k2 <- k2+1
  m2 <- conv2zipath(m2)
  s2 <- 1
  }
  whichCol <- match(m1y,at)
  #whichCol <- match(m1y,min(m1y):max(m1y))  ## which column, matrix of predicted probs
  
  m1p <- rep(NA,m1n)
  m2p <- rep(NA,m2n)
  if(class(m1)[1]!="zipath")
  p1 <- predprob(m1, which=s1, newdata=as.data.frame(m1$x), at=at)   ## likelihood contributions, model 1, cond on MLEs
  else
  p1 <- predprob(m1, which=s1, at=at)   ## likelihood contributions, model 1, cond on MLEs
  if(class(m2)[1]!="zipath")
  p2 <- predprob(m2, which=s2, newdata=as.data.frame(m2$x), at=at)   ## likelihood contributions, model 2
  else p2 <- predprob(m2, which=s2, at=at)   ## likelihood contributions, model 2
  for(i in 1:m1n){
    m1p[i] <- p1[i,whichCol[i]]  ## pick off correct column
    m2p[i] <- p2[i,whichCol[i]]
  }
  
  m <- log(m1p/m2p)  ## vector of likelihood ratios
  if(type=="AIC") m <- m-(k1-k2) 
  if(type=="BIC") m <- m-(k1-k2)/2 * log(m1n) 
  bad <- is.na(m) + is.nan(m) + (m==Inf) + (m==-Inf)
  if(any(bad)){
    cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
    cat("dropping these cases, but proceed with caution\n")
  }

  ## gather up degrees of freedom
#  k1 <- length(coef(m1))
#  k2 <- length(coef(m2))

  ## test statistic: Long (1997) p248
  mbar <- mean(m[!bad])
  s <- sd(m[!bad])
  v <- sqrt(sum(!bad))*mbar/s

  ## bundle up for output
  cat(paste("Vuong Non-Nested Hypothesis Test-Statistic:",
            signif(v,digits),
            "\n"))
  cat("(test-statistic is asymptotically distributed N(0,1) under the\n")
  cat(" null that the models are indistinguishible)\n")
  
  cat("in this case:\n")
  if(v>0)
    cat(paste("model1 > model2, with p-value",
              signif(1-pnorm(v),digits),
              "\n"))
  else
    cat(paste("model2 > model1, with p-value",
              signif(pnorm(v),digits),
              "\n"))
  
  invisible(NULL)
}

