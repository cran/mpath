C###Adapted from file MASS/R/glmregNB.R
C input: theta0 is useful if theta_fixed=1, in which case theta0 is a 
C  vector of length nlambda
      subroutine glmregnb_fortran(x, y, weights, n, m, offset, nlambda,
     +    lambda, penalty, alpha, gam, rescale, standardize, 
     +    penaltyfactor, thresh, maxit_theta, maxit, eps,epsbino,start, 
     +    etastart, mustart, thetastart, theta_fixed, theta0, trace,
     +    beta, b0, tht, yhat)
      implicit none
      integer k, i, j, n, m, trace, nlambda, penalty, maxit_theta, 
     +  maxit, rescale, standardize, iter, theta_fixed, satu
      double precision x(n, m), y(n), lambda(nlambda), weights(n),
     +  offset(n), start(m+1), etastart(n), mustart(n), thetastart, 
     +  del, epsbino, theta0(nlambda), d, tht(nlambda), eps, thresh, 
     +  alpha, gam, penaltyfactor, beta_1(m, 1), b0_1, 
     +  beta(m, nlambda), b0(nlambda), yhat(n)
      
      k = 1
C      convout <- twologlik <- rep(NA, nlambda)
      del = eps**0.25
10    if(k <= nlambda)then
        if(trace .EQ.1 )then
            call intpr("loop in lambda:", -1,  k, 1)
        endif
        if(theta_fixed .EQ. 1)then
            thetastart = theta0(k)
        endif
        iter = 0
        d = 10
30        if(d .GT. del .AND. iter <= maxit_theta)then
            call glmreg_fit_fortran(x, y, weights, n, m,start,
     +       etastart, mustart, offset, 1, lambda(k), alpha, gam, 
     +       rescale, standardize, penaltyfactor, thresh, epsbino,
     +       maxit,eps, thetastart, 4, penalty, trace, beta_1, 
     +       b0_1, yhat, satu)
            call dcopy(n, yhat, 1, mustart, 1)
            do i=1, n
            etastart(i) = dlog(mustart(i))
            enddo
            if(theta_fixed .EQ. 0)then
                call theta_ml(y, mustart, n, weights, 10,
     +           eps**0.25, thetastart, trace)
            endif
            d=0
            d=d+(start(1)-b0_1)**2
            start(1)=b0_1
            do j=1, m
            d=d+(start(j+1)-beta_1(j,1))**2
            start(j+1)=beta_1(j,1)
            enddo
C            del <- thetastart - th
C            Lm0 <- Lm
C            penval <- ifelse(standardize, n*fit$penval, fit$penval)
C            Lm <- loglik(n, th, mu, Y, w) - penval
C            fit$df.residual <- n - fit$df - 1
C            d1 <- sqrt(2 * max(1, fit$df.residual))
C            converged <- abs((Lm0 - Lm)/d1) < 1e-8
C            if(trace) {
C                Ls <- loglik(n, th, Y, Y, w)
C                Dev <- 2 * (Ls - Lm)
C                message("Theta(", iter, ") =", signif(th),
C                        ", 2(Ls - Lm) =", signif(Dev))
C            }
            iter=iter + 1
            goto 30
        endif
        tht(k) = thetastart
        call dcopy(m, beta_1, 1, beta(1:m, k), 1)
        b0(k) = b0_1
C        call loglikFor(n, y, mustart, thetastart, weights, 4, loglik(k))
C        nulldev(k) = nulldev0
C        resdev(k) = resdev0
C        if(trace) pll[,k] <- fit$pll
C        convout[k] <- converged
C        fitted[,k] <- fit$fitted.values
        k = k + 1
        goto 10
      endif
 
      return
      end
