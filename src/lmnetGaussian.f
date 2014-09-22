      subroutine preprocess(x,y,n,m,weights,family,standardize,normx,xd,
     +mu)
      implicit none
      double precision x(n,m),y(n),weights(n),normx(m), xd(m), 
     +wtnew(n),meanx(m),mu,wsum,ddot,xtmp(n,m)
      integer i,j,n,m, standardize,family

      do i=1, n
      do j=1,m
      xtmp(i,j) = x(i,j)
      enddo
      enddo
      wsum=0      
C compute weighted means sum(weights_i*y_i)
      mu = ddot(n, y, 1, weights, 1)
      do 20 i=1, n
      wsum = wsum + weights(i)
  20  continue
      mu = mu/wsum
      
C### Center x and y, and scale x, and save the means and sds 
C    for now, except for gaussian family, ignore weights
      do 30 i=1, n
      wtnew(i) = weights(i)/wsum
   30 continue
C compute weighted column averages meanx = x^(transpose) * wtnew

      call DGEMV('T',n, m, 1.0D0, x, n, wtnew, 1, 0.0D0, meanx, 1)

      if(standardize.EQ.1)then
          do 70 j=1, m
        do 60 i=1, n
        x(i,j) = x(i,j) - meanx(j)
C  deal the weights later
       x(i,j) = dsqrt(wtnew(i)) * x(i,j)
   60 continue
   70 continue
      if(family .EQ. 1)then
      do 80 i=1, n
        y(i) = y(i) - mu
      y(i)= dsqrt(wtnew(i)) * y(i)
   80 continue
      endif
      do 110 j=1, m
      normx(j)=0
      xd(j) =1
       do 100 i=1, n
       normx(j) = normx(j) + x(i,j) * x(i,j)
  100  continue
      normx(j) = dsqrt(normx(j))
  110 continue
      else
      do 120 j=1, m
      xd(j) = 0
      normx(j) = 1
       do 125 i=1, n
      xd(j) = xd(j) + weights(i) * x(i,j) * x(i,j)
  125 continue
  120 continue
      endif
      if(standardize.EQ.1 .AND. family .EQ. 1)then
        do 140 j=1, m
       do 130 i=1, n
      x(i,j) = x(i,j)/normx(j)
 130  continue
 140  continue
      endif
      if(standardize .EQ. 1 .AND. family .NE. 1)then
       do 150 j=1, m
       do 160 i=1, n
      x(i,j) = (xtmp(i,j)-meanx(j))/normx(j)
 160  continue
 150  continue
      endif
      
      return
      end

C                                        #coordinate descent algorithm
C                                        #ref: Regularization paths for generalized linear models via coordinate descent, Friedman et al., JSS, 2010, 33(1)
C  output:
C   beta
C   b0
C   mu
C   yhat
C   jj
C  input: others except the aboved mentioned output
      subroutine lmnetGaus(x, y, n, m, weights, w, lambda, alpha, 
     +gam, 
     +thresh, maxit, eps, standardize, family, 
     +trace,penalty,normx,xd,beta, b0, 
     +avg, yhat, jj, rescale, converged)

      implicit none
      double precision x(n, m), y(n), weights(n), lambda(m),
     +meanx(m), normx(m), xd(m),yhat(n),resid(n),
     +alpha, gam, thresh, eps, avg, wsum,w(n),
     + wtnew(n), b0, beta(m), ddot
      integer maxit, standardize, trace, penalty, n,m,i,j,jj,converged,
     +family, rescale 

      do i=1, n
      resid(i) = y(i)
      enddo
      if(standardize .EQ. 1)then
      b0 = 0
      else 
       b0 = avg
      endif
      call loop_gaussian(x,y, n,m,penalty,thresh,eps,maxit,standardize,
     +beta,b0,
     +resid,xd,lambda,alpha,gam,weights,avg,meanx,trace,jj,rescale, 
     +converged)
      jj = jj - 1
      if(standardize.EQ.1)then
      do 200 j=1, m
      beta(j) = beta(j)/normx(j)
      if(dabs(beta(j)) .LT. 1e-10)then
      beta(j) = 0
      endif
  200 continue
      endif
      b0 = 0
      do 210 j=1, m
      b0 = b0 + meanx(j) * beta(j)
  210 continue
      b0 = avg - b0

      return
      end

! coordinate descent algorithm

       subroutine enet(z, t, lone, ltwo, res)
       implicit none
       double precision z, t, lone, ltwo, res

           call soth(z, lone, res)
           res= res/(t + ltwo)
          return
          end
       
       subroutine mcp(z, t, lone, ltwo, gam, rescale, res)
       implicit none
       integer rescale
       double precision z, t, v, lone, ltwo, gam, res

       if(rescale .EQ. 1) then
       v = 1
       else
       v = t
       endif
                  if(dabs(z) .LE. v*gam * lone * (1+ltwo)) then
                    call soth(z, lone, res)                
                    if(rescale .EQ. 1) then
                     res =res/(t*(1 + ltwo - 1/gam))
                    else 
                     res =res/(t + ltwo - 1/gam)
                    endif
                  else 
                    if(rescale .EQ. 1) then
                     res = z/(t*(1 + ltwo))
                    else 
                     res = z/(t + ltwo)
                    endif
                  endif
                  return
                  end

       subroutine scad(z, t, lone, ltwo, gam, rescale, res)
       implicit none
       integer rescale
       double precision z, v, t, lone, ltwo, gam, res
       
       if(rescale .EQ. 1) then
       if(dabs(z) .LE. lone   + lone * (1+ltwo)) then
                    call soth(z, lone, res)                
C                    if(rescale .EQ. 1) then
                     res =res/(t*(1 + ltwo))
C                   else 
C                    res =res/(t + ltwo)
C                    endif
           else if(dabs(z) .LE. gam * lone * (1+ltwo)) then
                   call soth(z, gam*lone/(gam-1), res)
                     res=res/(t*(1-1/(gam-1)+ltwo))
                else 
                     res=z/(t*(1+ltwo))
        endif     
      else
      if(dabs(z) .LE. lone   + lone * (t+ltwo)) then
                    call soth(z, lone, res)                
                     res =res/(t + ltwo)
           else if(dabs(z) .LE. gam * lone * (t+ltwo)) then
                   call soth(z, gam*lone/(gam-1), res)
                     res=res/(t-1/(gam-1)+ltwo)
                else 
                     res=z/(t+ltwo)
                    endif
        endif     
       return
       end   

! input variables
! x
! y
! n
! m
! eps
! maxit
! threshold
! standardize
! beta
! b0
! xd
! lambda
! alpha
! gam
! 
! output variables
! beta
! jj       


      subroutine loop_gaussian(x, y, n, m, penalty, thresh, eps, maxit, 
     + standardize, beta, b0, resid,xd, lambda, alpha, gam,wtold,avg, 
     + meanx, trace, jj, rescale, converged)
      implicit none
      integer trace, n, m, maxit, standardize, jj, i,j,penalty,converged
      integer rescale
      double precision x(n, m), y(n), thresh, eps,beta(m),beta_old(m) 
      double precision lambda(m), alpha, gam, wtold(n), avg, meanx(m)
      double precision z, b0, xd(m), yhat(n), resid(n),b0_old
    
      jj = 1
      converged = 0
 1000    if(jj .LE. maxit .AND. converged .EQ.0)then
      do 100 j = 1, m
      beta_old(j) = beta(j)
  100 continue
      b0_old = b0
       do 10 i = 1, n
         yhat(i) = b0
         do 20 j = 1, m
         yhat(i) = yhat(i) + x(i,j) * beta(j)
   20 continue         
   10 continue
       do 30 i = 1, n
      resid(i) = y(i) - yhat(i)
   30 continue
! When cycling through the variables, we could restrict ourselves to the current 
! active set and visit all the variables e.g. every 10th iteration to update 
! the active set.
      do 40 j = 1, m
         z = 0.0D0
       if(standardize.EQ.1)then
         do 50 i = 1, n
         z = z + x(i,j) * resid(i)
   50 continue 
         z = z + beta(j)    
      else 
        do 60 i = 1, n
         z= z + wtold(i)*x(i,j) * (resid(i) + x(i,j) * beta(j))  
   60 continue
       endif
       if(penalty.EQ.1) then
       call enet(z, xd(j), lambda(j)*alpha,lambda(j)*(1-alpha), beta(j))
        else if(penalty.EQ.2) then
       call mcp(z, xd(j), lambda(j)*alpha, lambda(j)*(1-alpha), gam,
     + rescale, beta(j))
        else if(penalty.EQ.3) then
       call scad(z, xd(j), lambda(j)*alpha, 
     + lambda(j)*(1-alpha), gam, rescale, beta(j))
       endif
      if(dabs(beta(j) - beta_old(j)) .GT. eps)then
      do 70 i = 1, n
              resid(i) = resid(i) - x(i, j) * (beta(j) - beta_old(j))
   70 continue
      endif
   40 continue
      if(standardize.EQ.0)then
         b0 = 0.0D0
        do 80 j = 1, m
         b0 = b0 + meanx(j) * beta(j)
   80  continue
         b0 = avg - b0
      endif
      call checkConvergence(m, beta, beta_old, eps, thresh, converged)
      jj = jj + 1
      goto 1000
      endif

      return
      end

! soft-threshold=g operator
      subroutine soth(z, g, res)
      implicit none
      double precision res, z, g
      
      if(z .GT. g)then
              res=z-g
      else if (dabs(z) .LE. g) then
              res=0
           else if (z .LT. -g) then
              res=z+g
      endif
              return
              end

      subroutine checkConvergence(m, beta, betaold,eps,thresh,converged)
      implicit none
      integer m, j, converged
      double precision beta(m), betaold(m), eps, thresh

      converged=1
      do 10 j=1, m
      if(dabs(beta(j)) .GT. eps .AND. dabs(betaold(j)) .GT. eps) then
              if(dabs((beta(j)-betaold(j))/betaold(j)) .GT. thresh) then
                      converged=0
                      exit
              endif
      else if(dabs(beta(j)) .LE. eps .AND. dabs(betaold(j)) .GT. eps)
     +    then
              converged=0
              exit
      else if(dabs(beta(j)) .GT. eps .AND. dabs(betaold(j)) .LE. eps)
     + then
              converged=0
              exit
      endif
   10 continue

      return
      end

