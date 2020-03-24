C     outer loop: sequence of lambda
C     input:
C     start
C     mu: a vector
C     eta: a vector
C     output
C     b: beta
C     bz: b0
C     resdev: residual deviance
C     ypre: yhat
C     mu, eta: uppdated values

      subroutine outloop(x,y,weights, wt, n,m, penalty, nlambda, lam, 
     +     alpha,gam,theta,rescale,mu,eta, offset,
     #family,standardize, intercept, nulldev, 
     +     thresh, maxit, innermaxit, eps, trace, start, startv,b, bz,
     +     resdev,ypre, convout, satu, good, ep, outpll)

      implicit none
      integer n,m,i,j,k,kk,penalty, nlambda, family, standardize, maxit,
     +     innermaxit, trace,convmid,convout(nlambda), startv, rescale, 
     +     satu,good,convact,fullset(m),activeset(m),activesetold(m),
     +     jk,jkold,intercept
      double precision x(n,m), y(n), wt(n), lam(m, nlambda),alpha,
     +     gam, theta,mu(n),eta(n), nulldev,thresh, eps, b(m,nlambda),
     +     bz(nlambda),xold(n,m), yold(n), start(m+1), resdev(nlambda), 
     +     ypre(n,nlambda), lamk(m), beta(m), b0,dev,
     +     weights(n),yhat(n),ep, pll(maxit), outpll(maxit, nlambda),
     +     normx(m),xd(m),avg, offset(n) 
      if(family .NE. 1)then
         call preprocess(x, y, n, m, weights, family, standardize,
     +        normx, xd, avg)
      endif
C     keep a copy of x and y in case x and y will be changed in subroutine lmnet if standardize = 1
      do 100 j=1, m
         do 110 i=1, n
            xold(i,j)=x(i,j)
 110     continue
 100  continue
      call DCOPY(n, y, 1, yold, 1)
      if(startv .EQ. 0)then
         b0 = eta(1)
         do j=1, m
            beta(j) = 0
         enddo
      else 
         b0 = start(1)
         do j=1, m
            beta(j) = start(j+1)
         enddo
      endif

      k = 1
      satu = 0
C     good = nlambda
      do 401 j=1, m
         activeset(j)=j
         fullset(j)=j
 401  continue
      jk = m
 1000 if(k .LE. nlambda .AND. satu .EQ. 0)then
         if(trace.EQ.1)then
            call intpr("Outer loop: sequence of lambda", -1, 1, 1)
            call intpr("  lambda iteration", -1, k, 1)
            call dblepr("  lambda value", -1, lam(1,k), 1)
         endif
         do 10 j=1,m
            lamk(j) = lam(j,k)
 10      continue
C     if jk=0, it means an intercept-only model, thus, we update
C     activeset
C     if(jk .EQ. 0)then
C     do 18 j=1, m
C     activeset(j)=j
C     18         continue
C     jk = m
C     endif
C     Active set: begin with the 1st element of the sequence of lambda
C     values, i.e., k=1 in this subroutine. Now, the active set
C     contains all variables, cycle through all coefficients in the
C     active set with coordinate descent algorithm until convergence,
C     then cycle through the full set with all variables, but only ONCE.
C     This generates a new active set. Compare with the previous (old) 
C     active set. If no changes, then we are done with this 1st lambda. 
C     Otherwise, we repeat the above process until the updated active
C     set remains the same, or the number of iteration (convact) is met. 
C     Next, move to the 2nd element of the sequence of lambda values. We 
C     repeat the above process with the current active set.   
C     
C     70: if block --begin
         if(family .EQ. 1)then
C     For family=1, active set is implemented in midloop -> lmnetGaus 
C     -> loop_gaussian
            call midloop(n,m,x,y, xold,yold,weights,mu, eta, offset,
     +           family, penalty,lamk,alpha,gam,theta,rescale,
     +           standardize, intercept, eps,innermaxit, maxit, thresh,
     +           nulldev, wt, beta, b0, yhat, dev, trace, convmid, ep,
     +           normx,xd,avg, activeset, jk, fullset)
         else 
C     active set applies for family!=1. 
            convact=0
C     some loop, if no change to the active set, stop
            kk = 1
C     max.cycle: in case of nonconvergence the program looks for a cycle
C     that repeats itself, default is 2 (kk <= 2)
            do 2000 while (kk <= 2 .AND. convact == 0)
C     set maxit=1, and have a complete cycle through all the variables
               call midloopGLM(n,m,x,y,yold,weights,mu,eta,offset,
     +              family, penalty,lamk,alpha,gam,theta,rescale,
     +              standardize, intercept,
     +              eps,1,thresh,nulldev,wt,beta,b0,yhat, 
     +              dev,trace,convmid,satu,ep,pll,fullset,m)
               call find_activeset(m, beta, eps, activesetold, jkold)
C     it is possible, jk=0 if beta=0, like intercept-only model
               if(jkold .EQ. 0)then
                  convact=1
                  exit
               endif
C     now cycle through only the active set
               call midloopGLM(n,m,x,y,yold,weights,mu,eta,offset,
     +              family, penalty,lamk,alpha,gam,theta,rescale,
     +              standardize, intercept,
     +              eps,maxit,thresh,nulldev,wt,beta,b0,yhat,
     +              dev,trace,convmid,satu,ep,pll,
     +              activesetold, jkold)
C     set maxit=1, and have a complete cycle through all the variables
               call midloopGLM(n,m,x,y,yold,weights,mu,eta,offset,
     +              family, penalty,lamk,alpha,gam,theta,rescale,
     +              standardize, intercept,
     +              eps,1,thresh,nulldev,wt,beta,b0,yhat, 
     +              dev,trace,convmid,satu,ep,pll,fullset,m)
C     update the active set with only non-zero coefficients 
               call find_activeset(m, beta, eps, activeset, jk)
C     check if the active set was changed--begin
C     it is possible, jk=0 if beta=0, like intercept-only model
               if(jk .EQ. 0)then
                  convact=1
                  exit
               endif
               if(jkold .NE. jk)then
                  convact=0
               else
                  do 901 j=1, jk
                     if(activesetold(j) .NE. activeset(j))then
                        exit
                     endif
                     if(j .EQ. jk)then
                        convact = 1
                        exit
                     endif
 901              continue
               endif
C     check if the active set was changed--end
               kk=kk+1
 2000       continue
         endif
C     70: if block --end
         if(satu .EQ. 1)then
            good = k - 1
         endif
         convout(k)=convmid
         if(family .NE. 1)then
            do 15 i=1, maxit
               outpll(i,k) = pll(i)
 15         continue
         endif
         do 20 j=1, m
            b(j, k) = beta(j)
 20      continue
         bz(k) = b0
         resdev(k) = dev
         call linkinv(n, yhat, family, ypre(1:n,k))
C     call linkinv(n, yhat, family, v)
C     do 30 i=1, n
C     ypre(i,k) = v(i)
C     30      continue
         k = k + 1
         if(k .LE. nlambda .AND. satu .EQ. 0)then
            do 40 j=1,m
               b(j, k) = b(j, k-1)
 40         continue
         endif
         goto 1000
      endif
      
      return
      end  


