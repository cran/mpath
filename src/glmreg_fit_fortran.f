C similar to and derived from R/glmreg_fit
C input: maxit, start, etastart, mustart - are these changed in outloop
C subroutine?
c output: yhat is the updated mustart as output
      subroutine glmreg_fit_fortran(x, y, weights, n,m,start, etastart,
     +         mustart, offset, nlambda, lambda, alpha, gam, rescale, 
     +         standardize, intercept, penaltyfactor, thresh, epsbino, 
     +         maxit,eps,theta,family,penalty,trace,beta,b0, yhat, satu)
      implicit none
      integer n,m, i,j, penalty, nlambda, family, standardize, maxit,
     +     innermaxit, trace, rescale, good, intercept,
     +     satu, convout(nlambda), startv, tmpit
C      double precision, intent(in) :: etastart(n),mustart(n),start(m+1)
      double precision :: etastart(n),mustart(n),start(m+1)
      double precision x(n, m), xold(n,m), y(n), weights(n), 
Cstart(m+1),etastart(n), mustart(n), 
     +     offset(n), lambda(nlambda), lam(m, nlambda), 
     +     alpha, gam, thresh, epsbino, eps, theta, penaltyfactor(m),
     +     beta(m, nlambda), b0(nlambda), sumpen, outpll(maxit,nlambda),
     +     wt(n), meanx(m), normx(m),xd(m), nulldev,
     +     penfac(m), resdev(nlambda), yhat(n), mu(n), sumwt,
     +     crossprod_beta(nlambda), meany, meanoffset

      if(family.EQ.2)then
          do i=1, n
          if(y(i) < 0 .OR. y(i) > 1)then
          call rexit("y value should be between 0 and 1 in
     +     src/glmreg_fit_fortran")
      endif
          enddo
      endif
      call  deveval(n, y, mustart, theta, weights, family, nulldev)
      startv = 1
      sumpen = sum(penaltyfactor)
      do 110 i=1, m
      penfac(i) = penaltyfactor(i)/sumpen * m
110   continue 
C generate outer product of two vectors penfac and lambda, return lam
      call outprod(m, penfac, nlambda, lambda, lam)
      if(family.EQ.1)then 
        innermaxit = maxit
        tmpit = 1
      else 
              innermaxit = 1
              tmpit = maxit
      endif
                    
      sumwt = sum(weights)
      do 130 i=1, n
      wt(i) = weights(i)/sumwt
130   continue
C    wtnew <- weights/sum(weights)
C    meanx <- drop(wtnew %*% x)
C     compute weighted column averages meanx = x^(transpose) * wt
      call DGEMV('T',n, m, 1.0D0, x, n, wt, 1, 0.0D0, meanx, 1)

C    if(standardize){
C        xx <- scale(x, meanx, FALSE) # centers x
C        xx <- sqrt(wtnew) * xx
C        one <- rep(1,n)  
C        normx <- sqrt(one %*% xx^2)
C        xx <- scale(x, meanx, normx)
C    }
C choose family=2 (any interger except 1 to avoid centering/scaling y) below
      if(family .EQ. 1 .AND. standardize .EQ. 1)then
           do 1100 j=1, m
          do 1110 i=1, n
             xold(i,j)=x(i,j)
 1110     continue
 1100  continue
      endif
      if(standardize .EQ. 1)then
              call preprocess(x,y,n,m,weights, 2,standardize,normx,xd,
     +     mu)
      endif
C    need to check if satu is input or output
          satu = 0
          good = nlambda
      call outloop(x,y,weights,wt,n,m,penalty,nlambda,lam,alpha,gam,
     +            theta,rescale,mustart,etastart,offset,family,
     +            standardize,intercept, nulldev,thresh,tmpit,
     +            innermaxit,eps,trace,start,startv,beta,b0,resdev,yhat,
     +            convout, satu, good, epsbino,outpll)
      if (standardize .EQ. 1)then
         do 200 j=1, nlambda
         do 250 i=1, m
                beta(i, j)=beta(i, j)/normx(i)
 250      continue
 200      continue
C     compute crossproduct, like crossprod(meanx, beta) in R
C http://www.tat.physik.uni-tuebingen.de/~kley/lehre/ftn77/tutorial/blas.html
         call DGEMV('T',m, nlambda, 1.0D0, beta, m, meanx, 1, 0.0D0, 
     +             crossprod_beta, 1)
        if(intercept .EQ. 1)then
         do 300 j=1, nlambda
                b0(j)=b0(j) - crossprod_beta(j)
 300      continue
         if (family .EQ. 1)then
                meany=sum(y)/n
                meanoffset=sum(offset)/n
            do 400 j=1, nlambda
                b0(j)=b0(j) + meany - meanoffset
 400        continue
         endif
        endif
      endif
      if(family .EQ. 1 .AND. standardize .EQ. 1)then
           do 1200 j=1, m
          do 1210 i=1, n
             x(i,j)=xold(i,j)
 1210     continue
 1200  continue
         endif

      return
      end
