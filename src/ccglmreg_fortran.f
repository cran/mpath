C     used in ccglmreg.R ((adpated from nclreg_fortran.f)
C     output: nlambdacal: number of lambda actually computed including
C     repeated ones
      subroutine ccglmreg_fortran(x, y, weights, n,m,start, etastart,
     +     mustart, offset, iter, nlambda, lambda, alpha, gam, rescale, 
     +     standardize, intercept, penaltyfactor, maxit, eps, theta, 
     +     penalty, trace,del,cfun, dfun, s,thresh,decreasing, 
     +     active, beta, b0, yhat, los, pll, nlambdacal, delta,
     +     weights_cc)
      implicit none
      integer n,m,i,ii,j,jj,penalty,nlambda, standardize, maxit,
     +     trace, iter, jk, active, activeset(m), intercept,rescale, 
     +     m_act, nlambdacal, uturn, decreasing, cutpoint, cfun,
     +     dfun,dfunnew, 
     +     AllocateStatus, DeAllocateStatus, varsel(m), varsel_old(m)
      double precision x(n, m), y(n), weights(n), weights_update(n), 
     +     weights_cc(n, nlambda), start(m+1),etastart(n),
     +     mustart(n), offset(n), lambda(nlambda), alpha, gam, eps, 
     +     penaltyfactor(m), thresh, beta(m, nlambda),
     +     b0(nlambda), b0_1, yhat(n), del, lambda_i, s,sumw,wt(n),  
     +     fk(n), los(nlambda), pll(nlambda), penval, delta, theta
      double precision, dimension(:, :), allocatable :: x_act
      double precision, dimension(:), allocatable :: start_act, beta_1,
     +     penaltyfactor_act

      if(dfun .EQ. 8)then
         dfunnew= 3
         else if(dfun .EQ. 9)then
          dfunnew= 4
       endif
      m_act = m
      jk = m
      sumw = 0.0D0
       do i=1, n
          sumw=sumw+weights(i)
         enddo
         do i=1, n
          wt(i)=weights(i)/sumw
       enddo

      allocate(start_act(m_act+1), stat=AllocateStatus)
      if(AllocateStatus .NE. 0)then
         return
      endif
      allocate(penaltyfactor_act(m_act), stat=AllocateStatus)
      allocate(x_act(n, m), stat=AllocateStatus)
      if(AllocateStatus .NE. 0)then
         return
      endif
      call copymatrix(n, m, x, x_act)
      allocate(beta_1(m_act), stat=AllocateStatus)
      if(AllocateStatus .NE. 0)then
         return
      endif
      do 5 j=1, m
         beta_1(j)=0
         activeset(j)=j
 5    continue
      call dcopy(m+1, start, 1, start_act, 1)
      call dcopy(m, penaltyfactor, 1, penaltyfactor_act, 1)
      do 101 j=1, m
         varsel_old(j)=j
         varsel(j)=j
 101  continue
      i=1
      nlambdacal=0
      uturn=0
      cutpoint=1
 10   if(i .LE. nlambda)then
         if(cutpoint > 1 .AND. uturn==0 .AND. i < cutpoint)then
            i = cutpoint + 1
         endif
         if(trace .EQ. 1)then
            call intpr("ccglmreg_fortran lambda iteration", -1, i, 1)
            if(uturn==1)then
               call intpr("uturn=1", -1, 1, 1)
            endif
         endif
C unlike in nclreg_fortran, no need B: lambda_i=lambda(i)/B
         if(i .EQ. 101)then
             call intpr("i=", -1, i, 1)
         endif
         lambda_i=lambda(i)
         call ccglmreg_onelambda(x_act, y,weights, n,m_act,start_act,
     +        etastart, mustart, yhat, offset, lambda_i, alpha, gam, 
     +        rescale, standardize, intercept, penaltyfactor_act, 
     +        maxit, eps, theta, penalty, trace, iter, del, cfun, dfun,
     +        s, thresh, beta_1, b0_1, fk, delta, weights_update)
         nlambdacal=nlambdacal+1
          if(dfun .EQ. 1 .OR. dfun .EQ. 4 .OR. dfun .EQ. 5)then
              call loss2(n, y, fk, wt, cfun, dfun, s, delta, los(i))
          else
          call loss3(n,y,mustart,theta,wt,cfun,dfunnew,s,delta,
     +        los(i))
          endif
         call penGLM(beta_1, m_act, lambda_i*penaltyfactor_act, 
     +        alpha, gam, penalty, penval)
         if(standardize .EQ. 1)then
            pll(i)=los(i) + n*penval
         else 
            pll(i)=los(i) + penval
         endif
         b0(i) = b0_1
         do ii=1, n
            weights_cc(ii,i)=weights_update(ii)
         enddo
         if(jk .GT. 0)then
            do 200 ii = 1, m_act
               beta(varsel(ii), i) = beta_1(ii)
 200        continue
         endif
         if(decreasing==0 .AND. uturn==0 .AND. active .EQ. 1)then
            call find_activeset(m_act, beta_1, eps, activeset, jk)
C     this activeset is relative to the current x_act, but how about
C     relative to x instead? compute varsel for true index in x
            if(jk .NE. m_act .AND. jk .GT. 0)then
               deallocate(start_act, stat=DeAllocateStatus)
               allocate(start_act(jk+1), stat=AllocateStatus)
               deallocate(penaltyfactor_act, stat=DeAllocateStatus)
               allocate(penaltyfactor_act(jk), stat=AllocateStatus)
               start_act(1) = b0_1
               do 35 ii=1, jk
                  start_act(ii+1)=beta_1(activeset(ii))
                  varsel(ii)=varsel_old(activeset(ii))
                  varsel_old(ii)=varsel(ii)
                  penaltyfactor_act(ii)=penaltyfactor(varsel(ii))
 35            continue
               deallocate(beta_1, stat=DeAllocateStatus) 
               allocate(beta_1(jk), stat=AllocateStatus)
               do 37 ii=1, jk
               beta_1(ii)=0
 37            continue
               deallocate(x_act, stat=DeAllocateStatus)
               allocate(x_act(n, jk), stat=AllocateStatus)
C     update x_act matrix
               do 55 jj=1, n
                  do 45 ii=1, jk
                     x_act(jj, ii) = x(jj, varsel(ii))
 45               continue
 55            continue
               m_act = jk
            endif
         endif
C     redo (i-1)-lambda estimates with the current start_act for the
C     i-th lambda. Note, i=i-2 not i-1 will do this since i=i+1 is
C     computed before the end of the loop.
C      if(i > 1 .AND. uturn==0 .AND. cutpoint==1)then
C          if(abs(los(i)-los(i-1))/los(i) > epscycle)then
C            do j=1, m_act
C            start_act(j)=0
C            enddo
C            cutpoint = i
C            uturn=1
C          endif
C     endif
C      if(uturn==1)then
C         i = i - 2
C      endif
C      if(i==0 .AND. uturn==1)then
C         uturn=0
C      endif
      i = i + 1
      goto 10
      endif
      deallocate(beta_1, start_act, x_act, penaltyfactor_act)

      return
      end
