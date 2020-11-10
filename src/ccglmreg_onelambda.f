C     used in ccglmreg_fortran subroutine for one single lambda penalty
C     parameter (adapted from nclreg_onelambda.f)
C     input: start_act, etastart, mustart
C     output: start_act, etastart, mustart, beta_1, b0_1, fk
      subroutine ccglmreg_onelambda(x_act,y,weights,n,m_act,start_act, 
     +     etastart, mustart, yhat, offset, lambda_i,alpha,gam,rescale,
     +     standardize,intercept,penaltyfactor_act,maxit,eps,theta,
     +   penalty,trace,iter,del,cfun,dfun,s,thresh,beta_1,b0_1,fk,delta,
     +   weights_update)
      implicit none
      integer n,k,j,intercept,penalty,maxit,trace,iter,m_act,
     +        satu, cfun, dfun, dfunnew, i, ii, standardize, rescale
      double precision y(n), weights(n),etastart(n), mustart(n), 
     +     offset(n),lambda_i,alpha,gam,eps,
     +     thresh, b0_1, yhat(n), d, d1, del, fk_old(n), s, 
     +     fk(n), x_act(n, m_act), start_act(m_act+1),
     +     beta_1(m_act), penaltyfactor_act(m_act),
     +     weights_update(n), delta, ytmp(n), theta
      external update_wt, loss2, penGLM

      if(dfun==5)then
         do i=1, n
           ytmp(i)=(y(i)+1.0D0)/2.0D0
         enddo
       else 
         do i=1, n
           ytmp(i)=y(i)
         enddo
      endif
        k = 1
         d = 10
            do 10 ii=1, n
            weights_update(ii)=weights(ii)
 10         continue
 500     if(d .GT. del .AND. k .LE. iter)then
            if(trace .EQ. 1)then
               call intpr("  ccglmreg_onelambda iteration k=", -1, k, 1)
               call dblepr("     start_act", -1, start_act, m_act+1)
            endif
            call dcopy(n, yhat, 1, fk_old, 1)
CCCCCCCCCCCCCCCCCCC BEGIN
C instead of updating h values in nclreg_onelambda, update
C weights: call compute_h(rfamily, n, y, fk_old, s, B, h)
C     compute u
C            call compute_u(dfun, n, y, fk_old, u)
C     compute z=s(u)
C            call compute_z(dfun, n, u, z, s)
C     compute derivative of -g(z)
C            call compute_v(cfun, n, z, s, delta, v)
C     update weights
C            do 10 ii=1, n
C            weights_update(ii)=weights(ii)*v(ii)
C 10         continue
CCCCCCCCCCCCCCCCCCC END
C     check if h has NAN value, can be useful to debug
C            if(trace .EQ. 1)then
C            do 30 ii=1, n
C               if(h(ii) .NE. h(ii))then
C                  call intpr("  ccglmreg_onelambda iteration k=",-1,k,1)
C                  call intpr("    ii=", -1, ii, 1)
C                  call dblepr("     h(ii)", -1, h(ii), 1)
C                  call dblepr("     fk_old(ii)", -1, fk_old(ii), 1)
C                  exit
C               endif
C 30         continue
C            endif
C    unlike nclreg_onelambda, we update weights not h (here is y); also we may have different family depending on the dfun function
            if(dfun .EQ. 1 .OR. dfun .EQ. 4)then
                dfunnew = 1
            else if(dfun .EQ. 5)then 
                dfunnew= 2
              else if(dfun .EQ. 8)then 
                dfunnew= 3
              else if(dfun .EQ. 9)then 
                dfunnew= 4
              else
                  call rexit("not implemented yet")
            endif
            call glmreg_fit_fortran(x_act, ytmp,weights_update,n,m_act,
     +           start_act,etastart,mustart, offset, 1, lambda_i, alpha,
     +           gam,rescale,standardize,intercept, penaltyfactor_act, 
     +           thresh, 0.0D0, maxit, eps, theta, dfunnew,
     +           penalty, trace, beta_1, b0_1, yhat, satu)
            if(dfun .EQ. 1 .OR. dfun .EQ. 4 .OR. dfun .EQ. 5)then
                call update_wt(n,weights,y,etastart,cfun,dfun,s,delta,
     +          weights_update)
            else if(dfun .EQ. 8 .OR. dfun .EQ. 9)then
             call compute_wt3(n,y,mustart,weights,theta,cfun,dfunnew,s,
     +        delta,weights_update)
            endif
            call dcopy(n, yhat, 1, fk, 1)
            call dcopy(n, yhat, 1, mustart, 1)
            if(dfun .NE. 1 .OR. dfun .NE. 4)then
              call dcopy(n, yhat, 1, etastart, 1)
            endif
C     this is not needed for
C     the above glmreg_fit call (that calls zeval) if family=1.
        if(dfun .EQ. 1 .OR. dfun .EQ. 4)then
            start_act(1) = b0_1
            if(m_act .GT. 0)then
               do 100 j=1, m_act
                  start_act(j+1)=beta_1(j)
 100           continue
            endif
        endif
           d = 0.0D0
            do 120 ii=1, n
               d=d+(fk_old(ii) - fk(ii))**2.0D0
  120       continue
           d1 = 0.0D0
            do 130 ii=1, n
               d1=d1+fk(ii)**2.0D0
  130       continue
               d=d/d1
            if(trace .EQ. 1)then
               call dblepr("beta_1", -1, beta_1, m_act)
            endif
            k = k + 1
            goto 500
         endif

      return
      end
