C     used in nclreg_fortran subroutine for one single lambda penalty
C     parameter
C     input: start_act, etastart, mustart
C     output: start_act, etastart, mustart, beta_1, b0_1, fk
      subroutine nclreg_onelambda(x_act, y, weights, n,m_act,start_act, 
     +     etastart, mustart, yhat, offset, lambda_i, alpha, gam, 
     +     penaltyfactor_act, maxit, eps, penalty, trace, iter,
     +     del, rfamily, B, s, thresh, beta_1, b0_1, fk)
      implicit none
      integer n,k,j, penalty,maxit,trace,iter,rfamily,m_act,satu
      double precision y(n), weights(n),etastart(n), mustart(n), 
     +     offset(n),lambda_i,alpha,gam,eps,los,penval,pll,pll_old,
     +     thresh, b0_1, yhat(n), d, del, fk_old(n), s, B, 
     +     h(n), fk(n), x_act(n, m_act), start_act(m_act+1),
     +     beta_1(m_act), penaltyfactor_act(m_act)

         k = 1
         d = 10
         call loss(n, y, yhat, 0.5D0, rfamily, s, los)
         call penGLM(beta_1, m_act, lambda_i*penaltyfactor_act,
     +        alpha, gam, penalty, penval)
         pll_old=los + penval
 500     if(d .GT. del .AND. k .LE. iter)then
            if(trace .EQ. 1)then
               call intpr("  nclreg_onelambda iteration k=", -1, k, 1)
               call dblepr("     start_act", -1, start_act, m_act+1)
            endif
            call dcopy(n, yhat, 1, fk_old, 1)
            call compute_h(rfamily, n, y, fk_old, s, B, h)
C     check if h has NAN value, can be useful to debug
C            if(trace .EQ. 1)then
C            do 30 ii=1, n
C               if(h(ii) .NE. h(ii))then
C                  call intpr("  nclreg_onelambda iteration k=",-1,k,1)
C                  call intpr("    ii=", -1, ii, 1)
C                  call dblepr("     h(ii)", -1, h(ii), 1)
C                  call dblepr("     fk_old(ii)", -1, fk_old(ii), 1)
C                  exit
C               endif
C 30         continue
C            endif
            call glmreg_fit_fortran(x_act, h,weights,n,m_act,start_act, 
     +           etastart, mustart, offset, 1, lambda_i, alpha,
     +           gam, 0, 0, penaltyfactor_act, thresh,
     +           0.0D0, maxit, eps, 0.0D0, 1,
     +           penalty, trace, beta_1, b0_1, yhat, satu)
            call dcopy(n, yhat, 1, fk, 1)
            call dcopy(n, yhat, 1, mustart, 1)
C     call dcopy(n, yhat, 1, etastart, 1): this is not needed for
C     the above glmref_fit call (that calls zeval) since family=1.
C            d = 0
C            d=d+(start_act(1)-b0_1)**2
            start_act(1) = b0_1
            if(m_act .GT. 0)then
               do 100 j=1, m_act
C                  d=d+(start_act(j+1)-beta_1(j))**2
                  start_act(j+1)=beta_1(j)
 100           continue
            endif
         call loss(n, y, yhat, 0.5D0, rfamily, s, los)
         call penGLM(beta_1, m_act, lambda_i*penaltyfactor_act,
     +        alpha, gam, penalty, penval)
         pll=los + penval
         d=dabs((pll-pll_old)/pll_old)
         pll_old = pll
C           d = 0
C            do 120 ii=1, n
C               d=d+(fk_old(ii) - fk(ii))**2
C  120       continue
            if(trace .EQ. 1)then
               call dblepr("b0_1", -1, b0_1, 1)
               call dblepr("beta_1", -1, beta_1, m_act)
               call dblepr("updated start_act", -1, start_act, m_act+1)
               call dblepr("     d=", -1, d, 1)
            endif
            k = k + 1
            goto 500
         endif

      return
      end
