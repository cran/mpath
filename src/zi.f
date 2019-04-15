C     fit a zero inflated penalized regression for a single penalty
C     parameter. A building block used in zipath_fortran.f
C     inputs: family: 3 (poisson), 4 (negbin)
C       theta
C       m_count_act: number of count model variables of x having no intercept column
C       m_zero_act: number of zero model variables of z having no intercept column
C     outputs: betax, b0_x, betaz, b0z, theta
      subroutine zi(x_act, z_act, y, y1, probi, weights, n, m_count_act,
     +     m_zero_act, start_count_act, start_zero_act,mustart_count, 
     +     mustart_zero, offsetx, offsetz, lambda_count,
     +     lambda_zero, alpha_count, alpha_zero,  
     +     gam_count, gam_zero, standardize, penaltyfactor_count_act, 
     +     penaltyfactor_zero_act, maxit, eps, family,
     +     penalty, trace, yhat, iter, del, rescale, thresh, 
     +     epsbino, theta_fixed, maxit_theta, theta, 
     +     betax, b0_x, betaz, b0z)
      implicit none
      integer n,ii,k,j,penalty, family, 
     +     standardize, maxit, y1(n), trace, iter, 
     +     rescale, stopit,m_count_act, maxit_theta,
     +     m_zero_act, theta_fixed
      double precision weights(n), dpois, dnbinom, 
     +     etastart_count(n), etastart_zero(n),
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count, thetastart,
     +     lambda_zero, alpha_count, alpha_zero, gam_count, 
     +     gam_zero, eps, penaltyfactor_count_act(m_count_act), y(n),
     +     penaltyfactor_zero_act(m_zero_act), wt(n), probi(n), thresh, 
     +     epsbino, theta, b0_x, b0z,yhat(n), d, del, theta0,
     +     penval, x_act(n, m_count_act), 
     +     z_act(n, m_zero_act), 
     +     start_count_act(m_count_act+1), start_zero_act(m_zero_act+1),
     +     betax(m_count_act), betaz(m_zero_act)
      external :: dpois, dnbinom, gfunc

      b0_x=0
      b0z=0
      stopit = 0
         k = 1
         d = 10
 500     if(d .GT. del .AND. k .LE. iter)then
            if(trace .EQ. 1)then
               call intpr("  EM algorithm iteration k=", -1, k, 1)
               call dblepr("     d=", -1, d, 1)
               call dblepr("start_count_act=", -1,
     +                    start_count_act,m_count_act+1)
               call dblepr("start_zero_act=", -1,
     +                    start_zero_act,m_zero_act+1)
            endif
            do 30 ii=1, n
               wt(ii)=weights(ii)*(1-probi(ii)) 
 30         continue
            if(family .NE. 4 .OR. theta_fixed .EQ. 1)then
               call glmreg_fit_fortran(x_act,y,wt,n,m_count_act,
     +              start_count_act,etastart_count,mustart_count,
     +              offsetx,1, lambda_count, alpha_count, gam_count,
     +              rescale,0, penaltyfactor_count_act, thresh,
     +              epsbino, maxit, eps, theta, family,  
     +              penalty, 0, betax, b0_x, yhat)
            else
               thetastart = theta 
               call glmregnb_fortran(x_act,y,wt,n,m_count_act,offsetx,
     +              1, lambda_count, penalty,alpha_count, gam_count, 
     +              rescale, 0, penaltyfactor_count_act, thresh,
     +              maxit_theta, maxit, eps, epsbino, start_count_act, 
     +              etastart_count, mustart_count, thetastart, 0, 
     +              theta0, trace, betax, b0_x, theta, yhat)
            endif
C     yhat: the fitted mean values, obtained by transforming the
C     linear predictors by the inverse of the link function.
            call dcopy(n, yhat, 1, mustart_count, 1)
C     etastart_count: linear predictors, obtained by transforming the
C     mean values by the link function.
            call gfunc(mustart_count, n, family,epsbino,etastart_count)
            d = 0
            d=d+(start_count_act(1) - b0_x)**2
            start_count_act(1) = b0_x
            if(m_count_act .GT. 0)then
               do 100 j=1, m_count_act
                  d=d+(start_count_act(j+1) - betax(j))**2
                  start_count_act(j+1)=betax(j)
 100           continue
            endif
            do ii=1, n
               yhat(ii)=0
            enddo
            call glmreg_fit_fortran(z_act,probi,weights,n,m_zero_act,
     +           start_zero_act, etastart_zero, mustart_zero,offsetz,
     +           1, lambda_zero, alpha_zero, gam_zero, rescale,
     +           0, penaltyfactor_zero_act, thresh,
     +           epsbino, maxit, eps, theta, 2,  
     +           penalty, 0, betaz, b0z, yhat)
            call dcopy(n, yhat, 1, mustart_zero, 1)
            call gfunc(mustart_zero, n, 2, epsbino, etastart_zero)
            do ii=1, n
               if(y1(ii) .EQ. 0)then
                  probi(ii)=mustart_zero(ii) 
                  if(family .EQ. 3)then
                   probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +                    mustart_count(ii)))
                  else if(family .EQ. 4)then
                    probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*
     +                    dnbinom(0, theta, mustart_count(ii)))
                  endif
               endif
            enddo
            
            d=d+(start_zero_act(1) - b0z)**2
            start_zero_act(1) = b0z
            if(m_zero_act .GT. 0)then
               do 1100 j=1, m_zero_act
                  d=d+(start_zero_act(j+1) - betaz(j))**2
                  start_zero_act(j+1)=betaz(j)
 1100          continue
            endif
            if(trace .EQ. 1)then
               penval = 0.d0
               call penGLM(betax, m_count_act, 
     +              lambda_count*penaltyfactor_count_act, 
     +              alpha_count, gam_count, penalty, penval)
C               missing computing los
               if(standardize .EQ. 1)then
C                  pll(k)=los(k) + n*penval
C               else 
C                  pll(k)=los(k) + penval
               endif
            endif
            k = k + 1
            goto 500
            endif

            return
            end
