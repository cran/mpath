C     used in zipath.R
C     inputs: family: 3 (poisson), 4 (negbin)
C             theta
C     outputs: coefc, coefz, thetaout
      subroutine zipath_nonactive(x, z, y, y1, weights, n, kx, kz, 
     +     start_count, start_zero, mustart_count, mustart_zero, 
     +     offsetx, offsetz, nlambda, lambda_count,
     +     lambda_zero, alpha_count, alpha_zero,  
     +     gam_count, gam_zero, standardize, penaltyfactor_count, 
     +     penaltyfactor_zero, maxit, eps, family,
     +     penalty, trace, coefc, coefz, yhat, iter,
     +     del, rescale, thresh, epsbino, 
     +     theta_fixed, maxit_theta, theta, thetaout)
      implicit none
      integer n,i,ii,j,kx, kz, penalty,nlambda,family, 
     +     standardize, maxit, y1(n), trace, iter, 
     +     rescale, jk_count, jk_zero, stopit,m_count_act, maxit_theta,
     +     m_zero_act, AllocateStatus, satu, 
     +     varsel_count(kx), varsel_count_old(kx),
     +     varsel_zero(kz), varsel_zero_old(kz), theta_fixed
      double precision x(n, kx), z(n, kz), weights(n), 
     +     start_count(kx+1), dpois, dnbinom, start_zero(kz+1), 
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count(nlambda), thetaout(nlambda),
     +     lambda_zero(nlambda), alpha_count, alpha_zero, gam_count, 
     +     gam_zero, eps, penaltyfactor_count(kx), y(n),
     +     penaltyfactor_zero(kz), probi(n), thresh, epsbino, 
     +     theta, coefc(kx+1, nlambda), coefz(kz+1, nlambda), b0_x, b0z,
     +     yhat(n), del, start_count_act(kx+1), start_zero_act(kz+1), 
     +     x_act(n, kx), z_act(n, kz), betax(kx), betaz(kz),
     +     penaltyfactor_count_act(kx), penaltyfactor_zero_act(kz)
      external :: dpois, dnbinom, gfunc

      stopit = 0
      m_count_act = kx
      m_zero_act = kz
      jk_count = kx
      jk_zero = kz

      do ii=1, n
            if(y1(ii) .EQ. 1)then
               probi(ii)=0
            else
              probi(ii)=mustart_zero(ii) 
              if(family .EQ. 3)then
              probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +                  mustart_count(ii), 0))
              else if(family .EQ. 4)then
              probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dnbinom(0,
     +                  theta, mustart_count(ii), 0))
              endif
            endif
      enddo
      call copymatrix(n, kx, x, x_act)
      call copymatrix(n, kz, z, z_act)

      do 5 j=1, kx
         betax(j)=0
 5    continue
      do 105 j=1, kz
         betaz(j)=0
 105    continue
      call dcopy(kx+1, start_count, 1, start_count_act, 1)
      call dcopy(kz+1, start_zero, 1, start_zero_act, 1)
      call dcopy(kx, penaltyfactor_count, 1, penaltyfactor_count_act, 1)
      call dcopy(kz, penaltyfactor_zero, 1, penaltyfactor_zero_act, 1)
      do 101 j=1, kx
         varsel_count_old(j)=j
         varsel_count(j)=j
 101  continue
      do 102 j=1, kz
         varsel_zero_old(j)=j
         varsel_zero(j)=j
 102  continue

      i=1
 10   if(i .LE. nlambda)then
        if(trace .EQ. 1)then
            call intpr("Fortran lambda iteration i=", -1, i, 1)
        endif
        call zi_onelambda(x_act, z_act, y, y1, probi, weights, n, 
     +           m_count_act, m_zero_act, start_count_act, 
     +           start_zero_act, mustart_count, mustart_zero, offsetx, 
     +           offsetz, lambda_count(i), lambda_zero(i), alpha_count,
     +           alpha_zero, gam_count, gam_zero, standardize,
     +           penaltyfactor_count_act, penaltyfactor_zero_act,
     +           maxit, eps, family, penalty, trace, yhat, iter, del,
     +           rescale, thresh, epsbino, theta_fixed,
     +           maxit_theta, theta, betax, b0_x, betaz, b0z, satu)
         coefc(1, i) = b0_x
         if(jk_count .GT. 0)then
            do 200 ii = 1, m_count_act
               coefc(1+varsel_count(ii), i) = betax(ii)
 200        continue
         endif
         thetaout(i)=theta
         coefz(1, i) = b0z
         if(satu==1 .AND. i==1)then
             do ii=1, m_zero_act
             start_zero_act(ii)=0
             enddo
         else
         if(jk_zero .GT. 0)then
            do 210 ii = 1, m_zero_act
               coefz(1+varsel_zero(ii), i) = betaz(ii)
 210        continue
         endif
         endif
         i = i + 1
         goto 10
      endif

      return
      end
