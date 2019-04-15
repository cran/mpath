C     used in zipath.R
C     inputs: family: 3 (poisson), 4 (negbin)
C     theta
C     kx: number of variables of x having no intercept column
C     kz: number of variables of z having no intercept column
C     outputs: coefc, coefz, thetaout
      subroutine zipath_active(x, z, y, y1, weights, n, kx, kz, 
     +     start_count, start_zero, mustart_count, mustart_zero, 
     +     offsetx, offsetz, nlambda, lambda_count,
     +     lambda_zero, alpha_count, alpha_zero,  
     +     gam_count, gam_zero, standardize, penaltyfactor_count, 
     +     penaltyfactor_zero, maxit, eps, family,
     +     penalty, trace, coefc, coefz, yhat, iter,
     +     del, rescale, thresh, epsbino, 
     +     theta_fixed, maxit_theta, theta, thetaout)
      implicit none
      integer n,i,ii,j,jj,kx, kz, penalty,nlambda,family, 
     +     standardize, maxit, y1(n), trace, iter, 
     +     rescale, jk_count, jk_zero, activeset_count(kx), 
     +     activeset_count_old(kx), activeset_zero(kz),
     +     activeset_zero_old(kz), stopit,m_count_act, maxit_theta,
     +     m_zero_act, AllocateStatus, DeAllocateStatus, jc, jz, 
     +     nact, conv, theta_fixed
      double precision x(n, kx), z(n, kz), weights(n), 
     +     start_count(kx+1), dpois, dnbinom, b0_xall, b0zall, 
     +     start_zero(kz+1), etastart_count(n), etastart_zero(n),
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count(nlambda), thetaout(nlambda),
     +     lambda_zero(nlambda), alpha_count, alpha_zero, gam_count, 
     +     gam_zero, eps, penaltyfactor_count(kx), y(n),
     +     penaltyfactor_zero(kz), probi(n), thresh, epsbino, 
     +     theta, coefc(kx+1, nlambda), coefz(kz+1, nlambda), b0_x, b0z,
     +     yhat(n), del, betaxall(kx), betazall(kz)
      double precision, dimension(:, :), allocatable :: x_act, z_act
      double precision, dimension(:), allocatable :: start_count_act,
     +     start_zero_act, betax, betaz,
     +     penaltyfactor_count_act, penaltyfactor_zero_act
      external :: dpois, dnbinom, gfunc

      b0_x=0
      b0z=0
      stopit = 0
C     find current active set
      call find_activeset(kx, start_count(2: (kx+1)), eps, 
     +     activeset_count, jk_count)
      call find_activeset(kz, start_zero(2: (kz+1)), eps, 
     +     activeset_zero, jk_zero)
C     When all coef are zero except intercept, choose a predictor
      if(jk_count==0)then
         jk_count = 1
         activeset_count_old(1)=1
         activeset_count(1)=1
      endif
      if(jk_zero==0)then
         jk_zero = 1
         activeset_zero_old(1)=1
         activeset_zero(1)=1
      endif
      m_count_act = jk_count
      m_zero_act = jk_zero

      call gfunc(mustart_count, n, family, epsbino, etastart_count)
      call gfunc(mustart_zero, n, 2, epsbino, etastart_zero)
      do ii=1, n
         if(y1(ii) .EQ. 1)then
            probi(ii)=0
         else
            probi(ii)=mustart_zero(ii) 
            if(family .EQ. 3)then
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +              mustart_count(ii)))
            else if(family .EQ. 4)then
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dnbinom(0,
     +              theta, mustart_count(ii)))
            endif
         endif
      enddo
      allocate(start_count_act(jk_count+1), stat=AllocateStatus)
      allocate(penaltyfactor_count_act(jk_count),stat=AllocateStatus)
      allocate(betax(jk_count), stat=AllocateStatus)
      allocate(x_act(n, jk_count), stat=AllocateStatus)
      start_count_act(1)=start_count(1)
      do jj=1, n
         do ii=1, jk_count
            x_act(jj, ii)=x(jj, activeset_count(ii))
         enddo
      enddo
      do 5 j=1, jk_count
         betax(j)=0
         start_count_act(j+1)=start_count(1+activeset_count(j))
         penaltyfactor_count_act(j)=
     +        penaltyfactor_count(activeset_count(j))
 5    continue
      
      allocate(start_zero_act(jk_zero+1), stat=AllocateStatus)
      allocate(penaltyfactor_zero_act(jk_zero), stat=AllocateStatus)
      allocate(betaz(jk_zero), stat=AllocateStatus)
      allocate(z_act(n, jk_zero), stat=AllocateStatus)
      start_zero_act(1)=start_zero(1)
      do jj=1, n
         do ii=1, jk_zero
            z_act(jj, ii)=z(jj, activeset_zero(ii))
         enddo
      enddo
      do 105 j=1, jk_zero
         betaz(j)=0
         start_zero_act(j+1)=start_zero(1+activeset_zero(j))
         penaltyfactor_zero_act(j)=
     +        penaltyfactor_zero(activeset_zero(j))
 105  continue
      i=1
 10   if(i .LE. nlambda)then
         if(trace .EQ. 1)then
            call intpr("Fortran lambda iteration i=", -1, i, 1)
            call intpr("m_count_act", -1, m_count_act, 1)
            call intpr("m_zero_act", -1, m_zero_act, 1)
            call dblepr("start_count_act", -1, start_count_act,
     +           m_count_act+1)
            call dblepr("start_zero_act", -1, start_zero_act,
     +           m_zero_act+1)
            call dblepr("betax", -1, betax, m_count_act)
            call dblepr("betaz", -1, betaz, m_zero_act)
         endif
         j = 1
         nact = 2
         conv=0
13000    if(j <= nact .AND. conv==0)then
            if(trace==1)then
            call intpr("active set iteration, j=", -1, j, 1)
            call intpr("cycling through only active sets", -1, 1, 1)
            endif
            call zi(x_act, z_act, y, y1, probi, weights, n, m_count_act,
     +           m_zero_act, start_count_act, start_zero_act, 
     +           mustart_count, mustart_zero, offsetx, offsetz, 
     +           lambda_count(i), lambda_zero(i), alpha_count, 
     +           alpha_zero, gam_count, gam_zero, standardize, 
     +           penaltyfactor_count_act, penaltyfactor_zero_act,
     +           maxit, eps, family, penalty, trace, yhat, iter, del,
     +           rescale, thresh, epsbino, theta_fixed, 
     +           maxit_theta, theta, betax, b0_x, betaz, b0z)
C     update start_count with start_count_act, start_zero with
C     start_zero_act
               start_count(1)=b0_x
            do ii=1, jk_count
               start_count(activeset_count(ii)+1)=betax(ii)
            enddo
               start_zero(1)=b0z
            do ii=1, jk_zero
               start_zero(activeset_zero(ii)+1)=betaz(ii)
            enddo
            if(j .NE. nact)then
             call zi(x, z, y, y1, probi, weights, n, kx,
     +           kz, start_count, start_zero, mustart_count, 
     +           mustart_zero, offsetx, offsetz, lambda_count(i), 
     +           lambda_zero(i), alpha_count, alpha_zero, gam_count,
     +           gam_zero, standardize, penaltyfactor_count, 
     +           penaltyfactor_zero, maxit, eps, family, penalty, 
     +           trace, yhat, 1, del, rescale, thresh, 
     +           epsbino, theta_fixed, maxit_theta, theta, betaxall, 
     +           b0_xall, betazall, b0zall)
             call find_activeset(kx, betaxall, eps, activeset_count
     +           , jk_count)
             if(jk_count==0)then
               jk_count = 1
               activeset_count(1)=1
             endif
             call find_activeset(kz, betazall, eps, activeset_zero
     +           , jk_zero)
             if(jk_zero==0)then
               jk_zero = 1
               activeset_zero(1)=1
             endif
            endif
C     check if converged here!
            if(trace==1)then
            call intpr("activeset_count", -1, activeset_count, jk_count)
            call intpr("activeset_count_old",
     +           -1,activeset_count_old,m_count_act)
            call intpr("activeset_zero", -1, activeset_zero, jk_zero)
            call intpr("activeset_zero_old", -1,activeset_zero_old,
     +           m_zero_act)
            call intpr("jk_count=", -1, jk_count, 1)
            call intpr("m_count_act=", -1, m_count_act, 1)
            call intpr("jk_zero=", -1, jk_zero, 1)
            call intpr("m_zero_act=", -1, m_zero_act, 1)
            endif
            jc=0
            do ii=1, jk_count
               if(activeset_count(ii)==activeset_count_old(ii))then
                  jc=jc+1
               endif
            enddo
            jz=0
            do ii=1, jk_zero
               if(activeset_zero(ii)==activeset_zero_old(ii))then
                  jz=jz+1
               endif
            enddo
            if(jk_count==jc .AND. jk_zero==jz)then
               conv=1
            endif
            if(jk_count .NE. jc)then
               deallocate(start_count_act, stat=DeAllocateStatus)
               allocate(start_count_act(jk_count+1),stat=AllocateStatus)
               deallocate(penaltyfactor_count_act,stat=DeAllocateStatus)
               allocate(penaltyfactor_count_act(jk_count),stat=
     +              AllocateStatus)
               start_count_act(1) = b0_xall
               do 11135 ii=1, jk_count
                  start_count_act(ii+1)=betaxall(activeset_count(ii))
                  activeset_count_old(ii)=activeset_count(ii)
                  penaltyfactor_count_act(ii)=
     +                 penaltyfactor_count(activeset_count(ii))
11135          continue
               deallocate(betax, stat=DeAllocateStatus) 
               allocate(betax(jk_count), stat=AllocateStatus)
               deallocate(x_act, stat=DeAllocateStatus)
               allocate(x_act(n, jk_count), stat=AllocateStatus)
C     update x_act matrix
               do 11155 jj=1, n
                  do 11145 ii=1, jk_count
                     x_act(jj, ii) = x(jj, activeset_count(ii))
11145             continue
11155          continue
               m_count_act = jk_count
            endif
            if(jk_zero .NE. jz)then
               deallocate(start_zero_act, stat=DeAllocateStatus)
               allocate(start_zero_act(jk_zero+1),stat=AllocateStatus)
               deallocate(penaltyfactor_zero_act,stat=DeAllocateStatus)
               allocate(penaltyfactor_zero_act(jk_zero),stat=
     +              AllocateStatus)
               start_zero_act(1) = b0zall
               do 12135 ii=1, jk_zero
                  start_zero_act(ii+1)=betazall(activeset_zero(ii))
                  activeset_zero_old(ii)=activeset_zero(ii)
                  penaltyfactor_zero_act(ii)=
     +                 penaltyfactor_zero(activeset_zero(ii))
12135          continue
               deallocate(betaz, stat=DeAllocateStatus) 
               allocate(betaz(jk_zero), stat=AllocateStatus)
               deallocate(z_act, stat=DeAllocateStatus)
               allocate(z_act(n, jk_zero), stat=AllocateStatus)
               do 12155 jj=1, n
                  do 12145 ii=1, jk_zero
                     z_act(jj, ii) = z(jj, activeset_zero(ii))
12145             continue
12155          continue
               m_zero_act = jk_zero
            endif
            j=j+1
            goto 13000
         endif
         coefc(1, i) = b0_x
         if(jk_count .GT. 0)then
            do 200 ii = 1, m_count_act
               coefc(1+activeset_count(ii), i) = betax(ii)
 200        continue
         endif
         thetaout(i)=theta
         coefz(1, i) = b0z
         if(jk_zero .GT. 0)then
            do 210 ii = 1, m_zero_act
               coefz(1+activeset_zero(ii), i) = betaz(ii)
 210        continue
         endif
         i = i + 1
         goto 10
      endif
      deallocate(betax, start_count_act, x_act, 
     +     penaltyfactor_count_act)
      deallocate(betaz, start_zero_act, z_act, 
     +     penaltyfactor_zero_act)

      return
      end
