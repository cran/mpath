C     used in zipath.R
C     inputs: family: 3 (poisson), 4 (negbin)
C     theta
C     kx: number of variables of x having no intercept column
C     kz: number of variables of z having no intercept column
C     outputs: coefc, coefz, theta, thetaout
      subroutine zipath_active(x, z, y, y1, weights, n, kx, kz, 
     +     start_count, start_zero, mustart_count, mustart_zero, 
     +     offsetx, offsetz, intercept, nlambda, lambda_count,
     +     lambda_zero, alpha_count, alpha_zero,  
     +     gam_count, gam_zero, penaltyfactor_count, 
     +     penaltyfactor_zero, maxit, eps, family,
     +     penalty, trace, coefc, coefz, yhat, iter,
     +     del, rescale, thresh, epsbino, 
     +     theta_fixed, maxit_theta, theta, thetaout)
      implicit none
      integer n,i,ii,j,jj,kx, kz, penalty,nlambda,family, 
     +     maxit, y1(n), trace, iter, intercept, 
     +     rescale, jk_count, jk_zero, activeset_count(kx), 
     +     activeset_count_old(kx), activeset_zero(kz),
     +     activeset_zero_old(kz), m_count_act, maxit_theta,
     +     m_zero_act, AllocateStatus, jc, jz, 
     +     nact, conv, theta_fixed, fakec, fakez, satu
      double precision x(n, kx), z(n, kz), weights(n), 
     +     start_count(kx+1), dpois, dnbinom, b0_xall, b0zall, 
     +     start_zero(kz+1),
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count(nlambda), thetaout(nlambda),
     +     lambda_zero(nlambda), alpha_count, alpha_zero, gam_count, 
     +     gam_zero, eps, penaltyfactor_count(kx), y(n),
     +     penaltyfactor_zero(kz), probi(n), thresh, epsbino, 
     +     theta, thetaall, coefc(kx+1, nlambda), coefz(kz+1, nlambda),
     +     b0_x, b0z, yhat(n), del, betaxall(kx), betazall(kz)
      double precision, dimension(:, :), allocatable :: x_act, z_act
      double precision, dimension(:), allocatable :: start_count_act,
     +     start_zero_act, betax, betaz,
     +     penaltyfactor_count_act, penaltyfactor_zero_act
      external :: dpois, dnbinom, gfunc
      
      if(kx==0 .OR. kz==0)then
         return
      endif
      do ii=1, kx
         betaxall(ii)=0
         activeset_count(ii)=ii
         activeset_count_old(ii)=ii
      enddo
      do ii=1, kz
         betazall(ii)=0
         activeset_zero(ii)=ii
         activeset_zero_old(ii)=ii
      enddo
C     find current active set
      call find_activeset(kx, start_count(2: (kx+1)), eps, 
     +     activeset_count, jk_count)
      call find_activeset(kz, start_zero(2: (kz+1)), eps, 
     +     activeset_zero, jk_zero)
C     When all coef are zero except intercept, choose a predictor
      fakec=0
      fakez=0
      jz=0
      jc=0
      if(jk_count==0)then
         jk_count = 1
         activeset_count(1)=kx
         fakec=1
      endif
      if(jk_zero==0)then
         jk_zero = 1
         activeset_zero(1)=kz
         fakez=1
      endif
      do ii=1, jk_count
         activeset_count_old(ii)=activeset_count(ii)
      enddo
      do ii=1, jk_zero
         activeset_zero_old(ii)=activeset_zero(ii)
      enddo
      
      m_count_act = jk_count
      m_zero_act = jk_zero

      do ii=1, n
         if(y1(ii) .EQ. 1)then
            probi(ii)=0
         else
            probi(ii)=mustart_zero(ii) 
            if(family .EQ. 3)then
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +              mustart_count(ii), 0))
            else if(family .EQ. 4)then
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dnbinom(0,
     +              theta, mustart_count(ii), 0))
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
      nact = 2
      i=1
 10   if(i .LE. nlambda)then
         if(trace .EQ. 1)then
C            i is not array
C            call intpr("Fortran lambda iteration i=", -1, i, 1), 
            call intpr("kx=", -1, kx, 1)
            call intpr("kz=", -1, kz, 1)
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
         conv=0
13000    if(j <= nact .AND. conv==0)then
            if(trace==1)then
               call intpr("active set iteration, j=", -1, j, 1)
               call intpr("cycling through only active sets", -1, 1, 1)
            endif
            call zi_onelambda(x_act, z_act, y, y1, probi, weights, n, 
     +           m_count_act, m_zero_act, start_count_act, 
     +           start_zero_act, mustart_count, mustart_zero, offsetx, 
     +           offsetz, intercept, lambda_count(i), lambda_zero(i), 
     +           alpha_count, alpha_zero, gam_count, gam_zero, 
     +           penaltyfactor_count_act, penaltyfactor_zero_act,
     +           maxit, eps, family, penalty, trace, yhat, iter, del,
     +           rescale, thresh, epsbino, theta_fixed, maxit_theta, 
     +           theta, betax, b0_x, betaz, b0z, satu)
C     update start_count with start_count_act, start_zero with
C     start_zero_act
            start_count(1)=b0_x
            do ii=1, jk_count
               start_count(activeset_count(ii)+1)=betax(ii)
            enddo
            start_zero(1)=b0z
            if(satu==1 .AND. i==1)then
               do ii=1, jk_zero
                  start_zero_act(1+ii)=0
                  betaz(ii)=0
               enddo
            endif
            do ii=1, jk_zero
               start_zero(activeset_zero(ii)+1)=betaz(ii)
            enddo
            if(j .NE. nact)then
               thetaall = theta
               call zi_onelambda(x, z, y, y1, probi, weights, n, kx,
     +              kz, start_count, start_zero, mustart_count, 
     +              mustart_zero, offsetx, offsetz, intercept, 
     +              lambda_count(i), lambda_zero(i), alpha_count, 
     +              alpha_zero, gam_count,gam_zero,penaltyfactor_count, 
     +              penaltyfactor_zero, maxit, eps, family, penalty, 
     +              trace, yhat, 2, del, rescale, thresh, epsbino,
     +              theta_fixed, maxit_theta, thetaall, betaxall, 
     +              b0_xall, betazall, b0zall, satu)
               if(satu==1 .AND. i==1)then
                  do ii=1, kz
                     start_zero(1+ii)=0
                     betazall(ii)=0
                  enddo
               endif
               call find_activeset(kx, betaxall, eps, activeset_count
     +              , jk_count)
               if(jk_count==0)then
                  jk_count = 1
                  activeset_count(1)=kx
                  fakec=1
               endif
               call find_activeset(kz, betazall, eps, activeset_zero
     +              , jk_zero)
               if(jk_zero==0)then
                  jk_zero = 1
                  activeset_zero(1)=kz
                  fakez=1
               endif
            endif
C     check if converged here!
            if(trace==1)then
               call intpr("activeset_count",-1,activeset_count,jk_count)
               call intpr("activeset_count_old",
     +              -1,activeset_count_old,m_count_act)
               call intpr("activeset_zero", -1, activeset_zero, jk_zero)
               call intpr("activeset_zero_old", -1,activeset_zero_old,
     +              m_zero_act)
               call intpr("jk_count=", -1, jk_count, 1)
               call intpr("m_count_act=", -1, m_count_act, 1)
               call intpr("jk_zero=", -1, jk_zero, 1)
               call intpr("m_zero_act=", -1, m_zero_act, 1)
            endif
            if(jk_count .NE. m_count_act)then
               conv=0
            else if(jk_zero .NE. m_zero_act)then
               conv=0
            else
               jc=0
               do ii=1, max(jk_count, m_count_act)
                  if(activeset_count(ii)==activeset_count_old(ii))then
                     jc=jc+1
                  endif
               enddo
               jz=0
C     what if activeset_zero==activeset_zero_old after only one iteration as
C     in the above? Say, we only have betaz=0.2 with one element, but
C     jk_zero=kz now. Then we have a problem here to update the values as
C     stated below! Because in this case, the iteration stops.
               do ii=1, max(jk_zero, m_zero_act)
                  if(activeset_zero(ii)==activeset_zero_old(ii))then
                     jz=jz+1
                  endif
               enddo
               if(jk_count==jc .AND. jk_zero==jz)then
                  conv=1
               endif
            endif
            if(jk_count .NE. m_count_act .OR. jk_count .NE. jc)then
               theta = thetaall
               deallocate(betax, start_count_act, 
     +              penaltyfactor_count_act, x_act)
               allocate(betax(jk_count), stat=AllocateStatus)
               allocate(start_count_act(jk_count+1),stat=AllocateStatus)
               allocate(penaltyfactor_count_act(jk_count),stat=
     +              AllocateStatus)
               start_count_act(1) = b0_xall
               do 11135 ii=1, jk_count
                  betax(ii)=0
                  start_count_act(ii+1)=betaxall(activeset_count(ii))
                  activeset_count_old(ii)=activeset_count(ii)
                  penaltyfactor_count_act(ii)=
     +                 penaltyfactor_count(activeset_count(ii))
11135          continue
               allocate(x_act(n, jk_count), stat=AllocateStatus)
C     update x_act matrix
               do 11155 jj=1, n
                  do 11145 ii=1, jk_count
                     x_act(jj, ii) = x(jj, activeset_count(ii))
11145             continue
11155          continue
               m_count_act = jk_count
            else
               if(fakec==1)then
                  start_count_act(1)=b0_xall
                  start_count_act(2)=betaxall(1)
               endif
            endif
            if(jk_zero .NE. m_zero_act .OR. jk_zero .NE. jz)then
               deallocate(betaz, start_zero_act, 
     +              penaltyfactor_zero_act, z_act)
               allocate(betaz(jk_zero), stat=AllocateStatus)
               allocate(start_zero_act(jk_zero+1),stat=AllocateStatus)
               allocate(penaltyfactor_zero_act(jk_zero),stat=
     +              AllocateStatus)
               start_zero_act(1) = b0zall
               do 12135 ii=1, jk_zero
                  betaz(ii)=0
                  start_zero_act(ii+1)=betazall(activeset_zero(ii))
                  activeset_zero_old(ii)=activeset_zero(ii)
                  penaltyfactor_zero_act(ii)=
     +                 penaltyfactor_zero(activeset_zero(ii))
12135          continue
               allocate(z_act(n, jk_zero), stat=AllocateStatus)
               do 12155 jj=1, n
                  do 12145 ii=1, jk_zero
                     z_act(jj, ii) = z(jj, activeset_zero(ii))
12145             continue
12155          continue
               m_zero_act = jk_zero
            else
               if(fakez==1)then
                  start_zero_act(1)=b0zall
                  start_zero_act(2)=betazall(1)
               endif
            endif
            j=j+1
            goto 13000
         endif
         coefc(1, i) = b0_x
         if(jk_count .GT. 0)then
            do 200 ii = 1, jk_count
               coefc(1+activeset_count(ii), i) = betax(ii)
 200        continue
         endif
         thetaout(i)=theta
         coefz(1, i) = b0z
         if(jk_zero .GT. 0)then
            do 210 ii = 1, jk_zero
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
