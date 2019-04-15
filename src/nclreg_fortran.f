C     used in nclreg.R
      subroutine nclreg_fortran(x, y, weights, n,m,start, etastart,
     +     mustart, offset, nlambda, lambda, alpha, 
     +     gam, standardize, penaltyfactor, maxit, eps, family,
     +     penalty, trace, beta, b0, yhat, iter,
     +     del, rfamily, B, s, los, pll, rescale, thresh, epsbino, 
     +     theta, cost, active)
      implicit none
      integer n,m,i,ii,k,j,jj,penalty,nlambda,family,standardize, maxit,
     +     trace, iter, rfamily, rescale, jk, active, activeset(m), 
     +     stopit,m_act,
     +     AllocateStatus, DeAllocateStatus, varsel(m), varsel_old(m)
      double precision x(n, m), y(n), weights(n),start(m+1),etastart(n),
     +     mustart(n), offset(n), lambda(nlambda), alpha, gam, eps, 
     +     penaltyfactor(m), thresh, epsbino,  theta, beta(m, nlambda), 
     +     b0(nlambda), b0_1, yhat(n), d, del, lambda_i, fk_old(n), s, 
     +     B, h(n), fk(n), los(nlambda), pll(nlambda), cost, 
     +     penval
      double precision, dimension(:, :), allocatable :: x_act
      double precision, dimension(:), allocatable :: start_act, beta_1,
     +     penaltyfactor_act

      b0_1=0
      stopit = 0
      m_act = m
      jk = m
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
 10   if(i .LE. nlambda)then
         if(trace .EQ. 1)then
            call intpr("i=", -1, i, 1)
         endif
         k = 1
         d = 10
 500     if(d .GT. del .AND. k .LE. iter)then
            if(trace .EQ. 1)then
               call intpr("  k=", -1, k, 1)
               call intpr("  activeset", -1, activeset, jk)
               call dblepr("     d=", -1, d, 1)
            endif
            call dcopy(n, yhat, 1, fk_old, 1)
            call compute_h(rfamily, n, y, fk_old, s, B, h)
C     check if h has NAN value
            do 30 ii=1, n
               if(h(ii) .NE. h(ii))then
                  call intpr("i=", -1, i, 1)
                  call intpr("  k=", -1, k, 1)
                  call intpr("    ii=", -1, ii, 1)
                  call dblepr("     h(ii)", -1, h(ii), 1)
                  call dblepr("     fk_old(ii)", -1, fk_old(ii), 1)
                  stopit = 1
                  exit
               endif
 30         continue
            lambda_i = lambda(i)/B
            call glmreg_fit_fortran(x_act, h,weights,n,m_act,start_act, 
     +           etastart, mustart, offset, 1, lambda_i, alpha,
     +           gam, rescale, standardize, penaltyfactor_act, thresh,
     +           epsbino, maxit, eps, theta, family,  
     +           penalty, trace, beta_1, b0_1, yhat)
            call dcopy(n, yhat, 1, fk, 1)
            call dcopy(n, yhat, 1, mustart, 1)
C     call dcopy(n, yhat, 1, etastart, 1): this is not needed for
C     the above glmref_fit call (that calls zeval) since family=1.
            start_act(1) = b0_1
            d = 0
            d=d+(start_act(1)-b0_1)**2
            if(jk .GT. 0)then
               do 100 j=1, m_act
                  d=d+(start_act(j+1)-beta_1(j))**2
                  start_act(j+1)=beta_1(j)
 100           continue
            endif

            k = k + 1
            goto 500
         endif
         call loss(n, y, fk, cost, rfamily, s, los(i))
         penval = 0.d0
         call penGLM(beta_1, m_act, lambda_i*penaltyfactor_act, 
     +        alpha, gam, penalty, penval)
         if(standardize .EQ. 1)then
            pll(i)=los(i) + n*penval
         else 
            pll(i)=los(i) + penval
         endif
         b0(i) = b0_1
         if(jk .GT. 0)then
            do 200 ii = 1, m_act
               beta(varsel(ii), i) = beta_1(ii)
 200        continue
         endif
         if(active .EQ. 1)then
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
         i = i + 1
         goto 10
      endif
      deallocate(beta_1, start_act, x_act, penaltyfactor_act)

      return
      end
