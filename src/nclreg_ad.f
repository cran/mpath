C     used in nclreg.R for active=TRUE and decreasing=TRUE
C     output: nlambdacal: number of lambda actually computed including
C     repeated ones
C     
      subroutine nclreg_ad(x, y, weights, n,m,start, etastart,
     +     mustart, offset, iter, nlambda, lambda, alpha, gam, 
     +     standardize, penaltyfactor, maxit, eps, epscycle, 
     +     penalty, trace, del,rfamily, B, s, thresh, cost, 
     +     decreasing, beta, b0, yhat, los, pll, nlambdacal)
      implicit none
      integer n,m,i,ii,j,jj,penalty,nlambda, standardize, maxit,
     +     trace, iter, rfamily, jk, activeset(m), 
     +     activeset_old(m), m_act, nlambdacal, uturn, decreasing, 
     +     cutpoint, nact, conv, jc, fakejk, AllocateStatus, varsel(m) 
      double precision x(n, m), y(n), weights(n),start(m+1),etastart(n),
     +     mustart(n), offset(n), lambda(nlambda), alpha, gam, eps, 
     +     penaltyfactor(m), thresh, beta(m, nlambda), epscycle,
     +     betaall(m), b0all, b0(nlambda), b0_1, yhat(n), del, 
     +     lambda_i, s, B, fk(n), los(nlambda), 
     +     pll(nlambda), cost, penval
      double precision, dimension(:, :), allocatable :: x_act
      double precision, dimension(:), allocatable :: start_act, beta_1,
     +     penaltyfactor_act

      b0all=start(1)
      do ii=1, m
         activeset(ii)=ii
         activeset_old(ii)=ii
         betaall(ii)=start(ii+1)
      enddo
      call find_activeset(m, start(2:(m+1)), eps, activeset, jk)
      fakejk=0
      if(jk==0)then
         jk=1
         activeset(1)=1
         fakejk=1
      endif
      do ii=1, jk
         activeset_old(ii)=activeset(ii)
      enddo
      m_act = jk
      allocate(start_act(m_act+1), stat=AllocateStatus)
      if(AllocateStatus .NE. 0)then
         return
      endif
      allocate(penaltyfactor_act(m_act), stat=AllocateStatus)
      allocate(x_act(n, m_act), stat=AllocateStatus)
      if(AllocateStatus .NE. 0)then
         return
      endif
      do jj=1, n
         do ii=1, m_act
            x_act(jj, ii)=x(jj, activeset(ii))
         enddo
      enddo

      allocate(beta_1(m_act), stat=AllocateStatus)
      if(AllocateStatus .NE. 0)then
         return
      endif
      start_act(1)=start(1)
      do j=1, m_act
         beta_1(j)=0
         start_act(j+1)=start(1+activeset(j))
         penaltyfactor_act(j)=penaltyfactor(activeset(j))
      enddo
      i=1
      nlambdacal=0
      uturn=0
      cutpoint=1
 10   if(i .LE. nlambda)then
         if(trace .EQ. 1)then
            call intpr("nclreg_ad lambda iteration i=", -1, i, 1)
         endif
         lambda_i=lambda(i)/B
         j=1
         nact=2
         conv=0
13000    if(j <= nact .AND. conv==0)then
            if(trace .EQ. 1)then
               call intpr("begin activeset nclreg_onelambda", -1, 1, 1)
            endif
            call nclreg_onelambda(x_act, y,weights, n,m_act,start_act,
     +           etastart, mustart, yhat, offset, lambda_i, alpha, gam, 
     +           penaltyfactor_act, maxit, eps, penalty, trace, iter,
     +           del, rfamily, B, s, thresh, beta_1, b0_1, fk)
            start(1)=b0_1
            do ii=1, m_act
               start(activeset(ii)+1)=beta_1(ii)
            enddo
            if(j .NE. nact)then
               if(trace .EQ. 1)then
                  call dblepr("b0_1=", -1, b0_1, 1)
                  call dblepr("beta_1=", -1, beta_1, m_act)
                  call intpr("begin fullset nclreg_onelambda", -1, 1, 1)
                  call dblepr("start=", -1, start, m+1)
               endif
               call nclreg_onelambda(x, y,weights, n,m,start, etastart,
     +              mustart, yhat, offset, lambda_i, alpha, gam, 
     +              penaltyfactor, maxit, eps, penalty, trace, 1,
     +              del, rfamily, B, s, thresh, betaall, b0all, fk)
               call find_activeset(m, betaall, eps, activeset, jk)
               if(jk==0)then
                  jk=1
                  activeset(1)=1
                  fakejk=1
               endif
            endif
         if(jk .NE. m_act)then
             conv=0
         else
            jc=0
            do ii=1, max(jk, m_act)
               if(activeset(ii)==activeset_old(ii))then
                  jc=jc+1
               endif
            enddo
            if(jk==jc)then
               conv=1
            endif
        endif
            if(jk .NE. m_act .OR. jk .NE. jc)then
               deallocate(beta_1, start_act, penaltyfactor_act, x_act)
               allocate(beta_1(jk), stat=AllocateStatus)
               allocate(start_act(jk+1), stat=AllocateStatus)
               allocate(penaltyfactor_act(jk), stat=AllocateStatus)
               allocate(x_act(n, jk), stat=AllocateStatus)
               start_act(1)=b0all
               do 11135 ii=1, jk
                  beta_1(ii)=0
                  start_act(ii+1)=betaall(activeset(ii))
                  penaltyfactor_act(ii)=penaltyfactor(activeset(ii))
                  activeset_old(ii)=activeset(ii)
11135          continue
               do 11155 jj=1, n
                  do 11145 ii=1, jk
                     x_act(jj, ii)=x(jj, activeset(ii))
11145             continue
11155          continue
               m_act = jk
            else 
               if(fakejk==1)then
                  start_act(1)=b0all
                  start_act(2)=betaall(1)
               endif
            endif
            j=j+1
            goto 13000
         endif
         nlambdacal=nlambdacal+1
         b0(i)=b0_1
         if(jk .GT. 0)then
            do ii=1, jk
               beta(activeset(ii), i)=beta_1(ii)
            enddo
         endif
         call loss(n, y, fk, cost, rfamily, s, los(i))
         call penGLM(beta_1, m_act, lambda_i*penaltyfactor_act, 
     +        alpha, gam, penalty, penval)
         if(standardize .EQ. 1)then
            pll(i)=los(i) + n*penval
         else 
            pll(i)=los(i) + penval
         endif
C     redo (i-1)-lambda estimates with the current start_act for the
C     i-th lambda. Note, i=i-2 not i-1 will do this since i=i+1 is
C     computed before the end of the loop.
         if(decreasing==1)then
            if(i > 1)then
               if(abs(pll(i)-pll(i-1))/pll(i) > epscycle)then
                  if(cutpoint==1)then
                     cutpoint = i
                  endif
                  i = i - 2
                  uturn=1
               else 
                  uturn=0
               endif
            endif
         endif
         if(cutpoint > 1 .AND. uturn==0)then
            i = cutpoint
            cutpoint = 1
         endif
         i = i + 1
         goto 10
      endif
      deallocate(beta_1, start_act, x_act, penaltyfactor_act)
      return
      end
