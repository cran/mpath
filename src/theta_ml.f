C adapted from MASS/negbin theta.ml function but
C input n is different from theta.ml there
      subroutine theta_ml(y, mu, n, weights, limit, eps, t0, trace)

      implicit none
      integer i, n, it, limit, trace
      double precision a, t0, y(n), mu(n), weights(n), 
     +       eps, del, score, info
      external score, info
      
C     eps = .Machine$double.eps**0.25 
      eps = 0.0001220703 
      a=0
      do i=1, n
       a=a+weights(i)*(y(i)/mu(i) - 1.0D0)**2
      enddo
      t0 = sum(weights)/a
      it = 1
      del = 1
      if(trace .EQ. 1)then 
              call dblepr("initial theta=", -1, t0, 1)
      endif

30    if(it .LE. limit .AND. dabs(del) > eps)then
        t0 = dabs(t0)
        del=score(n, t0, mu, y, weights)/info(n,t0, mu,y,weights)
        t0 = t0 + del
        if(trace .EQ. 1)then 
              call intpr("theta iteration", -1, it, 1)
              call dblepr("        theta=", -1, t0, 1)
        endif
        it = it + 1
        goto 30
      endif
      if(t0 < 0)then
              t0 = 0
              call intpr("estimate truncated at zero", -1, 1, 1)
      endif
      if(it .EQ. limit .AND. trace .EQ. 1)then
              call intpr("iteration limit reached", -1, 1, 1)
      endif
      
      return
      end

      function score(n, th, mu, y, w)
      implicit none
      integer i, n
      double precision y(n), mu(n), th, w(n), rdigamma, 
     + score
      external :: rdigamma

      score = 0
      do i=1, n
       score=score+w(i)*(rdigamma(th + y(i)) - rdigamma(th) + dlog(th) +
     +          1.0D0 - dlog(th + mu(i)) - (y(i) + th)/(mu(i) + th))
      enddo
      
      return
      end function 

      function info(n, th, mu, y, w)
      implicit none
      integer i, n
      double precision y(n), mu(n), th, w(n), info, 
     + rtrigamma
      external :: rtrigamma 
      
      info = 0
      do i=1, n
       info=info+w(i)*(-rtrigamma(th + y(i)) + rtrigamma(th) - 
     +      1.0D0/th + 2.0D0/(mu(i) + th) - 
     +      (y(i) + th)/(mu(i) + th)**2.0D0)
      enddo
      
      return
      end function 

