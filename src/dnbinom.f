C compute density function of negative biomial distribution with
C parameter size and mu, in the same parameterization as in R dnbinom.
C rlgamma calls C function lgammafn: log gamma function, C Mathlib
C function
      function dnbinom(x, size_n, mu, log_true)
      integer :: x, j, Factorial, log_true
      double precision :: res, p, size_n, mu, dnbinom, rlgamma
      external :: Factorial, rlgamma

      if(size_n .LE. 0)then
        call dblepr("size should be strictly positive", -1, size_n, 1)
      endif

      p = size_n/(size_n+mu)
      if(log_true==0)then
      dnbinom=dexp(rlgamma(x*1.0D0+size_n)-rlgamma(size_n))/Factorial(x)
     +        *p**size_n*(1-p)**x
      else
        res=0 
        if(x > 0)then
            do j=1, x
            res=res+dlog(j*1.0D0)
            enddo
        endif
        dnbinom = rlgamma(x*1.0D0+size_n)-rlgamma(size_n) - res
     +  +size_n*dlog(p)+x*1.0D0*dlog(1-p) 
      endif

      return
      end function
