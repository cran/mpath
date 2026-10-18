C compute density function of Poisson distribution with parameter lambda
C if log_true=1, log density value is computed
      function dpois(x, lambda, log_true)
      integer :: x, j, log_true, Factorial
      double precision :: lambda, dpois, res
      external :: Factorial

      if(lambda .LT. 0)then
              call dblepr("lambda should be nonnegative", -1, lambda, 1)
      endif
      if(log_true==1)then
              res=0
          if(x > 0)then 
              do j=1, x
               res=res+dlog(j*1.0D0)
              enddo
          endif
        dpois = -lambda+x*1.0D0*dlog(lambda)- res
      else 
              dpois = dexp(-lambda)*lambda**x/Factorial(x)
      endif
      return
      end function
