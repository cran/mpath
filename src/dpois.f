C compute density function of Poisson distribution with parameter lambda
      function dpois(x, lambda)
      integer :: x, Factorial
      double precision :: lambda, dpois
      external :: Factorial

      if(lambda .LT. 0)then
              call intpr("lambda should be nonnegative")
      endif
      dpois = dexp(-lambda)*lambda**x/Factorial(x)

      return
      end function
