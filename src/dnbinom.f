C compute density function of negative biomial distribution with
C parameter size and mu, in the same parameterization as in R dnbinom.
C rlgamma calls C function lgammafn: log gamma function, C Mathlib
C function
      function dnbinom(x, size_n, mu)
      integer :: x, Factorial
      double precision :: p, size_n, mu, dnbinom, rlgamma
      external :: Factorial, rlgamma

      if(size_n .LE. 0)then
              call intpr("size should be strictly positive")
      endif

      p = size_n/(size_n+mu)
      dnbinom = dexp(rlgamma(x+size_n)-rlgamma(size_n))/Factorial(x)
     +*p**size_n*(1-p)**x

      return
      end function
