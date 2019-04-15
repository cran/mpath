! Lemmon and Schafer: Developing Statistical Software in Fortran 95,
! Springer 2005, page 39. No argument checking or overflow checking
!-----Factorial------------------------------------------------------
!
!  Function to calculate factorials resursively
!
!---------------------------------------------------------------------
      RECURSIVE INTEGER FUNCTION Factorial(n)  RESULT(Fact)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n

      IF (n == 0) THEN
         Fact = 1
      ELSE
         Fact = n * Factorial(n-1)
      END IF

      END FUNCTION Factorial
