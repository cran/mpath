C generate outer product between two vectors, similar to %o% in R
C input: m, x, n, y
C output: z
      subroutine outprod(m, x, n, y, z)
              implicit none
              integer m, n, i, j
              double precision x(m), y(n), z(m, n)

      do 10 i=1, m
         do 20 j=1, n
              z(i, j) = x(i)*y(j)
  20  continue
  10  continue
      
      return
      end
