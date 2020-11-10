C loss value in subroutine ccglmreg, cf: loss.f
C dfun=1,4,5
C input: f is linear predictors, not probability if dfun=5
C output: los is the average loss

      subroutine loss2(n, y, f, weights, cfun, dfun, s, delta, los)
       implicit none
       integer cfun, dfun, i, n
       double precision y(n),f(n),u(n),z(n),gval(n),s,los,weights(n),
     +   delta
       external compute_g, compute_z, compute_u 

      if(dfun .NE. 1 .AND. dfun .NE. 4 .AND. dfun .NE. 5)then
          call rexit("dfun not implmented in loss2")
      endif
      call compute_u(dfun, n, y, f, u)
C     s is not used in compute_z for dfun:1,4-7
      call compute_z(dfun, n, u, z, s)
      call compute_g(cfun, n, z, s, delta, gval)
      los = 0.d0
      do 20 i=1, n
      los = los + weights(i)*gval(i)
20    continue
C      los=los/n

      return
      end

