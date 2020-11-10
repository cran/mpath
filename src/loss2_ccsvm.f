C loss value in subroutine ccglmreg, cf: loss.f
C cost: Not implemented in the optimization problem!
C output: los is the average loss
C compare to loss2 subroutine, there is additional input eps for dfun= 2
C (eps-regression). cf: ../R/loss2_ccsvm function
      subroutine loss2_ccsvm(n, y, f, weights, cfun, dfun, s, eps,
     +      delta, los)
       implicit none
       integer cfun, dfun, i, n
       double precision y(n),f(n),u(n),z(n),gval(n),s,eps,los,delta,
     +   weights(n)
       external compute_g, compute_z, compute_u 

      call compute_u(dfun, n, y, f, u)
      call compute_z(dfun, n, u, z, eps)
      call compute_g(cfun, n, z, s, delta, gval)
      los = 0.d0
      do 20 i=1, n
      los = los + weights(i)*gval(i)
20    continue

      return
      end

