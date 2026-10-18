C update weights in R/ccsvm.R
C output: weights_update
      subroutine update_wt_ccsvm(n, weights, y, f, cfun, dfun, s, eps,
     +     delta, weights_update)
       implicit none
       integer cfun, dfun, i, n
       double precision weights(n),weights_update(n), y(n),f(n),
     + u(n),z(n), v(n), s, eps, delta
       external compute_v, compute_z, compute_u 
 
C     compute u
             call compute_u(dfun, n, y, f, u)
C     compute z=s(u)
             call compute_z(dfun, n, u, z, eps)
C     compute derivative of -g(z)
             call compute_v(cfun, n, z, s, delta, v)
C     update weights
             do 10 i=1, n
              weights_update(i)=-weights(i)*v(i)
10         continue
      return
      end
