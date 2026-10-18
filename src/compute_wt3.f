C ###compute updated weights for composite loss function with dfun=poisson and negbin
C update weights. cf subroutine update_wt for other dfun
C output: weights_update
      subroutine compute_wt3(n,y,mu,weights,theta,cfun,family,s,delta,
     +       weights_update)
       implicit none
       integer cfun, family, n, i
       double precision weights(n),weights_update(n),y(n),mu(n),z(n),
     + v(n), theta,s,delta
       external compute_v, loglikFor
 
       do i=1, n
       call loglikFor(1,y(i),mu(i),theta,1.0D0,family,z(i))
C convert MLE of log-likelihood to minimization of negative log-likelihood
       z(i)=-z(i)
       enddo
C     compute derivative of -g(z)
             call compute_v(cfun, n, z, s, delta, v)
C     update weights
             do 10 i=1, n
              weights_update(i)=-weights(i)*v(i)
10         continue
      return
      end
