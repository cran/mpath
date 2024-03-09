C ###compute composite loss values for dfun=poisson and negbin
C input: family=3,4. May work for 1 but not for 2 since y in irglmreg is +1/-1 not
C 1/0)
C ### output: z is negative log-likelihood value at the current
C parameter
C , i.e., s(u) value 
C ### output: los is the sum of g(s(u)) values
      subroutine loss3(n,y,mu,theta,weights,cfun,family,s,delta,los)
      implicit none
        integer cfun, family, n, i
        double precision weights(n),y(n),mu(n),z,gval,theta,s,delta,los
        double precision z_saturate, dfun_val
        external compute_g, loglikFor
 
       los=0.0D0
       do i=1, n
C ### z is log-likelihood value at the current
        call loglikFor(1,y(i),mu(i),theta,1.0D0,family,z)
C ### z_saturate is log-likelihood value in the saturate model
        call loglikFor(1,y(i),y(i),theta,1.0D0,family,z_saturate)
C  minimization problem, dfun s(x) must be convex, thus take negative value
C  of z, but -z can be negative, thus add z_saturate such that dfun_val
C  >=0
        dfun_val = z_saturate - z
        call compute_g(cfun,1,dfun_val,s,delta,gval)
        los=los+weights(i)*gval
       enddo

       return
       end
