C ###compute composite loss values for dfun=poisson and negbin
C input: family=3,4. May work for 1 but not for 2 since y in ccglmreg is +1/-1 not
C 1/0)
C ### output: z is maximum log-likelihood value, i.e., s(u) value 
C ### output: los is the sum of g(s(u)) values
      subroutine loss3(n,y,mu,theta,weights,cfun,family,s,delta,los)
      implicit none
        integer cfun, family, n, i
        double precision weights(n),y(n),mu(n),z,gval,theta,s,delta,los
        external compute_g, loglikFor
 
       los=0.0D0
       do i=1, n
        call loglikFor(1,y(i),mu(i),theta,1.0D0,family,z)
        call compute_g(cfun,1,z,s,delta,gval)
        los=los+weights(i)*gval
       enddo

       return
       end
