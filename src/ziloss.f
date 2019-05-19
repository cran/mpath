C     loss subroutine converted from R functions ziPoisson and ziNegBin
C     in package mpath/zipath.R
C     family=3 - ziPoisson; 4 - ziNegBin
C     output: los is the average loss

      subroutine ziloss(n, y, offsetx, offsetz, weights, fc, fz, 
     +           family, theta, los)
      implicit none
      integer family, i, n
      double precision x, y(n), offsetx(n), offsetz(n), weights(n), 
     +fc(n), fz(n), mu, phi, theta, los, dpois, dnbinom
      external dpois, dnbinom

      los = 0.d0
      do 20 i=1, n
         mu=dexp(fc(i)+offsetx(i))
         call linkinv(1, fz(i)+offsetz(i), 2, phi)
         if(family .EQ. 3)then
            if(y(i) < 1)then
               los=los+weights(i)*dlog(phi+dexp(dlog(1.0D0-phi)-mu))
            else
               los=los+weights(i)*(dlog(1.0D0-phi)
     +          +dpois(int(y(i)), mu, 1))
            endif
         else if(family .EQ. 4)then
            if(y(i) < 1)then
               los=los+weights(i)*dlog(phi+dexp(dlog(1.0D0-phi)+
     +                 dnbinom(0, theta, mu, 1)))
            else
               los=los+weights(i)*(dlog(1.0D0-phi)+
     +          dnbinom(int(y(i)),theta,mu,1))
            endif
         endif
 20   continue

      return
      end
