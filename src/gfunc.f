C input: mu, family, epsbino
C output: eta
      subroutine gfunc(mu, n, family, epsbino, eta)
      implicit none
      integer :: family, n, i
      double precision :: mu(n), epsbino, eta(n)

      do i=1, n
       if(family==1)then
         eta(i) = mu(i)
       else if(family==2)then
         if(1-mu(i) > epsbino .AND. mu(i) > epsbino)then
          eta(i)=dlog(mu(i)/(1-mu(i)))
         else if(mu(i) .LE. epsbino)then
           eta(i)=dlog(epsbino/(1-epsbino))
         else if(mu(i) .GE. 1-epsbino)then
          eta(i)=dlog(1-epsbino/epsbino)
         endif
       else if(family==3 .OR. family==4)then
          eta(i)=dlog(mu(i))
       endif
      enddo

      return
      end
