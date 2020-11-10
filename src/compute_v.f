C compute subdifferential of -g(z), step 4 in the COCO Algorithm
C used in ccglmreg_onelambda 
C output v
      subroutine compute_v(cfun, n, z, s, delta, v)
              implicit none
              integer cfun, i, n
              double precision PI, s, delta, s1, s2, z(n), v(n)

       if(cfun .EQ. 1)then
           s1 = s**2/2
       else if(cfun .EQ. 2)then
           PI=4.D0*DATAN(1.D0)
           s1 = s**2*PI**2/2
       else if(cfun .EQ. 3)then
           s1 = s**2
           s2 = 1.0/s1**2
       else if(cfun .EQ. 4)then
           s1 = s**2
       else if(cfun .EQ. 5)then
           s1 = dexp(-s)
       else if(cfun .EQ. 6)then
           s1=-delta**(s-1)/(delta+1)**(s+1)
       else if(cfun .EQ. 8)then
           PI=4.D0*DATAN(1.D0)
           s1 = 2*dexp(-delta/s)/dsqrt(PI*s*delta)
       endif

      do 10 i=1, n
      if(cfun .EQ. 1)then
          if(z(i) .LE. s1)then
              v(i)=-1
          else
              v(i)=-s*(2*z(i))**(-0.5)
          endif
      else if(cfun .EQ. 2)then
          if(z(i)==0)then
              v(i)=-1
          else if(z(i) > s1)then
              v(i)=0
          else
              v(i)=-s*dsin((2*z(i))**(0.5)/s)/(2*z(i))**(0.5)
          endif
      else if(cfun .EQ. 3)then
          if(z(i) .LE. s1/2)then
              v(i)=-s2*(2*z(i)-s1)**2
          else
              v(i)=0.D0
          endif
      else if(cfun .EQ. 4)then
              v(i)=-dexp(-z(i)/s1)
      else if(cfun .EQ. 5)then
          v(i)=-1/((z(i)+1)*(z(i)*s1+1))
      else if(cfun .EQ. 6)then
          if(z(i) .LE. delta)then
           v(i)=s1
          else
           v(i)=-z(i)**(s-1)/(z(i)+1)**(s+1)
          endif
      else if(cfun .EQ. 7)then
          if(z(i) .LE. s)then
              v(i)=-1
          else
              v(i)=0
          endif
      else if(cfun .EQ. 8)then
          if(z(i) .LE. delta)then
           v(i)=-s1
          else
           v(i)=-2*dexp(-z(i)/s)/dsqrt(PI*s*z(i))
          endif
C              if(v(i) .gt. 0)then
C               call   rexit("the v value should be less than 0")
      end if
 10   continue
             
            return
            end 
