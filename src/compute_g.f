C compute value of g(z) as a function of z
C input: z, cfun, s
C output: gval
      subroutine compute_g(cfun, n, z, s, delta, gval)
              implicit none
              integer cfun, i, n
              double precision PI, s, s1, s2, z(n), delta, gval(n)

      if(cfun .EQ. 1 .OR. cfun .EQ. 3)then
          s1 = s**2/2
      else if(cfun .EQ. 2)then
          PI=4.D0*DATAN(1.D0)
          s1 = s**2*PI**2/2
      else if(cfun .EQ. 4)then
          s1 = s**2
      else if(cfun .EQ. 5)then
          s1 = dexp(-s)
      else if(cfun .EQ. 6)then
          s1 = delta**(s-1)/(1+delta)**(s+1)
          s2 = -1.0D0/s*(delta/(1+delta))**s + delta**s/(1+delta)**(s+1)
      else if(cfun .EQ. 8)then
          PI=4.D0*DATAN(1.D0)
          s1 = 2*dexp(-delta/s)/dsqrt(PI*s*delta)
          s2 = s1*delta - erf(dsqrt(delta/s))
      endif
     
      do 10 i=1, n
      if(cfun .EQ. 1)then
        if(z(i) .LE. s1)then
          gval(i)=z(i)
        else 
          gval(i)=s*(2*z(i))**0.5-s1
        endif
      else if(cfun .EQ. 2)then
        if(z(i) .LE. s1)then
          gval(i)=s**2*(1-dcos((2*z(i))**0.5/s))
        else 
          gval(i)=2*s**2
        endif
      else if(cfun .EQ. 3)then
           if(z(i) .LE. s1)then
               gval(i)=1-(1-2*z(i)/s**2)**3
               gval(i)=gval(i)*s**2/6.0
           else
               gval(i)=s**2/6.0
           endif
      else if(cfun .EQ. 4)then
        gval(i)=1-dexp(-z(i)/s1)
        gval(i)=gval(i)*s1
      else if(cfun .EQ. 5)then
        gval(i)=dlog((1+z(i))/(1+z(i)*s1))
        gval(i)=gval(i)/(1-s1)
      else if(cfun .EQ. 6)then
          if(z(i) > delta)then
            gval(i)=(z(i)/(1+z(i)))**s/s + s2
          else
            gval(i)=s1*z(i)
          endif
      else if(cfun .EQ. 7)then
            if(z(i) .LE. s)then
                gval(i)=z(i)
            else
                gval(i)=s
            endif
      else if(cfun .EQ. 8)then
          if(z(i) > delta)then
            gval(i)=erf(dsqrt(z(i)/s)) + s2
          else
            gval(i)=s1*z(i)
          endif
      endif
 10   continue
             
            return
            end 
