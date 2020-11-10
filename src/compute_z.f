C compute z=s(u), step 3 in the COCO Algorithm
C used in ccglmreg_onelambda, 
C input dfun, n, u(n), s -- no changes
C       1: u^2/2
C output z(n)
      subroutine compute_z(dfun, n, u, z, eps)
              implicit none
              integer dfun, i, n
              double precision u(n), z(n), eps, s1

      if (dfun .EQ. 3)then
          s1=eps**2/2
      end if
      do 10 i=1, n
      if(dfun .EQ. 1)then
              z(i)=u(i)**2/2
      else if (dfun .EQ. 2)then
          if (dabs(u(i)) .LE. eps)then
              z(i)=0.0D0
          else
              z(i)=dabs(u(i))-eps
          endif
      else if (dfun .EQ. 3)then
          if (dabs(u(i)) .LE. eps)then
              z(i)=u(i)**2/2
          else
              z(i)=eps*dabs(u(i))-s1
          endif
      else if (dfun .EQ. 4)then
              z(i)=(1-u(i))**2/2
      else if (dfun .EQ. 5)then
          if(u(i) .GT. -10.0D0)then
              z(i)=dlog(1+dexp(-u(i)))
          else
              z(i)=-u(i)
          endif
      else if (dfun .EQ. 6)then
              z(i)=max(0.D0, 1-u(i))
      else if (dfun .EQ. 7)then
              z(i)=dexp(-u(i))
      end if
 10           continue
             
            return
            end 
