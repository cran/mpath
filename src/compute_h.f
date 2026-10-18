C converted from R mpath/R/compute.h function
C input rfamily, y, fk_old, s, B -- no changes
C       clossR: 11, closs: 12, gloss: 13, qloss: 14
C output h
      subroutine compute_h(rfamily, n, y, fk_old, s, B, h)
              implicit none
              integer rfamily, i, n
              double precision y(n), fk_old(n), s, B, u(n), g, h(n)

      if(rfamily .EQ. 11)then
              do 10 i=1, n
              u(i)=y(i)-fk_old(i)
              call gradient(rfamily, 1, u(i), s, g)
              h(i) = g/B + fk_old(i)
 10           continue
      else if (rfamily .EQ. 12 .OR. rfamily .EQ. 13 
     +         .OR. rfamily .EQ. 14) then
              do 30 i=1, n
              u(i)=y(i)*fk_old(i)
              call gradient(rfamily, 1, u(i), s, g)
              h(i) = -y(i)*g/B + fk_old(i)
  30          continue
      end if
             
            return
            end 
