C used in ccglmreg_onelambda, step 3 in the COCO Algorithm
C input dfun, if dfun=1-3, for regression; otherwise, if dfun=4-7 for
C classification, if dfun=8 or 9, Poisson or negative binomial
C output u
      subroutine compute_u(dfun, n, y, f, u)
              implicit none
              integer dfun, i, n
              double precision y(n), f(n), u(n)

      do 10 i=1, n
      if(dfun .LE. 3)then
              u(i)=y(i)-f(i)
          else if(dfun .GE. 4 .AND. dfun .LE. 7)then
              u(i)=y(i)*f(i)
          else 
              u(i)=f(i)
      end if
 10           continue
             
      return
      end 
