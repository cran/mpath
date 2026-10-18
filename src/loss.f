C loss subroutine converted from R package bst in bst.R function
C family = c("gaussian", "binom", "poisson"), s=-1, sh, fk=NULL
C      gausssian: 1, bino: 2, poisson:3, negative binomial: 4, clossR: 11, closs: 12, gloss:
C      13, qloss: 14
C sh and fk are not used below, but may needed in the future release
C cost: but not implemented in optimization problem, unlike in bst!
C output: los is the average loss

      subroutine loss(n, y, f, family, s, los)
       implicit none
       integer family, i, n
       double precision x, y(n),f(n),u, s, los
       external nonconvexloss

C      if(family .EQ. 2)then
C           do 10 i=1, n
C               if(y(i) .EQ. 1.0)then
C                 ly(i)=1-cost
C               else  
C                 ly(i)=1.0D0
C                 ly(i)=cost
C               endif
C10         continue 
C      endif
      los = 0.d0
      do 20 i=1, n
        if(family .EQ. 1)then
          los = los + 0.5D0*(y(i) - f(i))**2
        else if(family .EQ. 2)then
          los =  los + dlog(1+dexp(-y(i)*f(i))) 
C          los =  los + ly(i)*dlog(1+dexp(-y(i)*f(i))) 
        else if(family .EQ. 11)then
                u=y(i)-f(i)
          call nonconvexloss(family, u, s, x)
          los = los + x
        else if(family .EQ. 12 .OR. family .EQ. 13 .OR. 
     +          family .EQ. 14)then
                u=y(i)*f(i)
                 call nonconvexloss(family, u, s, x)
                 los = los + x
        endif
20    continue
        los=los/n

      return
      end

C ###clossR for C-loss regression with u= f - y
C ###closs  for C-loss classification with y={1, 1} and u=y * f
C family
C   clossR: 11, closs: 12, gloss: 13, qloss: 14
      subroutine nonconvexloss(family, u, s, los)
       implicit none
       integer family
       double precision cval, u, s, los
       external pnorm_fortran

       if(family .EQ. 11)then
        los=1-1/dexp(u**2/(2*s**2))
         else if(family .EQ. 12)then
           cval = 1/(1 - dexp(-1/(2*s**2)))
           los=cval*(1-1/dexp((1-u)**2/(2*s**2)))
         else if(family .EQ. 13)then
           los=2**s/(dexp(u)+1)**s 
         else if(family .EQ. 14)then
                 u=u/s
                 call pnorm_fortran(u)
                 los=2*(1-u)
       endif

        return
        end

C compute pnorm as in R, but we have the erf() function available in
C Fortran; it's a linear transformation of pnorm().  That is,
C erf(x) = 2*pnorm(x*sqrt(2)) - 1

        subroutine pnorm_fortran(x)
                implicit none
                double precision x
                
                  x=(erf(x/sqrt(2.0))+1)/2
                return
        end
