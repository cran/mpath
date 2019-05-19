C     see vignettes\kkt.pdf for how to find the lambda_max for both the
C     original loss and in EM algorithm 
C     For the EM algorith, uncomment related lines with respect to derz1 and
C     derx1
C     inputs: B with # of columns kx, G with # of columns kz are
C     count-model and zero-model design matrix with the first columns 1s
C     for intercept, this makes the subroutine valid even if there is an
C     intercept-only component in count-model or zero-model.
C     outputs: lmax_count, lmax_zero
      subroutine lmax_zipath(B, G, y, y1, weights, n, kx, kz, family, 
     +     coefc, coefz, alpha_count, alpha_zero, penaltyfactor_count, 
     +     penaltyfactor_zero, theta, lmax_count, lmax_zero)
      implicit none
      integer i, j, k, n, kx, kz, family, y1(n)
      double precision B(n, kx), G(n, kz), y(n), weights(n), coefc,
     +     coefz, alpha_count, alpha_zero, penaltyfactor_zero(kz), 
     +     penaltyfactor_count(kx), lmax_count, lmax_zero, derz(kz),
     +     derx(kx), penfac_zero(kz), penfac_count(kx), sumw, sumx, 
     +     sumz, p, derip, mu, a, c, theta, z(n)
C     +     , derz1(kz), derx1(kx)
      
      if(kz<1)then
           call rexit("zero model design matrix needs at least one
     +          column including intercept")
      endif
      if(kx<1)then
           call rexit("count model design matrix needs at least one
     +          column including intercept")
      endif
      do j=1, kz
         derz(j)=0
C     derz1(j)=0
      enddo
      do k=1, kx
         derx(k)=0
C     derx1(k)=0
      enddo
      sumw = sum(weights)
      do i=1, n
         weights(i)=weights(i)/sumw
      enddo
      sumx=sum(penaltyfactor_count)
      sumz=sum(penaltyfactor_zero)
      p=dexp(coefz)/(1+dexp(coefz))
      mu=dexp(coefc)
      if(family==3)then
         a=  dexp(coefz)/(dexp(coefz)+dexp(-dexp(coefc)))
         c=  dexp(-dexp(coefc))*
     +        (-dexp(coefc))/(dexp(coefz)+dexp(-dexp(coefc)))
         if(kz > 1)then
            do 10 j=2, kz
               do 20 i=1, n
                  if(y1(i)==0)then
                     derz(j)=derz(j)+weights(i)*a*G(i,j)
                  endif
                  derz(j)=derz(j)-weights(i)*p*G(i,j)
 20            enddo
               derz(j)=dabs(derz(j))
               penfac_zero(j)=penaltyfactor_zero(j)/sumz*kz
               derz(j)=derz(j)/(penfac_zero(j)*alpha_zero)
 10         enddo
         endif
         if(kx > 1)then
            do 30 j=2, kx
               do 40 i=1, n
                  if(y1(i)==0)then
                     derx(j)=derx(j)+weights(i)*c*B(i,j)
                  else
                     derx(j)=derx(j)+weights(i)*(y(i)*B(i,j)
     +                    -mu*B(i,j))
                  endif
 40            enddo
               derx(j)=dabs(derx(j))
               penfac_count(j)=penaltyfactor_count(j)/sumx*kx
               derx(j)=derx(j)/(penfac_count(j)*alpha_count)
 30         enddo
         endif
      else if(family==4)then
         a=(theta/(mu+theta))**theta
         do i=1, n
            if(y1(i)==0)then
               z(i)=1/(1+dexp(-coefz)*a)
            else 
               z(i)=0
            endif
         enddo
         if(kz > 1)then
            do 110 j=2, kz
               do 120 i=1, n
C     derivative of p
                  derip=G(i,j)*p**2/dexp(coefz)
                  if(y1(i)==0)then
                     derz(j)=derz(j)+weights(i)*derip*(1-a)/(p+(1-p)*a)
                  else
                     derz(j)=derz(j)-weights(i)*derip/(1-p)
                  endif
C     derz1(j)=derz1(j)+weights(i)*(z(i)*G(i,j)-p*G(i,j))
 120           enddo
               derz(j)=dabs(derz(j))
               penfac_zero(j)=penaltyfactor_zero(j)/sumz*kz
               derz(j)=derz(j)/(penfac_zero(j)*alpha_zero)
               
C     derz1(j)=dabs(derz1(j))
C     derz1(j)=derz1(j)/(penfac_zero(j)*alpha_zero)
 110        enddo
         endif
         if(kx > 1)then
            do 130 j=2, kx
               do 140 i=1, n
C     derivative of mu
                  derip=B(i,j)*mu
                  if(y1(i)==0)then
                     derx(j)=derx(j)-weights(i)*(1-p)*derip*theta*
     +                    a/(mu+theta)/(p+(1-p)*a)
                  else
                     derx(j)=derx(j)+weights(i)*derip*
     +                    (y(i)/mu-(y(i)+theta)/(mu+theta))
                  endif
C     derx1(j)=derx1(j)+weights(i)*(1-z(i))*derip*
C     +       (y(i)/mu-(y(i)+theta)/(mu+theta))
 140           enddo
               derx(j)=dabs(derx(j))
               penfac_count(j)=penaltyfactor_count(j)/sumx*kx
               derx(j)=derx(j)/(penfac_count(j)*alpha_count)
C     derx1(j)=dabs(derx1(j))
C     derx1(j)=derx1(j)/(penfac_count(j)*alpha_count)
 130        enddo
         endif
      endif

      lmax_count = maxval(derx)
      lmax_zero = maxval(derz)
C      call dblepr("lmax_count", -1, lmax_count, 1)
C     call dblepr("lmax_count with EM", -1, maxval(derx1), 1)
C      call dblepr("lmax_zero", -1, lmax_zero, 1)
C     call dblepr("lmax_zero with EM", -1, maxval(derz1), 1)
      return
      end
