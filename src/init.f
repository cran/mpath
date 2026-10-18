C### eta is the estimated beta_0 in the intercept-only model, derived
Cfrom ../R/init function
      subroutine init(n, wt, y, offset, family, mu, eta)
      implicit none
      integer i, n, family
      double precision wt(n), y(n), offset(n), mu1, mu(n), eta(n), ddot

C compute weighted means sum(wt_i * y_i)
      mu1 = ddot(n, y, 1, wt, 1)
      do 20 i= 1, n
           mu(i) = mu1 + offset(i)
 20   continue
      do 30 i= 1, n
      if(family.EQ.1)then
              eta(i)=mu(i)
        else if(family.EQ.2)then
              eta(i)=dlog(mu(i)/(1.0D0-mu(i)))
        else if(family.EQ.3 .OR. family.EQ.4)then
              eta(i)=dlog(DMAX1(1.0D0, mu(i)))
      endif
 30   continue
 
      return
      end
