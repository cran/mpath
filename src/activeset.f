C     update the active set with only non-zero coefficients
C     input: m, beta, eps
C     output: activeset, jk
      subroutine find_activeset(m, beta, eps, activeset, jk)

      implicit none
      integer j, jk, m, activeset(m)
      double precision beta(m), eps


      do 90 j=1, m
         activeset(j)=0
90    continue
      
C     jk: number of variables in active set
      jk = 0
C     it is possible, jk=0 if beta=0, like intercept-only
C     model
      do 60 j=1, m
      if(dabs(beta(j)) .GT. eps)then
C     used jk not j, to be indexed in loop_gaussian subroutine
         jk=jk+1
         activeset(jk)=j
      endif
60    continue
      
      return
      end
