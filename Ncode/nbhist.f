      SUBROUTINE NBHIST
*
*
*       Neighbour list histogram.
*       -------------------------
*
      INCLUDE 'common6.h'
      INTEGER  IHIST(32)
*
*
*       Initialize histogram counters.
      DO 10 J = 1,32
          IHIST(J) = 0
   10 CONTINUE
*
*       Loop over all single particles & c.m.
      JMAX = 0
      FAC = 1.0/LOG(1.9999999)
      DO 20 I = IFIRST,NTOT
          NNB = LIST(1,I)
          ZZ = NNB + 1
          J = 1 + LOG(ZZ)*FAC
          IHIST(J) = IHIST(J) + 1
          JMAX = MAX(J,JMAX)
   20 CONTINUE
*
*       Print histogram of neighbour list membership.
      JMAX = MIN(JMAX,12)
      WRITE (6,30)  (IHIST(J),J=1,JMAX)
   30 FORMAT (' #8  NBLIST    ',12I6)
*
      RETURN
*
      END
