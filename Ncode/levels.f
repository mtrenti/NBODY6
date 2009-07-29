      SUBROUTINE LEVELS
*
*
*       Histograms of block time-step levels.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      INTEGER  IHIST(32),IHISTR(32)
*
*
*       Initialize histogram counters.
      DO 10 J = 1,32
          IHIST(J) = 0
          IHISTR(J) = 0
   10 CONTINUE
*
*       Loop over all single particles & c.m.
      JMAX = 0
      JMAXR = 0
      FAC = 1.0/LOG(1.9999999)
      DO 20 I = IFIRST,NTOT
          IF (BODY(I).EQ.0.0D0) GO TO 20
          J = 1 - LOG(STEP(I))*FAC
          IHIST(J) = IHIST(J) + 1
          JMAX = MAX(J,JMAX)
          J = 1 - LOG(STEPR(I))*FAC
          IHISTR(J) = IHISTR(J) + 1
          JMAXR = MAX(J,JMAXR)
   20 CONTINUE
*
*       Print histograms of block-steps (STEPR with KZ(33) > 1).
      JMAX = MIN(JMAX,27)
      WRITE (6,30)  (IHIST(J),J=1,JMAX)
   30 FORMAT (' #6  STEP   ',12I6,5I5,10I4)
*
      IF (KZ(33).GT.1) THEN
          JMAXR = MIN(JMAXR,27)
          WRITE (6,35)  (IHISTR(J),J=1,JMAXR)
   35     FORMAT (' #7  STEPR  ',10I6,5I5,10I4)
      END IF
*
      RETURN
*
      END
