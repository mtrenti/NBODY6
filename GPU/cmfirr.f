      SUBROUTINE CMFIRR(I,I1,FIRR,FD)
*
*
*       Irregular force on perturbed binary from other c.m.
*       ---------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(6),FPD(6)
      INTEGER  JSAVE(LMAX),IPERT(LMAX)
*
*
*       Set second KS component and counters for current pair.
      I2 = I1 + 1
      NNB1 = LIST(1,I) + 1
      NNB2 = 0
      NNB = 1
*
*       Copy neighbours into single star approximation and nearby KS pairs.
      DO 10 L = 2,NNB1
*       Adopt convention of ILIST members starting in #2.
          J = LIST(L,I)
          IF (J.LE.N) THEN
              NNB = L
              JSAVE(L) = J
          ELSE
*       Use c.m. approximation to avoid closer unperturbed massive binary.
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2
     &                                   +  (X(3,I) - X(3,J))**2
              IF (CMSEP2*R(J-N)**2.LT.RIJ2) THEN
                  NNB = NNB + 1
                  JSAVE(NNB) = J
              ELSE
*       Note IPERT array should be thread-safe (not common) for nbintp.f.
                  NNB2 = NNB2 + 1
                  IPERT(NNB2) = J
              END IF
          END IF
   10 CONTINUE
*
*       Obtain c.m. contributions to FIRR for each component with C++.
      CALL CNBINT(I1,X,XDOT,BODY,NNB,JSAVE(2),FP,FPD)
      CALL CNBINT(I2,X,XDOT,BODY,NNB,JSAVE(2),FP(4),FPD(4))
*
*       Form irregular force components due to nearby perturbed KS pairs.
      DO 40 LL = 1,NNB2
          J = IPERT(LL)
          KDUM = 2*(J - N) - 1
          K = KDUM
*       Adopt c.m. of unperturbed binary for consistency.
          IF (LIST(1,KDUM).EQ.0) THEN
              K = J
          END IF
*       Evaluate perturbation on first component due to body #K.
   20     dr2 = 0.0
          drdv = 0.0
          DO 25 L = 1,3
              dx(L) = X(L,K) - X(L,I1)
              dv(L) = XDOT(L,K) - XDOT(L,I1)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   25     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 30 L = 1,3
              FP(L) = FP(L) + dx(L)*dr3i
              FPD(L) = FPD(L) + (dv(L) - dx(L)*drdv)*dr3i
   30     CONTINUE
*
*       Evaluate perturbation on second component due to body #K.
          dr2 = 0.0
          drdv = 0.0
          DO 32 L = 1,3
              dx(L) = X(L,K) - X(L,I2)
              dv(L) = XDOT(L,K) - XDOT(L,I2)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   32     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 35 L = 1,3
              FP(L+3) = FP(L+3) + dx(L)*dr3i
              FPD(L+3) = FPD(L+3) + (dv(L) - dx(L)*drdv)*dr3i
   35     CONTINUE
*
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 20
          END IF
   40 CONTINUE
*
*       Form mass-weighted perturbations.
      BODYIN = 1.0/BODY(I)
      DO 50 K = 1,3
          FIRR(K) = (BODY(I1)*FP(K) + BODY(I2)*FP(K+3))*BODYIN
          FD(K) = (BODY(I1)*FPD(K) + BODY(I2)*FPD(K+3))*BODYIN
   50 CONTINUE
*
*       Check force correction due to regularized chain.
      IF (NCH.GT.0) THEN
          CALL KCPERT(I,I1,FIRR,FD)
      END IF 
*
      RETURN
*
      END
