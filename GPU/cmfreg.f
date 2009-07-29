      SUBROUTINE CMFREG(I,RS2,NNB,FIRR,FREG,FD,FDR)
*
*
*       Regular & irregular force from J > N for perturbed KS.
*       ------------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(6),FPD(6),
     &        FREG(3),FDR(3),FREG2(3),FDR2(3)
      INTEGER JSAVE(KMAX)
*
*
*       Set individual KS components.
      I2 = 2*(I - N)
      I1 = I2 - 1
*
*       Save last single neighbour location and initialize c.m. counters.
      LS = NNB
      LJ = 0
      NNB2 = 0
      NRF = 1
*
*       Add any c.m. neighbours inside RS and form list of regular members.
      DO 10 J = N+1,NTOT
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          IF (RIJ2.LT.RS2) THEN
              IF (J.EQ.I) GO TO 10
              LJ = LJ + 1
              JSAVE(LJ) = J
*       Note that NNB count starts in ILIST(2) and NNB2 in JPERT(1).
              J1 = 2*(J - N) - 1
              IF (LIST(1,J1).EQ.0) THEN
                  NNB = NNB + 1
                  ILIST(NNB) = J
              ELSE
                  NNB2 = NNB2 + 1
                  JPERT(NNB2) = J
              END IF
          ELSE
              NRF = NRF + 1
              JLIST(NRF) = J
          END IF
   10 CONTINUE
*
*       Obtain c.m. contributions to FIRR & FREG with C++ routine.
      CALL CNBINT(I1,X,XDOT,BODY,NNB,ILIST(2),FP,FPD)
      CALL CNBINT(I2,X,XDOT,BODY,NNB,ILIST(2),FP(4),FPD(4))
      CALL CNBINT(I,X,XDOT,BODY,NRF,JLIST(2),FREG2,FDR2)
*
*       Add regular force contributions from any other J > N.
      DO 20 K = 1,3
          FREG(K) = FREG(K) + FREG2(K)
          FDR(K) = FDR(K) + FDR2(K)
   20 CONTINUE
*
*       Perform irregular force loop from nearby perturbed KS pairs.
      DO 40 LL = 1,NNB2
          JP = JPERT(LL)
*       Evaluate perturbation on first component due to body #K.
          KDUM = 2*(JP - N) - 1
          K = KDUM
   34     dr2 = 0.0
          drdv = 0.0
          DO 35 L = 1,3
              dx(L) = X(L,K) - X(L,I1)
              dv(L) = XDOT(L,K) - XDOT(L,I1)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   35     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 36 L = 1,3
              FP(L) = FP(L) + dx(L)*dr3i
              FPD(L) = FPD(L) + (dv(L) - dx(L)*drdv)*dr3i
   36     CONTINUE
*
*       Evaluate perturbation on second component due to body #K.
          dr2 = 0.0
          drdv = 0.0
          DO 37 L = 1,3
              dx(L) = X(L,K) - X(L,I2)
              dv(L) = XDOT(L,K) - XDOT(L,I2)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   37     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 38 L = 1,3
              FP(L+3) = FP(L+3) + dx(L)*dr3i
              FPD(L+3) = FPD(L+3) + (dv(L) - dx(L)*drdv)*dr3i
   38     CONTINUE
*
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 34
          END IF
   40 CONTINUE
*
*       Form mass-weighted c.m. force & derivative.
      BODYIN = 1.0/BODY(I)
      DO 50 K = 1,3
          FIRR(K) = (BODY(I1)*FP(K) + BODY(I2)*FP(K+3))*BODYIN
          FD(K) = (BODY(I1)*FPD(K) + BODY(I2)*FPD(K+3))*BODYIN
   50 CONTINUE
*
*       Copy any c.m. neighbours to current list and set membership.
      DO 60 L = 1,LJ
          ILIST(LS+L) = JSAVE(L)
   60 CONTINUE
      NNB = LS + LJ
*
*       Check force correction due to regularized chain.
      IF (NCH.GT.0) THEN
          CALL KCPERT(I,I1,FIRR,FD)
      END IF 
*
      RETURN
*
      END
