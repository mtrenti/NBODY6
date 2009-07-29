      SUBROUTINE CMFREG2(I,RS2,NNB,XI,XID,FIRR,FREG,FD,FDR)
*
*
*       Regular & irregular force from J > N (singles or unpert KS).
*       ------------------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XID(3),FIRR(3),FREG(3),DV(3),FD(3),FDR(3),
     &        FIRR2(3),FD2(3),FREG2(3),FDR2(3)
      INTEGER ISAVE(KMAX),JSAVE(KMAX)
*
*
*       Initialize counters (all neighbours, unpert & pert KS and regular).
      LJ = 0
      NI = 1
      NP = 0
      NNB1 = 1
*
*       Split all c.m. particles into regular and irregular parts.
      DO 5 J = N+1,NTOT
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          RIJ2 = A1**2 + A2**2 + A3**2
          IF (RIJ2.LT.RS2) THEN
              IF (J.EQ.I) GO TO 5
              LJ = LJ + 1
              JSAVE(LJ) = J
*       Distinguish between unperturbed and perturbed case.
              J1 = 2*(J - N) - 1
              IF (LIST(1,J1).EQ.0) THEN
                  NI = NI + 1
                  ISAVE(NI) = J
              ELSE
                  NP = NP + 1
                  JPERT(NP) = J
              END IF
          ELSE
*       Note that NNB1 count starts in JLIST(2) and NP in JPERT(1).
              NNB1 = NNB1 + 1
              JLIST(NNB1) = J
          END IF
    5 CONTINUE
*
*       Obtain irregular and regular force assuming c.m. approximation.
      CALL CNBINT(I,X,XDOT,BODY,NI,ISAVE(2),FIRR2,FD2)
      CALL CNBINT(I,X,XDOT,BODY,NNB1,JLIST(2),FREG2,FDR2)
*
*       Add any new contributions (zero otherwise).
      DO 10 K = 1,3
          FIRR(K) = FIRR(K) + FIRR2(K)
          FD(K) = FD(K) + FD2(K)
          FREG(K) = FREG(K) + FREG2(K)
          FDR(K) = FDR(K) + FDR2(K)
   10 CONTINUE
*
*       Perform differential force corrections for close c.m. neighbours.
      DO 20 L = 1,NP
          J = JPERT(L)
          J1 = 2*(J - N) - 1
*       Add the contributions from KS components of binary.
          K = J1
   15     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XID(1)
          DV(2) = XDOT(2,K) - XID(2)
          DV(3) = XDOT(3,K) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
*
          IF (K.EQ.J1) THEN
              K = K + 1
              GO TO 15
          END IF
   20 CONTINUE
*
*       Copy any c.m. neighbours to final list and set membership.
      DO 40 L = 1,LJ
          ILIST(NNB+L) = JSAVE(L)
   40 CONTINUE
      NNB = NNB + LJ
*
      RETURN
*
      END
