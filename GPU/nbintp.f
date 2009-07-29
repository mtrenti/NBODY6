      SUBROUTINE NBINTP(I,IR)
*
*
*       Parallel irregular integration.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDUM(3),DV(3)
*
*
*       Copy scalars and obtain total force & first derivative.
      DO 5 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
          FIRR(K) = 0.0D0
          FD(K) = 0.0D0
    5 CONTINUE
*
*       Assume small mass at centre for special case of no neighbours.
      NNB0 = LIST(1,I)
      IF (NNB0.EQ.0) THEN
          RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
          FIJ = 0.01*BODYM/(RI2*SQRT(RI2))
          RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                                 XI(3)*XIDOT(3))/RI2
          DO 10 K = 1,3
              FIRR(K) = -FIJ*XI(K)
              FD(K) = -(XIDOT(K) - RDOT*XI(K))*FIJ
   10     CONTINUE
          IF (I.GT.N) IPAIR = I - N
          GO TO 70
      END IF
*
*       Choose force loop for single particle or regularized c.m. body.
      IF (I.LE.N) GO TO 20
*
*       Set KS pair index.
      IPAIR = I - N
      I2 = 2*IPAIR
      I1 = I2 - 1
*
*       Adopt c.m. approximation for small total perturbation.
      IF (LIST(1,I1).GT.0) THEN
*       Obtain irregular force on perturbed c.m. body (including any chain).
          CALL CMFIRR(I,I1,FIRR,FD)
          GO TO 70
      END IF
*
*       Copy c.m. coordinates & velocities for rare unperturbed intruder.
      DO 15 K = 1,3
          X(K,I1) = XI(K)
          X(K,I2) = XI(K)
          XDOT(K,I1) = XIDOT(K)
          XDOT(K,I2) = XIDOT(K)
   15 CONTINUE
*
*       Set neighbour number & list index of the last single particle.
   20 NNB1 = NNB0 + 1
      NNB2 = NNB1
   25 IF (LIST(NNB2,I).LE.N) GO TO 30
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 25
*       Include special case of only c.m. neighbours.
      GO TO 40
*
*       Sum over single particles (unperturbed case included).
*  30 DO 35 L = 2,NNB2
*         K = LIST(L,I)
*         A1 = X(1,K) - XI(1)
*         A2 = X(2,K) - XI(2)
*         A3 = X(3,K) - XI(3)
*         DV(1) = XDOT(1,K) - XIDOT(1)
*         DV(2) = XDOT(2,K) - XIDOT(2)
*         DV(3) = XDOT(3,K) - XIDOT(3)
*         RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*         DR2I = 1.0/RIJ2
*         DR3I = BODY(K)*DR2I*SQRT(DR2I)
*         DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
*         FIRR(1) = FIRR(1) + A1*DR3I
*         FIRR(2) = FIRR(2) + A2*DR3I
*         FIRR(3) = FIRR(3) + A3*DR3I
*         FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
*         FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
*         FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
*  35 CONTINUE
*
   30 CALL CNBINT(I,X,XDOT,BODY,NNB2,LIST(2,I),FIRR,FD)
*
*       See whether any c.m. neighbours should be included.
      IF (NNB2.EQ.NNB1) GO TO 60
*
   40 NNB3 = NNB2 + 1
*       Set index for distinguishing c.m. or resolved components.
      KDUM = 0
*
*       Sum over regularized c.m. neighbours.
      DO 50 L = NNB3,NNB1
          J = LIST(L,I)
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          JP = J - N
          KDUM = 2*JP - 1
          K = KDUM
          IF (LIST(1,KDUM).EQ.0) THEN
              K = J
              GO TO 48
          END IF
*       Sum over individual components of pair #J.
   45     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
   48     DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
*
*       Adopt c.m. approximation outside the effective perturber sphere.
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
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 45
          END IF
   50 CONTINUE
*
*       Include differential force treatment for regularized subsystem.
   60 IF (NCH.GT.0) THEN
*       Distinguish between chain c.m. and any other particle.
          IF (NAME(I).EQ.0) THEN
              CALL CHFIRR(I,0,XI,XIDOT,FIRR,FD)
          ELSE
*       See if chain perturber list contains body #I.
              NP1 = LISTC(1) + 1
              DO 65 L = 2,NP1
                  J = LISTC(L)
                  IF (J.GT.I) GO TO 70
                  IF (J.EQ.I) THEN
                      CALL FCHAIN(I,0,XI,XIDOT,FIRR,FD)
                      GO TO 70
                  END IF
   65         CONTINUE
          END IF
      END IF 
*
*       Check option for external tidal field.
   70 IF (KZ(14).GT.0) THEN
          DO 75 K = 1,3
              FREG(K) = FR(K,I)
   75     CONTINUE
          CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDUM,0)
      END IF
*
*       Include the corrector and set new T0, F, FDOT, D1, D2 & D3.
      DT = TIME - T0(I)
      DTSQ = DT**2
      DT6 = 6.0D0/(DT*DTSQ)
      DT2 = 2.0D0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT
      T0(I) = TIME
*
      DO 80 K = 1,3
	  DF = FI(K,I) - FIRR(K)
	  FID = FIDOT(K,I)
	  SUM = FID + FD(K)
	  AT3 = 2.0D0*DF + DT*SUM
	  BT2 = -3.0D0*DF - DT*(SUM + FID)
*
	  X0(K,I) = XI(K) + (0.6D0*AT3 + BT2)*DTSQ12
	  X0DOT(K,I) = XIDOT(K) + (0.75D0*AT3 + BT2)*DT13
*
*         X0(K,I) = X(K,I)
*         X0DOT(K,I) = XDOT(K,I)
*
	  FI(K,I) = FIRR(K)
	  FIDOT(K,I) = FD(K)
*       Use total force for irregular step (cf. Makino & Aarseth PASJ, 1992).
          FDUM(K) = FIRR(K) + FR(K,I)
*
          D0(K,I) = FIRR(K)
          D1(K,I) = FD(K)
	  D2(K,I) = (3.0D0*AT3 + BT2)*DT2
	  D3(K,I) = AT3*DT6
*       NOTE: These are real derivatives!
   80 CONTINUE
*
*       Specify new time-step by standard criterion (STEPI version not good).
      TTMP = TSTEP(FDUM,FD,D2(1,I),D3(1,I),ETAI)
*
*     IF (I.GT.N) THEN
*       Check for hierarchical configuration but exclude small perturbations.
*         IF (H(IPAIR).LT.-ECLOSE.AND.KZ(36).GT.0) THEN
*             IF (GAMMA(IPAIR).GT.1.0E-04) THEN
*                 CALL KEPLER(I,TTMP)
*             END IF
*         END IF
*     END IF
*
*       Include convergence test for large step (cf. Makino, Ap.J. 369, 200).
      IF ((TTMP.LT.0.1*DTMIN.AND.TTMP.GT.0.01*DTMIN).OR.
     &     TTMP.GT.STEPJ) THEN
          DV2 = 0.0
          F2 = 0.0
          DO 85 K = 1,3
              DV2 = DV2 + (XIDOT(K) - X0DOT(K,I))**2
              F2 = F2 + FIRR(K)**2
   85     CONTINUE
*       Employ Jun's criterion to avoid over-shooting (cf. Book, 2.16).
          DTJ = STEP(I)*(1.0D-06*STEP(I)**2*F2/DV2)**0.1
          TTMP = MIN(TTMP,DTJ)
      END IF
      DT0 = TTMP
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEP(I)) THEN
          IF (DMOD(TIME,2.0*STEP(I)).EQ.0.0D0) THEN 
              TTMP = MIN(2.0*STEP(I),SMAX)
*       Include factor 4 increase for FPOLY initializations with small STEP.
              IF (DT0.GT.10.0*TTMP.AND.
     &            DMOD(TIME,2.0D0*TTMP).EQ.0.0D0) THEN
                  TTMP = MIN(2.0*TTMP,SMAX)
              END IF
          ELSE
              TTMP = STEP(I) 
          END IF
      ELSE IF (TTMP.LT.STEP(I)) THEN
          TTMP = 0.5*STEP(I)
          IF (TTMP.GT.DT0) THEN
              TTMP = 0.5*TTMP
          END IF
      ELSE
          TTMP = STEP(I)
      END IF
*
*       Set new block step and update next time.
      STEP(I) = TTMP
      TNEW(I) = STEP(I) + T0(I)
*
*       See whether total force & derivative needs updating.
      IF (IR.EQ.0) THEN
*       Extrapolate regular force & first derivatives to obtain F & FDOT.
          DTR = TIME - T0R(I)
          DO 90 K = 1,3
              F(K,I) = 0.5D0*(FRDOT(K,I)*DTR + FR(K,I) + FIRR(K))
              FDOT(K,I) = ONE6*(FRDOT(K,I) + FD(K))
   90     CONTINUE
      END IF
*
*       Increase step counter and count perturbed c.m. steps.
*     NSTEPI = NSTEPI + 1
      IF (I.GT.N) THEN
          IF (LIST(1,2*IPAIR-1).GT.0) NSTEPB = NSTEPB + 1
      END IF
*
      RETURN
*
      END
