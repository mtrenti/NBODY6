      SUBROUTINE REGINT2(I,XI,XIDOT)
*
*
*       Regular integration.
*       --------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),DV(3),FD(3),FDR(3)
      REAL*8  FRX(3),FDX(3)
*
*
*       Set neighbour number, time-step & choice of central distance.
      NNB0 = LIST(1,I)
      DTR = TIME - T0R(I)
      IRSKIP = 0
      IF (KZ(39).EQ.0) THEN
          RI2 = (XI(1) - RDENS(1))**2 + (XI(2) - RDENS(2))**2 +
     &                                  (XI(3) - RDENS(3))**2
          RH2 = RSCALE**2
      ELSE
          RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
          RH2 = 9.0*RSCALE**2
      END IF
*       Reduce neighbour radius to prevent future overflow.
      RS(I) = 0.9*RS(I)
*
*       Obtain irregular & regular force and determine current neighbours.
      RS2 = RS(I)**2
      RCRIT2 = 0.0
      VRFAC = 0.0
*       Take volume between inner and outer radius equal to basic sphere.
*   1 RCRIT2 = 1.59*RS2
*       Set radial velocity factor for the outer shell.
*     VRFAC = -0.1*RS2/DTR
*       Start count at 2 and subtract 1 at the end to avoid ILIST(NNB+1).
      NNB = 1
*
*       Initialize scalars for forces & derivatives.
      DO 5 K = 1,3
          FIRR(K) = 0.0D0
          FREG(K) = 0.0D0
          FD(K) = 0.0
          FDR(K) = 0.0
    5 CONTINUE
*
*       Perform fast force loop over single particles.
      DO 10 J = IFIRST,N
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
*       Predicted coordinates avoids spurious force differences.
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
*
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
*
*       First see whether the distance exceeds the outer shell radius.
*         IF (RIJ2.GT.RCRIT2) GO TO 8
*
*       Retain particles with small steps (large derivative corrections).
*         IF (RIJ2.GT.RS2.AND.STEP(J).GT.10.0*SMIN) THEN
*       Accept member if maximum penetration factor exceeds 8 per cent.
*             IF (DRDV.GT.VRFAC) GO TO 8
*         ELSE
*             IF (J.EQ.I) GO TO 10
*         END IF
          IF (RIJ2.GT.RS2) GO TO 8
          IF (J.EQ.I) GO TO 10
*
*       Increase neighbour counter and obtain current irregular force.
          NNB = NNB + 1
          ILIST(NNB) = J
          DR2I = 1.0/RIJ2
          DRDV = 3.0*DRDV*DR2I
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          GO TO 10
*
*       Obtain the regular force.
    8     DR2I = 1.0/RIJ2
          DRDV = 3.0*DRDV*DR2I
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
   10 CONTINUE
*
*       Choose appropriate force loop for single particle or c.m. body.
      IF (I.GT.N) THEN
*       Treat unperturbed KS in the single particle approximation.
          I1 = 2*(I - N) - 1
          IF (LIST(1,I1).GT.0) THEN
*       Obtain remaining force from J > N on c.m. particle.
             CALL CMFREG(I,RS2,NNB,FIRR,FREG,FD,FDR)
             GO TO 20
          END IF
      END IF
*
*       Add any contributions from c.m. particles (singles or unperturbed).
      IF (NPAIRS.GT.0) THEN
          CALL CMFREG2(I,RS2,NNB,XI,XIDOT,FIRR,FREG,FD,FDR)
      END IF
*
*       Include treatment for regularized subsystem.
      IF (NCH.GT.0) THEN
*       Distinguish between chain c.m. and any other particle.
          IF (NAME(I).EQ.0) THEN
              CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
          ELSE
*       Search the chain perturber list for #I.
              NP1 = LISTC(1) + 1
              DO 15 L = 2,NP1
                  J = LISTC(L)
                  IF (J.GT.I) GO TO 20
                  IF (J.EQ.I) THEN
                      CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
                      GO TO 20
                  END IF
   15         CONTINUE
          END IF
      END IF 
*
*       Check optional interstellar clouds.
   20 IF (KZ(13).GT.0) THEN
          CALL FCLOUD(I,FREG,FDR,1)
      END IF
*
*       See whether an external force should be added.
      IF (KZ(14).GT.0) THEN
*       Save current values for deriving work done by tides (#14 = 3).
          IF (KZ(14).EQ.3) THEN
              DO 22 K = 1,3
                  FRX(K) = FREG(K)
                  FDX(K) = FDR(K)
   22         CONTINUE
          END IF
*
*       Obtain the tidal perturbation (force and first derivative).
          CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,1)
*
*       Form rate of tidal energy change during last regular step.
          IF (KZ(14).EQ.3) THEN
              WDOT = 0.0
              W2DOT = 0.0
*             W3DOT = 0.0
              DO 24 K = 1,3
                  PX = FREG(K) - FRX(K)
                  DPX = FDR(K) - FDX(K)
                  WDOT = WDOT + XIDOT(K)*PX
                  W2DOT = W2DOT + (FREG(K) + FIRR(K))*PX + XIDOT(K)*DPX
*                 W3DOT = W3DOT + 2.0*(FREG(K) + FIRR(K))*DPX +
*    &                            (FDR(K) + FD(K))*PX
   24         CONTINUE
*       Note: second-order term derived by Douglas Heggie (Aug/03).
          END IF
      END IF
*
      NNB = NNB - 1
*       Check for zero neighbour number or distant particle with large STEP.
      IF (NNB.EQ.0.OR.(RI2.GT.100.0*RH2.AND.
     &    STEP(I).GT.200.0*DTMIN)) THEN
          IF (NNB.EQ.0) NBVOID = NBVOID + 1
*       Double the neighbour sphere and try again unless RI > 10*RSCALE.
          IF (RI2.GT.9.0*RH2.OR.LIST(1,I).EQ.0) THEN
              IRSKIP = 1
*       Assume small mass at centre for distant body or no neighbours.
              R2 = XI(1)**2 + XI(2)**2 + XI(3)**2
              FIJ = 0.01*BODYM/(R2*SQRT(R2))
              RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                                     XI(3)*XIDOT(3))/R2
              DO 25 K = 1,3
                  FIRR(K) = FIRR(K) - FIJ*XI(K)
                  FD(K) = FD(K) - (XIDOT(K) - RDOT*XI(K))*FIJ
   25         CONTINUE
*       Check maximum membership (note: NNB may be large).
              IF (NNB.GT.NNBMAX) THEN
                  RS2 = 0.9*RS2
*                 GO TO 1
              ELSE
*       Reduce neighbour sphere gradually but allow encounter detection.
                  RS(I) = MAX(0.75*RS(I),0.1*RSCALE)
                  GO TO 50
              END IF
          ELSE
              IRSKIP = 1
*       Assume small mass at centre for distant body or no neighbours.
              R2 = XI(1)**2 + XI(2)**2 + XI(3)**2
              FIJ = 0.01*BODYM/(R2*SQRT(R2))
              RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                                     XI(3)*XIDOT(3))/R2
              DO 28 K = 1,3
                  FIRR(K) = FIRR(K) - FIJ*XI(K)
                  FD(K) = FD(K) - (XIDOT(K) - RDOT*XI(K))*FIJ
   28         CONTINUE
              IF (NNB.EQ.0) THEN
                  RS(I) = MAX(0.75*RS(I),0.1*RSCALE)
                  GO TO 50
              END IF
              RS2 = 1.2*RS2
          END IF
          RS(I) = SQRT(RS2)
          IF (RS(I).GT.10.0*RSCALE.AND.KZ(39).EQ.0) IRSKIP = 1
*         GO TO 1
      END IF
*
*       Restrict neighbour number < NNBMAX to permit one normal addition.
      IF (NNB.LT.NNBMAX) GO TO 40
*
*       Reduce search radius by cube root of 90 % volume factor.
   30 NNB2 = 0.9*NNBMAX
      A1 = FLOAT(NNB2)/FLOAT(NNB)
      IF (RS(I).GT.5.0*RSCALE) THEN
          A1 = MIN(5.0*A1,0.9D0)
          IRSKIP = 1
      END IF
      RS2 = RS2*A1**0.66667
      RCRIT2 = 1.59*RS2
      RS(I) = SQRT(RS2)
      NNB1 = 0
*
      DO 35 L = 1,NNB
          J = ILIST(L+1)
          IF (L + NNB2.GT.NNB + NNB1) GO TO 32
*       Sum of neighbours (NNB1) & those left (NNB+1-L) set to NNB2.
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
*
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
          IF (RIJ2.GT.RCRIT2) GO TO 34
*       Retain outer shell members with negative radial velocity.
          IF (RIJ2.LT.RCRIT2.AND.DRDV.LT.0.0) GO TO 32
          IF (RIJ2.GT.RS2) GO TO 34
*
   32     NNB1 = NNB1 + 1
          JLIST(NNB1+1) = J
          GO TO 35
*
*       Subtract neighbour force included above and add to regular force.
   34     DRDV = 3.0*DRDV*DR2I
          FIRR(1) = FIRR(1) - A1*DR3I
          FIRR(2) = FIRR(2) - A2*DR3I
          FIRR(3) = FIRR(3) - A3*DR3I
          FD(1) = FD(1) - (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) - (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) - (DV(3) - A3*DRDV)*DR3I
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
   35 CONTINUE
*
      DO 38 L = 2,NNB1+1
          ILIST(L) = JLIST(L)
   38 CONTINUE
      NNB = NNB1
      NBFULL = NBFULL + 1
*       See whether to reduce NNB further.
      IF (NNB.GE.NNBMAX) GO TO 30
*
*       Stabilize NNB between ZNBMIN & ZNBMAX by square root of contrast.
   40 A3 = ALPHA*SQRT(FLOAT(NNB)*RS(I))/RS2
      NBP = A3
      A3 = MIN(A3,ZNBMAX) 
*       Reduce predicted membership slowly outside half-mass radius.
      IF (RI2.GT.RH2) THEN
          A3 = A3*RSCALE/SQRT(RI2)
      ELSE
          A3 = MAX(A3,ZNBMIN)
      END IF
      A4 = A3/FLOAT(NNB)
*       Include inertial factor to prevent resonance oscillations in RS.
      IF ((A3 - FLOAT(NNB0))*(A3 - FLOAT(NNB)).LT.0.0) A4 = SQRT(A4)
*
*       Modify volume ratio by radial velocity factor outside the core.
      IF (RI2.GT.RC2.AND.KZ(39).EQ.0.AND.RI2.LT.9.0*RH2) THEN
          RIDOT = (XI(1) - RDENS(1))*XIDOT(1) +
     &            (XI(2) - RDENS(2))*XIDOT(2) +
     &            (XI(3) - RDENS(3))*XIDOT(3)
          A4 = A4*(1.0 + RIDOT*DTR/RI2)
      END IF
*
*       See whether neighbour radius of c.m. body should be increased.
      IF (I.GT.N.AND.NNB.LT.ZNBMAX) THEN
*       Set perturber range (soft binaries & H > 0 have short duration).
          A2 = 100.0*BODY(I)/ABS(H(I-N))
*       Stabilize NNB on ZNBMAX if too few perturbers.
          IF (A2.GT.RS(I)) THEN
              A3 = MAX(1.0 - FLOAT(NNB)/ZNBMAX,0.0D0)
*       Modify volume ratio by approximate square root factor.
              A4 = 1.0 + 0.5*A3
          END IF
      END IF
*
*       Limit change of RS for small steps (RSFAC = MAX(25/TCR,1.5*VC/RSMIN).
      A5 = MIN(RSFAC*DTR,0.25D0)
*       Restrict volume ratio by inertial factor of 25 per cent either way.
      A4 = A4 - 1.0
      IF (A4.GT.A5) THEN
          A4 = A5
      ELSE
          A4 = MAX(A4,-A5)
      END IF
*
*       Modify neighbour sphere radius by volume factor.
      IF (IRSKIP.EQ.0) THEN
          A3 = ONE3*A4
          A1 = 1.0 + A3 - A3*A3
*       Second-order cube root expansion (maximum volume error < 0.3 %).
          IF (RS(I).GT.5.0*RSCALE) A1 = SQRT(A1)
*       Prevent reduction of small NNB if predicted value exceeds ZNBMIN.
          IF (NNB.LT.ZNBMIN.AND.NBP.GT.ZNBMIN) THEN
              IF (A1.LT.1.0) A1 = 1.05
          END IF
*       Skip modification for small changes (avoids oscillations in RS).
          IF (ABS(A1 - 1.0D0).GT.0.003) THEN
              RS(I) = A1*RS(I)
          END IF
          IF (NNB.LT.0.5*NBP) RS(I) = 1.05*RS(I)
      END IF
*
*       Calculate the radial velocity with respect to at most 3 neighbours.
      IF (NNB.LE.3.AND.RI2.LT.9.0*RH2) THEN
          A1 = 2.0*RS(I)
*
          DO 45 L = 1,NNB
              J = ILIST(L+1)
              RIJ = SQRT((XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                                         (XI(3) - X(3,J))**2)
              RSDOT = ((XI(1) - X(1,J))*(XIDOT(1) - XDOT(1,J)) +
     &                 (XI(2) - X(2,J))*(XIDOT(2) - XDOT(2,J)) +
     &                 (XI(3) - X(3,J))*(XIDOT(3) - XDOT(3,J)))/RIJ
*       Find smallest neighbour distance assuming constant regular step.
              A1 = MIN(A1,RIJ + RSDOT*DTR)
   45     CONTINUE
*
*       Increase neighbour sphere if all members are leaving inner region.
          RS(I) = MAX(A1,1.1*RS(I))
      END IF
*
*       Check minimum neighbour sphere since last output (skip NNB = 0).
      IF (LIST(1,I).GT.0) RSMIN = MIN(RSMIN,RS(I))
*
*       Check optional procedures for adding neighbours.
      IF ((KZ(37).EQ.1.AND.LISTV(1).GT.0).OR.KZ(37).GT.1) THEN
          CALL CHECKL(I,NNB,XI,XIDOT,RS2,FIRR,FREG,FD,FDR)
      END IF
*
*       Find loss or gain of neighbours at the same time.
   50 NBLOSS = 0
      NBGAIN = 0
*
*       Accumulate tidal energy change for general galactic potential.
      IF (KZ(14).EQ.3) THEN
*       Note: Taylor series at end of interval with negative argument.
*         ETIDE = ETIDE - BODY(I)*((ONE6*W3DOT*DTR - 0.5*W2DOT)*DTR +
*    &                                                   WDOT)*DTR
          ETIDE = ETIDE + BODY(I)*(0.5*W2DOT*DTR - WDOT)*DTR
*       Note: integral of Taylor series for V*P using final values.
      END IF
*
*       Check case of zero old membership (NBGAIN = NNB specifies net gain).
      IF (NNB0.EQ.0) THEN
          NBGAIN = NNB
          DO 55 L = 1,NNB
              JLIST(NNB0+L) = ILIST(L+1)
   55     CONTINUE
          GO TO 70
      END IF
*
      JMIN = 0
      L = 2
      LG = 2
*       Set termination value in ILIST(NNB+2) and save last list member.
      ILIST(NNB+2) = NTOT + 1
      ILIST(1) = LIST(NNB0+1,I)
*
*       Compare old and new list members in locations L & LG.
   56 IF (LIST(L,I).EQ.ILIST(LG)) GO TO 58
*
*       Now check whether inequality means gain or loss.
      IF (LIST(L,I).GE.ILIST(LG)) THEN
          NBGAIN = NBGAIN + 1
          JLIST(NNB0+NBGAIN) = ILIST(LG)
*       Number of neighbour losses can at most be NNB0.
          L = L - 1
*       The same location will be searched again after increasing L below.
      ELSE
          NBLOSS = NBLOSS + 1
          J = LIST(L,I)
          JLIST(NBLOSS) = J
*       Check SMIN step indicator (rare case permits fast skip below).
          IF (STEP(J).LT.SMIN) JMIN = J
          LG = LG - 1
      END IF
*
*       See whether the last location has been checked.
   58 IF (L.LE.NNB0) THEN
          L = L + 1
          LG = LG + 1
*       Last value of second search index is NNB + 2 which holds NTOT + 1.
          GO TO 56
      ELSE IF (LG.LE.NNB) THEN
          LG = LG + 1
          LIST(L,I) = NTOT + 1
*       Last location of list holds termination value (saved in ILIST(1)).
          GO TO 56
      END IF
*
*       See whether any old neighbour with small step should be retained.
      IF (JMIN.EQ.0) GO TO 70
*
      K = 1
   60 IF (NNB.GT.NNBMAX.OR.I.GT.N) GO TO 70
      J = JLIST(K)
*       A single regularized component will be replaced by the c.m.
      IF (STEP(J).GT.SMIN.OR.J.LT.IFIRST.OR.J.GT.N) GO TO 68
*       Retain old neighbour inside 2*RS to avoid large correction terms.
      RIJ2 = (XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                             (XI(3) - X(3,J))**2
      IF (RIJ2.GT.4.0*RS2) GO TO 68
*
      L = NNB + 1
   62 IF (ILIST(L).LT.J) GO TO 64
      ILIST(L+1) = ILIST(L)
      L = L - 1
      IF (L.GT.1) GO TO 62
*
*       Save index of body #J and update NNB & NBLOSS.
   64 ILIST(L+1) = J
      NNB = NNB + 1
      NBLOSS = NBLOSS - 1
*       Restore last old neighbour in case NBLOSS = 0 at end of search.
      LIST(NNB0+1,I) = ILIST(1)
      NBSMIN = NBSMIN + 1
*
*       Perform correction to irregular and regular force components.
      A1 = X(1,J) - XI(1)
      A2 = X(2,J) - XI(2)
      A3 = X(3,J) - XI(3)
      DV(1) = XDOT(1,J) - XIDOT(1)
      DV(2) = XDOT(2,J) - XIDOT(2)
      DV(3) = XDOT(3,J) - XIDOT(3)
*
      RIJ2 = A1*A1 + A2*A2 + A3*A3
      DR2I = 1.0/RIJ2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
      FIRR(1) = FIRR(1) + A1*DR3I
      FIRR(2) = FIRR(2) + A2*DR3I
      FIRR(3) = FIRR(3) + A3*DR3I
      FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
      FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
      FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
      FREG(1) = FREG(1) - A1*DR3I
      FREG(2) = FREG(2) - A2*DR3I
      FREG(3) = FREG(3) - A3*DR3I
      FDR(1) = FDR(1) - (DV(1) - A1*DRDV)*DR3I
      FDR(2) = FDR(2) - (DV(2) - A2*DRDV)*DR3I
      FDR(3) = FDR(3) - (DV(3) - A3*DRDV)*DR3I
*
*       Remove body #J from JLIST unless it is the last or only member.
      IF (K.GT.NBLOSS) GO TO 70
      DO 66 L = K,NBLOSS
          JLIST(L) = JLIST(L+1)
   66 CONTINUE
*       Index of last body to be moved up is L = NBLOSS + 1.
      K = K - 1
*       Check the same location again since a new body has appeared.
   68 K = K + 1
*       Last member to be considered is in location K = NBLOSS.
      IF (K.LE.NBLOSS) GO TO 60
*
*       Form time-step factors and update regular force time.
   70 DTSQ = DTR**2
      DT6 = 6.0D0/(DTR*DTSQ)
      DT2 = 2.0D0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DTR13 = ONE3*DTR
      T0R(I) = TIME
*
*       Suppress the corrector for large time-step ratios (experimental).
      IF (STEP(I).LT.5.0D0*DTMIN.AND.DTR.GT.50.0*STEP(I)) THEN
          DTR13 = 0.0
          DTSQ12 = 0.0
      END IF
*
*       Include the corrector and save F, FI, FR & polynomial derivatives.
      DO 75 K = 1,3
*       Subtract change of neighbour force to form actual first derivative.
          DFR = FR(K,I) - (FIRR(K) - FI(K,I)) - FREG(K)
          FDR0 = FDR(K) - (FIDOT(K,I) - FD(K))
*
          FRD = FRDOT(K,I)
	  SUM = FRD + FDR0
	  AT3 = 2.0D0*DFR + DTR*SUM
	  BT2 = -3.0D0*DFR - DTR*(SUM + FRD)
*
	  X0(K,I) = X0(K,I) + (0.6D0*AT3 + BT2)*DTSQ12
	  X0DOT(K,I) = X0DOT(K,I) + (0.75D0*AT3 + BT2)*DTR13
*
*         X0(K,I) = X(K,I)
*         X0DOT(K,I) = XDOT(K,I)
*
          FI(K,I) = FIRR(K)
	  FR(K,I) = FREG(K)
          F(K,I) = 0.5D0*(FREG(K) + FIRR(K))
          FIDOT(K,I) = FD(K)
	  FRDOT(K,I) = FDR(K)
          FDOT(K,I) = ONE6*(FDR(K) + FD(K))
*
          D0(K,I) = FIRR(K)
          D0R(K,I) = FREG(K)
*         D0R(K,I) = FREG(K) - (FI(K,I) - FIRR(K))
*         D1R(K,I) = FDR0
*       Use actual first derivatives (2nd derivs only consistent in FPCORR).
          D1(K,I) = FD(K)
          D1R(K,I) = FDR(K)
*       Set second & third derivatives based on old neighbours (cf. FPCORR).
	  D2R(K,I) = (3.0D0*AT3 + BT2)*DT2
	  D3R(K,I) = AT3*DT6
   75 CONTINUE
*
*       Correct force polynomials due to neighbour changes.
      IF (KZ(38).GT.0) THEN
          IF (KZ(38).EQ.1) THEN
              CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT)
          ELSE IF (KZ(38).EQ.2.AND.STEPR(I).GT.0.04) THEN
              CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT)
          END IF
      END IF
*
*       Copy new neighbour list if different from old list.
      IF (NBLOSS + NBGAIN.GT.0) THEN
          LIST(1,I) = NNB
          DO 80 L = 2,NNB+1
              LIST(L,I) = ILIST(L)
   80     CONTINUE
          NBFLUX = NBFLUX + NBLOSS + NBGAIN
      END IF
*
*       Obtain new regular integration step using composite expression.
*       STEPR = (ETAR*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
*
*       Determine new regular step.
      TTMP = TSTEP(FREG,FDR,D2R(1,I),D3R(1,I),ETAR)
*
*       Impose a smooth step reduction inside compact core.
      IF (NC.LT.50.AND.RI2.LT.RC2) THEN
          TTMP = TTMP*MIN(1.0D0,0.5D0*(1.0D0 + RI2*RC2IN))
      END IF
      DT0 = TTMP
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEPR(I)) THEN
          IF (DMOD(TIME,2.0*STEPR(I)).EQ.0.0D0) THEN 
              TTMP = MIN(2.0*STEPR(I),SMAX)
          ELSE
              TTMP = STEPR(I) 
          END IF
      ELSE IF (TTMP.LT.STEPR(I)) THEN
          TTMP = 0.5*STEPR(I)
*       Allow a second reduction to prevent spurious contributions.
          IF (TTMP.GT.DT0) THEN
              TTMP = 0.5*TTMP
          END IF
      ELSE
          TTMP = STEPR(I)
      END IF
*
*       Set new regular step and reduce irregular step if STEPR < STEP.
      STEPR(I) = TTMP
  110 IF (TTMP.LT.STEP(I)) THEN
          STEP(I) = 0.5D0*STEP(I)
          TNEW(I) = T0(I) + STEP(I)
          NICONV = NICONV + 1
          GO TO 110
      END IF
      NSTEPR = NSTEPR + 1
*
      RETURN
*
      END
