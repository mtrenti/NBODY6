      SUBROUTINE ENERGY2
*
*
*       Total energy from GPU potential.
*       --------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  GPUPHI(NMAX)
*
*
*       Sum the total energy of regularized pairs.
      EBIN = 0.0D0
      DO 1 IPAIR = 1,NPAIRS
*       Skip pairs with zero mass of c.m. particle (merged binary ghost).
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
*       Predict coordinates, velocities & binding energy.
              CALL RESOLV(IPAIR,1)
              EBIN = EBIN + BODY(2*IPAIR-1)*BODY(2*IPAIR)*HT/
     &                                                     BODY(N+IPAIR)
          END IF
    1 CONTINUE
*
*       Evaluate individual potentials on GPU (including all c.m.).
      VIR = 0.0
      POT = 0.0
      NN = NTOT - IFIRST + 1
      CALL GPUPOT(NN,BODY(IFIRST),X(1,IFIRST),GPUPHI)
*
*       Move the table entries down to give room for any KS components.
      I2 = 2*NPAIRS
      IF (NPAIRS.GT.0) THEN
          DO 5 I = NTOT,I2+1,-1
              GPUPHI(I) = GPUPHI(I-I2)
    5     CONTINUE
*
*       Copy c.m. potential to the components.
          DO 10 IPAIR = 1,NPAIRS
              I1 = 2*IPAIR - 1
              GPUPHI(I1) = GPUPHI(N+IPAIR)
              GPUPHI(I1+1) = GPUPHI(N+IPAIR)
   10     CONTINUE
      END IF
*
*       Sum individual contributions after differential correction.
      DO 20 I = IFIRST,NTOT
          CALL PHICOR(I,DPHI1,DPHI2)
          IF (I.LE.N) THEN
              GPUPHI(I) = GPUPHI(I) + DPHI1
              POT = POT + BODY(I)*GPUPHI(I)
          ELSE
              I1 = 2*(I - N) - 1
              GPUPHI(I1) = GPUPHI(I1) + DPHI1
              GPUPHI(I1+1) = GPUPHI(I1+1) + DPHI2
              POT = POT + BODY(I1)*GPUPHI(I1) + BODY(I1+1)*GPUPHI(I1+1)
          END IF
   20 CONTINUE
*
*       Take half the value because of double counting.
      POT = 0.5*POT
*
*     POTS = POT
*     CALL ENERGY
*     WRITE (6,40)  POTS, (POT - POTS)/POT
*  40 FORMAT (' POTENTIAL CHECK   POTS DP/P  ',F10.6,1P,E10.2)
*     POT = POTS
*       Sum the kinetic energy (include c.m. bodies but not components).
      ZKIN = 0.D00
      DO 50 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   50 CONTINUE
      ZKIN = 0.5D0*ZKIN
*
*       Obtain the tidal potential energy for linearized external field. 
      IF (KZ(14).EQ.0) THEN
*       Note: ETIDE holds accumulated tidal energy if KZ(14) = 3.
          ETIDE = 0.0D0
      ELSE
*       Employ general expression sum {m*r*F} for virial energy.
          CALL XTRNLV(1,N)
*       Form tidal energy with Plummer potential (note ETIDE use for #14=3).
          IF (KZ(14).EQ.4) THEN
              ETIDE = 0.0
              DO 60 I = 1,N
                  RI2 = AP2
                  DO 55 K = 1,3
                      RI2 = RI2 + X(K,I)**2
   55             CONTINUE
                  ETIDE = ETIDE - BODY(I)*MP/SQRT(RI2)
   60         CONTINUE
          END IF
      END IF
*
*       Check differential potential energy due to chain subsystem.
      IF (NCH.GT.0) THEN
          CALL CHPOT(DP)
          POT = POT + DP
      END IF
*
*       Total energy = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL.
*
      RETURN
*
      END
