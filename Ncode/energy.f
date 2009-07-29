      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common6.h'
*
*
*       Sum the total energy of regularized pairs.
      EBIN = 0.0D0
      DO 10 IPAIR = 1,NPAIRS
*       Skip pairs with zero mass of c.m. particle (merged binary ghost).
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
*       Predict coordinates, velocities & binding energy.
              CALL RESOLV(IPAIR,1)
              EBIN = EBIN + BODY(2*IPAIR-1)*BODY(2*IPAIR)*HT/
     &                                                     BODY(N+IPAIR)
          END IF
   10 CONTINUE
*
*       Calculate the potential energy.
      ZKIN = 0.D00
      POT = 0.0
      VIR = 0.0
      I = 1
   20 JMIN = I + 1
      IF (I.LE.2*NPAIRS) THEN
*       Binding energy of regularized pairs is included explicitly above.
          IPAIR = KVEC(I)
          JMIN = 2*IPAIR + 1
      END IF
      POTJ = 0.D00
*
      DO 30 J = JMIN,N
          IF (BODY(J).EQ.0.0D0) GO TO 30
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          POTJ = POTJ + BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3)
   30 CONTINUE
*
      POT = POT + BODY(I)*POTJ
      I = I + 1
      IF (I.LT.N) GO TO 20
*
*       Sum the kinetic energy (include c.m. bodies but not components).
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
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
              DO 50 I = 1,N
                  RI2 = AP2
                  DO 45 K = 1,3
                      RI2 = RI2 + X(K,I)**2
   45             CONTINUE
                  ETIDE = ETIDE - BODY(I)*MP/SQRT(RI2)
   50         CONTINUE
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
