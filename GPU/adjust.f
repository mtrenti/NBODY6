      SUBROUTINE ADJUST
*
*
*       Parameter adjustment and energy check.
*       --------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/ECHAIN/  ECH
      SAVE  DTOFF
      DATA  DTOFF /100.0D0/
*
*
*       Predict X & XDOT for all particles (except unperturbed pairs).
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time using GPU for potentials.
      CALL ENERGY2
*
*       Initialize c.m. terms.
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Obtain c.m. & angular momentum integrals and Z-moment of inertia.
      AZ = 0.0D0
      ZM = 0.0D0
      ZMASS = 0.0D0
      ISUN = 1
      DO 20 I = 1,N
          IF (NAME(I).EQ.1) ISUN = I
          ZMASS = ZMASS + BODY(I)
          DO 15 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   15     CONTINUE
*         RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
*    &                                   (X(3,I) - RDENS(3))**2
          AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
          ZM = ZM + BODY(I)*(X(1,I)**2 + X(2,I)**2)
   20 CONTINUE
*
*       Form c.m. coordinates & velocities (vectors & scalars).
      DO 25 K = 1,3
          CMR(K) = CMR(K)/ZMASS
          CMRDOT(K) = CMRDOT(K)/ZMASS
   25 CONTINUE
*
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
      CMRDOT(4) = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)
*
*       Subtract the kinetic energy of c.m. due to possible cloud effects.
      IF (KZ(13).GT.0) ZKIN = ZKIN - 0.5*ZMASS*CMRDOT(4)**2
*
*       Include non-zero virial energy for Plummer potential and/or 3D case.
      VIR = POT - VIR
*       Note angular momentum term is contained in virial energy (#14=1/2).
      Q = ZKIN/VIR
      E(3) = ZKIN - POT
*       Modify single particle energy by tidal energy (except pure 3D).
      IF (KZ(14).NE.3) THEN
          E(3) = E(3) + ETIDE
      END IF
*
*       Modify angular momentum integral using Chandrasekhar eq. (5.530).
      IF (KZ(14).EQ.1.OR.KZ(14).EQ.2) THEN
          AZ = AZ + 0.5*TIDAL(4)*ZM
      END IF
*
*       Define crossing time using single particle energy (cf. option 14).
      TCR = ZMASS**2.5/(2.0*ABS(E(3)))**1.5
*       Note: see below for better definition in Plummer or 3D orbit case.
      IF (Q.GT.1.0.AND.KZ(14).LT.3) THEN
          TCR = TCR*SQRT(2.0*Q)
      END IF
*       Form provisional total energy.
      ETOT = ZKIN - POT + ETIDE
*
*       Include KS pairs, triple, quad, mergers, collisions & chain.
      ETOT = ETOT + EBIN + ESUB + EMERGE + ECOLL + EMDOT + ECDOT
*       Add any contributions from chain regularization.
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Update energies and form the relative error (divide by ZKIN or E(3)).
      IF (TIME.LE.0.0D0) THEN
          DE = 0.0D0
          BE(1) = ETOT
          BE(3) = ETOT
          DELTA1 = 0.0D0
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = BE(3) - BE(2)
          DELTA1 = DE
          DETOT = DETOT + DE
          DE = DE/MAX(ZKIN,ABS(E(3)))
*       Save sum of relative energy error for main output and accumulate DE.
          ERROR = ERROR + DE
          ERRTOT = ERRTOT + DE
      END IF
*
*       Set provisional half-mass radius.
      RSCALE = 0.5*ZMASS**2/POT
*
*       Determine average neighbour number and smallest neighbour sphere.
      NNB = 0
      RS0 = RSCALE
      DO 30 I = IFIRST,NTOT
          NNB = NNB + LIST(1,I)
          IF (LIST(1,I).GT.0) RS0 = MIN(RS0,RS(I))
   30 CONTINUE
      NNB = NNB/(N - NPAIRS)
*
*       Use current value if minimum neighbour sphere not implemented.
      IF (RSMIN.EQ.0.0D0) RSMIN = RS0
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N-NPAIRS.GE.20.AND.KZ(29).EQ.0.AND.KZ(5).NE.3) THEN
          CALL CORE
      ELSE
          NC = N
          ZMC = ZMASS
          RHOD = 1.0
          RHOM = 1.0
          RC = RSCALE
          RC2 = RC**2
          RC2IN = 1.0/RC2
      END IF
*
*       Take the Sun as reference for plotting planetesimal disk members.
      IF (KZ(5).EQ.3) THEN
          DO 35 K = 1,3
              RDENS(K) = X(K,ISUN)
   35     CONTINUE
*
*       Determine the eccentricity dispersion and total energy of disk.
          DISP2 = 0.0
          EDISK = 0.0
          DO 40 I = 1,N
              IF (NAME(I).EQ.1.OR.NAME(I).EQ.NZERO) GO TO 40
              RI2 = (X(1,I) - X(1,ISUN))**2 + (X(2,I) - X(2,ISUN))**2
              VI2 = (XDOT(1,I) - XDOT(1,ISUN))**2 +
     &              (XDOT(2,I) - XDOT(2,ISUN))**2
              RRDOT = (X(1,I) - X(1,ISUN))*(XDOT(1,I) - XDOT(1,ISUN)) +
     &                (X(2,I) - X(2,ISUN))*(XDOT(2,I) - XDOT(2,ISUN))
              RI = SQRT(RI2)
              SEMI = 2.0/RI - VI2/(BODY(ISUN) + BODY(I))
              SEMI = 1.0/SEMI
              ECC2 = (1.0 - RI/SEMI)**2 +
     &                RRDOT**2/(SEMI*(BODY(I) + BODY(ISUN)))
              DISP2 = DISP2 + ECC2
              EDISK = EDISK - 0.5*BODY(I)/SEMI
   40     CONTINUE
          DISP = SQRT(DISP2/FLOAT(N-2))
          WRITE (35,41)  TTOT, DISP, EDISK
   41     FORMAT (' ',F8.1,1P,E10.2,E12.4)
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Scale average & maximum core density by the mean value.
      RHOD = 4.0*TWOPI*RHOD*RSCALE**3/(3.0*ZMASS)
      RHOM = 4.0*TWOPI*RHOM*RSCALE**3/(3.0*ZMASS)
*
*       Adopt density contrasts of unity for hot system.
      IF (KZ(29).GT.0.AND.ZKIN.GT.POT) THEN
          RHOD = 1.0
          RHOM = 1.0
      END IF
*
*       Check optional determination of regularization parameters.
      IF (KZ(16).GT.0) THEN
          RMIN0 = RMIN
*
*       Form close encounter distance from scale factor & density contrast.
          RMIN = 4.0*RSCALE/(FLOAT(N)*RHOD**0.3333)
*       Include alternative expression based on core radius (experimental).
          IF (KZ(16).GT.1.AND.NC.LT.0.01*N) THEN
              RMIN = 0.05*RC/FLOAT(NC)**0.3333
          END IF
*       Use harmonic mean to reduce fluctuations (avoid initial value).
          IF (TIME.GT.0.0D0) RMIN = SQRT(RMIN0*RMIN)
*       Impose maximum value for sufficient perturbers.
          RMIN = MIN(RMIN,RSMIN*GMIN**0.3333)
*       Define scaled DTMIN by RMIN & <M> and include ETAI for consistency.
          DTMIN = 0.04*SQRT(ETAI/0.02D0)*SQRT(RMIN**3/BODYM)
*       Specify binding energy per unit mass of hard binary (impose Q = 0.5).
          ECLOSE = 4.0*MAX(ZKIN,ABS(ZKIN - POT))/ZMASS
*       Adopt central velocity as upper limit (avoids large kick velocities).
          IF (2.0*ZKIN/ZMASS.GT.VC**2) ECLOSE = 2.0*VC**2
          IF (Q.GT.0.5) THEN
              ECLOSE = 0.5*ECLOSE/Q
          END IF
          ECLOSE = MIN(ECLOSE,1.0D0)
*       Initialize flag for increasing RMIN (used with #16 > 1).
          KSMAG = 0
      END IF
*
*       Check optional modification of DTMIN, ECLOSE & TCR for hot system.
      IF (KZ(29).GT.0.AND.Q.GT.1.0) THEN
          DTMIN = 0.04*SQRT(ETAI/0.02D0)*SQRT(RMIN**3/BODYM)
          SIGMA2 = 2.0*ZKIN/ZMASS
          VP2 = 4.0*BODYM/RMIN
          DTMIN = DTMIN*SQRT((VP2 + SIGMA2/Q)/(VP2 + 2.0D0*SIGMA2))
          ECLOSE = SIGMA2
          TCR = 2.0*RSCALE/SQRT(SIGMA2)
      END IF
*
*       Set useful scalars for the integrator.
      SMIN = 2.0*DTMIN
      RMIN2 = RMIN**2
      RMIN22 = 4.0*RMIN2
      EBH = -0.5*BODYM*ECLOSE
      IF (TIME.LE.0.0D0) THEN
          STEPJ = 0.01*(60000.0/FLOAT(N))**0.3333
          READ (5,*)  SMAX
          IF (DMOD(SMAX,DTK(10)).NE.0.0D0) THEN
              WRITE (6,42)  SMAX, SMAX/DTK(10)
   42         FORMAT (' FATAL ERROR!    SMAX SMAX/DTK(10) ',1P,2E10.2)
              STOP
          END IF
          WRITE (6,43)  SMAX
   43     FORMAT (/,12X,'MAXIMUM TIME-STEP ',F8.4)
*       Ensure no steps exceed maximum (large step could stay unchanged).
          DO 44 I = IFIRST,NTOT
              IF (STEPR(I).GT.SMAX) THEN
                  STEPR(I) = SMAX
                  STEP(I) = MIN(STEP(I),SMAX)
                  TNEW(I) = STEP(I)
              END IF
   44     CONTINUE
      END IF
*       Adopt 2*RSMIN for neighbour sphere volume factor in routine REGINT.
      RSFAC = MAX(25.0/TCR,3.0D0*VC/(2.0D0*RSMIN))
*
*       Update density contrast factor for neighbour sphere modification.
      IF (TIME.LE.0.0D0.OR.KZ(40).EQ.0) THEN
          ALPHA = FLOAT(NNBMAX)*SQRT(0.08D0*RSCALE**3/FLOAT(N-NPAIRS))
      END IF
*       Include optional stabilization to increase neighbour number.
      IF (KZ(40).EQ.1.AND.FLOAT(NNB).LT.0.25*NNBMAX) THEN
          FAC = 1.0 + (0.25*NNBMAX - NNB)/(0.25*FLOAT(NNBMAX))
          ALPHA = FAC*ALPHA
      END IF
*
*       Define tidal radius for isolated system (2*RTIDE used in ESCAPE).
      IF (KZ(14).EQ.0) RTIDE = 10.0*RSCALE
*       Redefine the crossing time for 3D cluster orbit or Plummer model.
      IF ((KZ(14).EQ.3.OR.KZ(14).EQ.4).AND.ZKIN.GT.0.0) THEN
          TCR = 2.0*RSCALE/SQRT(2.0*ZKIN/ZMASS)
      END IF
*
*       Print energy diagnostics & KS parameters.
      ICR = TTOT/TCR
*       Obtain elapsed wall clock time (hours, minutes & seconds).
      CALL WTIME(IHOUR,IMIN,ISEC)
      WRITE (6,45)  TTOT, Q, DE, BE(3),RMIN, DTMIN, ECLOSE, ICR, DELTA1,
     &              IHOUR, IMIN, ISEC
   45 FORMAT (/,' ADJUST:  TIME =',F8.2,'  Q =',F5.2,'  DE =',1P,E10.2,
     &          '  E =',0P,F10.6,'  RMIN =',1P,E8.1,'  DTMIN =',E8.1,
     &          '  ECLOSE =',0P,F5.2,'  TC =',I5,'  DELTA =',1P,E9.1,
     &          '  WTIME =',0P,I4,2I3)
      CALL FLUSH(6)
*
*       Perform automatic error control (RETURN on restart with KZ(2) > 1).
      CALL CHECK(DE)
      IF (ABS(DE).GT.5.0*QE) GO TO 70
*
*       Check for escaper removal.
      IF (KZ(23).GT.0) THEN
          CALL ESCAPE
      END IF
*
*       Include optional search of unstable triples.
      IF (KZ(18).GE.1) THEN
          CALL CHECK3
      END IF
*
*       Check correction for c.m. displacements.
      IF (KZ(31).GT.0) THEN
          CALL CMCORR
      END IF
*
*       Include diagnostics for massive binary (bound or unbound initially).
      IF (KZ(5).EQ.4) THEN
          IP = 0
          DO 50 IPAIR = 1,NPAIRS
              IF (NAME(2*IPAIR-1).LE.2.OR.NAME(2*IPAIR).LE.2) THEN
                  IP = IPAIR
              END IF
   50     CONTINUE
          IF (IP.GT.0) THEN
              I1 = 2*IP - 1
              I2 = I1 + 1
              SEMI = -0.5*BODY(N+IP)/H(IP)
              ECC2 = (1.0 - R(IP)/SEMI)**2 +
     &                                  TDOT2(IP)**2/(SEMI*BODY(N+IP))
              EB = BODY(I1)*BODY(I2)*H(IP)/BODY(N+IP)
              WRITE (35,52)  TTOT, SEMI, EB, E(3), SQRT(ECC2),
     &                       NAME(I1), NAME(I2)
   52         FORMAT (' ',F8.1,1P,3E12.4,0P,F7.3,2I5)
*  52         FORMAT (' T A E EB ECL NAME ',F8.1,1P,3E12.4,0P,F7.3,2I5)
              CALL FLUSH(35)
          END IF
      END IF
*
*       Allow a gradual decrease of NNBMAX due to escaper removal.
      IF (KZ(40).EQ.3) THEN
          NNBMAX = NBZERO*SQRT(FLOAT(N)/FLOAT(NZERO))
*       Note revised definition of ZNBMAX & ZNBMIN to reduce overflows.
          ZNBMAX = 0.8*NNBMAX
          ZNBMIN = 0.1*NNBMAX
      END IF
*
*       Include optional fine-tuning of neighbour number (#40 >= 2).
      IF (KZ(40).GE.2) THEN
          IF (NNB.LT.NNBMAX/5.0) THEN
              ALPHA = 1.1*ALPHA
          ELSE
              ALPHA = 0.9*ALPHA
          END IF
      END IF
*
*       See whether standard output is due.
      IOUT = 0
      IF (TIME.GE.TNEXT) THEN
          CALL OUTPUT
*       Check optional overflow diagnostics (#33 > 1: current & accumulated).
          IF (KZ(33).GT.1) THEN
              WRITE (6,55)  NOFL(1), NOFL(2), ALPHA, NNB, NNBMAX
   55         FORMAT (' #9  OVERFLOWS  ',I5,I9,'   ALPHA =',F6.2,
     &                                 '  <NNB> =',I4,'  NNBMAX =',I4)
              IOUT = 1
          END IF
      END IF
*
*       Include optional diagnostics for the hardest binary below ECLOSE.
      IF (KZ(33).GE.2.AND.IOUT.GT.0) THEN
          HP = 0.0
          IP = 0
          DO 60 IPAIR = 1,NPAIRS
*       Skip outer (ghost) binary of quadruple system.
              IF (H(IPAIR).LT.HP.AND.BODY(N+IPAIR).GT.0.0D0) THEN
                  HP = H(IPAIR)
                  IP = IPAIR
              END IF
   60     CONTINUE
          IF (IP.GT.0.AND.HP.LT.-ECLOSE) THEN
              I1 = 2*IP - 1
              I2 = I1 + 1
              SEMI = -0.5*BODY(N+IP)/H(IP)
              PB = DAYS*SEMI*SQRT(SEMI/BODY(N+IP))
              ECC2 = (1.0 - R(IP)/SEMI)**2 +
     &                                  TDOT2(IP)**2/(SEMI*BODY(N+IP))
              EB = BODY(I1)*BODY(I2)*H(IP)/BODY(N+IP)
              WRITE (39,62)  TTOT, NAME(I1), NAME(I2), KSTAR(N+IP),
     &                       LIST(1,I1), SQRT(ECC2), SEMI, PB, EB, E(3)
   62         FORMAT (' BINARY:   T NAME K* NP E A P EB E3 ',
     &                            F8.1,2I6,2I4,F7.3,1P,2E10.2,0P,2F9.4)
              CALL FLUSH(39)
          END IF
      END IF
*
*       Include optional KS reg of binaries outside standard criterion.
      IF (KZ(8).GT.0.AND.N.GE.5000) THEN
          DTCL = 30.0*DTMIN
          RCL = 10.0*RMIN
          CALL SWEEP(DTCL,RCL)
      END IF
*
*       Update time for next adjustment.
      TADJ = TADJ + DTADJ
*       Re-initialize marginal stability counter to avoid including old case.
      NMTRY = 0
*
*       Check optional truncation of time.
      IF (KZ(35).GT.0.AND.TIME.GE.DTOFF) THEN
          CALL OFFSET(DTOFF)
      END IF
*
*       Obtain elapsed CPU time and update total since last output/restart.
      CALL CPUTIM(TCOMP)
      CPUTOT = CPUTOT + TCOMP - CPU0
      CPU0 = TCOMP
*
*       Save COMMON after energy check (skip TRIPLE, QUAD, CHAIN).
      TDUMP = TIME
      IF (KZ(2).GE.1.AND.NSUB.EQ.0) CALL MYDUMP(1,2)
*       Check COMMON save on fort.1 at main output (#1 = 2).
      IF (KZ(1).EQ.2.AND.NSUB.EQ.0) THEN
          IF (IOUT.GT.0) CALL MYDUMP(1,1)
      END IF
*
*       Check termination criteria (TIME > TCRIT & N <= NCRIT).
      IF (TTOT.GE.TCRIT.OR.N.LE.NCRIT) THEN
*       Terminate after optional COMMON save.
          WRITE (6,65)  TTOT, CPUTOT/60.0, ERRTOT, DETOT
   65     FORMAT (//,9X,'END RUN',3X,'TIME =',F8.1,'  CPUTOT =',F7.1,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
          IF (KZ(1).GT.0.AND.NSUB.EQ.0) CALL MYDUMP(1,1)
          STOP
      END IF
*
   70 RETURN
*
      END
