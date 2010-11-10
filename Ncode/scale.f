      SUBROUTINE SCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common6.h'
      REAL*8  RSAVE(3,NMAX),VSAVE(3,NMAX),BSAVE(NMAX)
      real*8 G, M_sun, R_sun, pc, Km, Kmps
      real*8 mscale, lscale, vscale
*
*       Define physical constants (consistent with IAU & routine UNITS).
      G = 6.6743D-08
      M_sun = 1.9884D+33
      R_sun = 6.955D+10
      pc = 3.0856776D+18
      Km = 1.0D+05
      Kmps = 1.0D+05
*
*       Read virial ratio, rotation scaling factors & tidal radius.
      READ (5,*)  Q, VXROT, VZROT, RTIDE
      RSPH2 = RTIDE
      QVIR = Q
*
      ZMASS = 0.0D0
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Form total mass and centre of mass displacements.
      DO 30 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Adjust coordinates and velocities to c.m. rest frame.
      DO 40 I = 1,N
          DO 35 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   35     CONTINUE
   40 CONTINUE
*
*       Check optional generation of BH binary/single and other BHs.
      IF (KZ(11).GT.0) THEN
*       Read binary elements and component masses (SEMI = 0 for single BH).
          READ (5,*)  SEMI, ECC, BODY1, BODY2
*       Convert from M_sun to pre-scaled N-body units (ZMASS = N up to now).
          BODY(1) = (BODY1 + BODY2)/ZMBAR
          KSTAR(1) = 14
*       Place a stationary heavy binary/single BH at the centre.
          DO 41 K = 1,3
              X(K,1) = 0.0
              XDOT(K,1) = 0.0
   41     CONTINUE
          IF (SEMI.GT.0.0) THEN
*       Enforce initial KS for massive binary.
              NBIN0 = 1
              NBH0 = 2
              KSTAR(2) = 14
              ILAST = KZ(11) + 1
          ELSE
              NBH0 = 1
              ILAST = KZ(11)
          END IF
*      Add one or more optional free-floating BH.
          IF (KZ(11).GT.1) THEN
              DO 42 I = NBH0+1,ILAST
                  READ (5,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
                  BODY(I) = BODY(I)/ZMBAR
                  KSTAR(I) = 14
   42         CONTINUE
              NBH0 = ILAST
          END IF
      END IF
*
*       Include procedure for primordial binaries & singles from fort.10.
      IF ((KZ(22).EQ.4.OR.KZ(22).EQ.-1).AND.(NBIN0.GT.0)) THEN
*       Save the binaries for re-creating two-body motion after scaling.
          DO 45 I = 1,2*NBIN0
              BSAVE(I) = BODY(I)
              DO 44 K = 1,3
                  RSAVE(K,I) = X(K,I)
                  VSAVE(K,I) = XDOT(K,I)
   44         CONTINUE
   45     CONTINUE
*       Form the c.m. of each binary for standard scaling to E = -0.25.
          DO 47 I = 1,NBIN0
              I1 = 2*I - 1
              I2 = 2*I
*       Prevent over-writing BODY(I) when I = 1 and 2*I - 1 = 1.
              B1 = BODY(I1)
              B2 = BODY(I2)
              ZMB  = B1 + B2
              BODY(I) = ZMB
              DO 46 K = 1,3
                  X(K,I) = (B1*X(K,I1) + B2*X(K,I2))/ZMB
                  XDOT(K,I) = (B1*XDOT(K,I1) +
     &                         B2*XDOT(K,I2))/ZMB
   46         CONTINUE
   47     CONTINUE
*
*       Move any single stars up for scaling together with c.m. particles.
          NSING = N - 2*NBIN0
*       Note assumption ZMASS = 1 including any singles (skips DO 50 loop).
          DO 49 I = 1,NSING
             BODY(NBIN0+I) = BODY(2*NBIN0+I)
             DO 48 K = 1,3
                 X(K,NBIN0+I) = X(K,2*NBIN0+I)
                 XDOT(K,NBIN0+I) = XDOT(K,2*NBIN0+I)
   48        CONTINUE
   49     CONTINUE
*       Re-define N & NTOT for total energy calculation (NTOT = N for ZKIN).
          N = NBIN0 + NSING
          NTOT = N
*
          write (6,491)  NBIN0
  491     format (/,12X,'SCALE: ', I5,
     &                  ' Binaries converted to barycentres')
      END IF
*       Save total number of C.Ms.
      NCM = NTOT
*
*       Skip scaling of masses for unscaled upload or planetesimal disk.
      IF (KZ(22).GT.2.OR.KZ(22).EQ.-1.OR.KZ(5).EQ.3) GO TO 52
*
*       Scale masses to standard units of <M> = 1/N and set total mass.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
      ZMASS = 1.0
*
*       Obtain the total kinetic & potential energy.
   52 CALL ENERGY
*
      if (KZ(22).EQ.-1) then
*       Save K.E. and P.E. of C.Ms for TCR & TRH (units of M_sun, pc & km/s).
          ZCM = ZKIN
          PCM = POT
*
*       Evaluate total energy in physical units (M_sun pc2 s-2).
*       G (g-1 cm3 s-2) ==> G (M_sun-1 pc3 s-2)
          G = G/(pc**3/M_sun)
*       Form kinetic energy as KE (M_sun Km2 s-2) ==> KE (M_sun pc2 s-2).
          ZKIN = ZKIN*Km**2/pc**2
*       PE term = (M_sun-1 pc3 s-2)*(M_sun2 pc-1) = (M_sun pc2 s-2).
          EPH = ZKIN - G*POT
          IF (EPH.GT.0.0) EPH = -EPH
*       Set conversion factors (N-body to M_sun, pc, km/s; Heggie & Mathieu).
          mscale = ZMASS
          lscale = G*mscale**2*(-0.25/EPH)
          vscale = sqrt(G*mscale/lscale)*(pc/Km)
*       Specify virial radius and mean mass (pc & M_sun).
          RBAR = lscale
          ZMBAR = ZMASS/(N+NBIN0)
*
          write (6,53)  mscale, lscale, vscale
   53     format (/,12X,'SCALE: Conversion factors: MSCALE = ',F10.3,
     &                  ' M_sun;  LSCALE = ',F6.3,' pc;  VSCALE = ',
     &                    F6.3,' km/sec')
      end if
*
*       Use generalized virial theorem for external tidal field.
      IF (KZ(14).GT.0) THEN
          AZ = 0.0D0
          DO 55 I = 1,N
              AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
   55     CONTINUE
          IF (KZ(14).EQ.1) THEN
*       Use Chandrasekhar eq. (5.535) for virial ratio (rotating frame only).
              VIR = POT - 2.0*(ETIDE + 0.5*TIDAL(4)*AZ)
          ELSE
              VIR = POT - 2.0*ETIDE
          END IF
      ELSE
          VIR = POT
      END IF
*
*       Allow two optional ways of skipping standard velocity scaling.
      IF (KZ(22).EQ.3.OR.KZ(5).EQ.2.OR.KZ(5).EQ.3) THEN
          QV = SQRT(Q*VIR/ZKIN)
          E0 = ZKIN*QV**2 - POT + ETIDE
          SX = 1.0
*       Rescale velocities to new masses for two Plummer spheres.
          IF (KZ(5).EQ.2) THEN
              ZKIN = 0.0
              DO 57 I = 1,N
                  DO 56 K = 1,3
                      XDOT(K,I) = XDOT(K,I)*QV
                      ZKIN = ZKIN + 0.5*BODY(I)*XDOT(K,I)**2
   56             CONTINUE
   57         CONTINUE
              E0 = ZKIN - POT + ETIDE
              Q = ZKIN/POT
              WRITE (6,59)  E0, ZKIN/POT
   59         FORMAT (/,12X,'UNSCALED ENERGY    E =',F10.6,
     &                                       '  Q =',F6.2)
          ELSE
              IF (KZ(5).EQ.3) E0 = ZKIN - POT
              WRITE (6,54)  E0
   54         FORMAT (/,12X,'UNSCALED ENERGY    E =',F10.6)
          END IF
      ELSE IF(KZ(22).NE.-1) THEN
*       Scale non-zero velocities by virial theorem ratio.
          IF (ZKIN.GT.0.0D0) THEN
              QV = SQRT(Q*VIR/ZKIN)
              DO 60 I = 1,N
                  DO 58 K = 1,3
                      XDOT(K,I) = XDOT(K,I)*QV
   58             CONTINUE
   60         CONTINUE
          ELSE
              QV = 1.0
          END IF
*
*       Scale total energy to standard units (E = -0.25 for Q < 1).
          E0 = -0.25
*       Include case of hot system inside reflecting boundary.
          IF (KZ(29).GT.0.AND.Q.GT.1.0) E0 = 0.25
*         ETOT = (Q - 1.0)*POT
          ETOT = ZKIN*QV**2 - POT + ETIDE
*       Note that final ETOT will differ from -0.25 since ETIDE = 0.
          IF (Q.LT.1.0) THEN
              SX = E0/ETOT
          ELSE
              SX = 1.0
          END IF
*
          WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N)
   65     FORMAT (/,12X,'SCALING:    SX =',F6.2,'  E =',1PE10.2,
     &                  '  M(1) =',E9.2,'  M(N) =',E9.2,'  <M> =',E9.2)
*
*       Scale coordinates & velocities to the new units.
          DO 70 I = 1,N
              DO 68 K = 1,3
                  X(K,I) = X(K,I)/SX
                  XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   68         CONTINUE
   70     CONTINUE
      END IF
*
*       Perform second stage of the optional uploading procedure.
      IF ((KZ(22).EQ.4.OR.KZ(22).EQ.-1).AND.(NBIN0.GT.0)) THEN
*       Place any singles last and re-create the original binaries.
          DO 72 I = NSING,1,-1
              L1 = NBIN0
              L2 = 2*NBIN0
              BODY(L2+I) = BODY(L1+I)
              DO 71 K = 1,3
                  X(K,L2+I) = X(K,L1+I)
                  XDOT(K,L2+I) = XDOT(K,L1+I)
   71         CONTINUE
   72     CONTINUE
*       Split each c.m. into binary components using saved quantities.
          DO 74 I = 1,NBIN0
              IP = NBIN0 - I + 1
              BODY(2*IP-1) = BSAVE(2*IP-1)
              BODY(2*IP) = BSAVE(2*IP)
              ZMB = BODY(2*IP-1) + BODY(2*IP)
*       Use the original relative motion for unscaled two-body elements.
              DO 73 K = 1,3
                  XREL = RSAVE(K,2*IP-1) - RSAVE(K,2*IP)
                  VREL = VSAVE(K,2*IP-1) - VSAVE(K,2*IP)
*       Prevent over-writing of first location when IP = 1 and 2*IP-1 = 1.
                  X1 = X(K,IP)
                  V1 = XDOT(K,IP)
                  X(K,2*IP-1) = X(K,IP) + BSAVE(2*IP)*XREL/ZMB
                  X(K,2*IP) = X1 - BSAVE(2*IP-1)*XREL/ZMB
                  XDOT(K,2*IP-1) = XDOT(K,IP) + BSAVE(2*IP)*VREL/ZMB
                  XDOT(K,2*IP) = V1 - BSAVE(2*IP-1)*VREL/ZMB
   73         CONTINUE
   74     CONTINUE
*       Set final values of the particle number (increase from N + NBIN0).
          N = N + NBIN0
          NZERO = N
          NTOT = N
*
          write (6,741)  NBIN0
  741     format (/,12X,'SCALE: ', I5, ' Binaries restored')
      END IF
*
*       Scale masses, coordinates and velocities for input in physical units.
      if (KZ(22).EQ.-1) then
*       Re-calculate the scaled total mass.
          ZMASS = 0.0D0
          DO 743 I = 1,NTOT
             BODY(I) = BODY(I)/mscale
             ZMASS = ZMASS + BODY(I)
             DO 742 K = 1,3
                X(K,I) = X(K,I)/lscale
                XDOT(K,I) = XDOT(K,I)/vscale
  742        continue
  743     continue
*
*       Perform energy check after scaling.
          CALL ENERGY
          write (6,745)  ZKIN - POT
  745     format (/,12X,'SCALE: Scaled total energy: ',F8.5)
*       Set energy of singles & C.Ms only in N-body units for TCR & TRH.
*       escale = mscale*(vscale**2)
          SX = 1.0
          ZKIN = ZCM/(mscale*vscale**2)
          POT = PCM*lscale/mscale**2
          E0 = ZKIN - POT
          write (6,746)  E0
  746     format (/,12X,'SCALE: Scaled C.M. energy: ',F8.5)
      end if
*
*       Introduce optional BH binary (SEMI in scaled N-body units).
      IF (KZ(11).GT.0) THEN
          IF (SEMI.GT.0.0) THEN
              ZMB = BODY(1)
*       Split c.m. into components with given mass ratio.
              BODY(1) = BODY1/(BODY1 + BODY2)*ZMB
              BODY(2) = BODY2/(BODY1 + BODY2)*ZMB
              RAP = SEMI*(1.0 + ECC)
              VAP = SQRT(ZMB*(1.0 - ECC)/(SEMI*(1.0 + ECC)))
*       Specify two-body motion for SEMI & ECC in x-y plane.
              X(1,2) = X(1,1) - BODY(1)*RAP/ZMB
              X(1,1) = X(1,1) + BODY(2)*RAP/ZMB
              XDOT(2,2) = XDOT(2,1) - BODY(1)*VAP/ZMB
              XDOT(2,1) = XDOT(2,1) + BODY(2)*VAP/ZMB
              DO 75 I = 1,2
                  X(2,I) = 0.0
                  X(3,I) = 0.0
                  XDOT(1,I) = 0.0
                  XDOT(3,I) = 0.0
   75         CONTINUE
          END IF
          DO 77 I = 1,ILAST
              WRITE (6,76)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
   76         FORMAT (12X,'BH COMPONENTS    M X XDOT ',
     &                                      1P,E10.2,0P,2(2X,3F7.3))
   77     CONTINUE
      END IF
*
*       Check whether to include rotation (VXROT = 0 in standard case). 
      IF (VXROT.GT.0.0D0) THEN
*
*       Set angular velocity for retrograde motion (i.e. star clusters).
          OMEGA = -SX*SQRT(ZMASS*SX)
          WRITE (6,78)  VXROT, VZROT, OMEGA
   78     FORMAT (/,12X,'VXROT =',F6.2,'  VZROT =',F6.2,
     &                                                 '  OMEGA =',F7.2)
*
*       Add solid-body rotation about Z-axis (reduce random velocities).
          DO 80 I = 1,N
              XDOT(1,I) = XDOT(1,I)*VXROT - X(2,I)*OMEGA
              XDOT(2,I) = XDOT(2,I)*VXROT + X(1,I)*OMEGA
              XDOT(3,I) = XDOT(3,I)*VZROT
   80     CONTINUE
      END IF
*
*       Check option for writing the initial conditions on unit #10.
      IF (KZ(22).EQ.1) THEN
          DO 85 I = 1,N
              WRITE (10,84)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
   84         FORMAT (1P,7E14.6)
   85     CONTINUE
      END IF
*
*       Check for writing initial conditions (physical units) on unit #99.
      IF (KZ(22).EQ.-1.OR.KZ(22).EQ.5) THEN
          DO 86 I = 1,N
              write (99,84)  BODY(I)*mscale, (X(K,I)*lscale,K=1,3),
     &                        (XDOT(K,I)*vscale,K=1,3)
   86     CONTINUE
      END IF
*
*       Check option for reading initial subsystems.
      IF (KZ(24).GT.0) THEN
          K = KZ(24)
          DO 90 I = 1,K
              READ (5,*)  BODY(I), (X(J,I),J=1,3), (XDOT(J,I),J=1,3)
   90     CONTINUE
      END IF
*
*       Set initial crossing time in scaled units.
      TCR = ZMASS**2.5/(2.0D0*ABS(E0))**1.5
      TCR0 = TCR
*
*       Obtain approximate half-mass radius after scaling.
      RSCALE = 0.5*ZMASS**2/(SX*POT)
*       Set square radius of reflecting sphere (used with option 29).
      RSPH2 = (RSPH2*RSCALE)**2
*       Form equilibrium rms velocity (temporarily defined as VC).
      VC = SQRT(2.0D0*ABS(E0)/ZMASS)
*
*       Check for general binary search of initial condition.
      IF (KZ(4).GT.0) THEN
          CALL EVOLVE(0,0)
      END IF
*
*       Print half-mass relaxation time & equilibrium crossing time.
      A1 = FLOAT(NCM)
      TRH = 4.0*TWOPI/3.0*(VC*RSCALE)**3/(15.4*ZMASS**2*LOG(A1)/A1)
      WRITE (6,95)  TRH, TCR, 2.0*RSCALE/VC
   95 FORMAT (/,12X,'TIME SCALES:   TRH =',1PE8.1,'  TCR =',E8.1,
     &                                            '  2<R>/<V> =',E8.1,/)
*
      RETURN
*
      END


