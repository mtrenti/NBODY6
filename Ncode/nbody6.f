*
*             N B O D Y 6
*             ***********
*
*       Regularized N-body code with synthetic stellar evolution.
*       ---------------------------------------------------------
*
*       Hermite integration method with block-steps & AC scheme.
*       --------------------------------------------------------
*
*       Roche-lobe mass transfer & Mardling 2007 stability criterion.
*       -------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge (V 7.3 07/09).
*       ..........................................................
*
*       Binary and chain regularization: Seppo Mikkola, Turku.
*       ......................................................
*
*       Stellar evolution: Chris Tout & Jarrod Hurley (Swinburn).
*       .........................................................
*
*       GPU support: Keigo Nitadori, Riken (V 8.0.0 4/09).
*       ..................................................
*
      PROGRAM NBODY6
*
      INCLUDE 'common6.h'
      EXTERNAL MERGE
*
*
*       Initialize the timer.
      CALL CPUTIM(CPU0)
*
*       Read start/restart indicator & CPU time.
      READ (5,*)  KSTART, TCOMP
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          CALL ADJUST
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0
*       Set IPHASE = -1 for new time-step list in routine INTGRT.
          IPHASE = -1
*
*       Initialize evolution parameters which depend on metallicity.
          IF (KZ(19).GE.3) THEN
              CALL ZCNSTS(ZMET,ZPARS)
          END IF
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
      END IF
*
*       Advance solutions until next output or change of procedure.
    1 CALL INTGRT
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
          CALL KSREG
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
          CALL KSTERM
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          CALL ADJUST
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
          ISUB = 0
          CALL TRIPLE(ISUB)
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
          ISUB = 0
          CALL QUAD(ISUB)
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
          CALL MERGE
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
          CALL RESET
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
          ISUB = 0
          CALL CHAIN(ISUB)
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
