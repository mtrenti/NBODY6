      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
*
*
*       Make a formal call to define input parameters & counters.
      CALL DEFINE
*
*       Read & print the main input parameters.
      READ (5,*)  N, NFIX, NCRIT, NRAND, NNBMAX, NRUN
      READ (5,*)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT, QE, RBAR, ZMBAR
      READ (5,*)  (KZ(J),J=1,40)
      READ (5,*)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX
*
      WRITE (6,10)
   10 FORMAT (/////,15X,'N  NFIX  NCRIT  NRAND  NNBMAX  NRUN')
      WRITE (6,12)  N, NFIX, NCRIT, NRAND, NNBMAX, NRUN
   12 FORMAT (/,I16,I6,2I7,I8,I6)
      WRITE (6,15)
   15 FORMAT (//,12X,'ETAI      ETAR      RS0       DTADJ     DELTAT',
     &                        '    TCRIT     QE        RBAR      ZMBAR')
      WRITE (6,20)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT, QE, RBAR,
     &                                                             ZMBAR
   20 FORMAT (/,9X,1P,10E10.1)
      WRITE (6,22)
   22 FORMAT (//,12X,'OPTIONS')
      WRITE (6,24)  (J,J=1,40)
   24 FORMAT (/,9X,40I3)
      WRITE (6,26)  (KZ(J),J=1,40)
   26 FORMAT (/,9X,40I3)
      WRITE (6,28)
   28 FORMAT (//,12X,'DTMIN     RMIN      ETAU      ECLOSE    GMIN',
     &                                                     '      GMAX')
      WRITE (6,30)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX
   30 FORMAT (/,9X,1P,6E10.1)
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
      GPRINT(1) = 0.0
      DELTAS = 0.0
      IF (KZ(4).GT.0) THEN
*       Read parameters for binary evolution analysis.
          K = KZ(4)
          READ (5,*)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
          WRITE (6,40)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
   40     FORMAT (//,12X,'DELTAS =',F6.2,'  ORBITS(1) =',F6.2,
     &                                            '  GPRINT(J) =',9F7.3)
*       Modify binary output factor by perturbation at different levels.
          DO 50 L = 2,K
              ORBITS(L) = ORBITS(1)*(GPRINT(1)/GPRINT(L))**0.3333
   50     CONTINUE
      END IF
*
*       Set random number skip for routine DATA.
      IDUM1 = NRAND
*
*       Save square of c.m. approximation parameter (coupled to GMIN).
      CMSEP2 = GMIN**(-0.666666667)
*
*       Define total particle number & neighbour membership parameters.
      NTOT = N
      NZERO = N
      NBZERO = NNBMAX
      ZNBMIN = 0.2*FLOAT(NNBMAX)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
*       Save initial ETAI.
      ETA0 = ETAI
      RSMIN = RS0
      RC = RS0
*       Temporary save of initial neighbour sphere for OUTPUT & NBLIST.
*
      RETURN
*
      END
