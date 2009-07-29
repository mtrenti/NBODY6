      SUBROUTINE MODIFY(KSTART)
*
*
*       Parameter modification at restart.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
*
*
*       Read first, second or both lines (KSTART = 3, 4, 5).
      IF (KSTART.EQ.4) GO TO 10
*
*       Read new DTADJ, DELTAT, TADJ, TNEXT, TCRIT, QE & KZ(J) (if > 0).
      READ (5,*)  DTA, DT, TA, TN, TC, QE1, J, K
*
*       Copy old parameters if corresponding input is zero.
      IF (DTA.LE.0.0) THEN
          DTA = DTADJ
      END IF
*
      IF (DT.LE.0.0) THEN
          DT = DELTAT
      END IF
*
      IF (TA.LE.0.0) THEN
          TADJ = MAX(TADJ - DTADJ + DTA,TIME)
      ELSE
          TADJ = MAX(TA-TOFF,TIME)
      END IF
*
      IF (TN.LE.0.0) THEN
          TNEXT = MAX(TNEXT - DELTAT + DT,TIME)
      ELSE
          TNEXT = MAX(TN-TOFF,TIME)
      END IF
*
      DTADJ = DTA
      DELTAT = DT
      IF (TC.GT.0.0) TCRIT = TC
      IF (QE1.GT.0.0) QE = QE1
*
*       See whether any options should be changed.
      IF (J.GT.0) KZ(J) = K
*
      WRITE (6,5)  DTADJ, DELTAT, TCRIT, QE, J, K
    5 FORMAT (///,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,'  DELTAT =',
     &                            F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1,
     &                                            '  KZ(',I2,') =',I2,/)
*
*       Read new ETAI, ETAR, ETAU, DTMIN, RMIN, NCRIT (if > 0 & KSTART >= 4).
   10 IF (KSTART.GE.4) THEN
          READ (5,*)  ETA1, ETA2, ETA3, DTM, RM, NEWCR
*
*       Check modification of integration parameters.
          IF (ETA1.GT.0.0) ETAI = ETA1
          IF (ETA2.GT.0.0) ETAR = ETA2
          IF (ETA3.GT.0.0) ETAU = ETA3
          IF (DTM.GT.0.0) THEN
              DTMIN = DTM
              SMIN = 2.0*DTM
          END IF
          IF (RM.GT.0.0) THEN
              RMIN = RM
              RMIN2 = RM**2
              RMIN22 = 4.0*RMIN2
          END IF
          IF (NEWCR.GT.0) NCRIT = NEWCR
*
          WRITE (6,15)  ETAI, ETAR, ETAU, DTMIN, RMIN, NCRIT
   15     FORMAT (/,7X,'RESTART PARAMETERS:   ETAI =',F6.3,'  ETAR =',
     &                        F6.3,'  ETAU =',F6.3,'  DTMIN =',1P,E8.1,
     &                             '  RMIN =',E8.1,'  NCRIT =',0P,I5,/)
      END IF
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Save the new parameters on tape/disc unit #1 just in case.
      CALL MYDUMP(1,1)
*
      RETURN
*
      END
