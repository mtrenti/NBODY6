      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3)
      INTEGER  NXTLST(NMAX),IBL(LMAX),NBLIST(NMAX),LISTQ(NMAX),NL(20)
      LOGICAL LOOP,LSTEPM
      SAVE IQ,ICALL,NQ,LQ,LOOP,LSTEPM,STEPM,ISAVE,JSAVE
      DATA IQ,ICALL,LQ,LOOP,LSTEPM,STEPM /0,2,11,.TRUE.,.FALSE.,0.03125/
*
*
*       Enforce level search on return, except new and terminated KS.
      IF (IPHASE.NE.1.AND.IPHASE.NE.2) LOOP = .TRUE.
*
*       Update quantized value of STEPM for large N (first time only).
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
          K = (FLOAT(NZERO)/1024.0)**0.333333
          STEPM = 0.03125D0/2**(K-1)
          LSTEPM = .TRUE.
      END IF
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
      DTM = 1.0
      TPREV = TIME
*       Initialize end-point of integration times and set DTM.
      DO 1000 I = IFIRST,NTOT
          TNEW(I) = T0(I) + STEP(I)
          DTM = MIN(DTM,STEP(I))
 1000 CONTINUE
*
*       Determine level for the smallest step (ignore extreme values).
      LQS = 20
      DO 1001 L = 6,20
          IF (DTM.EQ.DTK(L)) THEN
              LQS = L
          END IF
 1001 CONTINUE
*
*       Specify upper level for optimized membership.
      LQB = LQS - 4
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
*       Enforce new block step search initially and on significant change.
      TLISTQ = TIME
*
*       Check updating new list of block steps with T0 + STEP =< TLISTQ.
    1 ICALL = ICALL + 1
*       Reset TMIN second & third time after change to catch new chain step.
      IF (TIME.GE.TLISTQ.OR.ICALL.LE.3) THEN
*       Update interval by optimization at major times (sqrt of N-NPAIRS).
          IF (DMOD(TLISTQ,2.0D0).EQ.0.0D0.OR.LOOP) THEN
              LOOP = .FALSE.
              DO 10 L = 1,20
                  NL(L) = 0
   10         CONTINUE
              DO 14 I = IFIRST,NTOT
*       Count steps at five different levels for the smallest values.
                  DO 12 L = LQB,LQS
                      IF (STEP(I).LT.DTK(L)) NL(L) = NL(L) + 1
   12             CONTINUE
   14         CONTINUE
              NLSUM = 0
*       Determine interval by summing smallest steps until near sqrt(N-N_b).
              NSQ = SQRT(FLOAT(N - NPAIRS))
              LQ = LQS
              DO 15 L = LQS,LQB,-1
                  NLSUM = NLSUM + NL(L)
                  IF (NLSUM.LE.NSQ) LQ = L
   15         CONTINUE
*             WRITE (6,16)  TIME+TOFF,NQ,NLSUM,LQ,(NL(K),K=LQB,LQS)
*  16         FORMAT (' LEVEL CHECK:    T NQ NLSUM LQ NL  ',
*    &                                  F9.3,3I5,2X,7I4)
          END IF
*
*       Increase interval by optimized value.
          NQ = 0
          TMIN = 1.0D+10
   18     TLISTQ = TLISTQ + DTK(LQ)
          DO 20 I = IFIRST,NTOT
              IF (TNEW(I).LE.TLISTQ) THEN
                  NQ = NQ + 1
                  LISTQ(NQ) = I
                  TMIN = MIN(TNEW(I),TMIN)
              END IF
   20     CONTINUE
*       Increase interval in rare case of zero membership.
          IF (NQ.EQ.0) GO TO 18
*       Make a slight adjustment for high levels and small membership.
          IF (LQ.LE.15.AND.NQ.LE.2) LQ = LQ - 1
      END IF
*
*       Find all particles in next block (TNEW = TMIN).
      CALL INEXT(NQ,LISTQ,TMIN,NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TIME = T0(I) + STEP(I)
      TBLOCK = TIME
      LI = 0
      IPRED = 0
*
*     ID = 0
*     IF (NSTEPI.EQ.1008936) ID = 1
*     IF (NSTEPI.GE.1008764) THEN
*     WRITE (6,22) NXTLEN, NSTEPI, I, NAME(I), STEP(I), STEPR(I)
* 22  FORMAT (' NEXT   LEN # I NM DT DTR  ',I5,I9,2I6,1P,2E10.2)
*     CALL FLUSH(6)
*     END IF
*       Re-determine list if current time exceeds boundary.
      IF (TIME.GT.TLISTQ) GO TO 1
*
*       Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
          END IF
      END IF
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
*         WRITE (6,25)  I, NAME(I), NSTEPI, TIME, STEP(I), TIME/STEP(I)
*  25     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I6,I11,F12.5,1P,E9.1,0P,F16.4)
*         STOP
*     END IF
*
*       Check for new regularization at end of block.
      IF (IKS.GT.0) THEN
          IPHASE = 1
*       Copy the saved component indices and time.
          ICOMP = ISAVE
          JCOMP = JSAVE
          TIME = TSAVE
          GO TO 100
      END IF
*
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance ARchain or KS at first new time.
      IF (TIME.GT.TPREV) THEN
          CALL SUBINT(IQ,I10)
          IF (IQ.LT.0) GO TO 999
      END IF
*
*       Check regular force condition for small block memberships.
      IR = 0
      IF (NXTLEN.LE.10) THEN
          DO 28 L = 1,NXTLEN
              J = NXTLST(L)
              IF (TNEW(J).GE.T0R(J) + STEPR(J)) THEN
                  IR = IR + 1
              END IF
   28     CONTINUE
      END IF
*
*       Decide between merging of neighbour lists or full N prediction.
      IF (NXTLEN.LE.10.AND.IR.EQ.0) THEN
*
*       Initialize pointers for neighbour lists.
          DO 30 L = 1,NXTLEN
              IBL(L) = NXTLST(L)
   30     CONTINUE
*
*       Merge all neighbour lists (with absent members of IBL added).
          CALL NBSORT(NXTLEN,IBL,NNB,NBLIST)
*
*       Predict coordinates & velocities of neighbours and #I to order FDOT.
          NBPRED = NBPRED + NNB
          DO 35 L = 1,NNB
              J = NBLIST(L)
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   35     CONTINUE
*       Ensure prediction of chain c.m. not on block-step (may be needed).
          IF (NCH.GT.0) THEN
              IF (TNEW(ICH).GT.TBLOCK) THEN
                  CALL XVPRED(ICH,0)
              END IF
          END IF
      ELSE
          IPRED = 1
          NNPRED = NNPRED + 1
          DO 40 J = IFIRST,NTOT
              IF (BODY(J).EQ.0.0D0) GO TO 40
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   40     CONTINUE
      END IF
*
*       Resolve perturbed KS pairs with c.m. prediction after NBSORT.
      JJ = -1
      DO 45 JPAIR = 1,NPAIRS
      JJ = JJ + 2
      IF (LIST(1,JJ).GT.0) THEN
*       Ignore c.m. prediction after full N loop (all active KS needed).
          IF (IPRED.EQ.0) THEN
              J = N + JPAIR
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
          END IF
          ZZ = 1.0
*       Distinguish between low and high-order prediction of U & UDOT.
          IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
          CALL KSRES2(JPAIR,J1,J2,ZZ)
      END IF
   45 CONTINUE
*
*       Save new time (output time at TIME > TADJ) and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
*
*       Advance the pointer (<= NXTLEN) and select next particle index.
   50 LI = LI + 1
      IF (LI.GT.NXTLEN) GO TO 1
      I = NXTLST(LI)
      TIME = T0(I) + STEP(I)
*
*       See whether the regular force needs to be updated (IR > 0).
      IF (T0R(I) + STEPR(I).LE.TIME) THEN
          IR = 1
      ELSE
          IR = 0
      END IF
*
*       Advance the irregular step.
      IKS0 = IKS
      CALL NBINT(I,IKS,IR,XI,XIDOT)
*
*       Save indices of first KS candidates in the block and TIME.
      IF (IKS0.EQ.0.AND.IKS.GT.0) THEN
          ISAVE = ICOMP
          JSAVE = JCOMP
          TSAVE = TIME
      END IF
*
*       See whether the regular step is due.
      IF (IR.GT.0) THEN
          CALL REGINT(I,XI,XIDOT)
      END IF
*
*       Determine next block time (note STEP may shrink in REGINT).
      TMIN = MIN(TNEW(I),TMIN)
* 
*       Copy current coordinates & velocities from corrected values.
      IF (LI.EQ.NXTLEN) THEN
          DO 60 L = 1,NXTLEN
              I = NXTLST(L)
              DO 55 K = 1,3
                  X(K,I) = X0(K,I)
                  XDOT(K,I) = X0DOT(K,I)
   55         CONTINUE
   60     CONTINUE
*
*       Check integration of tidal tail members.
          IF (NTAIL.GT.0) THEN
*       Allow large quantized interval with internal iteration.
              IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
                  DO 65 J = ITAIL0,NTTOT
                      IF (TNEW(J).LE.TIME) THEN
                          CALL NTINT(J)
                      END IF
   65             CONTINUE
              END IF
          END IF
      ELSE
*       Continue until last member has been done (improves reproducibility).
          GO TO 50
      END IF
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.NE.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
      IF (KZ(37).GT.0.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time at end of block-step.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 100-year step (new polynomials).
          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              IF (KZ(19).GE.3) THEN
                  CALL MDOT
              ELSE
                  CALL MLOSS
              END IF
              IF (IPHASE.LT.0) GO TO 999
          END IF
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.NMAX) GO TO 50

      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).GT.1) CALL MYDUMP(1,1)
      END IF
*
*       Check option for general binary search.
      IF (KZ(4).GT.0.AND.TIME - TLASTS.GT.DELTAS) THEN  
          CALL EVOLVE(0,0)
      END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUB.EQ.0)  WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 50
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 75 L = 1,NSUB
              STEPS(L) = 0.0D0
   75     CONTINUE
          NTIMER = NMAX
          GO TO 50
      END IF
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME+TOFF, TCOMP, CPUTOT/60.0, ERRTOT, DETOT
   80     FORMAT (/,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6)
      END IF
*
      STOP
*
*       Set current global time.
  100 TTOT = TIME + TOFF
*
      RETURN
*
      END
