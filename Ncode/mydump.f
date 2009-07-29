      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      INCLUDE 'params.h'
      PARAMETER  (NA=96*NMAX,NB=111*KMAX+7,
     &            NC=(LMAX+1)*NMAX+3*KMAX+MLR+MLD+MLV+65,
     &            NF=44*MMAX,NG=17*NMAX+390,NH=20*MCL+16,
     &            NO=24*NMAX,NP=2*NMAX+84,NM=32*NTMAX,NR=31*MMAX)
      REAL*4  A,B,C, E,F,G,H,O,P,Q,PL,R
      INTEGER  IC,ID,IR
*
*
      COMMON/NBODY/  A(NA)
      COMMON/PAIRS/  B(NB)
      COMMON/NAMES/  IC(NC)
      COMMON/COUNTS/ ID(70)
      COMMON/PARAMS/ E(334)
      COMMON/BINARY/ F(NF)
      COMMON/STARS/  G(NG)
      COMMON/MODES/  C(NM)
      COMMON/CLOUDS/ H(NH)
      COMMON/RAND2/  IR(99)
      COMMON/HERMIT/ O(NO)
      COMMON/BLOCKS/ P(NP)
      COMMON/GALAXY/ Q(40)
      COMMON/PLPOT/  PL(24)
      COMMON/RCHE/   R(NR)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, B, IC, ID, E, F, G, C, H, IR, O, P, Q, PL, R
      ELSE
          WRITE (J)  A, B, IC, ID, E, F, G, C, H, IR, O, P, Q, PL, R
          END FILE J
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
