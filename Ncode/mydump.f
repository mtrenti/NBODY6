      SUBROUTINE MYDUMP(II,J)
*
*
*       COMMON save or read.
*       --------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      INCLUDE 'params.h'
      PARAMETER  (NA=64,NB=70,NC=532,ND=390+MLR+MLD+MLV,NE=24,NM=40,
     &    NG=84+2*KMAX,NL=99,NO=20*MCL+16,NP=32*NTMAX,NQ=31*MMAX,
     &      NS=44*MMAX)
      REAL*4  A,B,C,D,E,G,L,M,O,P,Q,S
      INTEGER K,I,NTSAVE
*
      COMMON/NAMES/  NTOT,NPAIRS,NTTOT,A(NA)
      COMMON/COUNTS/ B(NB)
      COMMON/PARAMS/ C(NC)
      COMMON/STARS/  D(ND)
      COMMON/PLPOT/  E(NE)
      COMMON/BLOCKS/ G(NG)
      COMMON/RAND2/  L(NL)
      COMMON/GALAXY/ M(NM)
      COMMON/CLOUDS/ O(NO)
      COMMON/MODES/  P(NP)
      COMMON/RCHE/   Q(NQ)
      COMMON/BINARY/ S(NS)

      COMMON/NBODY/  X(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),F(3,NMAX),
     &               FDOT(3,NMAX),BODY(NMAX),RS(NMAX),XDOT(3,NMAX),
     &               FI(3,NMAX),D1(3,NMAX),D2(3,NMAX),D3(3,NMAX),
     &               FR(3,NMAX),D1R(3,NMAX),D2R(3,NMAX),D3R(3,NMAX),
     &               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX),
     &               TNEW(NMAX),RADIUS(NMAX),TEV(NMAX),TEV0(NMAX),
     &               BODY0(NMAX),EPOCH(NMAX),SPIN(NMAX),XSTAR(NMAX),
     &               ZLMSTY(NMAX),FIDOT(3,NMAX),D0(3,NMAX),
     &               FRDOT(3,NMAX),D0R(3,NMAX),KSTAR(NMAX)
*
      COMMON/PAIRS/  U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),FU(4,KMAX),
     &               FUDOT(4,KMAX),FUDOT2(4,KMAX),FUDOT3(4,KMAX),
     &               H(KMAX),HDOT(KMAX),HDOT2(KMAX),HDOT3(KMAX),
     &               HDOT4(KMAX),DTAU(KMAX),TDOT2(KMAX),TDOT3(KMAX),
     &               R(KMAX),R0(KMAX),GAMMA(KMAX),SF(7,KMAX),H0(KMAX),
     &               FP0(4,KMAX),FD0(4,KMAX),TBLIST,DTB,KBLIST(KMAX),
     &               KSLOW(KMAX),NAME(NMAX),LIST(LMAX,NMAX)
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (II.NE.0) THEN

*       Check expanding arrays to include possible tidal tails (up to NTTOT).
        NTSAVE = NTOT
        IF (NTTOT.GT.0) THEN
            NTOT = NTTOT
        END IF

        WRITE (J) ntot,npairs,nttot,a,b,c,d,e,g,l,m,o,p,q,s

        WRITE (J) ((x(k,i),k=1,3),i=1,ntot),((x0(k,i),k=1,3),i=1,ntot)
     *   ,((x0dot(k,i),k=1,3),i=1,ntot),((f(k,i),k=1,3),i=1,ntot),
     *   ((fdot(k,i),k=1,3),i=1,ntot),(body(i),i=1,ntot),
     *   (rs(i),i=1,ntot),((xdot(k,i),k=1,3),i=1,ntot),
     *   ((fi(k,i),k=1,3),i=1,ntot),((d1(k,i),k=1,3),i=1,ntot),
     *   ((d2(k,i),k=1,3),i=1,ntot),((d3(k,i),k=1,3),i=1,ntot),
     *   ((fr(k,i),k=1,3),i=1,ntot),((d1r(k,i),k=1,3),i=1,ntot),
     *   ((d2r(k,i),k=1,3),i=1,ntot),((d3r(k,i),k=1,3),i=1,ntot),
     *   (step(i),i=1,ntot),(t0(i),i=1,ntot),(stepr(i),i=1,ntot),
     *   (t0r(i),i=1,ntot),(tnew(i),i=1,ntot),(radius(i),i=1,ntot),
     *   (tev(i),i=1,ntot),
     *   (tev0(i),i=1,ntot),(body0(i),i=1,ntot),(epoch(i),i=1,ntot),
     *   (spin(i),i=1,ntot),(xstar(i),i=1,ntot),(zlmsty(i),i=1,ntot),
     *   ((fidot(k,i),k=1,3),i=1,ntot),((d0(k,i),k=1,3),i=1,ntot),
     *   ((frdot(k,i),k=1,3),i=1,ntot),((d0r(k,i),k=1,3),i=1,ntot),
     *   (kstar(i),i=1,ntot)

        write (J) ((u(k,i),k=1,4),i=1,npairs),((u0(k,i),k=1,4),i=1,
     *    npairs),((udot(k,i),k=1,4),i=1,npairs),((fu(k,i),k=1,4),i=1,
     *    npairs),((fudot(k,i),k=1,4),i=1,npairs),((fudot2(k,i),k=1,4),
     *    i=1,npairs),((fudot3(k,i),k=1,4),i=1,npairs),(h(i),i=1,
     *    npairs),(hdot(i),i=1,npairs),(hdot2(i),i=1,npairs),  
     *    (hdot3(i),i=1,npairs),(hdot4(i),i=1,npairs),(dtau(i),
     *    i=1,npairs),(tdot2(i),i=1,npairs),(tdot3(i),i=1,npairs),
     *    (r(i),i=1,npairs),(r0(i),i=1,npairs),(gamma(i),i=1,npairs),
     *    ((sf(k,i),k=1,7),i=1,npairs),(h0(i),i=1,npairs),((fp0(k,i), 
     *    k=1,4),i=1,npairs),((fd0(k,i),k=1,4),i=1,npairs),tblist,dtb,
     *    (kblist(i),i=1,kmax),(kslow(i),i=1,npairs),
     *       (name(i),i=1,ntot)

        write (J) ((list(k,i),k=1,list(1,i)+1),i=1,ntot)

        END FILE J
        CLOSE (UNIT=J)
*       Restore standard array pointer.
        NTOT = NTSAVE
      else
        READ (J) ntot,npairs,nttot,a,b,c,d,e,g,l,m,o,p,q,s

        if (ntot.gt.nmax) then
          write (*,*) "DANGER NTOT > NMAX !"
          stop
        end if

        if (npairs.gt.kmax) then
          write (*,*) "DANGER NPAIRS > KMAX !"
          stop
        end if 

        NTSAVE = NTOT
        IF (NTTOT.GT.0) THEN
            NTOT = NTTOT
        END IF
         
        read (J) ((x(k,i),k=1,3),i=1,ntot),((x0(k,i),k=1,3),i=1,ntot)
     *   ,((x0dot(k,i),k=1,3),i=1,ntot),((f(k,i),k=1,3),i=1,ntot),
     *   ((fdot(k,i),k=1,3),i=1,ntot),(body(i),i=1,ntot),
     *   (rs(i),i=1,ntot),((xdot(k,i),k=1,3),i=1,ntot),
     *   ((fi(k,i),k=1,3),i=1,ntot),((d1(k,i),k=1,3),i=1,ntot),
     *   ((d2(k,i),k=1,3),i=1,ntot),((d3(k,i),k=1,3),i=1,ntot),
     *   ((fr(k,i),k=1,3),i=1,ntot),((d1r(k,i),k=1,3),i=1,ntot),
     *   ((d2r(k,i),k=1,3),i=1,ntot),((d3r(k,i),k=1,3),i=1,ntot),
     *   (step(i),i=1,ntot),(t0(i),i=1,ntot),(stepr(i),i=1,ntot),
     *   (t0r(i),i=1,ntot),(tnew(i),i=1,ntot),(radius(i),i=1,ntot),
     *   (tev(i),i=1,ntot),
     *   (tev0(i),i=1,ntot),(body0(i),i=1,ntot),(epoch(i),i=1,ntot),
     *   (spin(i),i=1,ntot),(xstar(i),i=1,ntot),(zlmsty(i),i=1,ntot),
     *   ((fidot(k,i),k=1,3),i=1,ntot),((d0(k,i),k=1,3),i=1,ntot),
     *   ((frdot(k,i),k=1,3),i=1,ntot),((d0r(k,i),k=1,3),i=1,ntot),
     *   (kstar(i),i=1,ntot)

        read (J) ((u(k,i),k=1,4),i=1,npairs),((u0(k,i),k=1,4),i=1,
     *    npairs),((udot(k,i),k=1,4),i=1,npairs),((fu(k,i),k=1,4),i=1,
     *    npairs),((fudot(k,i),k=1,4),i=1,npairs),((fudot2(k,i),k=1,4),
     *    i=1,npairs),((fudot3(k,i),k=1,4),i=1,npairs),(h(i),i=1,
     *    npairs),(hdot(i),i=1,npairs),(hdot2(i),i=1,npairs),
     *    (hdot3(i),i=1,npairs),(hdot4(i),i=1,npairs),(dtau(i),
     *    i=1,npairs),(tdot2(i),i=1,npairs),(tdot3(i),i=1,npairs),
     *    (r(i),i=1,npairs),(r0(i),i=1,npairs),(gamma(i),i=1,npairs),
     *    ((sf(k,i),k=1,7),i=1,npairs),(h0(i),i=1,npairs),((fp0(k,i),
     *    k=1,4),i=1,npairs),((fd0(k,i),k=1,4),i=1,npairs),tblist,dtb,
     *    (kblist(i),i=1,kmax),(kslow(i),i=1,npairs),
     *          (name(i),i=1,ntot)

        read (J) (list(1,i),(list(k,i),k=2,list(1,i)+1),i=1,ntot)
        NTOT = NTSAVE
      END IF
*
      RETURN
*
      END
