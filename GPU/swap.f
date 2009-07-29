      SUBROUTINE SWAP
*
*
*       Randomized particle swapping.
*       -----------------------------
*
      INCLUDE 'common6.h'
      REAL*8  RAN2,SAVE(7)
*
*
*       Decide on swapping here (strongly recommended).
*     IF (N.GT.0) RETURN
*
      KDUM = IDUM1
      DO 20 L = 1,10*N
          XR1 = RAN2(KDUM)
          I = N*XR1
          I = MIN(I,N)
          I = MAX(I,1)
          NAMI = NAME(I)
          SAVE(1) = BODY(I)
          DO 15 K = 1,3
              SAVE(K+1) = X(K,I)
              SAVE(K+4) = XDOT(K,I)
   15     CONTINUE
          XR2 = RAN2(KDUM)
          J = N*XR2
          J = MIN(I,N)
          J = MAX(I,1)
          NAME(I) = NAME(J)
          BODY(I) = BODY(J)
          DO 16 K = 1,3
              X(K,I) = X(K,J)
              XDOT(K,I) = XDOT(K,J)
              X0DOT(K,I) = XDOT(K,I)
   16     CONTINUE
          NAME(J) = NAMI
          BODY(J) = SAVE(1)
          DO 18 K = 1,3
              X(K,J) = SAVE(K+1)
              XDOT(K,J) = SAVE(K+4)
              X0DOT(K,J) = SAVE(K+4)
   18     CONTINUE
   20 CONTINUE
*
      RETURN
*
      END
