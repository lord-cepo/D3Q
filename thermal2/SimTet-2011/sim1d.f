C file: s1d0onei.f
C date: 2011-03-20
C who:  S.Kaprzyk
C what: Complex linear form integral over standard 1-d simplex
C ----------------------------------------------------
C -         *                                        -
C -        *                         1               -
C - SIM0=  * dt1 ----------------------------------  -
C -        *     verm(2)+(verm(1)-verm(2))*t1        -
C -       *                                          -
C -    0<t1<1                                        -
C ----------------------------------------------------
      SUBROUTINE S1D0ONEI(SIM0, SIM0I, VERM)
      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX SIM0, SIM0I
      DOUBLE COMPLEX VERM(*)
C
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(2),LEPS(1)
c
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      INTEGER          I, K, N
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, EPS
      DATA             EPS/1.0D-6/
      DATA             ZERO/0.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      N = 2
      DO 100 I = 1, 1
       IF (DIMAG(VERM(N))*DIMAG(VERM(I)).LT.ZERO) THEN
         WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,2)
 9010    FORMAT (' ***s1d0onei: not signed ImgVERM()=',4(D13.6,1X))
c        STOP ' ***s1d0onei: '
       END IF
 100  CONTINUE
      DO 101 I = 1, 1
      IF (CDABS(VERM(I)).GT.CDABS(VERM(N))) N = I
 101  CONTINUE  
C Here are 2 cases, how w() are placed on complex-plane
C LONE(1..2) [w1<ONE]; [ONE<w1]
C LEPS(1) null
*-ASIS
      CALL S1D0LEPS(
     &  VERM, N, U, AS,
     &  LONE, EPS, LEPS,
     &  iuerr
     & )
c U(1) = W(1);
      IF (LONE(1)) THEN
        CALL S1D1ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
c U(1) = CONE/W(1)
      IF (LONE(2)) THEN
        CALL S1D2ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
C
      SIM0 =  SIM/VERM(N)
      SIM0I= CDLOG(VERM(N)) + SIMI - CONE
      RETURN
      END SUBROUTINE S1D0ONEI
C
C ----------------------------------------------------------------------
C file: s1d0twoi.f
C date: 2011-03-24
C who:  S.Kaprzyk
C what: Complex linear form integrals over standard 1d-simplex
C ------------------(i=1)--------------------------------------
C            *                                                -
C           *            t_i                                  -
C  VERL(i)= *dt1 ------------------------------------------   -
C           *        VERM(2)+[VERM(1)-VERM(2)]*t1             -
C          *                                                  -
C     0<t1<1                                                  -
C -------------------(i=2)-------------------------------------
C            *                                                -
C           *           1- t1                                 -
C  VERL(2)= *dt1 ------------------------------------------   -
C           *          VERM(2)+[VERM(1)-VERM(2)]*t1           -
C          *                                                  -
C     0<t1<1                                                  -
C -------------------------------------------------------------
      SUBROUTINE   S1D0TWOI(VERL, VERLI, VERM)
      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX VERL(*), VERLI(*), VERM(*)
C
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(2),LEPS(1)
c
      INTEGER          I, K, N
      DOUBLE PRECISION AS(1)
      DOUBLE COMPLEX   U(1)
      DOUBLE COMPLEX   CZERO
      DOUBLE PRECISION EPS
      DOUBLE PRECISION ZERO
      DATA             EPS/1.0D-6/
      DATA             ZERO/0.0D0/
C
      CZERO = (0.0D0,0.0D0)
      DO 50 I = 1, 1
        IF (DIMAG(VERM(2))*DIMAG(VERM(I)).LT.ZERO) THEN
         WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,2)
 9010    FORMAT (' ***s1d0twoi: not signed ImgVERM()=',2(D13.6,1X))
         STOP ' ***s1dt0woi: '
        END IF
 50    CONTINUE
      DO 100 N = 1, 2
      VERL(N) = CZERO
C Here are 2 cases, how w() are placed on complex-plane
C LONE(1..2) [w1<ONE]; [ONE<w1]
C LEPS(1) null
*-ASIS
       CALL S1D0LEPS(
     &  VERM, N, U, AS,
     &  LONE, EPS, LEPS,
     &  iuerr
     & )
c U(1) = W(1);
       IF (LONE(1)) THEN
        CALL S1D1TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
c U(1) = CONE/W(1);
       IF (LONE(2)) THEN
        CALL S1D2TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
c
       VERL(N)  = 0.50D0*SIM/VERM(N)
       VERLI(N) = 0.50D0*CDLOG(VERM(N)) + 0.25D0*SIMI -0.25D0
 100  CONTINUE ! end DO .. N = 1, 2
      RETURN
      END SUBROUTINE S1D0TWOI
C
C ----------------------------------------------------------------------
C file: s1d1onei.f
C date: 2010-09-22
C who:  S.Kaprzyk
C what: CASE a) : [w1<ONE]
C       U(1)=W(1)
C LEPS() null
      SUBROUTINE S1D1ONEI(U,AS,LEPS, SIM1, SIM1I)
      IMPLICIT       NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM1, SIM1I
C
      DOUBLE COMPLEX   SIM0UR0
      EXTERNAL         SIM0UR0
      DOUBLE COMPLEX   R0(1)
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
C ------------------------------------------------------------------
C      R0(U) = ARTH(U)/U
C      R0(U) = 1+U**2/3+U**6*X0(U)
C ------------------------------------------------------------------
      R0(1) = SIM0UR0(U(1))
      SIM1  = (CONE+U(1))*R0(1)
      SIM1I = (CONE-U(1))*R0(1)
      RETURN
      END SUBROUTINE S1D1ONEI
C
C ----------------------------------------------------------------------
C file: s1d1twoi.f
C date: 2011-03-25
C who:  S.Kaprzyk
C what: CASE a) : [w1<ONE]
C       U(1)=W(1)
C LEPS() null
      SUBROUTINE S1D1TWOI(U,AS,LEPS, SIM, SIMI)
      IMPLICIT   NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
C
      DOUBLE COMPLEX SIM0UX0
      EXTERNAL       SIM0UX0
c
      DOUBLE COMPLEX X0(1), R2(1), S2(1)
      DOUBLE COMPLEX AV
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
C -----------------------------------------------------------------
C -                                                               -
C -    R2(W)=1+(W-1)*(ARTH(W)/W-1)/W=[1+(W-1)*R0(W)]/W            -
C -    R2(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
C -                                                               -
C -----------------------------------------------------------------
       AV = U(1)
       X0(1) = SIM0UX0(AV)
*-ASIS
       R2(1) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 +
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(1)))))
       S2(1) = (CONE-U(1))/(CONE+U(1))*R2(1)
      SIM  = (CONE+U(1))*R2(1)
      SIMI = (CONE+U(1))*S2(1)
      RETURN
      END SUBROUTINE S1D1TWOI
C
C ----------------------------------------------------------------------
C file: s1d2onei.f
C date: 2011-03-20
C who:  S.Kaprzyk
C what: CASE b) : [ ONE<w1]
C       U(1) = CONE/W(1)
C LEPS(1) null
      SUBROUTINE  S1D2ONEI(U,AS,LEPS,SIM2,SIM2I)
      IMPLICIT       NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM2, SIM2I
c
      DOUBLE COMPLEX   SIM0UR0
      EXTERNAL         SIM0UR0
      DOUBLE COMPLEX   R0(1)
      DOUBLE COMPLEX   CZERO, CONE, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
C ------------------------------------------------------------------
C       R0(U)=ARTH(U)/U
C       R0(U)=1+U**2/3+U**6*X0(U)
C ------------------------------------------------------------------
      R0(1) = SIM0UR0(U(1))
      SIM2 = (U(1)+CONE)*(CIMAG*AS(1)+U(1)*R0(1))
      SIM2I= (U(1)-CONE)*(CIMAG*AS(1)+U(1)*R0(1))
c
      RETURN
      END SUBROUTINE S1D2ONEI
C
C ----------------------------------------------------------------------
C s1d2onei()
C file: s1d2twoi.f
C date: 2011-03-25
C who:  S.Kaprzyk
C what: CASE b) : [ ONE<w1]
C       U(1) = CONE/W(1)
C LEPS(1) null
      SUBROUTINE S1D2TWOI(U,AS,LEPS, SIM, SIMI)
      IMPLICIT   NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
C
      DOUBLE COMPLEX SIM0UX0
      EXTERNAL       SIM0UX0
c
      DOUBLE COMPLEX X0(2), R2(2), S2(2)
      DOUBLE COMPLEX AV
      INTEGER        I
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
C   -----------------------------------------------------------------
C   -                                                               -
C   -        R2(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
C   -    R2(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
C   -                                                               -
C   -----------------------------------------------------------------
       AV = U(1)
       X0(1) = SIM0UX0(AV)
*-ASIS
       R2(1) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 +
     &    AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(1)))))
C
      R2(1)= (CONE-U(1))*DCMPLX(ZERO,AS(1))+CONE+U(1)-U(1)*U(1)*R2(1)
      S2(1) = (U(1)-CONE)/(U(1)+CONE)*R2(1)
      SIM  = (U(1)+CONE)*R2(1)
      SIMI = (U(1)+CONE)*S2(1)
C
      RETURN
      END SUBROUTINE S1D2TWOI
C
C ----------------------------------------------------------------------
C file: s1d0leps.f
C date: 2011-04-02
C who:  S.Kaprzyk
C what: Return a proper case:
C       a) [w1<ONE]; b) [ONE<w1]
C INPUT:
C VERM(1,2) - two complex numbers
C N - integer number for W(I)=(VERM(I)-VERM(N))/(VERM(I)+VERM(N))
C EPS  - critical distance
C OUTPUT:
C LONE(2) [w1<ONE]; [ONE<w1]
C LEPS(1) null
      SUBROUTINE  S1D0LEPS(
     &  VERM, N, W, AS,
     &  LONE, EPS, LEPS,
     &  iuerr
     & )
      IMPLICIT       NONE
      DOUBLE COMPLEX   VERM(2)
      INTEGER          N
      DOUBLE COMPLEX   W(1)
      DOUBLE PRECISION AS(1)
      DOUBLE PRECISION EPS
      LOGICAL          LONE(2), LEPS(1)
      INTEGER          iuerr
C
c      DOUBLE PRECISION DLAMCH
c      EXTERNAl         DLAMCH
C
      LOGICAL          LWONE(1)
      DOUBLE PRECISION AL(1), AIW, OZERO
      INTEGER          I, II, K
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE, PI
      DATA             ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/
      DATA             OZERO/ 1.0D-13/
c      OZERO = DLAMCH('e')*10
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      IF ((N.GT.2).OR.(N.LT.1)) THEN
       WRITE (iuerr,9010) N
 9010  FORMAT (' ***s1d0leps: N =',I6,' must be 1, or 2')
       STOP ' ***s1d0leps: '
      END IF
C
      II = 0
      DO 200 I = 1, 2
       IF (I.EQ.N) GO TO 200
       II = II + 1
       IF (CDABS(VERM(N)-VERM(I)).LT.CDABS(VERM(N)+VERM(I))) THEN
        LWONE(II) = .TRUE.
        W(II) = (VERM(N)-VERM(I))/(VERM(N)+VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = 0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       ELSE
        LWONE(II) = .FALSE.
        W(II) = (VERM(N)+VERM(I))/(VERM(N)-VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = -0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       END IF
       AL(II) = CDABS(W(II))
 200  CONTINUE
C
      LONE(1) = .FALSE.
      LONE(2) = .FALSE.
      IF (LWONE(1)) THEN
       LONE(1) = .TRUE.
      END IF
      IF (.NOT.LWONE(1)) THEN
       LONE(2) = .TRUE.
      END IF
C here are cases, how close are w()
      LEPS(1) = .FALSE.
c
      IF ((.NOT.LONE(1)).AND.(.NOT.LONE(2))) THEN
       WRITE (iuerr,9030) (LONE(I),I=1,2)
 9030  FORMAT ('***s1d0leps: LONE()=',2(L1,1X),'not in order')
       STOP '***s1d0leps: '
      END IF
c
      RETURN
      END SUBROUTINE S1D0LEPS
C
