C file: s2d0onei.f
C date: 2011-04-01
C who:  S.Kaprzyk
C what: Complex linear form integral over standard 2-d simplex
C ----------------------------------------------------
C -         *                                        -
C -        *                         1               -
C - SIM0= 2* dt1dt2 ----------------------------------
C -        *     verm(3)+(verm(1)-verm(3))*t1+..     -
C -       *                                          -
C -    0<t1<1                                        -
C ----------------------------------------------------
      SUBROUTINE S2D0ONEI(SIM0, SIM0I, VERM)
      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX SIM0, SIM0I
      DOUBLE COMPLEX VERM (3)
C
      DOUBLE COMPLEX   VERM1D(2)
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(3),LEPS(1)
c
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2), AL(3)
      INTEGER          I, I1, I2, K, N
      DOUBLE PRECISION ZERO, EPS, SMALL
      DATA             EPS/1.0D-6/, SMALL/1.0D-5/
      DATA             ZERO/0.0D0/
C
      DO 100 I = 1, 2
       IF (DIMAG(VERM(3))*DIMAG(VERM(I)).LT.ZERO) THEN
        WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,3)
 9010   FORMAT (' ***s2d0onei: not signed ImgVERM()=',3(D13.6,1X))
c        STOP ' ***s2d0onei: '
       END IF
 100  CONTINUE
C
      DO 200 I = 1, 3
       AL(I) = CDABS(VERM(I))
 200  CONTINUE
      DO 300 N = 1, 3
       IF (AL(N).LT.SMALL) THEN
        I1 = MOD(N,3) + 1
        I2 = MOD(N+1,3) + 1
        VERM1D(1) = VERM(I1)
        VERM1D(2) = VERM(I2)
        CALL S1D0ONEI(SIM,SIMI,VERM1D)
        SIM0  = 2*SIM
        SIM0I = SIMI - 0.5D0
        RETURN
       END IF
 300  CONTINUE
C
      N = 3
      DO 400 I = 1, 2
       IF (AL(I).GT.AL(N)) N = I
 400  CONTINUE
C Here are 3 cases, how w() are placed on complex-plane
C LONE(1..3) [w1<w2<ONE];[w1<ONE<w2];[ONE<w1<w2]
C LEPS(1) [|w1-w2|<e
*-ASIS
      CALL S2D0LEPS(
     &  VERM, N, U, AS,
     &  LONE, EPS, LEPS,
     &  iuerr
     & )
c U(1) = W(1); U(2) = W(2);
      IF (LONE(1)) THEN
       CALL S2D1ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
c U(1) = W(1); U(2) = 1/W(2);
      IF (LONE(2)) THEN
       CALL S2D2ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
c U(1) = 1/W(1); U(2) = 1/W(2);
      IF (LONE(3)) THEN
       CALL S2D3ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
C
      SIM0 = SIM/VERM(N)
      SIM0I = CDLOG(VERM(N))+0.5D0*SIMI - 1.5D0
      RETURN
      END SUBROUTINE S2D0ONEI
C
C ----------------------------------------------------------------------
C file: s2d0twoi.f
C date: 2011-04-01
C who:  S.Kaprzyk
C what: Complex linear form integrals over standard 2-d simplex
C ------------------(i=1,2)----------------------------------------
C            *                                                    -
C           *            t_i                                      -
C  VERL(i)=2*dt1dt2 --------------------------------------------- -
C           *          VERM(3)+[VERM(1)-VERM(3)]*t1+...           -
C          *                                                      -
C     0<t1+t2<1                                                   -
C -------------------(i=3)-----------------------------------------
C            *                                                    -
C           *           1- t1 -t2
C  VERL(3)=2*dt1dt2 -----------------------------------------------
C           *          VERM(3)+[VERM(1)-VERM(3)]*t1+...           -
C          *                                                      -
C     0<t1+t2<1                                                   -
C -----------------------------------------------------------------
      SUBROUTINE   S2D0TWOI(VERL, VERLI, VERM)
      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX VERL(3), VERLI(3), VERM(3)
C
      DOUBLE COMPLEX  SIM, SIMI
      DOUBLE COMPLEX  VERL1D(2), VERL1DI(2), VERM1D(2)
      LOGICAL         LONE(3),LEPS(1)
c
      INTEGER          I, I1, I2, K, N
      DOUBLE PRECISION AS(2), AL(3)
      DOUBLE COMPLEX   U(2)
      DOUBLE COMPLEX   CZERO
      DOUBLE PRECISION EPS, SMALL
      DOUBLE PRECISION ZERO
      DATA             EPS/1.0D-6/, SMALL/1.0D-5/
      DATA             ZERO/0.0D0/
C
      CZERO = (0.0D0,0.0D0)
      DO 100 I = 1, 2
       IF (DIMAG(VERM(3))*DIMAG(VERM(I)).LT.ZERO) THEN
        WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,3)
 9010   FORMAT (' ***s2d0twoi: not signed ImgVERM()=',3(D13.6,1X))
        STOP ' ***s2d0twoi: '
       END IF
 100  CONTINUE
      DO 200 I = 1, 3
       AL(I) = CDABS(VERM(I))
 200  CONTINUE
C
      DO 300 N = 1, 3
       VERL(N) = CZERO
       VERLI(N) = CZERO
       I1 = MOD(N,3) + 1
       I2 = MOD(N+1,3) + 1
       IF (AL(I1).LT.SMALL) THEN
        VERM1D(1) = VERM(I2)
        VERM1D(2) = VERM(N)
        CALL S1D0TWOI(VERL1D,VERL1DI,VERM1D)
        VERL(N) = VERL1D(2)
        VERLI(N) = VERL1DI(2)*2/3
        GO TO 300
       END IF
       IF (AL(I2).LT.SMALL) THEN
        VERM1D(1) = VERM(I1)
        VERM1D(2) = VERM(N)
        CALL S1D0TWOI(VERL1D,VERL1DI,VERM1D)
        VERL(N) = VERL1D(2)
        VERLI(N) = VERL1DI(2)*2/3
        GO TO 300
       END IF
       IF (AL(N).LT.SMALL) THEN
        VERM1D(1) = VERM(I1)
        VERM1D(2) = VERM(I2)
        CALL S1D0ONEI(SIM,SIMI,VERM1D)
        VERL(N) = SIM
        VERLI(N) = SIMI/3 - 0.5D0
        GO TO 300
       END IF
C Here are 3 cases, how w() are placed on complex-plane
C LONE(1..3) [w1<w2<ONE];[w1<ONE<w2];[ONE<w1<w2]
C LEPS(1) [|w1-w2|<e
*-ASIS
       CALL S2D0LEPS(
     &  VERM, N, U, AS,
     &  LONE, EPS, LEPS,
     &  iuerr
     & )
c U(1) = W(1); U(2) = W(2)
       IF (LONE(1)) THEN
        CALL S2D1TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
c U(1) = W(1); U(2) = CONE/W(2)
       IF (LONE(2)) THEN
        CALL S2D2TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
c U(1) = CONE/W(1); U(2) = CONE/W(2)
       IF (LONE(3)) THEN
        CALL S2D3TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
c
       VERL(N) = SIM/VERM(N)/4
       VERLI(N) = (CDLOG(VERM(N))+SIMI/4)/3 - 5.0D0/18.0D0
 300  CONTINUE ! end DO ... N = 1, 3
      RETURN
      END SUBROUTINE S2D0TWOI
C
C ----------------------------------------------------------------------
C file: s2d1onei.f
C date: 2011-03-23
C who:  S.Kaprzyk
C what: CASE a) : [w1<w2<ONE]
C       U(1)=W(1); U(2)=W(2)
C LEPS(1)  [|w1-w2|<e]
      SUBROUTINE S2D1ONEI(U,AS,LEPS, SIM1, SIM1I)
      IMPLICIT       NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM1, SIM1I
C
      DOUBLE COMPLEX   SIM0UX0
      EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(2), DX0(4), R0(2), S0(2)
      DOUBLE COMPLEX   QX0(1), QR0(1), QS0(1)
      DOUBLE COMPLEX   AV, BV, CV, DU, G(5)
      INTEGER          I, K, N
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
      DO 100 I = 1, 2
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
C -----------------------------------------------------------------
C -   The DX0(I) contains  derivatives of the function X0(u)      -
C -----------------------------------------------------------------
c      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
      IF (LEPS(1)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
c      DO 200 I = 1, 3, 2
      I = 1
       DU = U(I) - U(2)
       IF (.NOT.LEPS(I)) THEN
        QX0(I) = (X0(I)-X0(2))/DU
       ELSE
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+
     &                     DX0(4)/24.0D0*DU)*DU)*DU
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
       QS0(I) = (-S0(2)-R0(2)+(CONE-U(I))*QR0(I))/(CONE+U(I))
C
      AV = (CONE+U(1))*(CONE+U(2))
      BV = (U(2)-CONE)
      SIM1  = AV*(R0(1) + QR0(1)*BV)
      SIM1I = AV*(S0(1) + QS0(1)*BV)
      RETURN
      END SUBROUTINE S2D1ONEI
C
C ----------------------------------------------------------------------
C file: s2d1twoi.f
C date: 2010-09-22
C who:  S.Kaprzyk
C what: CASE a) : [w1<w2<ONE]
C       U(1)=W(1); U(2)=W(2)
C LEPS(1)  [|w1-w2|<e]
      SUBROUTINE S2D1TWOI(U,AS,LEPS, SIM, SIMI)
      IMPLICIT   NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
C
      DOUBLE COMPLEX SIM0UX0
      EXTERNAL       SIM0UX0
c
      DOUBLE COMPLEX X0(2), DX0(4), R3(2), S3(2)
      DOUBLE COMPLEX QR3(1), QX0(1)
      DOUBLE COMPLEX QS3(1)
      DOUBLE COMPLEX AV, BV, UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CTWO  = (2.0D0,0.0D0)
C -----------------------------------------------------------------
C -                                                               -
C -        R3(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
C -    R3(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
C -                                                               -
C -----------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)
       X0(I) = SIM0UX0(AV)
*-ASIS
       R3(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 +
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
       S3(I) = (CONE-U(I))/(CONE+U(I))*R3(I)
 100  CONTINUE
C -----------------------------------------------------------------
C -   The DX0(I) Contains  Derivatives Of The Function X0(W)      -
C -----------------------------------------------------------------
      IF (LEPS(1)) THEN
        CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
C
c      DO 200 I = 1, 3, 2
      I = 1
       UV = U(I) - U(2)
       IF (LEPS(I)) THEN
       QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+
     &                    DX0(4)/24.0D0*UV)*UV)*UV
       ELSE
        QX0(I) = (X0(I)-X0(2))/UV
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR3(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 +
     &   G(3)/5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
       QS3(I) = (-S3(2)-R3(2)+(CONE-U(I))*QR3(I))/(CONE+U(I))
c 200  CONTINUE  ! DO 750 I = 1, 3, 2
c
      AV = (CONE+U(1))*(CONE+U(2))
      BV = U(2)-CONE
      SIM  = AV*(R3(1) + QR3(1)*BV)
      SIMI = AV*(S3(1) + QS3(1)*BV)
      RETURN
      END SUBROUTINE S2D1TWOI
C
C ----------------------------------------------------------------------
C file: s2d2onei.f
C date: 2011-03-24
C who:  S.Kaprzyk
C what: CASE b) : [w1<ONE<w2]
C       U(1)=W(1); U(2)= CONE/W(2)
C LEPS(1) [|w1-w2|<e]
      SUBROUTINE  S2D2ONEI(U,AS,LEPS, SIM2, SIM2I)
      IMPLICIT         NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM2, SIM2I
c
      DOUBLE COMPLEX   SIM0UX0
      EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(2), R0(2), S0(2), DX0(4)
c      DOUBLE COMPLEX   QX0(1), QR0(1), QS0(1)
      DOUBLE COMPLEX   AV, BV, CV, DV, G(5)
      DOUBLE COMPLEX   UD
      INTEGER          I, K, N
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CTWO  = (2.0D0,0.0D0)
C ------------------------------------------------------------------
C      R0(U) = ARTH(U)/U
C      R0(U) = 1+U**2/3+U**6*X0(U)
C ------------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
c
      AV = (CONE+U(1))*(U(2)+CONE)/(CONE-U(2)*U(1))
      BV = -(U(1)-CONE)
      CV = (CONE-U(2))
      SIM2 = AV*(BV*R0(1)+CV*(DCMPLX(ZERO,AS(2))+U(2)*R0(2)))
      CV = CV*(U(2)-CONE)/(U(2)+CONE)
      SIM2I= AV*(BV*S0(1)+CV*(DCMPLX(ZERO,AS(2))+U(2)*R0(2)))
c
      RETURN
      END SUBROUTINE S2D2ONEI
C
C ----------------------------------------------------------------------
C file: s2d2twoi.f
C date: 2011-03-26
C who:  S.Kaprzyk
C what: CASE b) : [w1<ONE<w2]
C       U(1)=W(1); U(2)= CONE/W(2)
C LEPS(1) [|w1-w2|<e]
      SUBROUTINE S2D2TWOI(U,AS,LEPS, SIM, SIMI)
      IMPLICIT   NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
C
      DOUBLE COMPLEX SIM0UX0
      EXTERNAL       SIM0UX0
c
      DOUBLE COMPLEX X0(2), R3(2), S3(2)
      DOUBLE COMPLEX AV, BV
      INTEGER          I
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
C   -----------------------------------------------------------------
C   -                                                               -
C   -        R3(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
C   -    R3(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
C   -                                                               -
C   -----------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)
       X0(I) = SIM0UX0(AV)
*-ASIS
       R3(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 +
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
 100  CONTINUE
      S3(1) = (CONE-U(1))/(CONE+U(1))*R3(1)
      AV = (CONE+U(1))*(U(2)+CONE)/(U(1)*U(2)-CONE)
      BV = (CONE-U(2))*DCMPLX(ZERO,AS(2))+CONE+U(2)-U(2)*U(2)*R3(2)
      SIM = AV*((U(1)-CONE)*R3(1) - (CONE-U(2))*BV)
      BV = (U(2)-CONE)/(U(2)+CONE)*BV
      SIMI= AV*((U(1)-CONE)*S3(1) - (CONE-U(2))*BV)
C
      RETURN
      END SUBROUTINE S2D2TWOI
C
C ----------------------------------------------------------------------
C file: s2d3onei.f
C date: 2011-04-01
C who:  S.Kaprzyk
C what: CASE c) : [ ONE<w1<w2]
C       U(1)=CONE/W(1); U(2)=CONE/W(2)
C LEPS(1) [|w1-w2|<e]
      SUBROUTINE  S2D3ONEI(U,AS,LEPS,SIM3,SIM3I)
      IMPLICIT       NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM3, SIM3I
c
      DOUBLE COMPLEX   SIM0UX0
      EXTERNAL         SIM0UX0
c
      DOUBLE COMPLEX   X0(2), QX0(1), DX0(4)
      DOUBLE COMPLEX   R0(2), UR0(2), QR0(1), QUR0(1)
      DOUBLE COMPLEX   US0(2), QUS0(1), QAS(1)
      DOUBLE COMPLEX   AV, UD,  G(5)
      INTEGER          I, K
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
      DO 100 I = 1, 2
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       UR0(I)= DCMPLX(ZERO,AS(I)) + U(I)*R0(I)
 100  CONTINUE
C -----------------------------------------------------------------
C -   THE DX0(I) CONTAINS  DERIVATIVES OF THE FUNCTION X0(U)
C -----------------------------------------------------------------
      IF (LEPS(1)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
c
      I = 1
       UD = U(I) - U(2)
       IF (.NOT.LEPS(I)) THEN
        QX0(I) = (X0(I)-X0(2))/UD
        QAS(I) = (AS(I)-AS(2))/UD
       ELSE
         QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0 +
     &            DX0(4)/24.0D0*UD)*UD)*UD
        QAS(I) = CZERO
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
c
       QUR0(1) = CIMAG*QAS(1) + R0(2) +U(1)*QR0(1)
       US0(2) = (U(2)-CONE)/(U(2)+CONE)*UR0(2)
       QUS0(1) = ( -US0(2)+UR0(2)+(U(1)-CONE)*QUR0(1))/(U(1)+CONE)
c
      AV = (U(1)+CONE)*(U(2)+CONE)
      SIM3  = AV*(UR0(2) + (U(1)-CONE)*QUR0(1))
      SIM3I = AV*(US0(2) + (U(1)-CONE)*QUS0(1))
      RETURN
      END SUBROUTINE S2D3ONEI
C
C ----------------------------------------------------------------------
C file: s2d3twoi.f
C date: 2011-03-28
C who:  S.Kaprzyk
C what: CASE c) : [ ONE<w1<w2]
C       U(1)=CONE/W(1); U(2)=CONE/W(2)
C LEPS(1) [|w1-w2|<e]
      SUBROUTINE S2D3TWOI(U,AS,LEPS,SIM,SIMI)
      IMPLICIT   NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
C
      DOUBLE COMPLEX SIM0UX0
      EXTERNAL       SIM0UX0
c
      DOUBLE COMPLEX X0(2), QX0(1), DX0(4)
      DOUBLE COMPLEX R3(2), UR3(2), QR3(1), QUR3(1)
      DOUBLE COMPLEX US3(2), QUS3(1)
      DOUBLE COMPLEX QAS(1), AV,  UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CTWO  = (2.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
C   -----------------------------------------------------------------
C   -                                                               -
C   -        R3(U)=1+(U-1)*(ARTH(U)/U-1)/U                          -
C   -    R3(U)=1-U/3+U**2/3-U**3/5+U**4/5+(U-1)*U**5*X0(U)          -
C   -                                                               -
C   -----------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)
       X0(I) = SIM0UX0(AV)
*-ASIS
       R3(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 +
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
       UR3(I)=(CONE-U(I))*DCMPLX(ZERO,AS(I))+CONE+U(I)-U(I)*U(I)*R3(I)
 100  CONTINUE  ! end DO ... I = 1, 2
C -----------------------------------------------------------------
C -   The DX0(I) Contains  Derivatives Of The Function X0(U)      -
C -----------------------------------------------------------------
      IF (LEPS(1)) THEN
        CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
C
       I = 1
       UV = U(I) - U(2)
       IF (LEPS(I)) THEN
        QX0(I) = DX0(1) +
     &     (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UV)*UV)*UV
        QAS(I) = CZERO
       ELSE
        QX0(I) = (X0(I)-X0(2))/UV
        QAS(I) = (AS(I)-AS(2))/UV
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR3(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 +
     &   G(3)/5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
       QUR3(1) = -DCMPLX(ZERO,AS(2))+CIMAG*(CONE-U(1))*QAS(1)+ CONE -
     &            (U(1)+U(2))*R3(2)-U(1)*U(1)*QR3(1)
C
       US3(2) = (U(2)-CONE)/(U(2)+CONE)*UR3(2)
       QUS3(1)=( -US3(2)+UR3(2)+(U(1)-CONE)*QUR3(1))/(U(1)+CONE)
c
      AV = (U(1)+CONE)*(U(2)+CONE)
      SIM = AV*(UR3(2) + (U(1)-CONE)*QUR3(1))
      SIMI= AV*(US3(2) + (U(1)-CONE)*QUS3(1))
c
      RETURN
      END SUBROUTINE S2D3TWOI
C
C ----------------------------------------------------------------------
C file: s2d0leps.f
C date: 2011-04-02
C who:  S.Kaprzyk
C what: Return a proper CASE:
C       a) [w1<w2<ONE]; b) [w1<ONE<w2]; c) [ONE<w1<w2]
C INPUT:
C VERM(3) - three complex numbers
C N - integer number for W(I)=(VERM(I)-VERM(N))/(VERM(I)+VERM(N))
C EPS  - critical distance |w1-w2| etc.
C OUTPUT:
C LONE(3) [w1<w2<ONE]; [w1<ONE<w2]; [ONE<w1<w2]
C LEPS(1) [|w1-w2|<eps]
      SUBROUTINE  S2D0LEPS(
     &  VERM, N, W, AS,
     &  LONE, EPS, LEPS,
     &  iuerr
     & )
      IMPLICIT       NONE
      DOUBLE COMPLEX   VERM(3)
      INTEGER          N
      DOUBLE COMPLEX   W(2)
      DOUBLE PRECISION AS(2)
      DOUBLE PRECISION EPS
      LOGICAL          LONE(3), LEPS(1)
      INTEGER          iuerr
C
c      DOUBLE PRECISION DLAMCH
c      EXTERNAl         DLAMCH
C
      LOGICAL          LWONE(2), LV
      DOUBLE COMPLEX   AV
      DOUBLE PRECISION AL(2), AA, AIW, OZERO
      INTEGER          I, II, K
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE, PI
      DATA             ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/
      DATA             OZERO/ 1.0D-13/
c      OZERO = DLAMCH('e')*10
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      IF ((N.GT.3).OR.(N.LT.1)) THEN
       WRITE (iuerr,9010) N
 9010  FORMAT (' ***s2d0leps: N =',I6,' must be 1, or 2')
       STOP ' ***s2d0leps: '
      END IF
      II = 0
      DO 200 I = 1, 3
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
      DO 300 I = 1, 2
       DO 250 K = I, 2
        IF (LWONE(I).AND.LWONE(K).AND.(AL(I).LT.AL(K))) GO TO 250
        IF (.NOT.LWONE(I).AND..NOT.LWONE(K).AND.(AL(K).LT.AL(I)))
     &   GO TO 250
        IF (.NOT.LWONE(I).AND.LWONE(K).AND.(ONE.LT.AL(I)*AL(K)))
     &   GO TO 250
        IF (LWONE(I).AND..NOT.LWONE(K).AND.(AL(I)*AL(K).LT.ONE))
     &   GO TO 250
        AA = AL(K)
        AL(K) = AL(I)
        AL(I) = AA
        AA = AS(K)
        AS(K) = AS(I)
        AS(I) = AA
        AV = W(K)
        W(K) = W(I)
        W(I) = AV
        LV = LWONE(K)
        LWONE(K) = LWONE(I)
        LWONE(I) = LV
 250   CONTINUE
 300  CONTINUE
C
      LONE(1) = .FALSE.
      LONE(2) = .FALSE.
      LONE(3) = .FALSE.
      IF (LWONE(1).AND.LWONE(2)) THEN
       LONE(1) = .TRUE.
      END IF
      IF (LWONE(1).AND..NOT.LWONE(2)) THEN
       LONE(2) = .TRUE.
      END IF
      IF (.NOT.LWONE(1).AND..NOT.LWONE(2)) THEN
       LONE(3) = .TRUE.
      END IF
C here are cases, how close are w()
      LEPS(1) = .FALSE.
      IF ( (LONE(1).OR.LONE(3)).AND.
     &   (CDABS(W(1)-W(2)).LT.EPS) ) LEPS(1)=.TRUE.
c
      IF ((.NOT.LONE(1)).AND.(.NOT.LONE(2)).AND.(.NOT.LONE(3))) THEN
        WRITE (iuerr,9030) (LONE(I),I=1,3)
 9030   FORMAT ('***s2d0leps: LONE()=',3(L1,1X),'not in order')
        STOP '***s2d0leps: '
      END IF
C
      RETURN
      END SUBROUTINE S2D0LEPS
C s2d0leps()
