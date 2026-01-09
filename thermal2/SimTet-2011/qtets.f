C file: qtets.f
C date: 2010-09-28
C who:  S.Kaprzyk
C what: This part of code splits the tetrahedron TET(3,4)
C       into small tetrahedrea and put its vertices into
C       QPTS(1..3, 1..NQPTS) - list. The mesh-points density
C       is controlled by the NMX-integer equal to the number
C       of points on each edge of input tetrahedron.
      SUBROUTINE QTETS
     & (
     & NMX, TET,
     & nqptsm, NQPTS, QPTS,
     & ntetsm, NTETS, N4TETS,
     & iuerr
     & )
      IMPLICIT   NONE
      INTEGER          NMX
      DOUBLE PRECISION TET(3,4)
      INTEGER          nqptsm, NQPTS
      DOUBLE PRECISION QPTS(3,nqptsm)
      INTEGER          ntetsm, NTETS, N4TETS(4,ntetsm)
      INTEGER          iuerr
C
      INTEGER  NTH(3,4,6), NABC
      INTEGER  I, J, K, L, M, IA(3)
      INTEGER  NA, NB, NBB, NC, NCC
      INTEGER  IPT, NMXPTS
      INTEGER  IW, IT, ITETS
      DOUBLE PRECISION DX, Q(3)
      DOUBLE PRECISION ONE
      DATA             ONE/1.0D0/
*-ASIS
      DATA (((NTH(I,J,K),I=1,3),J=1,4),K=1,3)/
     &  0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
     &  1, 0, 1,  1, 0, 0,  0, 1, 1,  0, 0, 1,
     &  0, 1, 0,  1, 0, 0,  0, 1, 1,  0, 0, 1/
*-ASIS
      DATA (((NTH(I,J,K),I=1,3),J=1,4),K=4,6)/
     &  0, 1, 0,  1, 0, 0,  0, 1, 1,  1, 1, 0,
     &  1, 0, 1,  1, 0, 0,  0, 1, 1,  1, 1, 0,
     &  1, 1, 1,  1, 0, 1,  0, 1, 1,  1, 1, 0/  
c
      NABC(NA,NB,NC,NMX)= NC + (3*(NB-1)*(2*NMX+4-2*NA-NB) +
     &  (NA-1)*(3*(NMX+1)*(NMX+2)-NA*(3*NMX+5-NA)))/6
C
      NMXPTS = NMX*(NMX+1)*(NMX+2)/6
      IF (NMXPTS.GT.nqptsm) THEN
       WRITE (iuerr,9010) NMX, NMXPTS, nqptsm
 9010  FORMAT (' ***qtets: NMX=',I3,' NMXPTS=',I8,' Increase nqptsm=',
     &  I8)
       STOP ' ***qtets '
      END IF
      NTETS = (NMX-1)**3
      IF (NTETS.GT.ntetsm) THEN
       WRITE (iuerr,9020) NMX, NTETS, ntetsm
 9020  FORMAT (' ***qtets: NMX=',I3,' NTETS=',I8,' Increase ntetsm=',I8)
       STOP ' ***qtets '
      END IF
c
      IPT = 0
      DX = ONE/(NMX-1)
      DO 100 NA = 1, NMX
       IA(1) = NA - 1
       NBB = NMX - NA + 1
       DO 50 NB = 1, NBB
        IA(2) = NB - 1
        NCC = NMX + 2 - NA - NB
        DO 40 NC = 1, NCC
        IA(3) = NC - 1
        DO 10 L = 1, 3
        Q(L) = TET(L,1)
        DO 5 M = 1, 3
        Q(L) = Q(L) + IA(M)*DX*(TET(L,M+1)-TET(L,1))
 5      CONTINUE
 10     CONTINUE
        IPT = IPT + 1
        IF (IPT.NE.NABC(NA,NB,NC,NMX)) THEN
        WRITE (*,9030) IPT, NABC(NA,NB,NC,NMX)
 9030   FORMAT ('***qtets:  IPT=',I8,' NABC=',I8)
        STOP ' ***qtets: must be equal '
        END IF
        DO 20 I = 1, 3
        QPTS(I,IPT) = Q(I)
 20     CONTINUE
 40     CONTINUE ! end DO ... NC = 1, NCC
 50    CONTINUE ! end DO ... NB = 1, NBB
 100  CONTINUE ! end DO ... NA = 1, NMX
C       NQPTS = IPT
      IF (NQPTS.NE.IPT) THEN
        WRITE (iuerr,9023) NMX, NQPTS, IPT
 9023   FORMAT (' ***qtets: NMX=',I3,' NQPTS=',I8,'must be  IPT=',I8)
       STOP ' ***qtets:  '
      END IF
C
      ITETS = 0
      DO 200 NA = 1, NMX
       NBB = NMX + 1 - NA
       DO 150 NB = 1, NBB
        NCC = NMX + 2 - NA - NB
        DO 120 NC = 1, NCC
        DO 110 IT = 1, 6
        IF ((IT.EQ.1).AND.(NA+NB+NC).GE.(NMX+2)) GO TO 120
        IF ((IT.GE.2.AND.IT.LE.5).AND.(NA+NB+NC).GE.(NMX+1)) GO TO 120
        IF ((IT.EQ.6).AND.(NA+NB+NC).GE.NMX) GO TO 120
        ITETS = ITETS + 1
        DO 105 IW = 1, 4
*-ASIS
          IPT = NABC( NA+NTH(1,IW,IT),
     &                NB+NTH(2,IW,IT),
     &                NC+NTH(3,IW,IT),NMX )
        N4TETS(IW,ITETS) = IPT
 105    CONTINUE
 110    CONTINUE
 120    CONTINUE  ! end DO ... NC = 1, NCC
 150   CONTINUE  ! end DO ... NB = 1, NBB
 200  CONTINUE  ! end DO ... NA = 1, NMX
c
      RETURN
      END SUBROUTINE QTETS
C qtets()
