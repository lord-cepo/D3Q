!
! This module is rewritten from the tetra.f90 in PW/src
!
MODULE thtetra
   !
   ! Tetrahedron method, linear and optimized. opt is better for all purpose
   ! weights for delta integration are obtained by multiplying gi with Iik in
   ! https://iopscience.iop.org/article/10.1088/0022-3719/12/15/008
   !
   ! another useful link, for delta integration where the integrand is 1 (DOS)
   ! http://staff.ustc.edu.cn/~zqj/posts/LinearTetrahedronMethod/#fn:tet_weight
   !
   ! optimized tetrahedron method, used in QE for theta integration in opt_tetra_weights
   ! https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.094515
   ! they multiply ni with Jik in the first article, then they transform (fit) through wlsm matrices
   USE kinds, ONLY: DP
   USE mpi_thermal, ONLY: my_id, num_procs
   !
   IMPLICIT NONE
   !
   PRIVATE
   SAVE
   !
   INTEGER :: ntetra = 0
   !! number of tetrahedra
   INTEGER :: nntetra
   !! k-points per tetrahedron used to compute weights.
   !! 4 for linear / 20 for optimized tetrahedron method
   INTEGER, ALLOCATABLE :: tetra(:,:)
   !! index of k-points in a given tetrahedron shape (nntetra,ntetra)
   REAL(DP), ALLOCATABLE :: wlsm(:,:)
   !! Weights for the optimized tetrahedron method
   LOGICAL :: opt_flag
   LOGICAL :: is_mpi_flag
   !
   PUBLIC :: tetra, ntetra, nntetra, wlsm
   PUBLIC :: tetra_init, tetra_weights_delta, deallocate_tetra, tetra_weights_theta
   !
CONTAINS
   !
   !--------------------------------------------------------------------------
   SUBROUTINE tetra_init(nq, bg, opt, is_mpi)
      !-----------------------------------------------------------------------------
      !! This rouotine sets the corners and additional points for each tetrahedron.
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nq(3)
      !! number of q-points in each direction
      !
      REAL(DP), INTENT(IN) :: bg(3,3)
      !! Reciplocal lattice vectors [2 pi / a]
      !
      LOGICAL, INTENT(IN) :: opt
      !! if .true., uses opt_tetra methods
      LOGICAL, INTENT(IN) :: is_mpi
      !! if .true., the grid is scattered
      REAL(DP), PARAMETER :: eps = 1e-5_dp
      !
      INTEGER :: i1, i2, i3, itet, itettot, ii, ik,  &
         ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3)
      !
      REAL(DP) :: l(4), bvec2(3,3), bvec3(3,4)


      IF(ntetra /= 0) RETURN
      opt_flag = opt
      is_mpi_flag = is_mpi
      !
      ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
      !
      bvec2(1:3,1) = bg(1:3,1) / REAL(nq(1), dp)
      bvec2(1:3,2) = bg(1:3,2) / REAL(nq(2), dp)
      bvec2(1:3,3) = bg(1:3,3) / REAL(nq(3), dp)
      !
      bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
      bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
      bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
      bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
      !
      DO ii = 1, 4
         l(ii) = DOT_PRODUCT(bvec3(1:3, ii), bvec3(1:3, ii))
      ENDDO
      !
      ii = MINLOC(l(1:4),1)
      !
      ivvec0(1:4) = (/ 0, 0, 0, 0 /)
      !
      divvec(1:4,1) = (/ 1, 0, 0, 0 /)
      divvec(1:4,2) = (/ 0, 1, 0, 0 /)
      divvec(1:4,3) = (/ 0, 0, 1, 0 /)
      divvec(1:4,4) = (/ 0, 0, 0, 1 /)
      !
      ivvec0(ii) = 1
      divvec(ii, ii) = - 1
      !
      ! Divide a subcell into 6 tetrahedra
      !
      itet = 0
      DO i1 = 1, 3
         DO i2 = 1, 3
            IF(i2 == i1) CYCLE
            DO i3 = 1, 3
               IF(i3 == i1 .OR. i3 == i2) CYCLE
               !
               itet = itet + 1
               !
               ivvec(1:3,1,itet) = ivvec0(1:3)
               ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
               ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
               ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
               !
            ENDDO
         ENDDO
      ENDDO
      !
      ! Additional points surrounding the tetrahedron
      !
      ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
      ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
      ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
      ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
      !
      ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
      ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
      ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
      ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
      !
      ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
      ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
      ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
      ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
      !
      ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
      ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
      ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
      ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
      !
      ! Set the weight for the each tetrahedron method
      !
      ntetra  = 6*nq(1)*nq(2)*nq(3)
      ! WRITE(stdout,*) "    [opt_tetra]  Optimized tetrahedron method is used."
      !
      IF (opt) THEN
         !
         nntetra = 20
         IF (.NOT. ALLOCATED(tetra)) ALLOCATE( tetra(nntetra,ntetra) )
         IF (.NOT. ALLOCATED(wlsm))  ALLOCATE( wlsm(4,nntetra) )
         !
         wlsm(1, 1: 4) = REAL((/1440,    0,   30,    0/), dp)
         wlsm(2, 1: 4) = REAL((/   0, 1440,    0,   30/), dp)
         wlsm(3, 1: 4) = REAL((/  30,    0, 1440,    0/), dp)
         wlsm(4, 1: 4) = REAL((/   0,   30,    0, 1440/), dp)
         !
         wlsm(1, 5: 8) = REAL((/ -38,    7,   17,  -28/), dp)
         wlsm(2, 5: 8) = REAL((/ -28,  -38,    7,   17/), dp)
         wlsm(3, 5: 8) = REAL((/  17,  -28,  -38,    7/), dp)
         wlsm(4, 5: 8) = REAL((/   7,   17,  -28,  -38/), dp)
         !
         wlsm(1, 9:12) = REAL((/ -56,    9,  -46,    9/), dp)
         wlsm(2, 9:12) = REAL((/   9,  -56,    9,  -46/), dp)
         wlsm(3, 9:12) = REAL((/ -46,    9,  -56,    9/), dp)
         wlsm(4, 9:12) = REAL((/   9,  -46,    9,  -56/), dp)
         !
         wlsm(1,13:16) = REAL((/ -38,  -28,   17,    7/), dp)
         wlsm(2,13:16) = REAL((/   7,  -38,  -28,   17/), dp)
         wlsm(3,13:16) = REAL((/  17,    7,  -38,  -28/), dp)
         wlsm(4,13:16) = REAL((/ -28,   17,    7,  -38/), dp)
         !
         wlsm(1,17:20) = REAL((/ -18,  -18,   12,  -18/), dp)
         wlsm(2,17:20) = REAL((/ -18,  -18,  -18,   12/), dp)
         wlsm(3,17:20) = REAL((/  12,  -18,  -18,  -18/), dp)
         wlsm(4,17:20) = REAL((/ -18,   12,  -18,  -18/), dp)
         !
         wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260.0_dp
         !
      ELSE
         !
         nntetra = 4
         IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra(nntetra,ntetra) )
         IF(.NOT. ALLOCATED(wlsm))  ALLOCATE ( wlsm(4,nntetra) )
         wlsm(:,:) = 0.0_dp
         !
         wlsm(1,1) = 1.0_dp
         wlsm(2,2) = 1.0_dp
         wlsm(3,3) = 1.0_dp
         wlsm(4,4) = 1.0_dp
         !
      ENDIF

      !
      !  locate k-points of the uniform grid in the list of irreducible k-points
      !  that was previously calculated
      !
      !  bring irreducible k-points to crystal axis
      !
      ! Construct tetrahedra
      !
      itettot = 0
      DO i1 = 1, nq(1)
         DO i2 = 1, nq(2)
            DO i3 = 1, nq(3)
               !
               DO itet = 1, 6
                  !
                  itettot = itettot + 1
                  !
                  DO ii = 1, nntetra
                     !
                     ikv(1:3) = (/i1, i2, i3/) - 1
                     ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
                     ikv(1:3) = MODULO(ikv(1:3), (/nq(1), nq(2), nq(3)/))
                     !
                     ik = ikv(3) + nq(3) * (ikv(2) + nq(2) * ikv(1)) + 1
                     !
                     tetra(ii, itettot) = ik
                     !
                  ENDDO ! ii
                  !
               ENDDO ! itet
               !
            ENDDO ! i3
         ENDDO ! i2
      ENDDO ! i1
      !
   END SUBROUTINE tetra_init
   !

   !--------------------------------------------------------------------------------------
   FUNCTION tetra_weights_delta( nqs, nbnd, et, ef, print) RESULT(wg)
      !-----------------------------------------------------------------------------------
      !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
      !! The resulting wg can be used as sum(Ak * wk)
      !-----------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nqs
      !! The total # of k in irr-BZ
      INTEGER, INTENT(IN) :: nbnd
      !! The # of bands
      REAL(DP), INTENT(IN) :: et(nbnd,nqs)
      !! Kohn Sham energy [Ry]
      REAL(DP) :: wg(nbnd,nqs)
      !! Integration weight of each k
      REAL(DP), INTENT(IN) :: ef
      !! The Fermi energy
      LOGICAL, INTENT(IN) :: print
      !! prints the 1->2 surface scattering

      ! ... local variables
      !
      INTEGER :: ik, nt, ibnd, i, ii, itetra(4), my_id_, num_procs_
      REAL(DP) :: e(4), wg0(4), C, a(4,4), print_surface(nqs)

      EXTERNAL hpsort
      !
      wg = 0._dp
      !
      IF(is_mpi_flag) THEN
         my_id_ = my_id
         num_procs_ = num_procs
      ELSE
         my_id_ = 0
         num_procs_ = 1
      ENDIF
      DO nt = 1+my_id_, ntetra, num_procs_
         !
         DO ibnd = 1, nbnd
            !
            e(1:4) = 0.0_dp
            DO ii = 1, nntetra
               !
               ik = tetra(ii, nt)
               IF(opt_flag) THEN
                  e(1:4) = e(1:4) + wlsm(1:4,ii) * et(ibnd,ik)
               ELSE
                  e(ii) = et(ibnd,ik)
               ENDIF
               !
            ENDDO
            !
            itetra(1) = 0
            CALL hpsort( 4, e, itetra )
            !
            IF( ef < e(4) .AND. ef > e(1) ) THEN
               IF(ibnd == 2 .and. print) THEN
                  DO ii = 1,4
                     print_surface(tetra(ii,nt)) = .true.
                  ENDDO
               ENDIF
               DO ii = 1, 4
                  DO i = 1, 4
                     IF ( ABS(e(i)-e(ii)) < 1.d-12 ) THEN
                        a(ii,i) = 0.0_dp
                     ELSE
                        a(ii,i) = ( ef - e(i) ) / (e(ii) - e(i) )
                     END IF
                  ENDDO
               ENDDO
               !
               IF( e(1) < ef .AND. ef < e(2) ) THEN
                  !
                  C = a(2,1) * a(3,1) * a(4,1) / (ef - e(1))
                  wg0(1) = a(1,2) + a(1,3) + a(1,4)
                  wg0(2:4) = a(2:4,1)

                  wg0 = wg0 * C
                  !
               ELSEIF( e(2) <= ef .AND. ef < e(3)) THEN
                  !
                  C = a(2,3) * a(3,1) + a(3,2) * a(2,4)
                  !
                  wg0(1) = a(1,4) * C + a(1,3) * a(3,1) * a(2,3)
                  wg0(2) = a(2,3) * C + a(2,4)**2 * a(3,2)
                  wg0(3) = a(3,2) * C + a(3,1)**2 * a(2,3)
                  wg0(4) = a(4,1) * C + a(4,2) * a(2,4) * a(3,2)

                  wg0 = wg0 / (e(4) - e(1))
                  !
               ELSEIF ( e(3) <= ef .AND. ef < e(4)) THEN
                  !
                  C = a(1,4) * a(2,4) * a(3,4) / (e(4) - ef)
                  !
                  wg0(1:3) = a(1:3,4)
                  wg0(4) = a(4,1) + a(4,2) + a(4,3)
                  !
                  wg0 = wg0 * C
                  !
               ENDIF
               !
               ! wg0(1:4) = wg0(1:4) / REAL(ntetra, dp)
               !
               DO ii = 1, nntetra
                  !
                  ik = tetra(ii, nt)
                  IF(opt_flag) THEN
                     wg(ibnd,ik) = wg(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(1:4),ii), wg0(1:4))
                  ELSE
                     wg(ibnd,ik) = wg(ibnd,ik) + wg0(ii)
                  ENDIF
               ENDDO
            ENDIF
            !
         ENDDO ! ibnd
         !
      ENDDO ! nt
      wg = wg / REAL(ntetra, dp)
      !
      !
      ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
      open(unit=10, file='surf.txt', status='replace', action='write')
      do i = 1, nqs
         write(10, *) print_surface(i)
      end do
      close(10)


   END FUNCTION tetra_weights_delta
   !

   !
   !--------------------------------------------------------------------
   SUBROUTINE deallocate_tetra( )
      !--------------------------------------------------------
      !! Deallocate tetra and wlsm
      !
      ntetra = 0
      nntetra = 0
      IF ( ALLOCATED(tetra) ) DEALLOCATE (tetra)
      IF ( ALLOCATED(wlsm ) ) DEALLOCATE (wlsm )
      !
   END SUBROUTINE deallocate_tetra
   !
   !
   FUNCTION tetra_weights_theta( nks, nbnd, et, ef) RESULT(wg)
      !--------------------------------------------------------------------
      !! Calculates weights with the tetrahedron method (P.E.Bloechl).
      !! Fermi energy has to be calculated in previous step.
      !! Generalization to noncollinear case courtesy of Iurii Timrov.
      !! @Note (P. Delugas 8/10/2019) Needs to be called only after initializations,
      !!       stops the program with an error call otherwise.
      !!
      !
      USE kinds
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nks
      !! Total # of k in irreducible BZ
      INTEGER, INTENT(IN) :: nbnd
      !! number of bands
      REAL(DP), INTENT(IN) :: et(nbnd,nks)
      !! eigenvalues of the hamiltonian
      REAL(DP) :: wg(nbnd,nks)
      !! the weight of each k point and band
      ! wg must be (inout) and not (out) because if is/=0 only terms for
      ! spin=is are initialized; the remaining terms should be kept, not lost.
      REAL(DP), INTENT(IN) :: ef
      !! Fermi energy

      EXTERNAL hpsort
      !
      ! ... local variables
      !
      REAL(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra(4), dosef
      INTEGER :: ibnd, nt, nk, i, kp1, kp2, kp3, kp4, itetra(4)
      !
      nk = 0
      wg = 0._dp
      !
      DO nt = 1, ntetra
         DO ibnd = 1, nbnd
            !
            ! etetra are the energies at the vertexes of the nt-th tetrahedron
            !
            DO i = 1, 4
               etetra(i) = et (ibnd, tetra(i,nt) + nk)
            ENDDO
            itetra (1) = 0
            CALL hpsort( 4, etetra, itetra )
            !
            ! ...sort in ascending order: e1 < e2 < e3 < e4
            !
            e1 = etetra(1)
            e2 = etetra(2)
            e3 = etetra(3)
            e4 = etetra(4)
            !
            ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
            !
            kp1 = tetra(itetra(1), nt) + nk
            kp2 = tetra(itetra(2), nt) + nk
            kp3 = tetra(itetra(3), nt) + nk
            kp4 = tetra(itetra(4), nt) + nk
            !
            ! calculate weights wg
            !
            IF (ef>=e4) THEN
               !
               wg(ibnd, kp1) = wg(ibnd, kp1) + 0.25d0 / ntetra
               wg(ibnd, kp2) = wg(ibnd, kp2) + 0.25d0 / ntetra
               wg(ibnd, kp3) = wg(ibnd, kp3) + 0.25d0 / ntetra
               wg(ibnd, kp4) = wg(ibnd, kp4) + 0.25d0 / ntetra
               !
            ELSEIF (ef<e4 .AND. ef>=e3) THEN
               !
               c4 = 0.25d0 / ntetra * (e4 - ef)**3 / (e4 - e1) / (e4 - e2) &
                  / (e4 - e3)
               dosef = 3.d0 / ntetra * (e4 - ef)**2 / (e4 - e1) / (e4 - e2) &
                  / (e4 - e3)
               wg(ibnd,kp1) = wg(ibnd,kp1) + 0.25d0 / ntetra - c4 * &
                  (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                  (ibnd, kp1) ) / 40.d0
               wg(ibnd,kp2) = wg(ibnd,kp2) + 0.25d0 / ntetra - c4 * &
                  (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                  (ibnd, kp2) ) / 40.d0
               wg(ibnd,kp3) = wg(ibnd,kp3) + 0.25d0 / ntetra - c4 * &
                  (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                  (ibnd, kp3) ) / 40.d0
               wg(ibnd,kp4) = wg(ibnd,kp4) + 0.25d0 / ntetra - c4 * &
                  (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
                  + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
                  et(ibnd,kp4) ) / 40.d0
               !
            ELSEIF (ef<e3 .AND. ef>=e2) THEN
               !
               c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
               c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
                  / (e4 - e1) / (e3 - e2) / (e3 - e1)
               c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
                  / (e3 - e2) / (e4 - e1)
               dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
                  (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
                  * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
               wg(ibnd, kp1) = wg(ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
                  / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
                  (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
               wg(ibnd, kp2) = wg(ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
                  * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
                  (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
               wg(ibnd, kp3) = wg(ibnd, kp3) + (c1 + c2) * (ef - e1) &
                  / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
                  (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
               wg(ibnd, kp4) = wg(ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
                  / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
                  e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
               !
            ELSEIF (ef<e2 .AND. ef>=e1) THEN
               !
               c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
                  / (e4 - e1)
               dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
                  / (e4 - e1)
               wg(ibnd, kp1) = wg(ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
                  * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
                  + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
               wg(ibnd, kp2) = wg(ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
                  + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
               wg(ibnd, kp3) = wg(ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
                  + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
               wg(ibnd, kp4) = wg(ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
                  + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0

               ! c4 = (ef-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)/4/ntetra
               ! wg(ibnd,kp1) = wg(ibnd,kp1) + (1 + (ef-e2)/(e1-e2) + (ef-e3)/(e1-e3) + (ef-e4)/(e1-e4))*c4
               ! wg(ibnd,kp2) = wg(ibnd,kp2) + (ef-e1)/(e2-e1)*c4
               ! wg(ibnd,kp3) = wg(ibnd,kp3) + (ef-e1)/(e3-e1)*c4
               ! wg(ibnd,kp4) = wg(ibnd,kp4) + (ef-e1)/(e4-e1)*c4
            ENDIF
            !
         ENDDO
      ENDDO
      !
   END FUNCTION tetra_weights_theta
END MODULE thtetra

