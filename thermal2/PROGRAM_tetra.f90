program test_tetra
  use defect
  use tetra_raja
  use fc3_interpolate,  ONLY : forceconst3
  use fc2_interpolate, only : mat2_diag
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use thutils, only: freq_in_grid
  use thtetra
  use constants, only: pi
  use parameters, only : npk
  use ktetra, only : old_init => tetra_init, &
    old_weights => tetra_weights_only
  USE symm_base, ONLY : set_sym_bl, set_sym, nsym, &
    s_symm_base => s, time_reversal, t_rev
  USE cell_base,  ONLY : at, bg


  IMPLICIT NONE
  !
  type(ph_system_info) :: S
  type(forceconst2_grid) :: fc2
  type(q_grid) :: grid, ingrid, symmgrid
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  !
  integer :: iw, i, nk3, iq, ibnd, jq, jbnd, kbnd, ik, j,k
  real(dp) :: omega, max_freq
  real(dp), allocatable :: freqs(:,:), brakets(:,:), infreqs(:,:), wgI(:,:), braketsr(:,:)
  complex(dp), allocatable :: wg(:,:), U(:,:,:), Uin(:,:,:), Ur(:,:,:)
  complex(dp), allocatable :: den(:,:), den1(:,:), wg_sym(:,:)
  INTEGER, ALLOCATABLE :: tetra(:,:), tetracount(:), tetramap(:,:,:)
  REAL(DP), ALLOCATABLE :: evals(:,:), tetraevals(:,:,:)
  integer :: SYMM, index
  real(dp) :: MEAN_ENERGY, wg1, q(3), e(4), ef


  CALL start_mpi()
  CALL READ_INPUT("TK", input, grid, S, fc2, fc3)
  at = S%at
  bg = S%bg

  call set_sym_bl()

  ! #################################################################
  ! initialize my tetra and some variables
  ! #################################################################
  nk3 = grid%nqtot
  ! call grid%symmetrize(S)
  allocate(freqs(S%nat3, grid%nqtot))
  allocate(wg(S%nat3, grid%nqtot))
  allocate(den(S%nat3, input%n_omega))
  allocate(den1(S%nat3, input%n_omega))
  allocate(U(S%nat3, S%nat3, nk3))
  allocate(Ur(S%nat3, S%nat3, nk3))
  CALL freq_in_grid(S, fc2, grid, freqs, U)
  max_freq = maxval(freqs)
  print*, max_freq
  ! CALL tetra_init_sym(grid, S, freqs, .true.)
  CALL tetra_init(grid%n, S%bg, freqs, .true.)

  ! #################################################################
  ! initialize elphbolt tetra
  ! #################################################################
  ALLOCATE(tetra(6*nk3, 4), tetracount(nk3), tetramap(2, nk3, 24))
  ALLOCATE(evals(S%nat3,nk3), tetraevals(6*nk3,S%nat3, 4))
  ! evals = TRANSPOSE(freqs)
  index = 0
  do k = 1, grid%n(3)
    do j = 1, grid%n(2)
      do i = 1, grid%n(1)
        q(1) = i / REAL(grid%n(1), dp)
        q(2) = j / REAL(grid%n(2), dp)
        q(3) = k / REAL(grid%n(3), dp)
        index = index + 1
        call cryst_to_cart(1, q, S%bg, 1)
        call freq_phq_safe(q, S, fc2, evals(:,index), Ur(:, :, index))
      enddo
    enddo
  enddo
  ! evals = freqs
  ! Ur = U
  CALL form_tetrahedra_3d(nk3, grid%n, tetra, tetracount, tetramap)
  CALL fill_tetrahedra_3d(nk3, S%nat3, tetra, evals, tetraevals)
  !
  den1 = 0._dp
  do iw = 1, input%n_omega
    omega = max_freq*iw/input%n_omega
    wg = tetra_weights_green(omega)
    DO iq = 1, nk3
      do ibnd = 1, S%nat3
        den1(ibnd,iw) = den1(ibnd,iw) + &
          CMPLX(real_tetra(omega, iq, ibnd, grid%n, tetramap, tetracount, tetraevals), &
          - pi * delta_fn_tetra(omega, iq, ibnd, grid%n, tetramap, tetracount, tetraevals), dp)
      END DO
    end do
  enddo

  den = den / nk3
  den1 = den1 / nk3
  open(unit=10, file='tetra.dat', status='unknown')
  open(unit=11, file='tetra1.dat', status='unknown')
  do iw = 1, input%n_omega
      write(10,'(13E13.4)') max_freq*iw/input%n_omega, (den(i, iw), i=1, S%nat3)
      ! write(11,'(13E13.4)') max_freq*iw/input%n_omega, (den1(i,iw), i=1, S%nat3)
  enddo
  close(10)
  close(11)


  ! ! #################################################################
  ! ! check if <q|q'> respects symmetries
  ! ! #################################################################

  ! CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
  !   input%nk_in(2), input%nk_in(3),&
  !   ingrid, scatter=.false., xq0=input%xk0_in)
  ! allocate(wgI(S%nat3, grid%nq))

  ! ! call grid%copy(symmgrid)
  ! ! call symmgrid%symmetrize(S)

  ! ! call old_init(nsym, s_symm_base, time_reversal, t_rev, S%at, S%bg, npk, &
  ! ! 0,0,0,symmgrid%n(1), symmgrid%n(2), symmgrid%n(3), symmgrid%nq, symmgrid%xq)

  ! do SYMM = 0, 1
  !   if (SYMM == 1) then
  !     call ingrid%symmetrize(S)
  !   endif
  !   allocate(infreqs(S%nat3, ingrid%nq))
  !   allocate(Uin(S%nat3, S%nat3,ingrid%nq))
  !   allocate(brakets(S%nat3, ingrid%nq))
  !   allocate(braketsr(S%nat3, ingrid%nq))

  !   do jq = 1, ingrid%nq
  !     call freq_phq_safe(ingrid%xq(:,jq), S, fc2, infreqs(:,jq), Uin(:,:,jq))
  !   enddo

  !   DO ik = 1, ingrid%nq
  !     DO ibnd = 1, S%nat3
  !       !
  !       wg1 = infreqs(ibnd,ik)
  !       !
  !       DO jbnd = ibnd + 1, S%nat3
  !         !
  !         IF (ABS(infreqs(ibnd,ik) - infreqs(jbnd,ik)) < 1e-4_dp * MEAN_ENERGY) THEN
  !           wg1 = wg1 + infreqs(jbnd,ik)
  !         ELSE
  !           !
  !           DO kbnd = ibnd, jbnd - 1
  !             infreqs(kbnd,ik) = wg1 / REAL(jbnd - ibnd, dp)
  !           ENDDO
  !           !
  !           EXIT
  !         ENDIF
  !         !
  !       ENDDO
  !       !
  !     ENDDO
  !   ENDDO

  !   ! call cryst_to_cart(ingrid%nq, ingrid%xq, S%at, -1)
  !   ! do iq = 1, ingrid%nq
  !   !   print"(I4,3F10.3)", iq, ingrid%xq(:,iq)
  !   ! enddo
  !   ! call cryst_to_cart(ingrid%nq, ingrid%xq, S%bg, 1)
  !   brakets = 0._dp
  !   braketsr = 0._dp
  !   ! print"(6E13.3)", infreqs(:,2)
  !   ! print"(6E13.3)", infreqs(:,3)
  !   do jq = 1, ingrid%nq
  !     do jbnd = 1, S%nat3
  !       call tetra_weights_delta(infreqs(jbnd,jq), wgI)
  !       ! call opt_tetra_weights_only(nirr, 1, S%nat3, freqs, spread(0, 1, nirr)&
  !       ! ef, wg, 0, spread(0, 1, nirr))
  !       ! if (SYMM == 1) &
  !       !   call old_weights(grid%nq, 1, 0, spread(0, 1, grid%nq), S%nat3, &
  !       !   -1._dp, freqs, infreqs(jbnd,jq), wgI)
  !       ! if (jq == 2 .or. jq == 3) then
  !       !   print"(6F12.4)", wgI(:,4)
  !       ! endif
  !       do iq = 1, grid%nq
  !         do ibnd = 1, S%nat3
  !           brakets(jbnd,jq) =  brakets(jbnd,jq) + ABS(dot_product(U(:,ibnd,iq), Uin(:,jbnd,jq)))**2 * &
  !             wgI(ibnd, iq)
  !           ! wgI(ibnd,iq) = delta_fn_tetra(infreqs(jbnd,jq), iq, ibnd, grid%n, tetramap, tetracount, tetraevals)
  !           ! braketsr(jbnd,jq) =  braketsr(jbnd,jq) + ABS(dot_product(Ur(:,ibnd,iq), Uin(:,jbnd,jq)))**2 * &
  !           !   wgI(ibnd, iq)
  !         enddo
  !       enddo
  !       ! if (SYMM == 0 .and. (jq == 2 .or. jq == 3)) then
  !       !   print*, jq, wgI(5,10)
  !       ! endif
  !     enddo
  !   enddo


  !   if (SYMM == 1) then
  !     open(unit=12, file='brakets-sym.dat', status='unknown')
  !     open(unit=13, file='braketsr-sym.dat', status='unknown')
  !   else
  !     open(unit=12, file='brakets.dat', status='unknown')
  !     open(unit=13, file='braketsr.dat', status='unknown')
  !   end if
  !   ! do iq = 1, grid%nq
  !   !   write(12,*) iq, grid%xq(:,iq), wgI(:,iq)
  !   ! enddo
  !   do jq = 1, ingrid%nq
  !     do jbnd = 1, S%nat3
  !       write(12,*) infreqs(jbnd,jq), brakets(jbnd,jq)
  !       write(13,*) infreqs(jbnd,jq), braketsr(jbnd,jq)
  !     enddo
  !   enddo
  !   close(12)
  !   close(13)
  !   deallocate(infreqs, Uin, brakets, braketsr)
  ! enddo

  ! ! ef = 1.0_dp
  ! ! e(1) = 0.9_dp
  ! ! e(2) = 0.9_dp
  ! ! e(3) = 1.1_dp
  ! ! e(4) = 1.2_dp

  ! ! call rm_degen_vertices(ef, e)
  ! ! print"(4F16.8)", e
  CALL stop_mpi()

end program
