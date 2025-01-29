program test_tetra
  use defect
  use tetra_raja
  use fc3_interpolate,  ONLY : forceconst3
  use fc2_interpolate, only : mat2_diag
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use thutils, only: freq_in_grid
  use thtetra, only: tetra_init
  use constants, only: pi
  IMPLICIT NONE
  !
  type(ph_system_info) :: S
  type(forceconst2_grid) :: fc2
  type(q_grid) :: grid
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  !
  integer :: iw, i, nk3, iq, ibnd
  real(dp) :: omega, max_freq
  real(dp), allocatable :: freqs(:,:)
  complex(dp), allocatable :: wg(:,:)
  complex(dp), allocatable :: den(:,:), den1(:,:)
  INTEGER, ALLOCATABLE :: tetra(:,:), tetracount(:), tetramap(:,:,:)
  REAL(DP), ALLOCATABLE :: evals(:,:), tetraevals(:,:,:)

  CALL start_mpi()
  CALL READ_INPUT("TK", input, grid, S, fc2, fc3)

  ! #################################################################
  ! initialize my tetra and some variables
  ! #################################################################
  nk3 = grid%nqtot
  allocate(freqs(S%nat3, nk3))
  allocate(wg(S%nat3, nk3))
  allocate(den(S%nat3, input%n_omega))
  allocate(den1(S%nat3, input%n_omega))
  CALL freq_in_grid(S, fc2, grid, freqs)
  max_freq = maxval(freqs)
  CALL tetra_init(grid%n, S%bg, freqs)
  !
  ! #################################################################
  ! initialize elphbolt tetra
  ! #################################################################
  ALLOCATE(tetra(6*nk3, 4), tetracount(nk3), tetramap(2, nk3, 24))
  ALLOCATE(evals(nk3,S%nat3), tetraevals(6*nk3,S%nat3, 4))
  evals = TRANSPOSE(freqs)
  CALL form_tetrahedra_3d(nk3, grid%n, tetra, tetracount, tetramap)
  CALL fill_tetrahedra_3d(nk3, S%nat3, tetra, evals, tetraevals)
  !
  den1 = 0._dp
  do iw = 1, input%n_omega
    omega = max_freq*iw/input%n_omega
    wg = tetra_weights_green(omega)
    den(:,iw) = SUM(wg, dim=2)
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
    write(11,'(13E13.4)') max_freq*iw/input%n_omega, (den1(i,iw), i=1, S%nat3)
  enddo
  close(10)
  close(11)

  CALL stop_mpi()

end program
