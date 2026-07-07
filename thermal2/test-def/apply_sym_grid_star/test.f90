program test_apply_sym_grid_star
  use kinds, only : dp
  use input_fc, only : ph_system_info
  use q_grids, only : q_grid, q_grid_copy, setup_simple_grid
  use thtetra, only : equiv_grid
  use symm_q_mat, only : apply_sym, get_symmetry_q_star, apply_sym_q_star
  use thutils, only : cryst2cart
  use test_print, only : test_log
  !
  implicit none
  !
  integer, parameter :: mesh(3) = [4,4,4]
  real(dp), parameter :: tol = 1.e-10_dp
  !
  type(ph_system_info) :: S
  type(q_grid) :: full_grid, sym_grid
  integer, allocatable :: equiv(:)
  complex(dp), allocatable :: dyn_raw(:,:,:), dyn_apply(:,:,:)
  real(dp), allocatable :: q_star(:,:)
  complex(dp), allocatable :: dyn_star(:,:,:)
  complex(dp), allocatable :: dyn_qstar(:,:)
  real(dp) :: max_diff, diff
  integer :: iq, istar, iq_star
  logical :: passed
  !
  call setup_body_centered_two_atom(S)
  call setup_simple_grid(S%bg, mesh(1), mesh(2), mesh(3), full_grid)
  call q_grid_copy(full_grid, sym_grid)
  call sym_grid%symmetrize(S)
  call equiv_grid(sym_grid, S, equiv)
  !
  allocate(dyn_raw(S%nat3,S%nat3,full_grid%nqtot))
  allocate(dyn_qstar(S%nat3,S%nat3))
  call fill_test_matrices(full_grid, dyn_raw)
  dyn_apply = dyn_raw
  call apply_sym(S, dyn_apply, equiv, full_grid%xq, .true.)
  !
  max_diff = 0._dp
  do iq = 1, full_grid%nqtot
    call get_symmetry_q_star(S, full_grid%xq(:,iq), q_star)
    allocate(dyn_star(S%nat3,S%nat3,size(q_star,2)))
    do istar = 1, size(q_star,2)
      call find_grid_q(full_grid, S, q_star(:,istar), iq_star)
      dyn_star(:,:,istar) = dyn_raw(:,:,iq_star)
    enddo
    call apply_sym_q_star(S, full_grid%xq(:,iq), q_star, dyn_star, dyn_qstar)
    diff = maxval(abs(dyn_apply(:,:,iq) - dyn_qstar))
    max_diff = max(max_diff, diff)
    deallocate(q_star, dyn_star)
  enddo
  !
  passed = max_diff < tol
  print"(A,ES18.8)", "max |apply_sym - apply_sym_q_star| = ", max_diff
  call test_log("apply_sym_grid_star", passed)
  !
contains
  !
  subroutine setup_body_centered_two_atom(S)
    type(ph_system_info), intent(out) :: S
    !
    S%ntyp = 2
    S%nat = 2
    S%nat3 = 3*S%nat
    S%nat32 = S%nat3**2
    S%nat33 = S%nat3**3
    S%ibrav = 1
    S%symm_type = 'cubic'
    S%celldm = 0._dp
    S%celldm(1) = 1._dp
    S%at = 0._dp
    S%at(1,1) = 1._dp
    S%at(2,2) = 1._dp
    S%at(3,3) = 1._dp
    S%bg = S%at
    S%omega = 1._dp
    S%alat = 1._dp
    S%tpiba = 1._dp
    S%atm = ''
    S%atm(1) = 'X'
    S%atm(2) = 'Y'
    S%amass = 1._dp
    S%amass(2) = 2._dp
    S%amass_variance = 0._dp
    S%epsil = 0._dp
    S%lrigid = .false.
    S%ldrigid = .false.
    S%nopbc = .false.
    allocate(S%tau(3,S%nat), S%ityp(S%nat))
    S%tau = 0._dp
    S%tau(:,2) = [0.5_dp, 0.5_dp, 0.5_dp]
    S%ityp = [1,2]
  end subroutine setup_body_centered_two_atom
  !
  subroutine fill_test_matrices(grid, dyn)
    type(q_grid), intent(in) :: grid
    complex(dp), intent(out) :: dyn(:,:,:)
    !
    integer :: iq, i, j
    real(dp) :: qc(3), arg
    !
    do iq = 1, grid%nqtot
      qc = cryst2cart(grid%xq(:,iq), S%at, -1)
      do j = 1, size(dyn,2)
        do i = 1, size(dyn,1)
          arg = real(7*i + 11*j + 13*iq, dp)
          dyn(i,j,iq) = cmplx( &
            sin(arg + 0.37_dp*qc(1) + 0.19_dp*qc(2)**2 + 0.11_dp*qc(3)), &
            cos(0.23_dp*arg + 0.41_dp*qc(1)*qc(2) + 0.17_dp*qc(3)), &
            kind=dp)
        enddo
      enddo
    enddo
  end subroutine fill_test_matrices
  !
  subroutine find_grid_q(grid, S, q, iq_match)
    type(q_grid), intent(in) :: grid
    type(ph_system_info), intent(in) :: S
    real(dp), intent(in) :: q(3)
    integer, intent(out) :: iq_match
    !
    integer :: iq
    real(dp) :: diff(3)
    !
    iq_match = 0
    do iq = 1, grid%nqtot
      diff = cryst2cart(q - grid%xq(:,iq), S%at, -1)
      if (norm2(diff - nint(diff)) < 1.e-5_dp) then
        iq_match = iq
        exit
      endif
    enddo
    if (iq_match == 0) call errore("find_grid_q", "q-star point not found in full grid", 1)
  end subroutine find_grid_q
  !
end program test_apply_sym_grid_star
