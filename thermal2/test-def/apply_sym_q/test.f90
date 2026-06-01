program test_apply_sym_q
  use kinds, only : dp
  use symm_q_mat, only : apply_sym_q
  use test_print, only : test_log
  use input_fc, only : ph_system_info
  implicit none
  !
  integer, parameter :: nat = 1
  integer, parameter :: nmat = 3
  real(dp), parameter :: tol = 1.0e-10_dp
  real(dp) :: at(3,3), bg(3,3), tau(3,nat)
  real(dp) :: q_gamma(3), q_axis(3)
  integer :: ityp(nat)
  complex(dp) :: mats(3,3,nmat), mat(3,3), mat_twice(3,3)
  logical :: test_passed
  integer :: imat
  type(ph_system_info) :: S
  !
  call setup_simple_cubic(at, bg, tau, ityp)
  call setup_matrices(mats)
  q_gamma = [0.0_dp, 0.0_dp, 0.0_dp]
  q_axis = [0.0_dp, 0.0_dp, 0.25_dp]
  test_passed = .true.
  S%at = at
  S%bg = bg
  S%nat = nat
  S%ityp = ityp
  S%tau = tau
  !
  do imat = 1, nmat
    mat = mats(:,:,imat)
    call apply_sym_q(S, q_gamma, mat)
    test_passed = test_passed .and. is_isotropic(mat, tol)
    !
    mat_twice = mat
    call apply_sym_q(S, q_gamma, mat_twice)
    test_passed = test_passed .and. is_close(mat, mat_twice, tol)
    !
    mat = mats(:,:,imat)
    call apply_sym_q(S, q_axis, mat)
    test_passed = test_passed .and. is_axis_symmetric(mat, tol)
    !
    mat_twice = mat
    call apply_sym_q(S, q_axis, mat_twice)
    test_passed = test_passed .and. is_close(mat, mat_twice, tol)
  enddo
  !
  call test_log("apply_sym_q", test_passed)
  if(.not. test_passed) error stop "apply_sym_q test failed"
  write(*, '(a)') "apply_sym_q: passed"
  !
contains
  !
  subroutine setup_simple_cubic(at, bg, tau, ityp)
    real(dp), intent(out) :: at(3,3), bg(3,3), tau(3,nat)
    integer, intent(out) :: ityp(nat)
    integer :: i
    !
    at = 0.0_dp
    bg = 0.0_dp
    do i = 1, 3
      at(i,i) = 1.0_dp
      bg(i,i) = 1.0_dp
    enddo
    tau(:,1) = 0.0_dp
    ityp(1) = 1
  end subroutine
  !
  subroutine setup_matrices(mats)
    complex(dp), intent(out) :: mats(3,3,nmat)
    !
    mats(:,:,1) = reshape([ &
      cmplx( 1.0_dp,  0.2_dp, dp), cmplx( 0.4_dp, -0.7_dp, dp), cmplx(-0.3_dp,  0.5_dp, dp), &
      cmplx(-1.2_dp,  0.1_dp, dp), cmplx( 2.5_dp, -0.4_dp, dp), cmplx( 0.8_dp,  0.9_dp, dp), &
      cmplx( 0.6_dp, -0.2_dp, dp), cmplx(-0.9_dp,  0.3_dp, dp), cmplx(-1.0_dp,  0.6_dp, dp) ], [3,3])
    !
    mats(:,:,2) = reshape([ &
      cmplx( 0.1_dp, -0.8_dp, dp), cmplx( 2.0_dp,  0.0_dp, dp), cmplx( 0.7_dp, -1.1_dp, dp), &
      cmplx( 2.0_dp,  0.0_dp, dp), cmplx(-0.5_dp,  0.4_dp, dp), cmplx(-0.2_dp,  0.3_dp, dp), &
      cmplx( 0.7_dp,  1.1_dp, dp), cmplx(-0.2_dp, -0.3_dp, dp), cmplx( 1.4_dp, -0.2_dp, dp) ], [3,3])
    !
    mats(:,:,3) = reshape([ &
      cmplx( 3.0_dp,  0.0_dp, dp), cmplx( 0.0_dp,  1.0_dp, dp), cmplx( 0.0_dp,  0.0_dp, dp), &
      cmplx( 0.0_dp, -2.0_dp, dp), cmplx(-1.0_dp,  0.0_dp, dp), cmplx( 4.0_dp,  0.0_dp, dp), &
      cmplx( 5.0_dp,  0.0_dp, dp), cmplx( 0.0_dp,  0.5_dp, dp), cmplx( 2.0_dp,  0.0_dp, dp) ], [3,3])
  end subroutine
  !
  logical function is_isotropic(mat, tol)
    complex(dp), intent(in) :: mat(3,3)
    real(dp), intent(in) :: tol
    complex(dp) :: expected(3,3), tr3
    !
    expected = cmplx(0.0_dp, 0.0_dp, dp)
    tr3 = (mat(1,1) + mat(2,2) + mat(3,3)) / 3.0_dp
    expected(1,1) = tr3
    expected(2,2) = tr3
    expected(3,3) = tr3
    is_isotropic = maxval(abs(mat - expected)) < tol
  end function
  !
  logical function is_axis_symmetric(mat, tol)
    complex(dp), intent(in) :: mat(3,3)
    real(dp), intent(in) :: tol
    complex(dp) :: expected(3,3)
    !
    expected = cmplx(0.0_dp, 0.0_dp, dp)
    expected(1,1) = 0.5_dp * (mat(1,1) + mat(2,2))
    expected(2,2) = expected(1,1)
    expected(3,3) = mat(3,3)
    is_axis_symmetric = maxval(abs(mat - expected)) < tol
  end function
  !
  logical function is_close(a, b, tol)
    complex(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(in) :: tol
    !
    is_close = maxval(abs(a - b)) < tol
  end function
  !
end program test_apply_sym_q
