program test_broyden_c2
  use kinds, only : dp
  use EPW_utilities, only : mix_broyden_full
  implicit none

  integer, parameter :: ndim = 2
  integer, parameter :: n_iter = 6
  integer, parameter :: max_iter = 40
  real(dp), parameter :: alphamix = 0.7_dp
  real(dp), parameter :: fixed_point_step = 0.35_dp
  real(dp), parameter :: tol = 1.0e-10_dp

  complex(dp) :: z(ndim)
  complex(dp) :: gz(ndim)
  complex(dp) :: fz(ndim)
  complex(dp) :: z_root(ndim)
  complex(dp) :: rhs(ndim)
  complex(dp) :: df(ndim, n_iter)
  complex(dp) :: dv(ndim, n_iter)
  real(dp) :: residual
  integer :: iter
  integer :: unit
  logical :: converged

  z_root(1) = cmplx( 0.40_dp, 0.20_dp, dp)
  z_root(2) = cmplx(-0.30_dp, 0.50_dp, dp)
  rhs = coupled_polynomial(z_root)
  print*, rhs

  z(1) = cmplx(-0.20_dp, 0.10_dp, dp)
  z(2) = cmplx( 0.20_dp, 0.00_dp, dp)
  df = cmplx(0.0_dp, 0.0_dp, dp)
  dv = cmplx(0.0_dp, 0.0_dp, dp)
  converged = .false.

  open(newunit=unit, file='broyden_c2_path.dat', status='replace', action='write')
  write(unit, '(a)') '# iter Re(z1) Im(z1) Re(z2) Im(z2) Re(F1) Im(F1) Re(F2) Im(F2) ||F||'
  write(*, '(a)') ' iter        Re(z1)        Im(z1)        Re(z2)        Im(z2)          ||F||'

  do iter = 1, max_iter
    fz = coupled_polynomial(z) - rhs
    residual = norm2_complex(fz)

    write(unit, '(i5, 9(1x, es24.16))') iter, real(z(1), dp), aimag(z(1)), &
      real(z(2), dp), aimag(z(2)), real(fz(1), dp), aimag(fz(1)), &
      real(fz(2), dp), aimag(fz(2)), residual
    write(*, '(i5, 5(1x, es13.5))') iter, real(z(1), dp), aimag(z(1)), &
      real(z(2), dp), aimag(z(2)), residual

    if (residual < tol) then
      converged = .true.
      exit
    endif

    gz = z - fixed_point_step * fz
    call mix_broyden_full(ndim, gz, z, alphamix, iter, n_iter, df, dv)
  enddo

  close(unit)

  if (.not. converged) error stop 'Broyden C2 test did not converge'
  if (norm2_complex(z - z_root) > 1.0e-7_dp) error stop 'Broyden C2 test converged to the wrong root'

  write(*, '(a)') 'broyden_c2: passed'

contains

  function coupled_polynomial(x) result(f)
    complex(dp), intent(in) :: x(ndim)
    complex(dp) :: f(ndim)

    f(1) = x(1)**2 + 0.25_dp * x(1) * x(2) + x(2)
    f(2) = x(2)**2 - 0.15_dp * x(1) * x(2) + x(1)
  end function coupled_polynomial

  function norm2_complex(x) result(norm)
    complex(dp), intent(in) :: x(:)
    real(dp) :: norm

    norm = sqrt(sum(abs(x)**2))
  end function norm2_complex

end program test_broyden_c2
