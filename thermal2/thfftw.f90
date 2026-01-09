#if defined(__FFTW3)

module thfftw
  ! use kinds, only: dp
! #include "fftw3.f03"
  USE, intrinsic :: iso_c_binding
  USE fft_param
#include "fftw3.f03"
  ! complex(dp) :: aq(2, 2, 64)
  ! complex(dp) :: ar(2, 2, 64)
  ! ! Initialize input with zeros and one delta
  ! ! aq = cmplx(0.0_c_double, 0.0_c_double, kind=C_DOUBLE_COMPLEX)
  ! ! aq(1,1,1,1,1) = cmplx(1.0_c_double,0.0_c_double, kind=C_DOUBLE_COMPLEX)
  ! aq = 0._dp
  ! aq(1,1,1) = (1.0_dp, 0.0_dp)
  ! !
  ! call fft_1d_3d(2, [4,4,4], aq, ar)
  ! !
  ! print*, "FFT result at (1,1,1,1,1): ", ar

contains
  subroutine fft_1d_3d(t, n, arr, arr_out, sign)
    integer, intent(in) :: t, n(3)
    integer, intent(in) :: sign
    complex(dp), intent(in) :: arr(t,t,product(n))
    complex(dp) :: arr_q(t,t,n(1), n(2), n(3))
    complex(dp) :: arr_r(t,t,n(1), n(2), n(3))
    !
    complex(dp), intent(out) :: arr_out(t,t,product(n))
    !
    arr_q = reshape(arr, [t, t, n(1), n(2), n(3)])
    !
    call fft3(t, n(1), n(2), n(3), arr_q, arr_r, sign)
    !
    arr_out = reshape(arr_r, [t, t, product(n)])
  end subroutine
  !
  subroutine fft3(t, nq1, nq2, nq3, inarr, outarr, sign)
    integer, intent(in) :: t, nq1, nq2, nq3
    complex(C_DOUBLE_COMPLEX), intent(in) :: inarr(t, t, nq1, nq2, nq3)
    complex(C_DOUBLE_COMPLEX), intent(out) :: outarr(t, t, nq1, nq2, nq3)
    integer(C_INTPTR_T) :: plan
    integer(C_INT) :: rank, n(3), howmany
    integer(C_INT) :: istride, ostride, idist, odist
    integer(C_INT) :: inembed(3), onembed(3)
    integer :: vol
    integer(C_INT) :: sign
    real(C_DOUBLE) :: scale
    integer :: i



    ! Setup FFT parameters
    rank = 3
    n = [nq1, nq2, nq3]
    howmany = t * t
    inembed = n
    onembed = n
    vol = nq1*nq2*nq3
    istride = howmany
    ostride = howmany
    idist = 1
    odist = 1
    ! sign = FFTW_BACKWARD  ! inverse FFT

    ! Create the FFTW plan for many 3D FFTs batched over t*t
    plan = 0
    call dfftw_plan_many_dft(plan, rank, n, howmany, inarr, inembed, istride, idist, &
      outarr, onembed, ostride, odist, sign, FFTW_ESTIMATE)
    if (plan == 0) then
      print*, "Plan creation failed!"
      stop
    endif

    ! Execute the FFT
    call dfftw_execute_dft(plan, inarr, outarr)

    ! Normalize output since FFTW does not
    ! scale = 1.0_c_double / real(vol, C_DOUBLE)
    ! outarr = outarr * scale

    ! Destroy the plan once done
    call dfftw_destroy_plan(plan)
  end subroutine

end module
#endif


