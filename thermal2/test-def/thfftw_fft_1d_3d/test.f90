program test_fft
  use kinds, only : dp
  use thfftw, only : fft_1d_3d
  use thutils, only : index2v_n
  use test_print, only : test_log
  !
  integer, parameter :: t = 2
  integer, parameter :: n(3) = [2,2,3]
  complex(dp) :: Aq(t,t,product(n))
  complex(dp) :: Aq1(t,t,product(n))
  complex(dp) :: Ar(t,t,product(n))
  integer :: iR
  integer :: R(product(n))
  integer :: Rn(n(1),n(2),n(3))
  !
  Aq = 0._dp
  Aq(1,1,5) = (1.0_dp, 0.0_dp)
  !
  call fft_1d_3d(t, n, Aq, Ar, 1)
  ! do ir = 1, product(n)
  !   R(ir) = iR
  !   ! print"(3I2,2F10.2)", index2v_n(ir,n), Ar(1,1,iR)
  ! enddo
  ! rn = reshape(R, n)
  ! print*, rn(1,1,1)
  call fft_1d_3d(t, n, Ar, Aq1, -1)
  !
  Aq1 = Aq1 / product(n)
  call test_log("thfftw_fft_1d_3d", all(abs(Aq - Aq1) <= 1e-14_dp))
end program

