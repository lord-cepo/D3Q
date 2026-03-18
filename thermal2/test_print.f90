module test_print
  use kinds, only: dp
  !
  interface allclose
    module procedure allclose_real_1d
    module procedure allclose_real_2d
    module procedure allclose_cplx_1d
    module procedure allclose_cplx_2d
  end interface allclose
  !
contains
  subroutine rnd_vec_cmplx(array, seed)
    implicit none

    ! Input/Output variables
    complex(dp), intent(out) :: array(:) ! Output: Complex array to initialize
    integer, intent(in), optional :: seed ! Optional: Seed for the random number generator

    ! Local variables
    integer :: i
    integer :: nrows
    real(dp) :: real_part, imag_part

    ! Get dimensions of the array
    nrows = size(array, 1)

    ! Initialize the random number generator with an optional seed
    if (present(seed)) then
      call random_seed(put=reshape([(seed + i, i = 0, 11)], shape=[12]))
    else
      call random_seed() ! Default seed
    end if

    ! Populate the array with random complex numbers
    do i = 1, nrows
      call random_number(real_part) ! Generate random number for the real part
      call random_number(imag_part) ! Generate random number for the imaginary part
      array(i) = cmplx(real_part, imag_part, dp)
    end do

  end subroutine
  !
  subroutine rnd_vec_real(array, seed)
    implicit none

    ! Input/Output variables
    real(dp), intent(out) :: array(:) ! Output: Complex array to initialize
    integer, intent(in), optional :: seed ! Optional: Seed for the random number generator

    ! Local variables
    integer :: i
    integer :: nrows
    real(dp) :: real_part

    ! Get dimensions of the array
    nrows = size(array, 1)

    ! Initialize the random number generator with an optional seed
    if (present(seed)) then
      call random_seed(put=reshape([(seed + i, i = 0, 11)], shape=[12]))
    else
      call random_seed() ! Default seed
    end if

    ! Populate the array with random complex numbers
    do i = 1, nrows
      call random_number(real_part) ! Generate random number for the real part
      array(i) = real_part
    end do

  end subroutine
  !
  subroutine rnd_mat_cmplx(array, seed)
    implicit none

    ! Input/Output variables
    complex(dp), intent(out) :: array(:,:) ! Output: Complex array to initialize
    integer, intent(in), optional :: seed ! Optional: Seed for the random number generator

    ! Local variables
    integer :: i, j
    integer :: nrows, ncols
    real(dp) :: real_part, imag_part

    ! Get dimensions of the array
    nrows = size(array, 1)
    ncols = size(array, 2)

    ! Initialize the random number generator with an optional seed
    if (present(seed)) then
      call random_seed(put=reshape([(seed + i, i = 0, 11)], shape=[12]))
    else
      call random_seed() ! Default seed
    end if

    ! Populate the array with random complex numbers
    do i = 1, nrows
      do j = 1, ncols
        call random_number(real_part) ! Generate random number for the real part
        call random_number(imag_part) ! Generate random number for the imaginary part
        array(i,j) = cmplx(real_part, imag_part, dp)
      end do
    end do

  end subroutine
  !
  subroutine rnd_mat_real(array, seed)
    implicit none

    ! Input/Output variables
    real(dp), intent(out) :: array(:,:) ! Output: Complex array to initialize
    integer, intent(in), optional :: seed ! Optional: Seed for the random number generator

    ! Local variables
    integer :: i, j
    integer :: nrows, ncols
    real(dp) :: real_part

    ! Get dimensions of the array
    nrows = size(array, 1)
    ncols = size(array, 2)

    ! Initialize the random number generator with an optional seed
    if (present(seed)) then
      call random_seed(put=reshape([(seed + i, i = 0, 11)], shape=[12]))
    else
      call random_seed() ! Default seed
    end if

    ! Populate the array with random complex numbers
    do i = 1, nrows
      do j = 1, ncols
        call random_number(real_part) ! Generate random number for the real part
        array(i,j) = real_part
      end do
    end do

  end subroutine
  !
  logical function allclose_real_1d(a, b, rtol)
    real(dp), intent(in) :: a(:), b(:)
    real(dp), intent(in) :: rtol
    allclose_real_1d = all(abs(a - b) <= rtol * max(abs(a), abs(b)))
  end function allclose_real_1d

  logical function allclose_real_2d(a, b, rtol)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(in) :: rtol
    allclose_real_2d = all(abs(a - b) <= rtol * max(abs(a), abs(b)))
  end function allclose_real_2d

  logical function allclose_cplx_1d(a, b, rtol)
    complex(dp), intent(in) :: a(:), b(:)
    real(dp), intent(in) :: rtol
    allclose_cplx_1d = all(abs(a - b) <= rtol * max(abs(a), abs(b)))
  end function allclose_cplx_1d

  logical function allclose_cplx_2d(a, b, rtol)
    complex(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(in) :: rtol
    allclose_cplx_2d = all(abs(a - b) <= rtol * max(abs(a), abs(b)))
  end function allclose_cplx_2d
  !
  subroutine orthonormalize(vectors)
    complex(dp), intent(inout) :: vectors(:,:)
    integer :: i, j, n
    complex(dp) :: proj
    real(dp) :: norm_vec

    n = size(vectors, 2)

    do i = 1, n
      ! Subtract projections onto previous vectors
      do j = 1, i-1
        proj = dot_product(vectors(:,j), vectors(:,i))
        vectors(:,i) = vectors(:,i) - proj * vectors(:,j)
      end do
      ! Normalize the vector
      norm_vec = SQRT(SUM(ABS(vectors(:,i))**2))
      if (norm_vec > 0.0_dp) then
        vectors(:,i) = vectors(:,i) / norm_vec
      else
        ! Handle the case where vector is zero after projections
        vectors(:,i) = 0.0_dp
      end if
    end do

  end subroutine orthonormalize
  !
  subroutine test_log(message, test_passed)
    character(len=*), intent(in) :: message
    logical, intent(in) :: test_passed
    !
    open(10, file='../../test.log', action='write', position='append')
    if(test_passed) then
      write(10,*) "Test passed: ", trim(message)
    else
      write(10,*) "----------------> Test FAILED: ", trim(message)
    end if
    close(10)
  end subroutine test_log
  !
end module
