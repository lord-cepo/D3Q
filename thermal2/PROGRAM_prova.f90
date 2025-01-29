program prova
  use kinds, only: dp
  use functions, only: invzmat
  use thutils, only: v2index, e_iqr, braket
  !
  implicit none
  integer, parameter :: nat3 = 4, N = 3, N3 = N**3
  integer, parameter :: mesh(3) = [N,N,N]
  integer :: i, j, k, v(3), iq, jq, rn, mrn
  complex(dp) :: VKK(nat3*N3, nat3*N3), VK(nat3, nat3)
  real(dp) :: VRR(nat3*N3, nat3*N3)
  complex(dp) :: eR(nat3*N3), phases(N3), e(nat3)
  real(dp) :: lambda(3), rand, q(3,N3), R(3,N3), mlambda(3)
  complex(dp) :: full_braket, reduced_braket
  !
  do i = 0, N-1
    do j = 0, N-1
      do k = 0, N-1
        v = [i, j, k]
        R(:,v2index(v, mesh)) = real(v, dp)
        q(:,v2index(v, mesh)) = real(v, dp) / N
      enddo
    enddo
  enddo
  !
  call random_seed()
  call random_number(rand)
  rn = NINT(rand*N**3)
  lambda = q(:, rn)
  mlambda = - lambda
  do i = 1, 3
    if ( mlambda(i) < -1e-10 ) mlambda(i) = mlambda(i) + 1
  enddo
  mrn = v2index(NINT(mlambda*mesh), mesh)
  !
  VKK = 0._dp
  call rnd_mat_real(VRR)
  call rnd_vec_cmplx(e)
  do iq = 1, N3
    do jq = 1, N3
      do i = 1, N3
        eR((i-1)*nat3+1:i*nat3) = e * e_iqr(lambda, R(:,i))
        do j = 1, N3
          VKK((iq-1)*nat3+1:iq*nat3, (jq-1)*nat3+1:jq*nat3) = &
            VKK((iq-1)*nat3+1:iq*nat3, (jq-1)*nat3+1:jq*nat3) + &
            VRR((i-1)*nat3+1:i*nat3, (j-1)*nat3+1:j*nat3) * &
            e_iqr(-q(:,iq), R(:, i)) * e_iqr(-q(:,jq), R(:,j))
        enddo
      enddo
    enddo
  enddo
  !

  full_braket = braket(eR, CMPLX(VRR, 0.0_dp, dp))
  reduced_braket = braket(e, VKK( &
    (rn-1)*nat3+1:rn*nat3, (mrn-1)*nat3+1:mrn*nat3))
  if (ABS(full_braket - reduced_braket) < 1e-8) then
    print *, 'Test passed ----------------------'
  else
    print *, 'Test failed'
    print'(A,F7.2,1X,F7.2)', 'Full braket:    ', full_braket
    print'(A,F7.2,1X,F7.2)', 'Reduced braket: ', reduced_braket
  endif
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
  pure function outer_product(A) result(AA)
    COMPLEX(dp), intent(in) :: A(:)
    COMPLEX(dp), allocatable :: AA(:,:)
    COMPLEX(dp) :: r
    integer :: nA, i, j
    nA=size(A)
    allocate(AA(nA,nA))

    DO i = 1, nA
      DO j = i+1, nA
        r = A(i) * CONJG(A(j))
        AA(i,j) = r
        AA(j,i) = CONJG(r)
      END DO
      AA(i,i) = A(i) * CONJG(A(i))
    END DO
  end function
end program prova
