program prova
  use kinds, only: dp
  use functions, only: invzmat
  use thutils, only: v2index, e_iqr, braket, index2v, outer_product, outer_product2
  !
  implicit none
  integer, parameter :: nat3 = 2, N = 3, N3 = N**3
  integer, parameter :: mesh(3) = [N,N,N]
  integer :: i, j, k, v(3), iq, jq, rn, mrn, Rij(3), inner, na
  integer :: na1, na2, miq, mjq
  complex(dp), dimension(nat3, nat3, N3, N3) :: VKK, TRR, GRR, temp
  real(dp) :: VRR(nat3, nat3, N3, N3), mq(3)
  ! complex(dp) :: eR(nat3,N3), phases(N3)
  real(dp) :: lambda(3), rand, q(3,N3), R(3,N3), mlambda(3)
  complex(dp) :: full_braket, reduced_braket, flat_braket
  complex(dp) :: U(nat3, N3), UR_lambda(nat3, N3), UR(nat3, N3, N3)
  complex(dp), dimension(nat3*N3, nat3*N3) :: G_flat, T_flat
  complex(dp) :: U_lambda_flat(nat3*N3), UR_flat(nat3*N3, N3)
  complex(dp), dimension(nat3*N3, nat3*N3) :: V_flat
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
  if (rn == 0) rn = 1
  mrn = minus_ind(q(:,rn), mesh)
  !
  !> builds VRR
  ! do i = 1, N3
  !   do j = 1, N3!_q grid for UR
  !     call rnd_mat_real(VRR(:,:,i,j))
  !   enddo
  ! enddo

  !> symmetrize VRR
  do na1 = 1, nat3
    do na2 = na1, nat3
      do i = 1, N3
        do j = i, N3
          CALL random_number(VRR(na1,na2,i,j))
          VRR(na2,na1,j,i) = VRR(na1,na2,i,j)
        enddo
      enddo
    enddo
  enddo

  !> builds U
  do iq = 1, N3
    call rnd_vec_cmplx(U(:,iq))
  enddo

  !> builds UR, it's better to separate them (?)
  do iq = 1, N3
    do i = 1, N3
      UR(:,i,iq) = U(:,iq) * e_iqr(q(:,iq),R(:,i))
    enddo
  enddo

  !
  !> builds VKK from back FFT2
  VKK = 0._dp
  do iq = 1, N3
    do jq = 1, N3
      do i = 1, N3
        do j = 1, N3
          ! if (ABS(e_iqr(-q(:,iq), R(:, i)) - e_iqr(-q(:,jq), R(:, j))) > 1e-10) print*, "no"
          VKK(:,:,iq,jq) = VKK(:,:,iq,jq) + VRR(:,:,i,j) * &
            e_iqr(-q(:,iq), R(:, i)) * e_iqr(-q(:,jq), R(:,j))
        enddo
      enddo
    enddo
  enddo
  !
  !> flatten UR
  do iq = 1, N3
    do i = 1, N3
      do na = 1, nat3
        UR_flat(na + (i-1)*nat3,iq) = UR(na,i,iq)
      enddo
    enddo
  enddo
  !
  !> builds GRR with FFT
  GRR = 0._dp
  do i = 1, N3
    do j = 1, N3
      Rij = index2v(i, mesh) - index2v(j, mesh)
      do iq = 1, N3
        GRR(:,:,i,j) = GRR(:,:,i,j) + outer_product2(UR(:,i,iq),UR(:,j,iq))
      enddo
    enddo
  enddo

  !
  V_flat = reshape_RR_real(VRR, nat3, N3)
  !
  G_flat = reshape_RR_cmplx(GRR, nat3, N3)
  T_flat = matmul(matmul(V_flat, G_flat), V_flat)

  ! !> first part of TRR product
  ! temp = 0._dp
  ! do i = 1, N3
  !   do j = 1, N3
  !     do inner = 1, N3
  !       temp(:,:,i,j) = temp(:,:,i,j) + matmul( &
  !         VRR(:,:,i,inner), GRR(:,:,inner,j))
  !     enddo
  !   enddo
  ! enddo

  ! !> second part of TRR product
  ! TRR = 0._dp
  ! do i = 1, N3
  !   do j = 1, N3
  !     do inner = 1, N3
  !       TRR(:,:,i,j) = TRR(:,:,i,j) + &
  !         matmul(temp(:,:,i,inner), VRR(:,:,inner,j))
  !     enddo
  !   enddo
  ! enddo

  !> RR braket
  ! full_braket = 0.0_dp
  ! do i = 1, N3
  !   do j = 1, N3
  !     full_braket = full_braket + &
  !       braket(UR(:,i,rn), TRR(:,:,i,j), UR(:,j,rn))
  !   enddo
  ! enddo

  ! !>
  full_braket = 0.0_dp
  do iq = 1, N3
    full_braket = full_braket + &
      braket(UR_flat(:,rn), V_flat, UR_flat(:,iq)) * &
      braket(UR_flat(:,iq), V_flat, UR_flat(:,rn))
  enddo

  !> flattened braket
  flat_braket = braket(UR_flat(:,rn), T_flat)

  !> reduced KK braket
  reduced_braket = 0.0_dp
  do iq = 1, N3
    miq = minus_ind(q(:,iq), mesh)
    reduced_braket = reduced_braket + &
      ABS(braket(U(:,rn), VKK(:,:,rn,miq), U(:,iq)))**2
  enddo

  ! do iq = 1, N3
  !   do jq = 1, N3
  !     do na1 = 1, nat3
  !       do na2 = 1, nat3
  !         print"(4I4)", iq, jq
  !         print"(4F9.3)", VKK(na1,na2,iq,jq), VKK(na2,na1,jq,iq)
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! reduced_braket = braket(U(:,rn), VKK(:,:,rn,mrn))
  if (ABS(full_braket - reduced_braket) < 1e-8) then
    print *, 'Test passed ----------------------'
  else
    print *, 'Test failed'
  endif
    print'(A,F9.1,1X,F9.1)', 'flattened:      ', flat_braket
    print'(A,F9.1,1X,F9.1)', 'Full braket:    ', full_braket
    print'(A,F9.1,1X,F9.1)', 'Reduced braket: ', reduced_braket
  !
contains
  function minus_ind(q, mesh)
    real(dp), intent(in) :: q(3)
    integer, intent(in) :: mesh(3)
    integer :: minus_ind
    real(dp) :: mq(3)
    !
    mq = -q
    where(mq < 0) mq = mq + 1
    minus_ind = v2index(NINT(mq*mesh), mesh)
  end function
  !
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
  function reshape_RR_cmplx(mat, nat3, N3) result(mat_flat)
    complex(dp) :: mat(nat3,nat3,N3,N3)
    integer :: nat3, N3
    !
    complex(dp), dimension(nat3*N3, nat3*N3) :: mat_flat
    integer :: i, j, na1, na2, n1, n2
    !
    do n1 = 1, N3
      do n2 = 1, N3
        do na1 = 1, nat3
          do na2 = 1, nat3
            i = na1 + (n1-1)*nat3
            j = na2 + (n2-1)*nat3
            mat_flat(i,j) = mat(na1,na2,n1,n2)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function reshape_RR_real(mat, nat3, N3) result(mat_flat)
    real(dp) :: mat(nat3,nat3,N3,N3)
    integer :: nat3, N3
    !
    real(dp), dimension(nat3*N3, nat3*N3) :: mat_flat
    integer :: i, j, na1, na2, n1, n2
    !
    do n1 = 1, N3
      do n2 = 1, N3
        do na1 = 1, nat3
          do na2 = 1, nat3
            i = na1 + (n1-1)*nat3
            j = na2 + (n2-1)*nat3
            mat_flat(i,j) = mat(na1,na2,n1,n2)
          enddo
        enddo
      enddo
    enddo
  end function
end program prova
