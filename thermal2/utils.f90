module thutils
  use kinds, only: dp
  use input_fc, only: ph_system_info
  use fc2_interpolate, only: freq_phq_safe, forceconst2_grid
  use q_grids, only: q_grid
  use mpi_thermal, only: mpi_bsum
  !
contains
  !
  subroutine freq_in_grid(S, fc2, grid, freqs, Us)
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    type(forceconst2_grid), intent(in) :: fc2
    !
    real(dp), intent(out):: freqs(S%nat3, grid%nqtot)
    complex(dp), intent(out), optional :: Us(S%nat3, S%nat3, grid%nqtot)
    !
    integer :: iq, iqp
    !
    freqs = 0.0_dp
    if (PRESENT(Us)) Us = 0.0_dp
    DO iq = 1, grid%nq
      iqp = iq + grid%iq0
      if (present(Us)) then
        CALL freq_phq_safe(grid%xq(:,iq), S, fc2, freqs(:,iqp), Us(:,:,iqp))
      else
        CALL freq_phq_safe(grid%xq(:,iq), S, fc2, freqs(:,iqp))
      end if
    END DO
    if (grid%scattered) CALL mpi_bsum(S%nat3, grid%nqtot, freqs)
    IF (grid%scattered .and. PRESENT(Us)) CALL mpi_bsum(S%nat3, S%nat3, grid%nqtot, Us)
  end subroutine
  !
  pure function outer_product(A) result(AA)
    COMPLEX(dp), intent(in) :: A(:)
    COMPLEX(dp), allocatable :: AA(:,:)
    COMPLEX(dp) :: r
    integer :: nA, i, j
    nA=size(A)
    allocate(AA(nA,nA))

    AA = 0.0_dp
    DO i = 1, nA
      DO j = i+1, nA
        r = A(i) * CONJG(A(j))
        AA(i,j) = r
        AA(j,i) = CONJG(r)
      END DO
      AA(i,i) = A(i) * CONJG(A(i))
    END DO
  end function
  !
  pure function outer_product2(A,B) result(AA)
    COMPLEX(dp), intent(in) :: A(:), B(:)
    COMPLEX(dp), allocatable :: AA(:,:)
    COMPLEX(dp) :: r
    integer :: nA, i, j,nB
    nA=size(A)
    nB=size(B)
    allocate(AA(nA,nB))

    AA = 0.0_dp
    DO i = 1, nA
      DO j = 1, nB
        AA(i,j) = A(i) * CONJG(B(j))
      END DO
    END DO
  end function
  !
  function interp1_matrix(matrices, x)
    !> interpolates linearly a 2D matrix on a 1D grid,
    !> the last dimension of matrices is the grid index.
    !> The matrices should be calculated from 0 to N included
    complex(dp), intent(in) :: matrices(:,:,0:)
    real(dp), intent(in) :: x
    complex(dp), allocatable :: interp1_matrix(:,:)
    real(dp) :: dx
    integer :: x0
    !
    if(x > size(matrices, 3)) call errore('interp1_matrix', ': x out of range', INT(x))
    allocate(interp1_matrix(size(matrices, 1), size(matrices, 2)))
    !
    x0 = INT(x)
    dx = x - x0

    interp1_matrix = (1.0_dp - dx) * matrices(:,:,x0) + &
      dx * matrices(:,:,x0+1)
    !
  end function
  !
  function interp1_tns4(matrices, x)
    !> interpolates linearly a 2D matrix on a 1D grid,
    !> the last dimension of matrices is the grid index.
    !> The matrices should be calculated from 0 to N included
    complex(dp), intent(in) :: matrices(:,:,:,:,0:)
    real(dp), intent(in) :: x
    complex(dp), allocatable :: interp1_tns4(:,:,:,:)
    real(dp) :: dx
    integer :: x0
    !
    if(x > size(matrices, 5)) call errore('interp1_tns4', ': x out of range', INT(x))
    allocate(interp1_tns4, source=matrices(:,:,:,:,0))
    !
    x0 = INT(x)
    dx = x - x0

    interp1_tns4 = (1.0_dp - dx) * matrices(:,:,:,:,x0) + &
      dx * matrices(:,:,:,:,x0+1)
    !
  end function
  !
  function id_mat(n)
    integer, intent(in) :: n
    real(dp) :: id_mat(n,n)
    integer :: i
    !
    id_mat = 0.0_dp
    do i = 1, n
      id_mat(i,i) = 1.0_dp
    enddo
  end function
  !
  pure function braket(left, matrix, right)
    complex(dp), intent(in) :: left(:)
    complex(dp), intent(in) :: matrix(:,:)
    complex(dp), intent(in), optional :: right(:)
    complex(dp) :: braket
    !
    if(PRESENT(right)) then
      braket = dot_product(left, matmul(matrix, right))
    else
      braket = dot_product(left, matmul(matrix, left))
    end if
  end function
  !
  pure function e_iqr(q, R)
    use constants, only: tpi
    real(DP), intent(in) :: q(3), R(3)
    complex(DP) :: e_iqr
    !
    real(DP) :: arg
    arg = tpi * dot_product(q, R)
    e_iqr = CMPLX(COS(arg), SIN(arg), kind=DP)
  end function
  !
  pure function v2index(v, mesh)
    !! Multiplex index of a single wave vector.
    !! Output is always 1-based.
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base states whether v has 0- or 1-based indexing.

    integer, intent(in) :: v(3), mesh(3)
    integer :: v2index

    v2index = (v(1)*mesh(2) + v(2))*mesh(3) + v(3) + 1
  end function
  !
  pure function index2v(i, mesh)
    !! Demultiplex index of a single wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base chooses whether v has 0- or 1-based indexing.

    integer, intent(in) :: i, mesh(3)
    integer :: index2v(3)
    integer :: aux

    call int_div(i - 1, mesh(3), aux, index2v(3))
    call int_div(aux, mesh(2), index2v(1), index2v(2))
  end function
  !
  pure subroutine int_div(num, denom, q, r)
    !! Quotient(q) and remainder(r) of the integer division num/denom.
    integer, intent(in) :: num, denom
    integer, intent(out) :: q, r

    q = num/denom
    r = mod(num, denom)
  end subroutine int_div
  !
  function grid_vec(mesh)
    !! Generates a grid of wave vectors.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.

    integer, intent(in) :: mesh(3)
    integer, allocatable :: grid_vec(:,:)
    integer :: i, j, k, n

    n = product(mesh)
    allocate(grid_vec(3,n))

    n = 0
    do i = 0, mesh(1)-1
      do j = 0, mesh(2)-1
        do k = 0, mesh(3)-1
          n = n + 1
          grid_vec(:,n) = [i, j, k]
        end do
      end do
    end do
  end function
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
  !
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
  subroutine print_message(message)
    use mpi_thermal, only: ionode
    character(len=*), intent(in) :: message
    !
    character(len(message)) :: dashes
    !
    dashes = repeat('-', len(message))
    if (ionode) print*, " "
    if (ionode) print*, dashes
    if (ionode) print*, message
    if (ionode) print*, dashes
    if (ionode) print*, " "
    !
  end subroutine
  !
end module
