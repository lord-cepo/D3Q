module thutils
  use kinds, only: dp
  use input_fc, only: ph_system_info
  use fc2_interpolate, only: freq_phq_safe, forceconst2_grid, freq_phq
  use q_grids, only: q_grid
  use mpi_thermal, only: mpi_bsum
  !
  interface cryst2cart
    module procedure cryst2cart_more
    module procedure cryst2cart_one
  end interface
!
contains
  !
  function contain(str, sub_str)
    character(len=*), intent(in) :: str, sub_str
    logical :: contain
    !
    contain = index(str, sub_str) > 0
  end function
  !
  subroutine freq_in_grid(S, fc2, grid, freqs, Us)
    use merge_degenerate, only: merge_degen
    !
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
      call merge_degen(S%nat3, freqs(:,iq), freqs(:,iq))
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
    complex(dp), intent(in) :: matrices(:,:,:)
    real(dp), intent(in) :: x
    complex(dp), allocatable :: interp1_matrix(:,:)
    real(dp) :: dx
    integer :: x0
    !
    if(x > size(matrices, 3)) call errore('interp1_matrix', ': x out of range', INT(x))
    allocate(interp1_matrix(size(matrices, 1), size(matrices, 2)))
    !
    x0 = INT(x+1)
    dx = x+1 - x0

    interp1_matrix = (1.0_dp - dx) * matrices(:,:,x0) + &
      dx * matrices(:,:,x0+1)
    !
  end function
  !
  function interp1_scl(scl, x)
    !> interpolates linearly a 2D matrix on a 1D grid,
    !> the last dimension of scalars is the grid index.
    !> The scalars should be calculated from 0 to N included
    complex(dp), intent(in) :: scl(:)
    real(dp), intent(in) :: x
    complex(dp) :: interp1_scl
    real(dp) :: dx
    integer :: x0
    !
    if(x > size(scl)) call errore('interp1_scl', ': x out of range', INT(x))
    !
    x0 = INT(x) + 1
    dx = x+1 - x0

    interp1_scl = (1.0_dp - dx) * scl(x0) + &
      dx * scl(x0+1)
    !
  end function
  !
  function interp1_vector(vectors, x)
    !> interpolates linearly a 2D matrix on a 1D grid,
    !> the last dimension of vectors is the grid index.
    !> The vectors should be calculated from 0 to N included
    complex(dp), intent(in) :: vectors(:,:)
    real(dp), intent(in) :: x
    complex(dp), allocatable :: interp1_vector(:)
    real(dp) :: dx
    integer :: x0
    !
    if(x > size(vectors, 2)) call errore('interp1_vector', ': x out of range', INT(x))
    allocate(interp1_vector(size(vectors, 1)))
    !
    x0 = INT(x) + 1
    dx = x+1 - x0

    interp1_vector = (1.0_dp - dx) * vectors(:,x0) + &
      dx * vectors(:,x0+1)
    !
  end function
  !
  function interp1_tns4(matrices, x)
    !> interpolates linearly a 2D matrix on a 1D grid,
    !> the last dimension of matrices is the grid index.
    !> The matrices should be calculated from 0 to N included
    complex(dp), intent(in) :: matrices(:,:,:,:,:)
    real(dp), intent(in) :: x
    complex(dp), allocatable :: interp1_tns4(:,:,:,:)
    real(dp) :: dx
    integer :: x0
    !
    if(x > size(matrices, 5)) call errore('interp1_tns4', ': x out of range', INT(x))
    allocate(interp1_tns4, source=matrices(:,:,:,:,1))
    !
    x0 = INT(x) + 1
    dx = x+1 - x0

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
  pure function v2index_n(v, mesh)
    !! Multiplex index of a single wave vector.
    !! Output is always 1-based.
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base states whether v has 0- or 1-based indexing.

    integer, intent(in) :: v(3), mesh(3)
    integer :: v2index_n

    v2index_n = (v(3)*mesh(2) + v(2))*mesh(1) + v(1) + 1
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
  pure function index2v_n(i, mesh)
    !! Demultiplex index of a single wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base chooses whether v has 0- or 1-based indexing.

    integer, intent(in) :: i, mesh(3)
    integer :: index2v_n(3)
    integer :: aux

    call int_div(i - 1, mesh(1), aux, index2v_n(1))
    call int_div(aux, mesh(2), index2v_n(3), index2v_n(2))
  end function
  !
  pure function index2v_cart(i, mesh, at)
    !! Demultiplex index of a single wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base chooses whether v has 0- or 1-based indexing.
    integer, intent(in) :: i, mesh(3)
    real(dp), intent(in) :: at(3,3)
    real(dp) :: index2v_cart(3)
    !
    index2v_cart = cryst2cart(real(index2v(i, mesh), dp), at, 1)
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
  function grid_vec_cryst(mesh, center, natural)
    !! Generates a grid of wave vectors.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.

    integer, intent(in) :: mesh(3)
    integer,allocatable :: grid_vec_cryst(:,:)
    logical, intent(in), optional :: center, natural
    integer :: i, j, k, n
    integer :: center_
    logical :: natural_

    !
    natural_ = .false.
    if(present(natural)) natural_ = natural
    !
    allocate(grid_vec_cryst(3,product(mesh)))
    if(present(center)) then
      if(center) then
        center_ = 1
      else
        center_ = 0
      end if
    else
      center_ = 0
    end if

    n = 0
    do i = 0, mesh(1)-1
      do j = 0, mesh(2)-1
        do k = 0, mesh(3)-1
          n = n + 1
          if(natural_) then
            grid_vec_cryst(:,n) = [k, j, i] - center_ * mesh / 2
          else
            grid_vec_cryst(:,n) = [i, j, k] - center_ * mesh / 2
          endif
        end do
      end do
    end do
  end function
  !
  function grid_vec_cart(mesh, at, center, divide, natural)
    !! Generates a grid of wave vectors.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.

    integer, intent(in) :: mesh(3)
    real(dp) :: at(3,3)
    logical, intent(in), optional :: center, divide, natural
    real(dp), allocatable :: grid_vec_cart(:,:)
    logical :: center_, divide_, natural_

    integer :: i
    integer, allocatable :: grid_vec_cryst_(:,:)
    !
    if(present(center)) then
      center_ = center
    else
      center_ = .false.
    end if
    !
    if(present(natural)) then
      natural_ = natural
    else
      natural_ = .false.
    endif
    !
    if(present(divide)) then
      divide_ = divide
    else
      divide_ = .false.
    end if
    !
    allocate(grid_vec_cryst_(3,product(mesh)))
    allocate(grid_vec_cart(3,product(mesh)))
    grid_vec_cryst_ = grid_vec_cryst(mesh, center_, natural_)
    grid_vec_cart = REAL(grid_vec_cryst_, DP)
    deallocate(grid_vec_cryst_)
    if (divide_) then
      do i = 1, size(grid_vec_cart,2)
        grid_vec_cart(:,i) = grid_vec_cart(:,i) / REAL(mesh, DP)
      end do
    endif
    call cryst_to_cart(product(mesh), grid_vec_cart, at, 1)
    !
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
  pure function bz2simple(R, grid)
    integer, intent(in) :: R(3), grid(3)
    integer :: bz2simple(3)
    integer :: i, j
    !
    bz2simple = mod(R, grid)
    do i = 1, 3
      do while (bz2simple(i) < 0)
        bz2simple(i) = bz2simple(i) + grid(i)
      enddo
      do while (bz2simple(i) >= grid(i))
        bz2simple(i) = bz2simple(i) - grid(i)
      enddo
      ! if(bz2simple(i) < 0 .or. bz2simple(i) >= grid(i)) &
      ! call errore('bz2simple', ': R out of range', 1)
      ! stop 1
    end do
  end function
  !
  function near(val1, val2, thr)
    use constants, only: eps6
    real(dp), intent(in) :: val1
    real(dp), intent(in), optional :: val2
    real(dp), intent(in), optional :: thr
    !
    logical :: near
    real(dp) :: threshold, val2_
    !
    if (PRESENT(thr)) then
      threshold = thr
    else
      threshold = eps6
    end if
    !
    if (present(val2)) then
      val2_ = val2
    else
      val2_ = 0.0_dp
    end if
    !
    near = 2*abs(val1 - val2_)/(val1 + val2_) < threshold
  end function
  !
  subroutine qr_gauge_fix(U)
    implicit none
    complex(dp), intent(inout) :: U(:,:)

    integer :: M, N, INFO, i, LWORK
    complex(dp), allocatable :: TAU(:), WORK(:)
    complex(dp) :: R(size(U, 2)), WORK_TEST(1)

    M = size(U, 1)
    N = size(U, 2)
    allocate(TAU(min(M,N)))

    ! === QUERY optimal LWORK for ZGEQRF ===
    LWORK = -1
    call ZGEQRF(M, N, U, M, TAU, WORK_TEST, LWORK, INFO)
    LWORK = int(real(WORK_TEST(1)))
    allocate(WORK(LWORK))

    ! === ACTUAL ZGEQRF ===
    call ZGEQRF(M, N, U, M, TAU, WORK, LWORK, INFO)
    if (INFO /= 0) call errore('qr_gauge_fix', ': ZGEQRF failed', INFO)

    ! === Normalize R diagonals ===
    do i = 1, N
      R(i) = U(i,i) / ABS(U(i,i))
    end do

    ! === QUERY optimal LWORK for ZUNGQR ===
    LWORK = -1
    call ZUNGQR(M, N, N, U, M, TAU, WORK_TEST, LWORK, INFO)
    LWORK = int(real(WORK_TEST(1)))
    deallocate(WORK)
    allocate(WORK(LWORK))

    ! === ACTUAL ZUNGQR ===
    call ZUNGQR(M, N, N, U, M, TAU, WORK, LWORK, INFO)
    if (INFO /= 0) call errore('qr_gauge_fix', ': ZUNGQR failed', INFO)

    ! === Gauge fix: absorb R phase into Q ===
    do i = 1, N
      U(:,i) = U(:,i) * R(i)
    end do

    deallocate(WORK, TAU)
  end subroutine
  !
  subroutine polar_gauge_fix(U)
    use test_print, only : allclose
    use fc2_interpolate, only: mat2_diag
    use iso_fortran_env, only: dp => real64
    implicit none
    complex(dp), intent(inout) :: U(:,:)
    integer :: M, N, i
    complex(dp), allocatable :: H(:,:)
    real(dp), allocatable :: s(:)

    M = size(U, 1)
    N = size(U, 2)

    ! Compute Gram matrix: H = U^† U
    allocate(H(N,N))
    H = matmul(transpose(conjg(U)), U)

    ! Eigen-decomposition: H = V * diag(s) * V^†
    allocate(s(N))
    call mat2_diag(N, H, s)
    ! H now contains eigenvectors (columns), s contains eigenvalues

    ! Build H^{-1/2}
    do i = 1, N
      if (s(i) > 1.0e-12_dp) then
        s(i) = 1.0_dp / sqrt(s(i))
      else
        s(i) = 0.0_dp
      end if
    end do

    ! Reconstruct H^{-1/2}
    ! H := V * diag(1/sqrt(s)) * V^†
    ! call scale_columns(H, s) ! Scale columns of eigenvector matrix
    H = matmul(H, matmul(diag(s), transpose(conjg(H))))

    ! Now H is H^{-1/2}

    ! Final orthonormalization: U := U * H^{-1}
    U = matmul(U, H)

    deallocate(H, s)
  end subroutine polar_gauge_fix
  !
  function diag(s)
    real(dp), intent(in) :: s(:)
    real(dp) :: diag(size(s), size(s))
    integer :: i
    !
    diag = 0._dp
    do i = 1, size(s)
      diag(i,i) = s(i)
    enddo
  end function
  !
  function diag_cmplx(s)
    complex(dp), intent(in) :: s(:)
    complex(dp) :: diag_cmplx(size(s), size(s))
    integer :: i
    !
    diag_cmplx = 0._dp
    do i = 1, size(s)
      diag_cmplx(i,i) = s(i)
    enddo
  end function
  !
  subroutine scale_columns(A, s)
    use iso_fortran_env, only: dp => real64
    implicit none
    complex(dp), intent(inout) :: A(:,:)
    real(dp), intent(in) :: s(:)
    integer :: i, N, M

    M = size(A, 1)
    N = size(A, 2)

    do i = 1, N
      A(:,i) = A(:,i) * s(i)
    end do
  end subroutine scale_columns
  !
  subroutine proj_gauge_fix(U)
    use fc2_interpolate, only: mat2_diag
    use test_print, only : orthonormalize
    !
    complex(dp), intent(inout) :: U(:,:)
    complex(dp) :: P(size(U, 1), size(U, 1))
    real(dp) :: s(size(U, 1))
    integer :: i, j
    !
    P = matmul(U, transpose(conjg(U)))
    do i = 1, size(U, 2)
      P = P + 1e-3_dp * outer_product(U(:,i))
    end do
    call mat2_diag(size(U, 1), P, s)
    !
    j = 1
    do i = 1, size(U, 1)
      if (s(i) > 1.0e-10_dp) then
        U(:,j) = P(:,i)
        j = j + 1
      endif
    end do
  end subroutine
  !
  pure function cryst2cart_more (vec, trmat, iflag)
    !-----------------------------------------------------------------------
    !
    !     This routine transforms the atomic positions or the k-point
    !     components from crystallographic to cartesian coordinates
    !     ( iflag=1 ) and viceversa ( iflag=-1 ).
    !     Output cartesian coordinates are stored in the input ('vec') array
    !
    !
    implicit none
    !
    integer, intent(in) :: iflag
    ! nvec:  number of vectors (atomic positions or k-points)
    !        to be transformed from crystal to cartesian and vice versa
    ! iflag: gives the direction of the transformation
    real(dp), intent(in) :: trmat (3, 3)
    ! trmat: transformation matrix
    ! if iflag=1:
    !    trmat = at ,  basis of the real-space lattice,       for atoms   or
    !          = bg ,  basis of the reciprocal-space lattice, for k-points
    ! if iflag=-1: the opposite
    real(dp), intent(in) :: vec (:,:)
    ! coordinates of the vector (atomic positions or k-points) to be
    ! transformed - overwritten on output
    !
    !    local variables
    !
    integer :: nv, kpol
    ! counter on vectors
    ! counter on polarizations
    real(dp):: vau (3)
    real(dp) :: cryst2cart_more(3, size(vec,2))
    ! workspace
    !
    !     Compute the cartesian coordinates of each vectors
    !     (atomic positions or k-points components)
    !
    do nv = 1, size(vec,2)
      if (iflag.eq.1) then
        do kpol = 1, 3
          vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
            * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
        end do
      else
        do kpol = 1, 3
          vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
            * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
        end do
      endif
      do kpol = 1, 3
        cryst2cart_more(kpol, nv) = vau(kpol)
      end do
    end do
  end function
  !
  pure function cryst2cart_one (vec, trmat, iflag)
    !-----------------------------------------------------------------------
    !
    !     This routine transforms the atomic positions or the k-point
    !     components from crystallographic to cartesian coordinates
    !     ( iflag=1 ) and viceversa ( iflag=-1 ).
    !     Output cartesian coordinates are stored in the input ('vec') array
    !
    !
    implicit none
    !
    integer, intent(in) :: iflag
    ! nvec:  number of vectors (atomic positions or k-points)
    !        to be transformed from crystal to cartesian and vice versa
    ! iflag: gives the direction of the transformation
    real(dp), intent(in) :: trmat (3, 3)
    ! trmat: transformation matrix
    ! if iflag=1:
    !    trmat = at ,  basis of the real-space lattice,       for atoms   or
    !          = bg ,  basis of the reciprocal-space lattice, for k-points
    ! if iflag=-1: the opposite
    real(dp), intent(in) :: vec (3)
    ! coordinates of the vector (atomic positions or k-points) to be
    ! transformed - overwritten on output
    !
    !    local variables
    !
    integer :: kpol
    ! counter on vectors
    ! counter on polarizations
    real(dp):: vau (3)
    real(dp) :: cryst2cart_one(3)
    ! workspace
    !
    !     Compute the cartesian coordinates of each vectors
    !     (atomic positions or k-points components)
    !
    if (iflag.eq.1) then
      do kpol = 1, 3
        vau (kpol) = trmat (kpol, 1) * vec(1) + trmat (kpol, 2) &
          * vec (2) + trmat (kpol, 3) * vec (3)
      end do
    else
      do kpol = 1, 3
        vau (kpol) = trmat (1, kpol) * vec (1) + trmat (2, kpol) &
          * vec (2) + trmat (3, kpol) * vec (3)
      end do
    endif
    do kpol = 1, 3
      cryst2cart_one(kpol) = vau(kpol)
    end do
  end function
  !
end module
