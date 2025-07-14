program test_gauge_fix
  use test_print
  use thutils, only : diag, qr_gauge_fix
  implicit none
  !
  complex(dp) :: Udeg_old(6,2), Udeg1_old(6,2)
  complex(dp) :: Udeg(6,2), Udeg1(6,2), P(6,6)
  !
  call rnd_mat_cmplx(Udeg)
  call orthonormalize(Udeg)
  P = matmul(Udeg, transpose(conjg(Udeg)))
  Udeg = P(:,1:2)
  ! call orthonormalize(Udeg)
  !
  Udeg_old = Udeg
  Udeg1(:,1) = Udeg(:,1) + 0.7 * Udeg(:,2)
  Udeg1(:,2) = Udeg(:,2) + 0.5 * Udeg(:,1)
  call orthonormalize(Udeg1)
  ! P = matmul(Udeg1, transpose(conjg(Udeg1)))
  ! Udeg1 = P(:,1:2)
  ! call orthonormalize(Udeg1)
  call proj_gauge_fix(Udeg)
  call proj_gauge_fix(Udeg1)
  !
  print*, allclose(Udeg, Udeg1, 1e-8_dp)
contains
  !
  subroutine align_subspaces(v1, v1_prime, subspace)
    use thutils, only : outer_product
    implicit none
    complex(dp), intent(in)    :: v1(:), v1_prime(:)  ! Input vectors (must be normalized)
    complex(dp), intent(inout) :: subspace(:,:)       ! Subspace to rotate (each column is a vector)
    complex(dp), allocatable   :: u(:), H(:,:)
    real(dp)                   :: norm_u
    integer                    :: N, i

    N = size(v1)

    ! --- Check if v1 and v1_prime are already aligned ---
    if (all(abs(v1 - v1_prime) < 1.0e-12_dp)) return

    ! --- Compute Householder reflector ---
    allocate(u(N), H(N,N))
    u = v1 - v1_prime
    norm_u = sqrt(real(dot_product(u, u), kind=dp))

    if (norm_u > 1.0e-12_dp) then
      u = u / norm_u
      H = -2.0_dp * outer_product(u)  ! H = I - 2uu†
      forall (i=1:N) H(i,i) = H(i,i) + 1.0_dp  ! Add identity
    else
      H = 0.0_dp
      forall (i=1:N) H(i,i) = 1.0_dp  ! No rotation if v1 ≈ v1'
    endif

    ! --- Apply rotation to the subspace ---
    subspace = matmul(H, subspace)

    deallocate(u, H)
  end subroutine
  !
  subroutine polar_gauge_fix(U)
    use functions, only : invzmat
    use test_print, only : allclose
    use fc2_interpolate, only: mat2_diag
    use iso_fortran_env, only: dp => real64
    implicit none
    complex(dp), intent(inout) :: U(:,:)
    integer :: M, N, i
    complex(dp), allocatable :: H(:,:), H2(:,:)
    real(dp), allocatable :: s(:)

    M = size(U, 1)
    N = size(U, 2)

    ! Compute Gram matrix: H = U^† U
    allocate(H(N,N), H2(N,N))
    H = matmul(transpose(conjg(U)), U)

    ! Eigen-decomposition: H = V * diag(s) * V^†
    allocate(s(N))
    H2 = H
    call invzmat(N, H2)
    call mat2_diag(N, H, s)
    ! H now contains eigenvectors (columns), s contains eigenvalues

    ! Build H^{-1/2}
    do i = 1, N
      if (s(i) > 1.0e-12_dp) then
        s(i) = 1/sqrt(s(i))
      else
        s(i) = 0.0_dp
      end if
    end do

    H = matmul(H, matmul(diag(s), transpose(conjg(H))))
    print*, allclose(H2, matmul(H,H), 1e-6_dp)
    ! Reconstruct H^{-1/2}
    ! H := V * diag(1/sqrt(s)) * V^†
    ! call scale_columns(H, s) ! Scale columns of eigenvector matrix

    ! Now H is H^{-1/2}

    ! Final orthonormalization: U := U * H^{-1}
    U = matmul(U, H)

    deallocate(H, s)
  end subroutine polar_gauge_fix
!
  ! subroutine proj_gauge_fix(U)
  !   use constants, only: eps8
  !   use fc2_interpolate, only: mat2_diag
  !   use thutils, only: near
  !   implicit none
  !   complex(dp), intent(inout) :: U(:,:)  ! Input: N x M matrix (M vectors in N-dim space)
  !   complex(dp), allocatable :: Q(:,:), R(:,:)
  !   integer :: i, j, k, N, M

  !   N = size(U, 1)  ! Dimension of the full space
  !   M = size(U, 2)  ! Number of vectors in the subspace

  !   ! --- Step 1: Orthonormalize the vectors (e.g., via QR decomposition) ---
  !   allocate(Q(N, M), R(M, M))
  !   call qr_decomposition(U, Q, R)  ! Replace with your QR routine if needed

  !   ! --- Step 2: Fix the gauge (make first non-zero element real and positive) ---
  !   do j = 1, M
  !     do i = 1, N
  !       if (abs(Q(i,j)) > eps8) then
  !         ! Normalize phase so the first significant element is real and positive
  !         Q(:,j) = Q(:,j) * conjg(Q(i,j)) / abs(Q(i,j))
  !         exit
  !       endif
  !     enddo
  !   enddo

  !   ! --- Step 3: Overwrite U with the gauge-fixed basis ---
  !   U(:,1:M) = Q(:,1:M)

  !   deallocate(Q, R)
  ! end subroutine
  !
  ! Example QR decomposition (replace with LAPACK's ZGEQRF if available)
  subroutine qr_decomposition(A, Q, R)
    complex(dp), intent(in) :: A(:,:)
    complex(dp), intent(out) :: Q(:,:), R(:,:)
    integer :: i, j, k
    complex(dp) :: v(size(A,1)), u(size(A,1))
    Q = A
    R = 0.0_dp
    do k = 1, size(A,2)
      v = Q(:,k)
      do j = 1, k-1
        R(j,k) = dot_product(Q(:,j), v)
        v = v - R(j,k) * Q(:,j)
      enddo
      R(k,k) = sqrt(real(dot_product(v, v), kind=dp))
      Q(:,k) = v / R(k,k)
    enddo
  end subroutine
  !
  subroutine proj_gauge_fix(U)
    use fc2_interpolate, only: mat2_diag
    use thutils, only : near, outer_product
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
end program
