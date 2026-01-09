module defect_eig
  use kinds, only: dp
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, mat2_diag
  !   fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  ! only: outer_product, freq_in_grid, interp1_matrix, &
  ! id_mat, braket, index2v, v2index, e_iqr, interp1_tns4, grid_vec
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use quter_defect, only : forceconst2_sc, flatten_RR_cmplx
contains
  function D_or_V(S, fc2, fc2_sc, grid, which)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid
    character(1) :: which
    !
    integer :: iq, jq, ibnd, jbnd, compi, compj
    complex(dp) :: D_or_V(grid%nqtot*S%nat3,grid%nqtot*S%nat3)
    complex(dp) :: U(S%nat3,S%nat3,grid%nqtot)
    complex(dp) :: V(S%nat3,S%nat3)
    real(dp) :: w0(S%nat3*grid%nqtot), w2(S%nat3*grid%nqtot)
    !
    D_or_V = 0._dp
    do iq = 1, grid%nqtot
      do ibnd = 1, S%nat3
        compi = ibnd + S%nat3*(iq-1)
        call freq_phq_safe(grid%xq(:,iq), S, fc2, w0((iq-1)*S%nat3+1:iq*S%nat3), U(:,:,iq))
        if(which == "D") &
          D_or_V(compi,compi) = w0(compi)**2
      enddo
    enddo
    !
    do iq = 1, grid%nqtot
      call fc2_sc%r2q(grid%xq(:,iq))
      do jq = 1, grid%nqtot
        call fc2_sc%r2q(grid%xq(:,jq), V)
        V = matmul(conjg(transpose(U(:,:,iq))), matmul(V, U(:,:,jq)))
        do ibnd = 1, S%nat3
          do jbnd = 1, S%nat3
            compi = ibnd + S%nat3*(iq-1)
            compj = jbnd + S%nat3*(jq-1)
            D_or_V(compi,compj) = D_or_V(compi,compj) + V(ibnd,jbnd) / grid%nqtot
          enddo
        enddo
      enddo
    enddo
    !
  end function
  !
  subroutine full_diag(S, fc2, fc2_sc, grid, E)
    use constants, only : RY_TO_CMM1
    use thutils, only : freq_in_grid
    !
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid
    real(dp), intent(in) :: E
    !
    integer :: iq, jq, ibnd, compi
    complex(dp) :: D(grid%nqtot*S%nat3,grid%nqtot*S%nat3)
    real(dp) :: w0(S%nat3*grid%nqtot), w2(S%nat3*grid%nqtot)
    complex(dp) :: U_dagger(S%nat3, S%nat3)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    real(dp) :: freqs(S%nat3, grid%nqtot)
    !
    call freq_in_grid(S, fc2, grid, freqs, Us)
    !
    D = D_or_V(S, fc2, fc2_sc, grid, "D")
    do iq = 1, grid%nqtot
      U_dagger = conjg(transpose(Us(:,:,iq)))
      U_dagger(:,4:) = 0._dp
      do jq = 1, grid%nqtot
        D((iq-1)*S%nat3+1:iq*S%nat3, (jq-1)*S%nat3+1:jq*S%nat3) = &
          D((iq-1)*S%nat3+1:iq*S%nat3, (jq-1)*S%nat3+1:jq*S%nat3) + &
          matmul(U_dagger, Us(:,:,jq)) * fc2_sc%mass_ratios(1)*E**2 / grid%nqtot
      end do
    end do
    !
    call mat2_diag(grid%nqtot*S%nat3, D, w2)
    !
    do iq = 1, grid%nqtot
      do ibnd = 1, S%nat3
        call freq_phq_safe(grid%xq(:,iq), S, fc2, w0((iq-1)*S%nat3+1:iq*S%nat3))
      enddo
    enddo
    !
    open(10, file="eig.dat", status='replace', action='write')
    do compi = 1, size(w0)
      write(10, "(2E20.8)") w0(compi)*RY_TO_CMM1, SQRT(w2(compi))*RY_TO_CMM1
    enddo
    close(10)
    !
  end subroutine
  !
  function largest_eig(n, A, i)
    integer , intent(in) :: n, i
    complex(dp), intent(in) :: A(n,n)
    !
    real(dp) :: W(n)
    integer :: info, M, ISUPPZ(2*n)
    complex(dp) :: Z(n,1)
    real(dp), allocatable :: RWORK(:)
    complex(dp), allocatable :: WORK(:)
    integer, allocatable :: IWORK(:)
    real(dp) :: rwork_query, work_query
    integer :: iwork_query
    real(dp) :: ABSTOL
    !
    real(dp) :: largest_eig(i)
    !
    ! Workspace query
    allocate(WORK(1), RWORK(1), IWORK(1))
    ABSTOL = 0.0_dp
    call zheevr('N','I','U', n, A, n, 0.0_dp, 0.0_dp, n-i+1, n, ABSTOL, M, W, Z, n, ISUPPZ, &
      WORK, -1, RWORK, -1, IWORK, -1, info)
    if (info /= 0) stop "Workspace query failed"
    ! Extract recommended workspace sizes
    work_query  = max(1.0_dp, real(WORK(1)))
    rwork_query = max(1.0_dp, RWORK(1))
    iwork_query = max(1, IWORK(1))

    ! Allocate with those sizes (convert to integer safely)
    deallocate(WORK, RWORK, IWORK)
    allocate(WORK(int(work_query)))
    allocate(RWORK(int(rwork_query)))
    allocate(IWORK(int(iwork_query)))

    call zheevr('N','I','U', n, A, n, 0.0_dp, 0.0_dp, n-i+1, n, ABSTOL, M, W, Z, n, ISUPPZ, &
      WORK, size(WORK), RWORK, size(RWORK), IWORK, size(IWORK), info)

    if (info /= 0) stop "ZHEEVR failed"

    largest_eig = W(1:i)
  end function
  !
  function lowest_eig(n, A, i)
    integer , intent(in) :: n, i
    complex(dp), intent(in) :: A(n,n)
    !
    real(dp) :: W(n)
    integer :: info, M, ISUPPZ(2*n)
    complex(dp) :: Z(n,1)
    real(dp), allocatable :: RWORK(:)
    complex(dp), allocatable :: WORK(:)
    integer, allocatable :: IWORK(:)
    real(dp) :: rwork_query, work_query
    integer :: iwork_query
    real(dp) :: ABSTOL
    !
    real(dp) :: lowest_eig(i)
    !
    ! Workspace query
    allocate(WORK(1), RWORK(1), IWORK(1))
    ABSTOL = 0.0_dp
    call zheevr('N','I','U', n, A, n, 0.0_dp, 0.0_dp, 1, i, ABSTOL, M, W, Z, n, ISUPPZ, &
      WORK, -1, RWORK, -1, IWORK, -1, info)
    if (info /= 0) stop "Workspace query failed"
    ! Extract recommended workspace sizes
    work_query  = max(1.0_dp, real(WORK(1)))
    rwork_query = max(1.0_dp, RWORK(1))
    iwork_query = max(1, IWORK(1))

    ! Allocate with those sizes (convert to integer safely)
    deallocate(WORK, RWORK, IWORK)
    allocate(WORK(int(work_query)))
    allocate(RWORK(int(rwork_query)))
    allocate(IWORK(int(iwork_query)))

    call zheevr('N','I','U', n, A, n, 0.0_dp, 0.0_dp, 1, i, ABSTOL, M, W, Z, n, ISUPPZ, &
      WORK, size(WORK), RWORK, size(RWORK), IWORK, size(IWORK), info)

    if (info /= 0) stop "ZHEEVR failed"

    lowest_eig = W(1:i)
  end function
  !
  subroutine norm_gV(S, fc2, fc2_sc, grid, n_omega, eta)
    use thtetra, only : set_wg, tetra_output
    use constants, only : RY_TO_CMM1
    use thutils, only : freq_in_grid
    use simtet, only : tetra_init_sym_cmplx, tetra_weights_green_cmplx
    use mpi_thermal, only : num_procs, my_id
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid
    integer, intent(in) :: n_omega
    real(dp), intent(in) :: eta(:)
    !
    real(dp), dimension(grid%nqtot*S%nat3,grid%nqtot*S%nat3) :: id
    complex(dp), dimension(grid%nqtot*S%nat3,grid%nqtot*S%nat3) :: VK, A, one_minus_gV, V
    integer :: iq, jq, ibnd, jbnd, compi, compj, iw, N, ieta
    type(tetra_output) :: wg
    real(dp) :: norm_VK, norm_gV_, w2(grid%nqtot*S%nat3)
    real(dp) :: freqs(S%nat3,grid%nqtot)
    complex(dp), allocatable :: wg_cmplx(:,:)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: U_dagger(S%nat3, S%nat3)
    !
    V = D_or_V(S, fc2, fc2_sc, grid, "V")
    call set_wg(S, fc2, grid, n_omega, wg)
    allocate(wg_cmplx(size((wg%w),1), size(wg%w,2)))
    call freq_in_grid(S, fc2, grid, freqs, Us)
    !
    ! norm_VK = frobenius_norm(VK)
    ! print*, "norm V matrix", norm_VK
    N = grid%nqtot * S%nat3
    id = id_mat(N)
    !
    open(10, file="norm.dat")
    do ieta = 1, size(eta)
      print*, "eta =", eta(ieta)
      call tetra_init_sym_cmplx(grid, S, cmplx(freqs**2, eta(ieta), dp))
      do iw = my_id+1, size(wg%en), num_procs
        ! do concurrent(iq=1:grid%nqtot, jq=1:grid%nqtot, ibnd=1:S%nat3, jbnd=1:S%nat3)
        !   compi = (iq-1)*S%nat3 + ibnd
        !   compj = (jq-1)*S%nat3 + jbnd
        !   A(compi,compj) = V(compi,compj) * wg%w(ibnd,wg%e(iq),iw)
        ! enddo
        ! V = VK
        wg_cmplx = tetra_weights_green_cmplx(wg%en(iw)**2)
        ! do iq = 1, grid%nqtot
        !   U_dagger = conjg(transpose(Us(:,:,iq)))
        !   U_dagger(:,4:) = 0._dp
        !   do jq = 1, grid%nqtot
        !     V((iq-1)*S%nat3+1:iq*S%nat3, (jq-1)*S%nat3+1:jq*S%nat3) = &
        !       V((iq-1)*S%nat3+1:iq*S%nat3, (jq-1)*S%nat3+1:jq*S%nat3) + &
        !       matmul(U_dagger, Us(:,:,jq)) * fc2_sc%mass_ratios(1)*E(iw)**2 / grid%nqtot
        !   end do
        ! end do
        !
        do concurrent(jq=1:grid%nqtot, jbnd=1:S%nat3)
          compj = (jq-1)*S%nat3 + jbnd
          one_minus_gV(:,compj) = V(:,compj) / ((wg%en(iw) - eta(ieta))**2 - freqs(jbnd,jq)**2)
        enddo
        !
        one_minus_gV = matmul(conjg(transpose(id - one_minus_gV)), id - one_minus_gV)
        ! norm_gV_ = frobenius_norm(gVVg)
        ! print*, "iw =", iw, " norm gVg =", norm_gV_
        write(10, "(20E25.8)") wg%en(iw) * RY_TO_CMM1, &
          eta(ieta) * RY_TO_CMM1, SQRT(lowest_eig(N, one_minus_gV, 1))
      enddo
    enddo
    close(10)
    !
  contains
    function frobenius_norm(M)
      complex(dp), intent(in) :: M(:,:)
      real(dp) :: frobenius_norm
      !
      frobenius_norm = sqrt(sum(abs(M)**2))
      !
    end function
    !
    function frobenius_norm_vec(wg, iw)
      integer, intent(in) :: iw
      type(tetra_output), intent(in) :: wg
      real(dp) :: frobenius_norm_vec
      integer :: i, Nq
      !
      ! N = size(wg%w, 1) * size(wg%w, 2)
      Nq = size(wg%w, 2)
      frobenius_norm_vec = 0._dp
      do iq = 1, size(wg%w, 2)
        do ibnd = 1, size(wg%w, 1)
          frobenius_norm_vec = frobenius_norm_vec + Nq * wg%qw(iq) * abs(wg%w(ibnd,iq,iw))**2
        enddo
      enddo
      frobenius_norm_vec = sqrt(frobenius_norm_vec)
      !
    end function
  end subroutine
  !
  subroutine norm_V(S, fc2, fc2_sc)
    use q_grids, only : setup_simple_grid
    use fc2_interpolate, only : mat2_diag
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    !
    integer :: i, N, j, n_eigs
    complex(dp) :: trace
    complex(dp), allocatable :: V(:,:), V1(:,:)
    type(q_grid) :: grid
    real(dp) :: eigs(7**3*S%nat3,8)
    real(dp), allocatable :: eig_v(:)
    !
    n_eigs = size(eigs,1)
    eigs = 0._dp
    open(10, file="eig-V.dat", status='replace', action='write')
    do i = 3, 7
      call setup_simple_grid(S%bg, i, i, i, grid)
      N = grid%nqtot * S%nat3
      allocate(V(N,N), V1(N,N), eig_v(N))
      V = D_or_V(S, fc2, fc2_sc, grid, "V")
      !
      V1 = V
      trace = 0._dp
      do j = 1, N
        trace = trace + V(j,j)
      enddo
      print*, "Grid ", i, "x", i, "x", i, " trace V: ", trace
      ! eigs(:n_eigs/2,i) = lowest_eig(N, V, n_eigs/2)
      ! eigs(n_eigs/2+1:,i) = largest_eig(N, V1, n_eigs/2)
      call mat2_diag(N, V, eig_V)
      eigs(:N/2,i) = eig_V(1:N/2)
      eigs(n_eigs-N/2+1:,i) = eig_V(N/2+1:)
      call grid%destroy()
      deallocate(V, V1, eig_V)
    enddo
    !
    do i = 1, n_eigs
      write(10, "(100E20.8)") eigs(i,:)
    enddo
    close(10)
  end subroutine
  !
  ! subroutine green_inversion(S, fc2, fc2_sc, out_grid)
  !   type(ph_system_info), intent(in) :: S
  !   type(forceconst2_grid), intent(in) :: fc2
  !   type(forceconst2_sc), intent(inout) :: fc2_sc
  !   type(q_grid), intent(in) :: out_grid
  !   ! type(code_input_type), intent(in) :: input
  !   !
  !   type(q_grid) :: grid_sym
  !   integer :: nR, iq, jq, ibnd, k, iR, ieta, nq, i, j, jbnd, compi, compj
  !   integer, parameter :: SMEARINGS = 5
  !   real(dp), allocatable :: Eoo(:), E0(:,:), R(:,:), lws(:,:,:), E0_sym(:,:)
  !   complex(dp), allocatable :: U0(:,:,:), Uoo(:,:), braks(:,:,:), &
  !     G_mat(:,:), S_mat(:,:), V(:,:,:,:), V_sum(:,:)
  !   complex(dp) :: G(SMEARINGS), brak, D(S%nat3, S%nat3)
  !   character(20) :: filename
  !   character(20) :: file_energies
  !   logical :: file_exists
  !   real(dp) :: eta(SMEARINGS)
  !   !
  !   nR = product(fc2%nq)
  !   nq = out_grid%nqtot
  !   allocate(Eoo(S%nat3*nq), E0(S%nat3, nq))
  !   allocate(R(3, nR))
  !   allocate(lws(SMEARINGS, S%nat3, nq))
  !   allocate(U0(S%nat3, S%nat3, nq))
  !   allocate(Uoo(S%nat3*nq, S%nat3*nq))
  !   allocate(braks(S%nat3*nR, S%nat3, nq))
  !   allocate(G_mat(S%nat3*nq, S%nat3*nq))
  !   allocate(S_mat(S%nat3*nq, S%nat3*nq))
  !   allocate(V(S%nat3, S%nat3, nq, nq))
  !   allocate(V_sum(nq, nq))
  !   !
  !   call freq_in_grid(S, fc2, out_grid, E0, U0)
  !   E0 = E0**2
  !   print*, "mean E0",  sum(E0)  / real(S%nat3*nq, dp)
  !   do iq = 1, nq
  !     print"(6E15.3)", E0(:,iq)
  !   enddo
  !   write(file_energies, '(A2,I2.2,A4)') "E-", out_grid%n(1), '.dat'
  !   inquire(file=file_energies, exist=file_exists)

  !   file_exists = .false.
  !   print*, file_exists
  !   if (file_exists) then
  !     open(10, file=file_energies, status='old')
  !     ! do j = 1, S%nat3*nq
  !     !   do i = 1, S%nat3*nq
  !     !     read(10, "(2F15.6)") Uoo(i,j)
  !     !   enddo
  !     ! enddo
  !     !
  !     do j = 1, S%nat3*nq
  !       read(10, "(2E15.5)") Eoo(j)
  !     enddo
  !     close(10)
  !     call print_message("read energies from file")
  !   else
  !     ! do iq = 1, out_grid%nq
  !     !   call fc2_sc%r2q(out_grid%xq(:,iq))
  !     !   call fc2_sc%r2q(out_grid%xq(:,iq), D)
  !     !   call mat2_diag(S%nat3, D, E0(:,iq))
  !     !   print"(6E15.3)", E0(:,iq)
  !     ! enddo
  !     !
  !     call full_V(S, fc2, out_grid, fc2_sc, V, U0)
  !     do iq = 1, nq
  !       ! ! call mat2_diag(S%nat3, V(:,:,iq,iq), E0(:,iq))
  !       call fc2_sc%r2q(out_grid%xq(:,iq))
  !       ! do ibnd = 1, S%nat3
  !       !   E0(ibnd,iq) = braket(U0(:,ibnd,iq), D)
  !       ! enddo
  !       ! print"(6E15.3)", E0(:,iq)
  !       do jq = 1, nq
  !         call fc2_sc%r2q(out_grid%xq(:,jq), D)
  !         V_sum(jq,iq) = sum(D)
  !       enddo
  !     enddo
  !     Uoo = flatten_RR_cmplx(V)
  !     !
  !     call print_message("end of interpolation")
  !     ! Uoo = fc2d%fc(:,:,1)
  !     call mat2_diag(S%nat3*nq, Uoo, Eoo)
  !     call print_message("end of frequency calculation")
  !     !
  !     open(10, file=file_energies, status='replace', action='write')
  !     ! do j = 1, S%nat3*nq
  !     !   do i = 1, S%nat3*nq
  !     !     write(10, "(2F15.6)") Uoo(i,j)
  !     !   enddo
  !     ! enddo
  !     !
  !     do j = 1, S%nat3*nq
  !       write(10, "(E15.5)") Eoo(j)
  !     enddo
  !     close(10)
  !   endif

  !   print*, "mean Eoo", sum(Eoo) / real(S%nat3*nq, dp)
  !   ! call q_grid_copy(out_grid, grid_sym)
  !   ! call grid_sym%symmetrize(S)
  !   ! allocate(E0_sym(S%nat3, grid_sym%nqtot))
  !   ! call freq_in_grid(S, fc2, grid_sym, E0_sym)
  !   ! E0_sym = E0_sym**2
  !   ! call tetra_init_sym(grid_sym, S, E0_sym, .false.)

  !   eta = [1e-4_dp, 1e-6_dp, 1e-8_dp, 1e-10_dp, 1e-12_dp]
  !   ! do iR = 1, nR
  !   !   R(:,iR) = REAL(index2v(iR, fc2%nq), dp)
  !   ! enddo
  !   ! call cryst_to_cart(nR, R, S%at, 1)
  !   !
  !   ! do iq = 1, nq
  !   !   do ibnd = 1, S%nat3
  !   !     compi = ibnd + (iq-1)*S%nat3
  !   !     do k = 1, S%nat3 * nq
  !   !       S_mat(k,compi) = dot_product(Uoo((iq-1)*S%nat3+1:iq*S%nat3,k), U0(:,ibnd,iq))
  !   !     enddo
  !   !   enddo
  !   ! enddo
  !   ! !
  !   ! call print_message("end of S matrix calculation")
  !   ! G_mat = 0._dp
  !   ! do iq = 1, nq
  !   !   do ibnd = 1, S%nat3
  !   !     compi = ibnd + (iq-1)*S%nat3
  !   !     do jq = 1, nq
  !   !       do jbnd = 1, S%nat3
  !   !         compj = jbnd + (jq-1)*S%nat3
  !   !         G = 0._dp
  !   !         do k = 1, S%nat3 * nq
  !   !           ! brak = 0._dp
  !   !           ! do iR = 1, nR
  !   !           !   brak = brak + &
  !   !           !     dot_product(U0(:,ibnd,iq), Uoo((iR-1)*S%nat3+1:iR*S%nat3, k)) * &
  !   !           !     e_iqr(out_grid%xq(:,iq), R(:,iR))
  !   !           ! enddo
  !   !           ! brak = dot_product(U0(:,ibnd,iq), Uoo((iq-1)*S%nat3+1:iq*S%nat3,k))
  !   !           ! do ieta = 1, SMEARINGS
  !   !           !   G(ieta) = G(ieta) + abs(brak)**2 / cmplx(1.14e-5_dp-Eoo(k), eta(ieta), dp)
  !   !           ! enddo
  !   !           G_mat(compi, compj) = G_mat(compi, compj) + &
  !   !             conjg(S_mat(k, compi)) * S_mat(k, compj) / cmplx(1.14e-5_dp-Eoo(k), eta(3), dp)
  !   !         enddo
  !   !       enddo
  !   !     enddo
  !   !   enddo
  !   ! enddo

  !   ! ! write(filename, '(A2,I2.2,A4)') "G-", out_grid%n(1), '.dat'
  !   ! ! open(10, file=filename, status='unknown')
  !   ! ! do iq = 1, nq
  !   ! !   do ibnd = 1, S%nat3
  !   ! !     write(10, "(e14.5,5e14.5,I3,3F7.2)") E0(ibnd,iq), lws(:,ibnd,iq), ibnd, out_grid%xq(:,iq)
  !   ! !   enddo
  !   ! ! enddo
  !   ! ! close(10)
  !   ! write(filename, '(A2,I2.2,A4)') "M-", out_grid%n(1), '.dat'
  !   ! open(10, file=filename, status='unknown')
  !   ! do i = 1, S%nat3*nq
  !   !   write(10, "(10000E14.5)") G_mat(:,i)
  !   ! enddo
  !   ! close(10)

  !   ! write(filename, '(A3,I2.2,A4)') "E0-", out_grid%n(1), '.dat'
  !   ! open(10, file=filename, status='unknown')
  !   ! do iq = 1, nq
  !   !   do ibnd = 1, S%nat3
  !   !     write(10, "(E15.5,I2,3F7.2)") E0(ibnd,iq), ibnd, out_grid%xq(:,iq)
  !   !   enddo
  !   ! enddo
  !   ! close(10)

  !   write(filename, '(A2,I2.2,A4)') "V-", out_grid%n(1), '.dat'
  !   open(10, file=filename, status='unknown')
  !   do i = 1, nq
  !     write(10, "(10000E14.5)") V_sum(:,i)
  !   enddo
  !   close(10)
  ! end subroutine
  !
end module
