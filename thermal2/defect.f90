module defect
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi, pi
  use quter_defect
  use functions, only: f_gauss
  use fc3_interpolate, only: forceconst3, sparse, d3_mixed, sum_R3
  use merge_degenerate, only: merge_degen
  use test_print, only: allclose
  use simtet, only: tetra_init_sym_cmplx, tetra_weights_green_cmplx
  use ph_velocity, only : velocity
  use constants, only : RY_TO_CMM1
  use mpi_thermal, only : my_id, num_procs, mpi_bsum, ionode
  use symm_q_mat, only : apply_sym_q
  !
  implicit none
  !
contains
  !
  subroutine tetra_from_self_cart(S, grid, D_plus_Sigma, en2, weights, mpi)
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    complex(dp), intent(in) :: D_plus_Sigma(S%nat3, S%nat3, grid%nqtot)
    real(dp), intent(in) :: en2
    !
    complex(dp), intent(out) :: weights(S%nat3, S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3, grid%nqtot) :: den_weights
    complex(dp), dimension(S%nat3,S%nat3,grid%nqtot) :: den_UL, den_UR
    complex(dp), dimension(S%nat3, grid%nqtot) :: overlap
    logical, optional, intent(in) :: mpi
    !
    complex(dp), dimension(S%nat3, S%nat3) :: DL, DR, temp_U_idx
    integer :: iq, i, idx(S%nat3), iqp, j, k
    complex(dp) :: den_eig(S%nat3, grid%nqtot)
    real(dp) :: real_eig(S%nat3)
    logical :: mpi_
    !
    if(present(mpi)) then
      mpi_ = mpi
    else
      mpi_ = .false.
    endif
    !
    den_eig = 0._dp
    den_UL = 0._dp
    den_UR = 0._dp
    overlap = 0._dp
    !
    do iq = 1, grid%nq
      iqp = iq + grid%iq0
      DL = D_plus_Sigma(:,:,iqp)
      call mat2_diag(S%nat3, DL, DR, den_eig(:,iqp))
      real_eig = real(den_eig(:,iqp),dp)
      idx = 0
      CALL hpsort(S%nat3, real_eig, idx)
      den_eig(:,iqp) = den_eig(idx,iqp)
      den_UL(:,:,iqp) = DL
      den_UR(:,:,iqp) = DR
      temp_U_idx = den_UR(:, idx, iqp)
      den_UR(:, :, iqp) = temp_U_idx
      temp_U_idx = den_UL(:, idx, iqp)
      den_UL(:, :, iqp) = temp_U_idx
      do i = 1, S%nat3
        overlap(i,iqp) = dot_product(den_UL(:,i,iqp), den_UR(:,i,iqp))
      enddo
      ! call merge_degen(S%nat3, den_eig(:,iq), den_eig(:,iq))
    enddo
    if (grid%scattered) then
      call mpi_bsum(S%nat3, grid%nqtot, den_eig)
      call mpi_bsum(S%nat3, S%nat3, grid%nqtot, den_UL)
      call mpi_bsum(S%nat3, S%nat3, grid%nqtot, den_UR)
      call mpi_bsum(S%nat3, grid%nqtot, overlap)
    endif
    where(aimag(den_eig)> 0._dp) den_eig = conjg(den_eig)
    !
    !{ weights of diagonal tetra
    call tetra_init_sym_cmplx(grid, S, den_eig, mpi=mpi_)
    den_weights = tetra_weights_green_cmplx(en2)
    where(isnan(abs(den_weights))) den_weights = 0._dp
    !
    weights = 0._dp
    do iq = 1, grid%nq
      iqp = iq + grid%iq0
      do k = 1, S%nat3
        do j = 1, S%nat3
          do i = 1, S%nat3
            weights(i,j,iqp) = weights(i,j,iqp) + &
              den_weights(k,iqp) * conjg(den_UL(j,k,iqp)) * den_UR(i,k,iqp) / overlap(k,iqp)
          enddo
        enddo
      enddo
    enddo
    call mpi_bsum(S%nat3, S%nat3, grid%nqtot, weights)
    !
  end subroutine
  !
  subroutine tetra_from_self(S, grid, freqs, self, en2, den_weights, den_UL, den_UR, overlap, mpi)
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    real(dp), intent(in) :: freqs(S%nat3, grid%nqtot)
    complex(dp), intent(in) :: self(S%nat3, S%nat3, grid%nqtot)
    real(dp), intent(in) :: en2
    !
    complex(dp), intent(out), dimension(S%nat3, grid%nqtot) :: den_weights
    complex(dp), optional, intent(out), dimension(S%nat3,S%nat3,grid%nqtot) :: den_UL, den_UR
    complex(dp), optional, intent(out), dimension(S%nat3, grid%nqtot) :: overlap
    logical, optional, intent(in) :: mpi
    !
    complex(dp), dimension(S%nat3, S%nat3) :: DL, DR, temp_U_idx
    integer :: iq, i, idx(S%nat3), iqp
    complex(dp) :: den_eig(S%nat3, grid%nqtot)
    real(dp) :: real_eig(S%nat3)
    logical :: mpi_
    !
    if(present(mpi)) then
      mpi_ = mpi
    else
      mpi_ = .false.
    endif
    !
    den_eig = 0._dp
    if(present(den_UL)) den_UL = 0._dp
    if(present(den_UR)) den_UR = 0._dp
    if(present(overlap)) overlap = 0._dp
    !
    do iq = 1, grid%nq
      iqp = iq + grid%iq0
      DL = diag(freqs(:,iqp)**2) + self(:,:,iqp)
      call mat2_diag(S%nat3, DL, DR, den_eig(:,iqp))
      real_eig = real(den_eig(:,iqp),dp)
      idx = 0
      CALL hpsort(S%nat3, real_eig, idx)
      den_eig(:,iqp) = den_eig(idx,iqp)
      if(present(den_UL)) then
        den_UL(:,:,iqp) = DL
        den_UR(:,:,iqp) = DR
        temp_U_idx = den_UR(:, idx, iqp)
        den_UR(:, :, iqp) = temp_U_idx
        temp_U_idx = den_UL(:, idx, iqp)
        den_UL(:, :, iqp) = temp_U_idx
        do i = 1, S%nat3
          overlap(i,iqp) = dot_product(den_UL(:,i,iqp), den_UR(:,i,iqp))
        enddo
      endif
      ! call merge_degen(S%nat3, den_eig(:,iq), den_eig(:,iq))
    enddo
    if (grid%scattered) then
      call mpi_bsum(S%nat3, grid%nqtot, den_eig)
      if(present(den_UL)) call mpi_bsum(S%nat3, S%nat3, grid%nqtot, den_UL)
      if(present(den_UR)) call mpi_bsum(S%nat3, S%nat3, grid%nqtot, den_UR)
      if(present(overlap)) call mpi_bsum(S%nat3, grid%nqtot, overlap)
    endif
    where(aimag(den_eig)> 0._dp) den_eig = conjg(den_eig)
    !
    !{ weights of diagonal tetra
    call tetra_init_sym_cmplx(grid, S, den_eig, mpi=mpi_)
    den_weights = tetra_weights_green_cmplx(en2)
    where(isnan(abs(den_weights))) den_weights = 0._dp
  end subroutine
  !
  subroutine tetra_from_self_diag(S, grid, w2_plus_sigma, en2, den_weights, mpi)
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    complex(dp), intent(in) :: w2_plus_sigma(S%nat3, grid%nqtot)
    real(dp), intent(in) :: en2
    !
    complex(dp), intent(out), dimension(S%nat3, grid%nqtot) :: den_weights
    logical, optional, intent(in) :: mpi
    !
    complex(dp), dimension(S%nat3, S%nat3) :: DL, DR, temp_U_idx
    integer :: iq, i, idx(S%nat3), iqp
    complex(dp) :: den_eig(S%nat3, grid%nqtot)
    real(dp) :: real_eig(S%nat3)
    logical :: mpi_
    !
    if(present(mpi)) then
      mpi_ = mpi
    else
      mpi_ = .false.
    endif
    !
    den_eig = w2_plus_sigma
    where(aimag(den_eig)> 0._dp) den_eig = conjg(den_eig)
    !
    !{ weights of diagonal tetra
    call tetra_init_sym_cmplx(grid, S, den_eig, mpi=mpi_)
    den_weights = tetra_weights_green_cmplx(en2)
    where(isnan(abs(den_weights))) den_weights = 0._dp
  end subroutine
  !
  subroutine create_diffs(fc2_sc, diffs)
    type(forceconst2_sc), intent(in) :: fc2_sc
    integer, intent(out), pointer :: diffs(:,:)
    !
    integer :: ir1, ir2, iR, nR, R(3)
    !
    allocate(diffs(3,1))
    diffs(:,1) = fc2_sc%yR1(:,1,1) - fc2_sc%yR2(:,1)
    nR = 1
    !
    do ir2 = 1, fc2_sc%n_R2
      do ir1 = 1, fc2_sc%n_R1(ir2)
        R = fc2_sc%yR1(:,ir1,ir2)-fc2_sc%yR2(:,ir2)
        call find_where(R, diffs, iR)
        if (iR == -1) call enlarge_R(R, diffs, nR)
      enddo
    enddo
    !
  end subroutine
  !
  subroutine enlarge_R(R, list_of_R, nR)
    integer, intent(in) :: R(3)
    integer, intent(inout) :: nR
    integer, intent(inout), pointer :: list_of_R(:,:)
    integer, pointer :: R_new(:,:)
    !
    allocate(R_new(3, nR+1))
    R_new(:,1:nR) = list_of_R(:,1:nR)
    R_new(:,nR+1) = R
    nR = nR + 1
    deallocate(list_of_R)
    list_of_R => R_new
  end subroutine
  !
  subroutine build_csr_sym(A)!, row_ptr, col_ind, values)
    complex(dp), intent(in) :: A(:,:)
    ! integer, pointer, intent(out) :: row_ptr(:), col_ind(:)
    ! complex(dp), pointer, intent(out) :: values(:)
    !
    integer :: nrows, ncols, i, j, count, N2
    integer, allocatable :: col_ind_(:)
    complex(dp), allocatable :: values_(:)
    !
    nrows = size(A, 1)
    ncols = size(A, 2)
    N2 = nrows * ncols
    !
    ! allocate(row_ptr(nrows + 1))
    allocate(col_ind_(N2))
    allocate(values_(N2))
    !
    count = 0
    ! row_ptr(1) = 1
    do i = 1, nrows
      do j = i, ncols
        if (A(i,j) /= 0._dp) then
          count = count + 1
          col_ind_(count) = j
          values_(count) = A(i,j)
        end if
      end do
      ! row_ptr(i + 1) = count + 1
    end do
    !
    ! Resize col_ind and values arrays to actual number of non-zero elements
    ! allocate(col_ind(count))
    ! allocate(values(count))
    ! col_ind = col_ind_(1:count)
    ! values = values_(1:count)
    deallocate(col_ind_)
    deallocate(values_)
    print*, "Non-zero elements in CSR matrix: ", count / real(N2, dp) * 100.0_dp, "%"
    !
  end subroutine
  !
  subroutine full_born_center(S, input, fc2, fc2_sc, grid, sym_grid, out_grid)
    use quter_defect, only: fc_sc2RR, allocate_fc2_sc
    type(ph_system_info), intent(in) :: S
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid, sym_grid, out_grid
    !
    integer :: iR1, iR2, iR, nR, nR_large
    integer :: i, j, iq, ibnd, N
    real(dp) :: c
    integer, pointer :: R_list(:,:), diff_list(:,:)
    integer :: R(3)
    real(dp), allocatable :: Rij_cart(:,:,:)
    real(dp), allocatable :: diffs(:,:)
    complex(dp), allocatable :: g0(:,:,:)
    complex(dp), allocatable :: T(:,:,:)
    complex(dp) :: Tq(S%nat3, S%nat3,out_grid%nqtot, input%n_omega)
    complex(dp) :: self_energy(S%nat3, out_grid%nqtot, input%n_omega)
    integer :: iw, jn, defat
    type(tetra_output) :: wg, wg_out
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    real(dp) :: freqs(S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    complex(dp) :: out_Us_c(S%nat3, S%nat3, out_grid%nqtot)
    real(dp) :: out_freqs(S%nat3, out_grid%nqtot)
    integer, dimension(maxval(fc2_sc%n_R1),fc2_sc%n_R2):: iR_diff
    integer, pointer :: diff_list_large(:,:)
    real(dp), allocatable :: diff_large(:,:)
    integer, allocatable :: iR_large(:,:)
    complex(dp), allocatable :: phases_out(:,:)
    complex(dp), allocatable :: TRR__(:,:,:)
    !
    complex(dp), allocatable :: VK(:,:), gV__(:,:), g0__(:,:), &
      T__(:,:), I_gV__(:,:), VM(:,:), V__(:,:)
    complex(dp), allocatable :: S__(:,:), S1__(:,:), Sg__(:,:), &
      I_Sg__(:,:), Gm__(:,:), GVS__(:,:), I_GVS__(:,:), G__(:,:)
    integer, allocatable :: g0_iR(:)
    real(dp), parameter :: alpha = 1.0_dp
    complex(dp) :: w_self(S%nat3, out_grid%nqtot)
    real(dp) :: dos(input%n_omega)

    !
    call set_wg(S, fc2, sym_grid, input, wg)
    if(input%calculation == 'dos' .or. input%calculation == 'test') &
      call set_wg(S, fc2, out_grid, input, wg_out)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)

    c = input%conc * size(fc2_sc%defects,2)
    do iq = 1, out_grid%nqtot
      out_Us_c(:,:,iq) = conjg(transpose(out_Us(:,:,iq)))
    enddo

    allocate(R_list(3,1))
    R_list(:,1) = fc2_sc%yR2(:,1)
    nR = 1

    do iR2 = 1, fc2_sc%n_R2
      call find_where(fc2_sc%yR2(:,iR2), R_list, iR)
      if(iR == -1) call enlarge_R(fc2_sc%yR2(:,iR2), R_list, nR)
      do iR1 = 1, fc2_sc%n_R1(iR2)
        call find_where(fc2_sc%yR1(:,iR1,iR2), R_list, iR)
        if(iR == -1) call enlarge_R(fc2_sc%yR1(:,iR1,iR2), R_list, nR)
      enddo
    enddo
    if(ionode) print*, "Number of unique R vectors: ", nR
    !
    call create_diffs(fc2_sc, diff_list)
    allocate(diffs(3,size(diff_list,2)))
    diffs = cryst2cart(real(diff_list,dp), S%at, 1)

    nR_large = 1
    allocate(diff_list_large(3,1))
    diff_list_large(:,1) = [0,0,0]
    do i = 1, nR
      do j = 1, nR
        R = R_list(:,i) - R_list(:,j)
        call find_where(R, diff_list_large, iR)
        if (iR == -1) call enlarge_R(R, diff_list_large, nR_large)
      enddo
    enddo
    !
    allocate(Rij_cart(3,nR,nR))
    allocate(iR_large(nR,nR))
    allocate(phases_out(nR_large,out_grid%nqtot))
    allocate(TRR__(S%nat3,S%nat3,nR_large))
    allocate(g0_iR(nR_large))
    do iq = 1, out_grid%nqtot
      do j = 1, nR
        do i = 1, nR
          if (iq == 1) then
            Rij_cart(:,i,j) = cryst2cart(real(R_list(:,i) - R_list(:,j),dp), S%at, 1)
            call find_where(R_list(:,i) - R_list(:,j), diff_list_large, iR)
            iR_large(i,j) = iR
            R = diff_list_large(:,iR)
            g0_iR(iR) = v2index_n(bz2simple(R, grid%n), grid%n)
          endif
          phases_out(iR_large(i,j),iq) = e_iqr(out_grid%xq(:,iq), -Rij_cart(:,i,j))
        enddo
      enddo
    enddo
    !
    allocate(diff_large(3,nR_large))
    diff_large = cryst2cart(real(diff_list_large,dp), S%at, 1)
    if(ionode) print*, "Number of unique large diff vectors: ", nR_large

    N = S%nat3 * nR
    allocate(VK(N,N))
    allocate(VM(N,N))
    allocate(V__(N,N))
    allocate(g0__(N,N))
    allocate(Gm__, gV__, T__, I_gV__, S__, S1__, &
      Sg__, I_Sg__, GVS__, I_GVS__, G__, source=g0__)
    !
    gV__ = 0._dp
    T__ = 0._dp
    Tq = 0._dp
    dos = 0._dp
    self_energy = 0._dp
    do iR2 = 1, fc2_sc%n_R2
      do iR1 = 1, fc2_sc%n_R1(iR2)
        call find_where(fc2_sc%yR1(:,iR1,iR2)-fc2_sc%yR2(:,iR2), diff_list, iR_diff(iR1,iR2))
      enddo
    enddo

    allocate(T(S%nat3,S%nat3,size(diff_list,2)))
    !
    VM = 0._dp
    VK = 0._dp
    defat = fc2_sc%defects(1,1)
    do iR2 = 1, fc2_sc%n_R2
      call find_where(fc2_sc%yR2(:,iR2), R_list, j)
      do iR1 = 1, fc2_sc%n_R1(iR2)
        call find_where(fc2_sc%yR1(:,iR1,iR2), R_list, i)
        if(all(fc2_sc%yR1(:,iR1,iR2) == fc2_sc%yR2(:,iR2)) .and. &
          all(fc2_sc%yR2(:,iR2) == index2v(fc2_sc%defects(2,1), fc2_sc%nq))) then
          do jn = 1, 3
            VM( (i-1)*S%nat3+(defat-1)*3+jn, (j-1)*S%nat3+(defat-1)*3+jn ) = fc2_sc%eps
          enddo
        endif
        VK( (i-1)*S%nat3+1:i*S%nat3, (j-1)*S%nat3+1:j*S%nat3 ) = &
          cmplx( fc2_sc%fc(:,:,iR1,iR2), 0._dp, dp )
      enddo
    enddo
    !

    do iw = 1+my_id, input%n_omega, num_procs
      print*, "Frequency index: ", iw
      g0 = green_0_c(iw, wg, S, grid, Us, diff_large, g0_iR)
      g0__ = 0._dp
      do i = 1, nR
        do j = 1, nR
          g0__( (i-1)*S%nat3+1:i*S%nat3, (j-1)*S%nat3+1:j*S%nat3 ) = &
            g0(:,:,iR_large(i,j))
        enddo
      enddo
      !
      V__ = VK + wg%en(iw)**2 * VM
      call zgemm_N(N, g0__, V__, gV__)
      !
      I_gV__ = id_mat(N) - gV__ * (1-c)
      call invzmat(N, I_gV__)
      !
      call zgemm_N(N, c*V__, I_gV__, T__)
      !
      TRR__ = 0._dp
      do j = 1, nR
        do i = 1, nR
          TRR__(:,:,iR_large(i,j)) = TRR__(:,:,iR_large(i,j)) + &
            T__((i-1)*S%nat3+1:i*S%nat3, (j-1)*S%nat3+1:j*S%nat3)
        enddo
      enddo
      !
      do iq = 1, out_grid%nqtot
        do iR = 1, nR_large
          Tq(:,:,iq,iw) = Tq(:,:,iq,iw) + &
            TRR__(:,:,iR) * &
            phases_out(iR,iq)
        enddo
        call apply_sym_q(S, out_grid%xq(:,iq), Tq(:,:,iq,iw))
        Tq(:,:,iq,iw) = matmul(out_Us_c(:,:,iq), matmul(Tq(:,:,iq,iw), out_Us(:,:,iq)))
        do ibnd = 1, S%nat3
          self_energy(ibnd,iq,iw) = Tq(ibnd,ibnd,iq,iw)
        enddo
        ! call merge_degen(S%nat3, self_energy(:,iq,iw), out_freqs(:,iq))
      enddo
      if(input%calculation == 'dos' .or. input%calculation == 'test') then
        call tetra_from_self(S, out_grid, out_freqs, Tq(:,:,:,iw), wg_out%en(iw)**2, w_self)
        dos(iw) = sum(matmul(AIMAG(w_self), wg_out%qw)) * product(out_grid%n)
      endif
    enddo !iw
    call mpi_bsum(input%n_omega, dos)
    call mpi_bsum( S%nat3, out_grid%nqtot, input%n_omega, self_energy)
    call mpi_bsum( S%nat3, S%nat3, out_grid%nqtot, input%n_omega, Tq)

    !
    ! do iq = 1, out_grid%nqtot
    !   do ibnd = 1, S%nat3
    !     x = out_freqs(ibnd,iq)*input%n_omega/wg%max_f + 1
    !     x0 = INT(x)
    !     dx = x - x0
    !     lws(ibnd,iq) = (1.0_dp - dx) * Tq(ibnd,ibnd,iq,x0) + dx * Tq(ibnd,ibnd,iq,x0+1)
    !   enddo
    !   call merge_degen(S%nat3, lws(:,iq), out_freqs(:,iq))
    ! enddo
    ! !
    ! call write_file(out_freqs, lws, "spectral-full.dat", out_grid%type)
    select case(input%calculation)
     case('self')
      call write_self("self-fb.dat", wg%en, self_energy)
     case('spf-def')
      call write_spf_ndiag('spf-fb-ndiag.dat', wg%en, Tq, out_freqs)
      call write_spf('spf-fb.dat', wg%en, self_energy, out_freqs)
      call write_self("self-fb.dat", wg%en, self_energy)
     case('dos')
      call write_dos("dos-fb.dat", wg%en, dos)
     case('test')
      call write_spf_ndiag('spf-fb-ndiag-test.dat', wg%en, Tq(:,:,:1,:), out_freqs(:,:1))
      call write_spf('spf-fb-test.dat', wg%en, self_energy(:,:1,:), out_freqs(:,:1))
      call write_self("self-fb-test.dat", wg%en, self_energy(:,:1,:))
      call write_dos("dos-fb-test.dat", wg%en, dos)
     case default
      if (ionode) print*, "Unknown calculation type for fb: ", input%calculation
    end select
  contains
    subroutine t_symmetrize(TRR, iR_list, nR, TR)
      complex(dp), intent(in) :: TRR(:,:)
      integer, intent(in) :: iR_list(:,:)
      integer, intent(in) :: nR
      complex(dp), intent(out) :: TR(:,:,:)
      !
      TR = 0._dp
      do j = 1, nR
        do i = 1, nR
          TR(:,:,ir_list(i,j)) = TR(:,:,ir_list(i,j)) + &
            TRR((i-1)*S%nat3+1:i*S%nat3, (j-1)*S%nat3+1:j*S%nat3)
        enddo
      enddo
    end subroutine
    !
    function green_0_c(iw, wg, S, grid, U, diffs, g0_iR) result(g0)
      use thfftw, only: fft_1d_3d
      integer, intent(in) :: iw
      type(tetra_output), intent(in) :: wg
      type(ph_system_info), intent(in) :: S
      type(q_grid), intent(in) :: grid
      complex(dp), intent(in) :: U(S%nat3,S%nat3,grid%nqtot)
      real(dp), intent(in) :: diffs(:,:)
      integer, intent(in) :: g0_iR(:)
      !
      complex(dp), allocatable :: g0(:,:,:)
      complex(dp) :: g0_(S%nat3,S%nat3,size(diffs,2))
      complex(dp) :: g0_R(S%nat3, S%nat3, grid%nqtot)
      complex(dp) :: g0_q(S%nat3, S%nat3, grid%nqtot)
      integer :: iR, iq, iq_tetra, ibnd
      !
      allocate(g0(S%nat3, S%nat3, size(diffs,2)))
      g0_q = 0._dp
      do iq = 1, grid%nqtot
        ! FFT buffers follow natural Fortran order; tetrahedron weights use
        ! the grid order expected by the tetrahedron machinery.
        iq_tetra = v2index_n(index2v(iq, grid%n), grid%n)
        do ibnd = 1, S%nat3
          g0_q(:,:,iq_tetra) = g0_q(:,:,iq_tetra) + &
            wg%w(ibnd, wg%e(iq), iw) * outer_product(U(:,ibnd,iq))
        enddo
      enddo
      !
      call fft_1d_3d(S%nat3, grid%n, g0_q, g0_R, 1)
      !
      do iR = 1, size(diffs,2)
        g0(:,:,iR) = g0_R(:,:,g0_iR(iR))
      enddo
      !
    end function
    !
  end subroutine
  !
  subroutine full_born_periodic(S, input, fc2, fc2_sc, grid, sym_grid, out_grid)
    use defutils, only : flatten_RR_cmplx, flatten_RR_real
    type(ph_system_info), intent(in) :: S
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid, sym_grid, out_grid
    !
    integer :: iR1, iR2, iR, nR, nR_large, jR
    integer :: i, j, iq, ibnd, N
    real(dp) :: c
    integer :: R(3)
    complex(dp), allocatable :: g0(:,:,:,:)
    complex(dp), allocatable :: T(:,:,:)
    complex(dp) :: Tq(S%nat3, S%nat3,out_grid%nqtot, input%n_omega)
    complex(dp) :: self_energy(S%nat3, out_grid%nqtot, input%n_omega)
    integer :: iw
    type(tetra_output) :: wg, wg_out
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    real(dp) :: freqs(S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    complex(dp) :: out_Us_c(S%nat3, S%nat3, out_grid%nqtot)
    real(dp) :: out_freqs(S%nat3, out_grid%nqtot)
    !
    complex(dp), allocatable :: V__(:,:), gV__(:,:), g0__(:,:), T__(:,:), I_gV__(:,:)
    complex(dp), allocatable :: S__(:,:), S1__(:,:), Sg__(:,:), &
      I_Sg__(:,:), Gm__(:,:), GVS__(:,:), I_GVS__(:,:), G__(:,:)
    complex(dp) :: w_self(S%nat3, out_grid%nqtot)
    real(dp) :: dos(input%n_omega)
    !
    call set_wg(S, fc2, sym_grid, input, wg)
    if(input%calculation == 'dos' .or. input%calculation == 'test') &
      call set_wg(S, fc2, out_grid, input, wg_out)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)

    c = input%conc * size(fc2_sc%defects,2)
    do iq = 1, out_grid%nqtot
      out_Us_c(:,:,iq) = conjg(transpose(out_Us(:,:,iq)))
    enddo

    nR = product(fc2_sc%nq)

    N = S%nat3 * nR
    allocate(V__(N,N))
    allocate(g0__(N,N))
    allocate(Gm__, gV__, T__, I_gV__, S__, S1__, &
      Sg__, I_Sg__, GVS__, I_GVS__, G__, source=g0__)
    allocate(g0(S%nat3, S%nat3, nR, nR))
    !
    V__ = flatten_RR_real(fc2_sc%fc)
    gV__ = 0._dp
    T__ = 0._dp
    Tq = 0._dp
    dos = 0._dp
    self_energy = 0._dp
    !
    do iw = 3+my_id, input%n_omega, num_procs
      print*, "Periodic frequency index: ", iw
      g0 = 0._dp
      do jR = 1, nR
        do iR = 1, nR
          do iq = 1, grid%nqtot
            do ibnd = 1, S%nat3
              g0(:,:,iR,jR) = g0(:,:,iR,jR) + &
                wg%w(ibnd, wg%e(iq), iw) * outer_product(Us(:,ibnd,iq)) * &
                e_iqr(grid%xq(:,iq), fc2_sc%xR1(:,iR,jR) - fc2_sc%xR2(:,jR))
            enddo
          enddo
        enddo
      enddo
      g0__ = flatten_RR_cmplx(g0)
      !
      call zgemm_N(N, g0__, V__, gV__)
      !
      I_gV__ = id_mat(N) - gV__ * (1-c)
      call invzmat(N, I_gV__)
      !
      call zgemm_N(N, c*V__, I_gV__, T__)
      !
      do iq = 1, out_grid%nqtot
        do jR = 1, nR
          do iR = 1, nR
            Tq(:,:,iq,iw) = Tq(:,:,iq,iw) + &
              T__((iR-1)*S%nat3+1:iR*S%nat3, (jR-1)*S%nat3+1:jR*S%nat3) * &
              e_iqr(out_grid%xq(:,iq), fc2_sc%xR2(:,jR) - fc2_sc%xR1(:,iR,jR))
          enddo
        enddo
        call apply_sym_q(S, out_grid%xq(:,iq), Tq(:,:,iq,iw))
        Tq(:,:,iq,iw) = matmul(out_Us_c(:,:,iq), matmul(Tq(:,:,iq,iw), out_Us(:,:,iq)))
        do ibnd = 1, S%nat3
          self_energy(ibnd,iq,iw) = Tq(ibnd,ibnd,iq,iw)
        enddo
      enddo
      !
      if(input%calculation == 'dos' .or. input%calculation == 'test') then
        call tetra_from_self(S, out_grid, out_freqs, Tq(:,:,:,iw), wg_out%en(iw)**2, w_self)
        dos(iw) = sum(matmul(AIMAG(w_self), wg_out%qw)) * product(out_grid%n)
      endif
    enddo !iw
    call mpi_bsum(input%n_omega, dos)
    call mpi_bsum( S%nat3, out_grid%nqtot, input%n_omega, self_energy)
    call mpi_bsum( S%nat3, S%nat3, out_grid%nqtot, input%n_omega, Tq)

    select case(input%calculation)
     case('self')
      call write_self("self-fbp.dat", wg%en, self_energy)
     case('spf-def')
      call write_spf_ndiag('spf-fbp-ndiag.dat', wg%en, Tq, out_freqs)
      call write_spf('spf-fbp.dat', wg%en, self_energy, out_freqs)
      call write_self("self-fbp.dat", wg%en, self_energy)
     case('dos')
      call write_dos("dos-fbp.dat", wg%en, dos)
     case('test')
      call write_spf_ndiag('spf-fbp-ndiag-test.dat', wg%en, Tq(:,:,:1,:), out_freqs(:,:1))
      call write_spf('spf-fbp-test.dat', wg%en, self_energy(:,:1,:), out_freqs(:,:1))
      call write_self("self-fbp-test.dat", wg%en, self_energy(:,:1,:))
      call write_dos("dos-fbp-test.dat", wg%en, dos)
     case default
      if (ionode) print*, "Unknown calculation type for fbp: ", input%calculation
    end select
  contains
    function green_0_sum_c(iw, wg, S, grid, U, diffs) result(g0)
      integer, intent(in) :: iw
      type(tetra_output), intent(in) :: wg
      type(ph_system_info), intent(in) :: S
      type(q_grid), intent(in) :: grid
      complex(dp), intent(in) :: U(S%nat3,S%nat3,grid%nqtot)
      real(dp), intent(in) :: diffs(:,:)
      !
      complex(dp), allocatable :: g0(:,:,:)
      complex(dp) :: g0_q(S%nat3,S%nat3)
      integer :: iR, iq, ibnd
      !
      allocate(g0(S%nat3, S%nat3, size(diffs,2)))
      g0 = 0._dp
      do iq = 1, grid%nqtot
        g0_q = 0._dp
        do ibnd = 1, S%nat3
          g0_q = g0_q + &
            wg%w(ibnd, wg%e(iq), iw) * outer_product(U(:,ibnd,iq))
        enddo
        do iR = 1, size(diffs,2)
          g0(:,:,iR) = g0(:,:,iR) + g0_q * e_iqr(grid%xq(:,iq), diffs(:,iR))
        enddo
      enddo
    end function
    !
  end subroutine
  !
  subroutine find_where(R, list_of_R, idx) !res(idx)
    integer, intent(in) :: R(3)
    integer, intent(in) :: list_of_R(:,:)
    integer, intent(out) :: idx
    integer :: iR
    !
    idx = -1
    do iR = 1, size(list_of_R,2)
      if (all(R - list_of_R(:,iR) == 0)) then
        idx = iR
        exit
      endif
    enddo
  end subroutine
  !
  subroutine main_defect(S, fc2, fc2_sc, grid, sym_grid, out_grid, input)
    use constants, only: BOHR_RADIUS_CM, RY_TO_CMM1
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(forceconst2_grid), intent(in) :: fc2
    type(q_grid), intent(in) :: grid, sym_grid, out_grid
    type(code_input_type), intent(in) :: input
    type(tetra_output) :: w_in, w_out
    !
    real(dp), dimension(S%nat3, grid%nqtot) :: freqs
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs
    complex(dp), dimension(S%nat3, out_grid%nqtot) :: lws_out
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3, S%nat3, out_grid%nqtot) :: out_Us!, out_us_copy
    complex(dp), allocatable :: interp(:,:)
    complex(dp), allocatable, dimension(:,:,:,:) :: V
    real(dp) :: omega
    integer :: iq, iw, ibnd, nR, jq, jbnd
    complex(dp), allocatable, dimension(:,:,:) :: Vqqs_out, Vms
    character(15) :: filename
    real(dp) :: spectral_function(input%n_omega), c
    real(dp), dimension(input%n_omega) :: dos, dos0
    complex(dp) :: self_energy(S%nat3, out_grid%nqtot, input%n_omega)
    complex(dp) :: w_self(S%nat3, out_grid%nqtot)
    complex(dp) :: self_ndiag(S%nat3, S%nat3, out_grid%nqtot, input%n_omega)
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    complex(dp) :: phase_factor(grid%nqtot)
    ! complex(dp) :: vk_plus_vm
    complex(dp), dimension(S%nat3,S%nat3) :: U_dagger
    !
    !> full born quantities in real space
    !
    c = input%conc * size(fc2_sc%defects,2)
    nR   = product(fc2%nq)
    allocate(Vms(S%nat3,S%nat3,grid%nqtot))
    allocate(Vqqs_out(S%nat3,S%nat3,out_grid%nqtot))
    allocate(V(S%nat3,S%nat3,grid%nqtot,out_grid%nqtot))
    !
    !> input files are periodic, it can be changed
    !> SC: grid type   (big_nat3,big_nat3,1)
    !> UC: grid type   (nat3,nat3,nR)
    !> RR: new SC type (nat3,nat3,nR,nR)

    call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    !
    ! call freq_in_grid_degen(S, fc2, fc2_sc, out_grid, out_freqs, out_Us, out_freqs1)
    ! call freq_in_grid_degen(S, fc2, fc2_sc, grid, freqs, Us, freqs1)

    !
    call set_wg(S, fc2, sym_grid, input, w_in)
    call set_wg(S, fc2, out_grid, input, w_out)
    call print_message("end of tetra initialization")
    !
    do iq = 1, out_grid%nq
      call fc2_sc%r2q(out_grid%xq(:,iq))
      call fc2_sc%r2q(out_grid%xq(:,iq), Vqqs_out(:,:,iq))
    enddo
    !
    do iq = 1, out_grid%nqtot
      call fc2_sc%r2q(out_grid%xq(:,iq))
      U_dagger = conjg(transpose(out_Us(:,:,iq)))
      do jq = 1, grid%nqtot
        call fc2_sc%r2q(grid%xq(:,jq), V(:,:,jq,iq))
        V(:,:,jq,iq) = matmul(U_dagger, matmul(V(:,:,jq,iq), Us(:,:,jq)))
        ! pixel(jq,iq) = sum(V(:,:,jq))
      enddo
    enddo
    !
    self_energy = 0._dp
    self_ndiag = 0._dp
    w_self = 0._dp
    dos = 0._dp
    dos0 = 0._dp
    allocate(interp(S%nat3, sym_grid%nqtot))
    do iw = 1+my_id, input%n_omega, num_procs
      !
      do iq = 1, out_grid%nqtot
        ! if(mod(iw, input%n_omega/10) == 0) &
        ! print"(A,A,I3,A)", input%calculation, " progress ", NINT(100 * REAL(iw,DP) / input%n_omega), "%"
        omega = w_in%en(iw)
        if (trim(input%calculation) /= 'lw') interp = w_in%w(:,:,iw)
        if (iw > 1 .and. trim(input%calculation) == 'lw') exit
        lws_out(:,iq) = 0._dp
        do ibnd = 1, S%nat3
          lws_out(ibnd,iq) = lws_out(ibnd,iq) + REAL(braket(out_us(:,ibnd,iq), Vqqs_out(:,:,iq)), dp)
          if (trim(input%calculation) == 'lw') then
            omega = out_freqs(ibnd,iq)
            interp = interp1_matrix(w_in%w, omega*input%n_omega/w_in%max_f) !/ omega
          endif
          do jq = 1, grid%nqtot
            do jbnd = 1, S%nat3
              lws_out(ibnd,iq) = lws_out(ibnd,iq) + abs(V(ibnd,jbnd,jq,iq))**2 * interp(jbnd,w_in%e(jq))
            enddo
          enddo
        enddo
        ! call merge_degen(S%nat3, lws_out(:,iq), out_freqs(:,iq))
        !
        self_energy(:,iq,iw) = c*lws_out(:,iq)
        do ibnd = 1, S%nat3
          self_ndiag(ibnd,ibnd,iq,iw) = c * lws_out(ibnd,iq)
        enddo
      enddo ! iq
      if(input%calculation == 'dos' .or. input%calculation == 'test') then
        call tetra_from_self(S, out_grid, out_freqs, self_ndiag(:,:,:,iw), w_out%en(iw)**2, w_self)
        dos(iw) = sum(matmul(AIMAG(w_self), w_out%qw)) * product(out_grid%n)
        dos0(iw) = sum(matmul(AIMAG(w_out%w(:,:,iw)), w_out%qw)) * product(out_grid%n)
      end if
    enddo ! iw
    call mpi_bsum(S%nat3, out_grid%nqtot, input%n_omega, self_energy)
    call mpi_bsum(S%nat3, S%nat3, out_grid%nqtot, input%n_omega, self_ndiag)
    call mpi_bsum(input%n_omega, dos)
    call mpi_bsum(input%n_omega, dos0)
    !
    select case(trim(input%calculation))
     case("lw")
      if(trim(out_grid%type) == 'path') then
        write(filename, '(A4,I2.2,A4)') 'path', grid%n(1), '.dat'
      else
        write(filename, '(I2.2,I2.2,A4)') out_grid%n(1), grid%n(1), '.dat'
      endif
      call write_file(out_freqs, lws_out, filename, out_grid%type)
     case("spf-def")
      call write_spf('spf-1b.dat', w_in%en, self_energy, freqs)
      call write_self("self-1b.dat", w_in%en, self_energy)
     case("self")
      call write_self("self-1b.dat", w_in%en, self_energy)
     case("dos")
      if(trim(out_grid%type) /= 'grid' .and. trim(out_grid%type) /= 'simple' .and. ionode) &
        print*, "WARNING: DOS calculation in a non-grid q-point set: ", out_grid%type
      call write_dos("dos-1b.dat", w_out%en, dos)
      call write_dos("dos-0.dat", w_out%en, dos0)
     case("test")
      call write_spf('spf-1b-test.dat', w_in%en, self_energy(:,:1,:), freqs(:,:1))
      call write_self("self-1b-test.dat", w_in%en, self_energy(:,:1,:))
      call write_dos("dos-1b-test.dat", w_out%en, dos)
      call write_dos("dos-0-test.dat", w_out%en, dos0)
     case default
      if (ionode) print*, "Unknown calculation type for 1B: ", input%calculation
    end select
    call print_message("end of 1B calculation")
    !
  end subroutine
  !
  subroutine write_file(freqs_, lws_, filename_, type)
    real(dp), intent(in) :: freqs_(:,:)
    complex(dp), intent(in) :: lws_(:,:)
    character(*), intent(in) :: filename_
    character(*), intent(in) :: type
    !
    integer :: iq_, ibnd_
    !
    open(10, file=filename_)
    if (trim(type) == 'path') then
      do iq_ = 1, size(freqs_, 2)
        write(10, "(1000e14.5)") freqs_(:,iq_), -AIMAG(lws_(:,iq_))
      enddo
    else
      do iq_ = 1, size(freqs_, 2)
        do ibnd_ = 1, size(freqs_, 1)
          write(10, "(3E15.5,I2)") freqs_(ibnd_,iq_), lws_(ibnd_,iq_), ibnd_
        enddo
      enddo
    endif
    close(10)
  end subroutine
  !
  subroutine read_lw(filename, nq, nbnd, lws)
    character(*), intent(in) :: filename
    integer, intent(in) :: nq, nbnd
    real(dp), intent(out) :: lws(nbnd, nq)
    real(dp) :: data(nbnd*4+5,nq)
    !
    integer :: i
    character(100) :: line
    !

    open(278, file=filename, status='old', action='read')
    !
    read(278, '(A)') line
    read(278, '(A)') line
    !
    do i = 1, nq
      read(278, *) data(:,i)
    end do
    close(278)
    lws = data(5+nbnd+1:5+nbnd*2, :)
    !
  end subroutine
  !
  subroutine write_self(filename, en, self_energy)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: en(:)
    complex(dp), intent(in) :: self_energy(:,:,:)
    !
    integer :: iw, iq
    !
    if(.not. ionode) return
    open(10, file=filename)
    do iw = 1, size(self_energy,3)
      do iq = 1, size(self_energy,2)
        write(10, "(E20.8,X,I3,100E20.8)") &
          en(iw)*RY_TO_CMM1, iq, self_energy(:,iq,iw) * RY_TO_CMM1**2
      enddo
    enddo
    close(10)
  end subroutine
  !
  subroutine write_self_discrepancy(filename, en, self_energy, self_ref)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: en(:)
    complex(dp), intent(in) :: self_energy(:,:,:), self_ref(:,:,:)
    !
    integer :: iw
    real(dp) :: diff_norm, ref_norm, max_abs, max_ref
    !
    if(.not. ionode) return
    open(10, file=filename)
    do iw = 1, size(self_energy,3)
      max_abs = maxval(abs(self_energy(:,:,iw) - self_ref(:,:,iw)))
      max_ref = maxval(abs(self_ref(:,:,iw)))
      diff_norm = sqrt(sum(abs(self_energy(:,:,iw) - self_ref(:,:,iw))**2))
      ref_norm = sqrt(sum(abs(self_ref(:,:,iw))**2))
      write(10, "(5E20.8)") en(iw)*RY_TO_CMM1, &
        max_abs * RY_TO_CMM1**2, max_abs / max(max_ref, tiny(1._dp)), &
        diff_norm * RY_TO_CMM1**2, diff_norm / max(ref_norm, tiny(1._dp))
    enddo
    close(10)
  end subroutine
  !
  subroutine write_freq(filename, xq, freqs)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: xq(:, :)
    real(dp), intent(in) :: freqs(:, :)
    !
    integer :: iq
    !
    if(.not. ionode) return
    open(10, file=filename)
    do iq = 1, size(xq, 2)
      write(10, "(100E20.8)") xq(:, iq), freqs(:, iq)
    enddo
    close(10)
  end subroutine
  !
  subroutine write_spf_weights(filename, en, weights)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: en(:)
    complex(dp), intent(in) :: weights(:,:,:)
    !
    integer :: iw, iq
    !
    if(.not. ionode) return
    open(10, file=filename)
    do iw = 2, size(weights,3)
      do iq = 1, size(weights,2)
        write(10, "(E20.8,X,I3,100E20.8)") en(iw)*RY_TO_CMM1, iq, &
          -2 * en(iw) / pi * AIMAG(weights(:,iq,iw)) / RY_TO_CMM1
      enddo
    enddo
    close(10)
  end subroutine
  !
  subroutine write_spf(filename, en, self_energy, freqs)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: en(:)
    complex(dp), intent(in) :: self_energy(:,:,:)
    real(dp), intent(in) :: freqs(:,:)
    !
    real(dp), dimension(size(self_energy,1)) :: r,s
    integer :: iw, iq
    !
    if(.not. ionode) return
    open(10, file=filename)
    do iw = 2, size(self_energy,3)
      do iq = 1, size(self_energy,2)
        r = real(self_energy(:,iq,iw), dp)
        s = AIMAG(self_energy(:,iq,iw))
        write(10, "(E20.8,X,I3,100E20.8)") en(iw)*RY_TO_CMM1, iq, -2 / pi * en(iw) * s / &
          ((en(iw)**2 - freqs(:,iq)**2 - r)**2 + s**2) / RY_TO_CMM1
      enddo
    enddo
    close(10)
  end subroutine
  !
  subroutine write_spf_ndiag(filename, en, self_energy, freqs)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: en(:)
    complex(dp), intent(in) :: self_energy(:,:,:,:)
    real(dp), intent(in) :: freqs(:,:)
    !
    integer :: iq, i, iw
    complex(dp) :: M(size(self_energy,1), size(self_energy,2))
    real(dp) :: spf(size(self_energy,1))
    !
    if(.not. ionode) return
    open(10, file=filename)
    do iw = 2, size(en)
      do iq = 1, size(self_energy,3)
        M = - self_energy(:,:,iq,iw)
        do i = 1, size(self_energy,1)
          M(i,i) = M(i,i) + en(iw)**2 - freqs(i,iq)**2
        enddo
        call invzmat(size(M,1), M)
        do i = 1, size(self_energy,1)
          spf(i) = -2 / pi * en(iw) * aimag(M(i,i))
        enddo
        write(10, "(E20.8,X,I3,100E20.8)") en(iw)*RY_TO_CMM1, iq, spf / RY_TO_CMM1
      enddo
    enddo
    close(10)
  end subroutine
  !
  subroutine write_dos(filename, en, dos)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: en(:)
    real(dp), intent(in) :: dos(:)
    !
    integer :: iw
    !
    if(.not. ionode) return
    open(17, file=filename)
    do iw = 1, size(en)
      if(ionode) WRITE(17, "(2E20.8)") en(iw) * RY_TO_CMM1, &
        - dos(iw) / pi * 2 * en(iw) / RY_TO_CMM1
    enddo
    close(17)
  end subroutine
end module
