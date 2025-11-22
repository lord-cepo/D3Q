module defect
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  ! only: outer_product, freq_in_grid, interp1_matrix, &
  ! id_mat, braket, index2v, v2index, e_iqr, interp1_tns4, grid_vec
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi
  use quter_defect
  use functions, only: f_gauss
  use fc3_interpolate, only: forceconst3, sparse, d3_mixed, sum_R3
  use merge_degenerate, only: merge_degen
  use test_print, only: allclose
  use simtet, only: tetra_init_sym_cmplx, tetra_weights_green_cmplx
  use ph_velocity, only : velocity
  ! use tetra_raja
  !
  implicit none
  !
contains
  !
  ! subroutine frobenius_triangle(fc2_sc, S, filename)
  !   type(forceconst2_sc), intent(in) :: fc2_sc
  !   type(ph_system_info), intent(in) :: S
  !   character(len=*), intent(in) :: filename
  !   !
  !   integer :: ir1, ir2, na1, na2
  !   real(dp), dimension(3) :: R1, R2
  !   real(dp) :: perimeter
  !   !
  !   open(unit=10, file=filename, status="replace")
  !   do ir1 = 1, fc2_sc%n_R1
  !     do na1 = 1, S%nat
  !       R1 = fc2_sc%xR1(:,ir1) + S%tau(:,na1)
  !       do ir2 = 1, fc2_sc%n_R2(ir1)
  !         do na2 = 1, S%nat
  !           R2 = fc2_sc%xR2(:,ir2,ir1) + S%tau(:,na2)
  !           perimeter = norm2(R2 - R1) + &
  !             norm2(R1 - fc2_sc%taudef) + &
  !             norm2(R2 - fc2_sc%taudef)
  !           write(10,*) perimeter, fc2_sc%fc((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,ir2,ir1)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   close(10)

  ! end subroutine
  !
  subroutine full_V(S, out_grid, fc2_sc, V_flat, out_Us, grid, Us)
    type(ph_system_info), intent(in) :: S
    type(q_grid), target, intent(in) :: out_grid
    type(forceconst2_sc), intent(inout) :: fc2_sc
    complex(dp), allocatable, intent(out) :: V_flat(:,:)
    complex(dp), allocatable :: V(:,:,:,:)
    complex(dp), intent(out), target, optional :: out_us(:,:,:)
    type(q_grid), intent(in), target, optional :: grid
    complex(dp), intent(in), target, optional :: Us(:,:,:)
    !
    integer :: iq, jq
    complex(dp), dimension(S%nat3, S%nat3) :: D, D1
    type(q_grid), pointer :: grid_
    complex(dp), pointer :: us_(:,:,:)
    !
    if (present(grid)) then
      grid_ => grid
      if(present(out_Us)) us_ => us
    else
      grid_ => out_grid
      if (present(out_Us)) us_ => out_Us
    endif
    allocate(V(S%nat3, S%nat3, grid_%nqtot, out_grid%nqtot))
    !
    do iq = 1, grid_%nqtot
      call fc2_sc%r2q(grid_%xq(:,iq))
      do jq = 1, out_grid%nq
        call fc2_sc%r2q(out_grid%xq(:,jq), D)
        ! if (ALL(ABS(out_grid%xq(:,jq) - grid_%xq(:,iq)) < 1e-10_dp)) then
        !   call fftinterp_mat2(out_grid%xq(:,jq), S, fc2, D1)
        !   D = D - D1
        ! endif
        if (present(out_Us)) &
          D = matmul(conjg(transpose(Us_(:,:,iq))), matmul(D, out_Us(:,:,jq)))
        V(:,:,iq,jq) = D
      enddo
    enddo
    !
    call mpi_bsum(S%nat3, S%nat3, out_grid%nqtot, grid_%nqtot, V)
    V_flat = flatten_RR_cmplx(V)
  end subroutine
  !
  ! subroutine check_derivative_swap(fc, nat3)
  !   type(forceconst2_sc), intent(in) :: fc
  !   integer, intent(in) :: nat3
  !   !
  !   real(dp) :: mean, max_diff_fc2d, max_diff_rel
  !   integer :: ibnd, jbnd, ir1, ir2, ir1p, ir2p
  !   integer(dp) :: R(3)
  !   logical :: ok
  !   !
  !   max_diff_fc2d = 0._dp
  !   max_diff_rel = 0._dp
  !   do ir1 = 1, fc%n_R1
  !     do ir2 = 1, fc%n_R2(ir1)
  !       R = fc%yR2(:,ir2,ir1)
  !       ok = .false.
  !       do ir1p = 1, fc%n_R1
  !         if (ALL(fc%yR1(:,ir1p) == R)) then
  !           ok = .true.
  !           exit
  !         endif
  !       enddo
  !       if (.not. ok) then
  !         do ir1p = 1, fc%n_R1
  !           print"(3I2)", fc%yR1(:,ir1p)
  !         enddo
  !         print*, R
  !         call errore("check_derivative_swap", "R1 not found in R2", 1)
  !       endif
  !       R = fc%yR1(:,ir1)
  !       ok = .false.
  !       do ir2p = 1, fc%n_R2(ir1p)
  !         if (ALL(fc%yR2(:,ir2p,ir1p) == R)) then
  !           ok = .true.
  !           exit
  !         endif
  !       enddo
  !       if (.not. ok) &
  !         call errore("check_derivative_swap", "R2 not found in R1", 1)
  !       do ibnd = 1, nat3
  !         do jbnd = ibnd+1, nat3
  !           mean = (fc%fc(ibnd,jbnd,ir2,ir1) + fc%fc(jbnd,ibnd,ir2p,ir1p)) / 2.0_dp
  !           if (ABS(fc%fc(ibnd,jbnd,ir2,ir1) - fc%fc(jbnd,ibnd,ir2p,ir1p)) > max_diff_fc2d) then
  !             max_diff_fc2d = ABS(fc%fc(ibnd,jbnd,ir2,ir1) - fc%fc(jbnd,ibnd,ir2p,ir1p))
  !             max_diff_rel = max_diff_fc2d / ABS(mean)
  !           endif
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   if (ionode) print"(A,E14.4)", "max_diff_fc2d", max_diff_fc2d
  !   if (ionode) print"(A,E14.4)", "max_diff_rel", max_diff_rel
  ! end subroutine
  !
  !
  subroutine create_diffs(fc2_sc, diffs)
    type(forceconst2_sc), intent(in) :: fc2_sc
    integer, intent(out), pointer :: diffs(:,:)
    !
    integer :: ir1, ir2, iR, nR
    logical :: is_new
    !
    allocate(diffs(3,1))
    diffs(:,1) = fc2_sc%yR1(:,1,1) - fc2_sc%yR2(:,1)
    nR = 1
    !
    do ir2 = 1, fc2_sc%n_R2
      ! call find_where(fc2_sc%yR2(:,ir2), R_list, j)
      do ir1 = 1, fc2_sc%n_R1(ir2)
        ! call find_where(fc2_sc%yR1(:,ir1,ir2), R_list, i)
        ! V((i-1)*S%nat3+1:i*S%nat3,(j-1)*S%nat3+1:j*S%nat3) = &
        ! fc2_sc%fc(:,:,ir1,ir2)
        is_new = .true.
        do iR = 1, nR
          if (all( (fc2_sc%yR1(:,ir1,ir2))-fc2_sc%yR2(:,ir2) - diffs(:,iR) == 0)) then
            is_new = .false.
            exit
          endif
        enddo
        if(is_new) call enlarge_R( (fc2_sc%yR1(:,ir1,ir2))-fc2_sc%yR2(:,ir2), diffs, nR)
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
  subroutine full_born_center(S, input, fc2, fc2_sc, grid, out_grid)
    use quter_defect, only: fc_sc2RR, allocate_fc2_sc, minimal_image
    type(ph_system_info), intent(in) :: S
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(in) :: fc2_sc
    type(q_grid), intent(in) :: grid, out_grid
    !
    integer :: iR1, iR2, iR, nR, Ri, Rj, Rk, nR_large
    integer :: i, j, iq, na1, na2, x0, ibnd, ii
    real(dp) :: x, dx
    integer, pointer :: R_list(:,:), diff_list(:,:)
    integer :: R(3)
    real(dp), allocatable :: Rij_cart(:,:,:)
    real(dp), allocatable :: diffs(:,:)
    logical :: is_new
    ! real(dp), allocatable :: V(:,:)
    complex(dp), allocatable :: g0(:,:,:)
    ! complex(dp), allocatable :: g0_RR(:,:)
    complex(dp) :: T_RR(S%nat3,S%nat3,maxval(fc2_sc%n_R1),fc2_sc%n_R2)
    complex(dp) :: g0_RR(S%nat3,S%nat3,maxval(fc2_sc%n_R1),fc2_sc%n_R2)
    complex(dp) :: gV(S%nat3,S%nat3,maxval(fc2_sc%n_R1),fc2_sc%n_R2)
    complex(dp), allocatable :: T(:,:,:)
    complex(dp) :: Tq(S%nat3, S%nat3, out_grid%nqtot, input%n_omega)
    integer :: iw
    type(tetra_output) :: wg
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    real(dp) :: freqs(S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    complex(dp) :: out_Us_c(S%nat3, S%nat3, out_grid%nqtot)
    real(dp) :: out_freqs(S%nat3, out_grid%nqtot)
    complex(dp) :: lws(S%nat3, out_grid%nqtot), lws1(S%nat3)
    integer, dimension(maxval(fc2_sc%n_R1),fc2_sc%n_R2):: iR_diff
    integer, dimension(fc2_sc%n_R2,fc2_sc%n_R2):: iR21
    integer, dimension(maxval(fc2_sc%n_R1),maxval(fc2_sc%n_R1)):: iR12
    complex(dp) :: mat(S%nat3,S%nat3)
    complex(dp) :: V(S%nat3,S%nat3,maxval(fc2_sc%n_R1),fc2_sc%n_R2)
    integer, pointer :: diff_list_large(:,:)
    real(dp), allocatable :: diff_large(:,:)
    integer, allocatable :: iR_large(:,:)
    complex(dp), allocatable :: phases_out(:,:)
    complex(dp), allocatable :: TRR__(:,:,:)
    !
    complex(dp), allocatable :: V__(:,:), gV__(:,:), g__(:,:), T__(:,:), I_gV__(:,:)
    integer :: i_list(maxval(fc2_sc%n_R1),fc2_sc%n_R2)
    integer :: j_list(fc2_sc%n_R2)
    character(10) :: filename
    complex(dp) :: Tij
    integer, allocatable :: g0_iR(:)
    !
    call set_wg(S, fc2, grid, input%n_omega, 1.5_dp, wg)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)
    !
    do iq = 1, out_grid%nqtot
      out_Us_c(:,:,iq) = conjg(transpose(out_Us(:,:,iq)))
    enddo
    !
    allocate(R_list(3,1))
    R_list(:,1) = fc2_sc%yR2(:,1)
    nR = 1

    do iR2 = 1, fc2_sc%n_R2
      is_new = .true.
      do iR = 1, nR
        if (all(fc2_sc%yR2(:,iR2) - R_list(:,iR) == 0)) then
          is_new = .false.
          exit
        endif
      enddo
      if(is_new) call enlarge_R(fc2_sc%yR2(:,iR2), R_list, nR)
      do iR1 = 1, fc2_sc%n_R1(iR2)
        is_new = .true.
        do iR = 1, nR
          if (all(fc2_sc%yR1(:,iR1,iR2) - R_list(:,iR) == 0)) then
            is_new = .false.
            exit
          endif
        enddo
        if(is_new) call enlarge_R(fc2_sc%yR1(:,iR1,iR2), R_list, nR)
      enddo
    enddo
    print*, "Number of unique R vectors: ", nR
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
          endif
          R = diff_list_large(:,iR)
          do ii = 1, 3
            if (R(ii) < 0) R(ii) = R(ii) + grid%n(ii)
          enddo
          ! if(any(R < 0 .or. R > grid%n)) stop diff_list_large1
          g0_iR(iR) = v2index_n(R, grid%n)
          phases_out(iR_large(i,j),iq) = e_iqr(out_grid%xq(:,iq), -Rij_cart(:,i,j))
        enddo
      enddo
    enddo
    !
    allocate(diff_large(3,nR_large))
    diff_large = cryst2cart(real(diff_list_large,dp), S%at, 1)
    print*, "Number of unique large diff vectors: ", nR_large
    !

    allocate(V__(S%nat3*nR,S%nat3*nR))
    allocate(gV__(S%nat3*nR,S%nat3*nR))
    allocate(g__(S%nat3*nR,S%nat3*nR))
    allocate(T__(S%nat3*nR,S%nat3*nR))
    allocate(I_gV__(S%nat3*nR,S%nat3*nR))
    !
    V__ = 0._dp
    gV__ = 0._dp
    T__ = 0._dp
    Tq = 0._dp
    do iR2 = 1, fc2_sc%n_R2
      do iR1 = 1, fc2_sc%n_R1(iR2)
        call find_where(fc2_sc%yR1(:,iR1,iR2)-fc2_sc%yR2(:,iR2), diff_list, iR_diff(iR1,iR2))
      enddo
    enddo

    allocate(T(S%nat3,S%nat3,size(diff_list,2)))
    !
    do iR2 = 1, fc2_sc%n_R2
      call find_where(fc2_sc%yR2(:,iR2), R_list, j)
      j_list(iR2) = j
      do iR1 = 1, fc2_sc%n_R1(iR2)
        call find_where(fc2_sc%yR1(:,iR1,iR2), R_list, i)
        i_list(iR1,iR2) = i
        V__( (i-1)*S%nat3+1:i*S%nat3, (j-1)*S%nat3+1:j*S%nat3 ) = &
          cmplx( fc2_sc%fc(:,:,iR1,iR2), 0._dp, dp )
      enddo
    enddo
    !
    do iw = 1, input%n_omega
      print*, "Frequency index: ", iw
      g0 = green_0_c(iw, wg, S, grid, Us, diff_large, g0_iR)
      g__ = 0._dp
      do i = 1, nR
        do j = 1, nR
          g__( (i-1)*S%nat3+1:i*S%nat3, (j-1)*S%nat3+1:j*S%nat3 ) = &
          ! g0(:,:,iR_diff(iR1,iR2))
            g0(:,:,iR_large(i,j))
        enddo
      enddo
      !
      gV__ = 0._dp
      call zgemm_N(S%nat3*nR, g__, V__, gV__)
      !
      ! I_gV__ = id_mat(S%nat3*nR) - gV__
      ! call invzmat(S%nat3*nR, I_gV__)
      T__ = 0._dp
      call zgemm_N(S%nat3*nR, V__, gV__, T__)
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
        Tq(:,:,iq,iw) = matmul(out_Us_c(:,:,iq), matmul(Tq(:,:,iq,iw), out_Us(:,:,iq)))
      enddo
    enddo
    ! open(10, file="spectral-full.dat")
    ! open(11, file="weights.dat")
    do iq = 1, out_grid%nqtot
      do ibnd = 1, S%nat3
        x = out_freqs(ibnd,iq)*input%n_omega/wg%max_f + 1
        x0 = INT(x)
        dx = x - x0
        lws(ibnd,iq) = (1.0_dp - dx) * Tq(ibnd,ibnd,iq,x0) + dx * Tq(ibnd,ibnd,iq,x0+1)
        ! lws1(ibnd) = (1.0_dp - dx) * wg%w(ibnd, iq, x0) + dx * wg%w(ibnd, iq, x0+1)
        ! do iw = 1, input%n_omega
        !   write(11, "(3E20.8)") wg%en(iw), wg%w(ibnd, iq, iw)
        ! enddo
      enddo
      call merge_degen(S%nat3, lws(:,iq), out_freqs(:,iq))
    enddo
    !
    call write_file(out_freqs, lws, "spectral-full.dat", out_grid%type)
    open(103, file="spectral-full-re.dat")
    do iw = 1, input%n_omega
      do iq = 1, out_grid%nqtot
        do ibnd = 1, 1
          lws(ibnd,iq) = Tq(ibnd,ibnd,iq,iw)
        enddo
      enddo
      write(103, "(1200E20.8)") aimag(lws(1,:))
    enddo
    close(103)
    ! close(10)
    ! close(11)
    !
  contains
    subroutine zgemm_N(N, A, B, C)
      integer, intent(in) :: N
      complex(dp), intent(in) :: A(N,N), B(N,N)
      complex(dp), intent(inout) :: C(N,N)
      !
      complex(dp) :: one
      one = (1.0_dp, 0.0_dp)
      call zgemm('N','N', N, N, N, one, A, N, B, N, one, C, N)
    end subroutine
    !
    subroutine t_symmetrize(fc2_sc, diff_list, T_RR, T_R)
      type(forceconst2_sc), intent(in) :: fc2_sc
      integer, intent(in) :: diff_list(:,:)
      complex(dp), intent(in) :: T_RR(:,:,:,:)
      complex(dp), intent(out) :: T_R(size(T_RR,1),size(T_RR,1),size(diff_list,2))
      !
      integer :: iR1, iR2, iR, niR(size(diff_list,2))
      !
      T_R = 0._dp
      niR = 0
      do iR2 = 1, fc2_sc%n_R2
        do iR1 = 1, fc2_sc%n_R1(iR2)
          call find_where(fc2_sc%yR1(:,iR1,iR2)-fc2_sc%yR2(:,iR2), diff_list, iR)
          niR(iR) = niR(iR) + 1
          T_R(:,:,iR) = T_R(:,:,iR) + T_RR(:,:,iR1,iR2)
        enddo
      enddo
      !
      do iR = 1, size(diff_list,2)
        T_R(:,:,iR) = T_R(:,:,iR) / real(niR(iR), dp)
      enddo
    end subroutine
    !
    subroutine t_symmetrize_1(fc2_sc, R_list, diff_list, T_RR, nat3, T_R)
      type(forceconst2_sc), intent(in) :: fc2_sc
      integer, intent(in) :: R_list(:,:)
      integer, intent(in) :: diff_list(:,:)
      complex(dp), intent(in) :: T_RR(:,:)
      integer, intent(in) :: nat3
      complex(dp), intent(out) :: T_R(nat3,nat3,size(diff_list,2))
      !
      integer :: iR1, iR2, iR
      !
      do iR2 = 1, size(R_list,2)
        do iR1 = 1, size(R_list,1)
          call find_where(R_list(:,iR1)-R_list(:,iR2), diff_list, iR)
          if(iR == -1) cycle
          T_R(:,:,iR) = T_R(:,:,iR) + T_RR((iR1-1)*nat3+1:iR1*nat3,(iR2-1)*nat3+1:iR2*nat3)
        enddo
      enddo
    end subroutine
    !
    function green_0_c(iw, wg, S, grid, U, diffs, g0_iR) result(g0)
      use fftw_prova, only: fft_1d_3d
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
      integer :: iR, iq, ibnd
      !
      allocate(g0(S%nat3, S%nat3, size(diffs,2)))
      g0_q = 0._dp
      do iq = 1, grid%nqtot
        do ibnd = 1, S%nat3
          g0_q(:,:,iq) = g0_q(:,:,iq) + &
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
      ! g0_ = 0._dp
      ! do iR = 1, size(diffs,2)
      !   do iq = 1, grid%nqtot
      !     g0_(:,:,iR) = g0_(:,:,iR) + &
      !       g0_q(:,:,iq)  * e_iqr(grid%xq(:,iq), diffs(:,iR))
      !   enddo
      !   ! if(any(abs(g0_(:,:,iR) - g0(:,:,iR)) > 1e-5)) then
      !   !   print*, g0(:,:,iR)
      !   !   print*, "----"
      !   !   print*, g0_(:,:,iR)
      !   !   stop 1
      !   ! endif
      ! enddo
      ! !
      ! do iR = 1, size(g0_R,3)
      !   if(abs(sum(g0_R(:,:,iR) - g0_(:,:,1))) < 1e-1_dp) &
      !     print*, iR
      ! enddo
    end function
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
  end subroutine
  !
  function green_0(sc_grid, iw, wg, S, grid, U, diffs, n_weights, weights)
    integer, intent(in) :: sc_grid(3)
    integer, intent(in) :: iw
    type(tetra_output), intent(in) :: wg
    type(ph_system_info), intent(in) :: S
    real(dp), allocatable, intent(in) :: diffs(:,:,:,:,:), weights(:,:,:,:)
    integer, allocatable, intent(in) :: n_weights(:,:,:)
    type(q_grid), intent(in) :: grid
    complex(dp), intent(in) :: U(S%nat3,S%nat3,grid%nqtot)
    !<grid%S%nat3>
    complex(dp) :: green_0(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: proj(S%nat3,S%nat3,product(sc_grid))
    integer :: iR, iq, ibnd, nR, na1, na2, j, iqf
    real(dp), allocatable :: R_grid(:,:), q_grid(:,:)

    !
    nR = product(sc_grid)
    R_grid = grid_vec_cart(sc_grid, S%at)
    q_grid = grid_vec_cart(sc_grid, S%bg) / sc_grid(1)
    !
    ! do concurrent(Ri=1:nR, Rj=1:nR, iq=1:grid%nqtot, ibnd=1:S%nat3)
    proj = 0._dp
    do iq = 1, grid%nqtot
      do ibnd = 1, S%nat3
        proj(:,:,iq) = outer_product(U(:,ibnd,iq)) * wg%w(ibnd, wg%e(iq), iw)
      enddo
    enddo
    !
    green_0 = 0._dp
    ! do concurrent(na1=1:S%nat, na2=1:S%nat, ir=1:nR)
    !   do j = 1, n_weights(iR,na1,na2)
    !     do iq = 1, grid%nqtot
    !       green_0(j,(na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iR) = &
    !         green_0(j,(na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iR) + &
    !         proj((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq) * &
    !         e_iqr(grid%xq(:,iq), -diffs(:,j,iR,na1,na2)) * &
    !         weights(j,iR,na1,na2)
    !     enddo
    !   enddo
    ! enddo
    !
    do iR = 1, nR
      do iq = 1, grid%nqtot
        green_0(:, :, iR) = green_0(:, :, iR) + &
          proj(:, :, iq) * e_iqr(grid%xq(:,iq), R_grid(:,iR)) / nR
      enddo
    enddo
  end function
  !
  subroutine full_born_real(S, Sd, fc2, input, sym_grid, grid, V_sc)
    use quter_defect, only: fc_sc2RR, allocate_fc2_sc, minimal_image
    type(ph_system_info), intent(in) :: S, Sd
    type(forceconst2_grid), intent(in) :: fc2
    type(code_input_type), intent(in) :: input
    type(q_grid), intent(in) :: grid, sym_grid
    real(dp), intent(in) :: V_sc(:,:)
    !
    type(tetra_output) :: wg
    real(dp) :: freqs(S%nat3, grid%nqtot), sym_freqs(S%nat3, sym_grid%nqtot)
    complex(dp), dimension(S%nat3, S%nat3, grid%nqtot) :: Us
    complex(dp), dimension(S%nat3, S%nat3, sym_grid%nqtot) :: sym_Us, sym_Us_c
    complex(dp) :: G_R(S%nat3,S%nat3,product(input%sc_grid))
    complex(dp) :: G_q(S%nat3,S%nat3,sym_grid%nqtot,input%n_omega)
    complex(dp), dimension(size(V_sc,1), size(V_sc,2)) :: G, I_VG, gV, g0_RR
    complex(dp) :: g0(S%nat3,S%nat3,product(input%sc_grid))
    integer :: iw, iq, iR, na1, na2, i, j, k, ibnd, x0, jR, kR, jbnd, kbnd
    real(dp) :: x, dx
    integer :: iR_diff(product(input%sc_grid),product(input%sc_grid)), nR
    real(dp) :: A(input%n_omega)
    complex(dp) :: lws(S%nat3)
    complex(dp) :: dos
    !
    real(dp), allocatable :: diffs(:,:,:,:,:), weights(:,:,:,:)
    integer, allocatable:: n_weights(:,:,:)
    complex(dp) :: w2(size(V_sc,1))
    real(dp), allocatable :: R_grid(:,:)
    !
    call set_wg(S, fc2, sym_grid, input%n_omega, 1.5_dp, wg)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    call freq_in_grid(S, fc2, sym_grid, sym_freqs, sym_Us)
    open(10, file="partial_dos.dat")
    do iq = 1, sym_grid%nqtot
      do ibnd = 1, S%nat3
        x = sym_freqs(ibnd,iq)*input%n_omega/wg%max_f + 1
        x0 = INT(x)
        dx = x - x0
        write(10, "(3E20.8)") sym_freqs(ibnd,iq), wg%w(ibnd, iq,x0)*(1.0_dp-dx) + &
          wg%w(ibnd, iq,x0+1)*dx
      enddo
    enddo
    close(10)
    open(10, file="dos.dat")
    do iw = 1, input%n_omega
      write(10, "(3E20.8)") wg%en(iw), sum(matmul(wg%w(:, :, iw), wg%qw(:)))
    enddo
    close(10)
    !
    R_grid = grid_vec_cart(input%sc_grid, S%at)
    nR = product(input%sc_grid)
    A = 0._dp
    !
    do concurrent(iR=1:nR, jR=1:nR)
      iR_diff(iR,jR) = v2index(bz2simple(index2v(iR, input%sc_grid)-index2v(jR, &
        input%sc_grid), input%sc_grid), input%sc_grid)
    enddo

    !
    call minimal_image_2(S, input%sc_grid, diffs, n_weights, weights)
    do iw = 1, input%n_omega
      print*, iw
      g0 = green_0(input%sc_grid, iw, wg, S, grid, Us, diffs, n_weights, weights)
      ! print*, any(aimag(g0) > 1e-10_dp)
      ! print*, "g0"
      ! G = g0 + matmul(g0, matmul(matmul(V_sc, I_VG), g0))
      ! print*, "G"
      !
      gV = 0._dp
      do concurrent(iR=1:nR, jR=1:nR)
        g0_RR((iR-1)*S%nat3+1:iR*S%nat3,(jR-1)*S%nat3+1:jR*S%nat3) = &
          g0(:,:,iR_diff(iR,jR))
      enddo
      !
      ! I_VG = id_mat(size(V_sc,1)) - gV
      ! call invzmat(size(I_VG, 1), I_VG)
      !
      G = matmul(V_sc, matmul(g0_RR, V_sc))
      ! do concurrent(iR=1:nR, jR=1:nR)
      !   gV((iR-1)*S%nat3+1:iR*S%nat3,(jR-1)*S%nat3+1:jR*S%nat3) = &
      !     g0(:,:,iR_diff(iR,jR))
      ! enddo
      ! !
      ! G = gV!matmul(V_sc, gV)
      ! call mat2_diag(size(G,1), G, w2)
      ! ! print*, aimag(w2)
      ! print*, any(aimag(w2) > 1e-8_dp)
      !
      call t_symmetrize(input%sc_grid, fc_sc2RR_cmplx(input%sc_grid, S, Sd, G), G_R)
      do iq = 1, sym_grid%nqtot
        ! do iR = 1, nR
        !   do na1 = 1, S%nat
        !     do na2 = 1, S%nat
        !       do j = 1, n_weights(iR,na1,na2)
        !         G_q(  (na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq,iw) = &
        !           G_q((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq,iw) + &
        !           g0( (na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iR) + &
        !           e_iqr(sym_grid%xq(:,iq), diffs(:,j,iR,na1,na2)) * weights(j,iR,na1,na2)
        !       enddo
        !     enddo
        !   enddo
        ! enddo
        do iR = 1, nR
          G_q(:,:,iq,iw) = G_q(:,:,iq,iw) + &
            G_R(:,:,iR) * e_iqr(sym_grid%xq(:,iq), -R_grid(:,iR))
        enddo
        if(iw == 1) sym_Us_c(:,:,iq) = conjg(transpose(sym_Us(:,:,iq)))
        G_q(:,:,iq,iw) = matmul(sym_Us_c(:,:,iq), matmul(G_q(:,:,iq,iw), sym_Us(:,:,iq)))
        ! do j = 1, S%nat3
        !   A(iw) = A(iw) + AIMAG(G_q(j,j,iq,iw))
        ! enddo
      enddo
    enddo
    !
    open(10, file="spectral-full.dat")
    ! do iw = 1, input%n_omega
    !   write(10, "(2E20.8)") wg%en(iw), A(iw)
    ! enddo
    ! do iw = 1, input%n_omega
    !   dos = 0._dp
    !   do iq = 1, sym_grid%nqtot
    !     do ibnd = 1, S%nat3
    !       dos = dos + G_q(ibnd,ibnd,iq,iw) * wg%qw(iq)
    !     enddo
    !   enddo
    !   write(10, "(3E20.8)") wg%en(iw), dos
    ! enddo
    do iq = 1, sym_grid%nqtot
      do ibnd = 1, S%nat3
        x = sym_freqs(ibnd,iq)*input%n_omega/wg%max_f + 1
        x0 = INT(x)
        dx = x - x0
        lws(ibnd) = (1.0_dp - dx) * G_q(ibnd,ibnd,iq,x0) + dx * G_q(ibnd,ibnd,iq,x0+1)
      enddo
      call merge_degen(S%nat3, lws, sym_freqs(:,iq))
      do ibnd = 1, S%nat3
        write(10, "(4E208.8)") sym_freqs(ibnd,iq), lws(ibnd)/sym_freqs(ibnd,iq)
      enddo
    enddo
    close(10)
    !
  contains
    subroutine t_symmetrize(sc_grid, T_RR, T_R)
      complex(dp), intent(in) :: T_RR(:,:,:,:)
      complex(dp), intent(out) :: T_R(size(T_RR,1), size(T_RR,2), size(T_RR,3))
      integer, intent(in) :: sc_grid(3)
      !
      integer :: R1, T, R_sum(3)
      T_R = 0._dp
      do R1 = 1, size(T_RR,4)
        do T = 1, size(T_RR,3)
          R_sum = bz2simple(index2v(T, sc_grid) + index2v(R1,sc_grid), sc_grid)
          T_R(:,:,R1) = T_R(:,:,R1) + T_RR(:,:,v2index(R_sum, sc_grid),T)
        enddo
      enddo
      T_R = T_R / size(T_RR,4)
    end subroutine
  end subroutine
  !
  subroutine full_born(NR, S, fc2, fc2_sc, input, grid, out_grid)
    use quter_defect, only: fc_sc2RR, allocate_fc2_sc, minimal_image
    use q_grids, only : setup_simple_grid
    integer, intent(in) :: NR
    type(ph_system_info), intent(in) :: S
    ! complex(dp), intent(in) :: V_sc(:,:)
    type(forceconst2_sc), intent(in) :: fc2_sc
    type(forceconst2_grid), intent(in) :: fc2
    type(code_input_type), intent(in) :: input
    type(q_grid), intent(in) :: grid, out_grid
    !
    type(q_grid) :: q_grid_sc
    complex(dp), dimension(S%nat3,S%nat3,out_grid%nqtot,input%n_omega) :: Tq
    !
    type(tetra_output) :: wg
    real(dp) :: freqs(S%nat3, grid%nqtot)
    real(dp) :: freqs_sc(S%nat3, nR)
    complex(dp), dimension(S%nat3, S%nat3, nR) :: Us_sc, Us_sc_C
    real(dp) :: out_freqs(S%nat3, out_grid%nqtot)
    complex(dp), dimension(S%nat3,S%nat3,out_grid%nqtot) :: out_Us, out_Us_C
    real(dp),allocatable :: R(:,:), q(:,:)
    complex(dp), dimension(S%nat3, S%nat3, grid%nqtot) :: Us, Us_C
    complex(dp), dimension(NR*S%nat3, NR*S%nat3) :: I_VG, gV, Vt, T_MAT
    complex(dp), dimension(NR*S%nat3) :: G
    complex(dp) :: g0(S%nat3,S%nat3,nR)
    complex(dp) :: g0_R(S%nat3,S%nat3,grid%nqtot)
    complex(dp), allocatable :: TR(:,:,:)
    integer :: iw, iq, j, ibnd, jbnd, iR, na1, na2, x0, jR, jq
    real(dp) :: dx, x
    complex(dp) :: lws(S%nat3), lws1(S%nat3)
    real(dp) :: A(input%n_omega)
    integer :: R_diff(3)
    complex(dp) :: D(S%nat3,S%nat3)
    complex(dp) :: dos
    !
    real(dp), allocatable :: diffs(:,:)
    integer, pointer :: diffs_int(:,:)
    ! real(dp), allocatable :: diffs(:,:,:,:,:), weights(:,:,:,:)
    ! integer, allocatable:: n_weights(:,:,:)
    !
    call setup_simple_grid(S%bg, input%sc_grid(1), input%sc_grid(2), &
      input%sc_grid(3), q_grid_sc)
    call set_wg(S, fc2, grid, input%n_omega, 1.5_dp, wg)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    call freq_in_grid(S, fc2, q_grid_sc, freqs_sc, Us_sc)
    call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)
    !
    call create_diffs(fc2_sc, diffs_int)
    allocate(diffs(3, size(diffs_int,2)))
    allocate(TR(S%nat3,S%nat3,size(diffs_int,2)))
    diffs = cryst2cart(real(diffs_int,dp), S%at, 1)
    !
    ! call full_V(S, q_grid_sc, fc2_sc, V, Us)
    ! V = (V + transpose(conjg(V))) / 2
    ! Vt = transpose(conjg(V_sc))
    do iq = 1, grid%nqtot
      Us_C(:,:,iq) = transpose(conjg(Us(:,:,iq)))
    enddo
    do iq = 1, out_grid%nqtot
      out_Us_C(:,:,iq) = transpose(conjg(out_Us(:,:,iq)))
    enddo
    do iq = 1, nR
      Us_sc_C(:,:,iq) = transpose(conjg(Us_sc(:,:,iq)))
    enddo
    !
    open(10, file="dos.dat")
    do iw = 1, input%n_omega
      write(10, "(3E20.8)") wg%en(iw), sum(matmul(wg%w(:, :, iw), wg%qw(:)))
    enddo
    close(10)
    !
    ! T = 0._dp
    Tq = 0._dp
    R = grid_vec_cart(input%sc_grid, S%at)
    q = grid_vec_cart(input%sc_grid, S%bg)
    do iq = 1, nR
      q(:,iq) = q(:,iq) / input%sc_grid
    enddo
    !
    ! call minimal_image_2(S, input%sc_grid, diffs, n_weights, weights)
    ! open(155, file="fc.dat")
    do iw = 1, input%n_omega
      ! g0 = 0._dp
      print*, iw
      !
      g0 = green_0_q(input%sc_grid, iw, wg, S, grid, Us, diffs)
      ! g0_R = green_0(grid%n, iw, wg, S, grid, Us, diffs, n_weights, weights)
      ! do concurrent(na1=1:S%nat, na2=1:S%nat, iR=1:grid%nqtot)
      !   distance = norm2(cryst2cart(, trmat, iflag) + S%tau(:,na1) - S%tau(:,na2))
      !     g0_R((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq) = &
      !       g0_R((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq) + &
      !       g0( (na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq) * &
      !       e_iqr(grid%xq(:,iq), -diffs(:,j,iq,na1,na2)) * &
      !       weights(j,iq,na1,na2)
      !   enddo

      !
      ! do iq = 1, nR
      !   do jq = 1, nR
      !     gV((iq-1)*S%nat3+1:iq*S%nat3,(jq-1)*S%nat3+1:jq*S%nat3) = &
      !       matmul(g0(:,:,iq), V_sc((iq-1)*S%nat3+1:iq*S%nat3,(jq-1)*S%nat3+1:jq*S%nat3))
      !   enddo
      ! enddo
      ! print*, any(aimag(matmul(V,gV)) > 1e-10_dp)
      ! I_VG = id_mat(N) - gV
      ! call invzmat(N, I_VG)
      ! T_mat = matmul(V, gV)
      ! do concurrent(ibnd = 1:S%nat3, jbnd = 1:S%nat3, iq = 1:nR)
      !   T(ibnd,jbnd,iq) = &
      !     dot_product(Vt(:,(iq-1)*S%nat3+ibnd), gV(:,(iq-1)*S%nat3+jbnd))
      !   ! V((iq-1)*S%nat3+ibnd,(iq-1)*S%nat3+jbnd)
      ! enddo
      !
      ! do iq = 1, nR
      !   g0(:,:,iq) = matmul(Us_sc_C(:,:,iq), matmul(g0(:,:,iq), Us_sc(:,:,iq)))
      ! enddo
      !
      TR = 0._dp
      do iR = 1, size(diffs,2)
        do iq = 1, nR
          TR(:,:,iR) = TR(:,:,iR) + g0(:,:,iq) * &
            e_iqr(q(:,iq), diffs(:,iR)) / nR
        enddo
      enddo
      !
      ! do concurrent(na1=1:S%nat, na2=1:S%nat, iR=1:grid%nqtot)
      !   write(155, "(2E20.8)") &
      !   norm2(- diffs(:,1,na1,na2,ir) + S%tau(:,na1) - S%tau(:,na2)), &
      !   sum(abs(TR((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iR)))
      ! enddo
      !
      do iq = 1, out_grid%nqtot
        ! do concurrent(iR = 1:nR, na1=1:S%nat, na2=1:S%nat)
        !   do j = 1, n_weights(iR,na1,na2)
        !     Tq((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq,iw) = &
        !       Tq((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq,iw) + &
        !       TR((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iR) * &
        !       e_iqr(out_grid%xq(:,iq), diffs(:,j,iR,na1,na2)) * weights(j,iR,na1,na2)
        !   enddo
        ! enddo
        do iR = 1, size(diffs,2)
          Tq(:,:,iq,iw) = Tq(:,:,iq,iw) + &
            TR(:,:,iR) * e_iqr(out_grid%xq(:,iq), -diffs(:,iR))
        enddo
        Tq(:,:,iq,iw) = matmul(out_Us_C(:,:,iq), matmul(Tq(:,:,iq,iw), out_Us(:,:,iq)))
      enddo
      !
    enddo
    ! close(155)
    !

    open(111, file="fb-lw.dat")
    do iq = 1, out_grid%nqtot
      do ibnd = 1, S%nat3
        x = out_freqs(ibnd,iq)*input%n_omega/wg%max_f + 1
        x0 = INT(x)
        dx = x - x0
        lws(ibnd) = (1.0_dp - dx) * Tq(ibnd,ibnd,iq,x0) + dx * Tq(ibnd,ibnd,iq,x0+1)
        lws1(ibnd) = (1.0_dp - dx) * wg%w(ibnd,iq,x0) + dx * wg%w(ibnd,iq,x0+1)
      enddo
      call merge_degen(S%nat3, lws, out_freqs(:,iq))
      do ibnd = 1, S%nat3
        write(111, "(5E20.8)") out_freqs(ibnd,iq), lws(ibnd), lws1(ibnd)
      enddo
    enddo
    close(111)
    !
  contains
    function green_0_q(sc_grid, iw, wg, S, grid, U, diffs) result(g0)
      integer, intent(in) :: sc_grid(3)
      integer, intent(in) :: iw
      type(tetra_output), intent(in) :: wg
      type(ph_system_info), intent(in) :: S
      type(q_grid), intent(in) :: grid
      complex(dp), intent(in) :: U(S%nat3,S%nat3,grid%nqtot)
      real(dp), intent(in) :: diffs(:,:)
      !
      complex(dp) :: g0(S%nat3, S%nat3, product(sc_grid))
      complex(dp) :: g0_R(S%nat3, S%nat3, size(diffs,2))
      complex(dp) :: g0_R1(S%nat3, S%nat3, size(diffs,2))
      complex(dp) :: g0_q(S%nat3, S%nat3, grid%nqtot)
      integer :: iq, ibnd, nR, iqf
      real(dp) :: small_q(3, product(sc_grid))
      !
      small_q = grid_vec_cart(sc_grid, S%bg)
      do iq = 1, nR
        small_q(:,iq) = small_q(:,iq) / sc_grid
      enddo
      !
      nR = product(sc_grid)
      g0_q = 0._dp
      do iq = 1, grid%nqtot
        do ibnd = 1, S%nat3
          g0_q(:,:,iq) = g0_q(:,:,iq) + &
            wg%w(ibnd, wg%e(iq), iw) * outer_product(U(:,ibnd,iq))
        enddo
      enddo
      !
      g0_R = 0._dp
      do iR = 1, size(diffs,2)
        do iq = 1, grid%nqtot
          ! iqf = v2index(mod(index2v(iq, grid%n), sc_grid), sc_grid)
          g0_R(:,:,iR) = g0_R(:,:,iR) + &
          ! wg%w(ibnd, wg%e(iq), iw) * outer_product(U(:,ibnd,iq)) * &
            g0_q(:,:,iq) * e_iqr(grid%xq(:,iq), diffs(:,iR))
        enddo
      enddo
      !
      g0 = 0._dp
      do iq = 1, nR
        do iR = 1, size(diffs,2)
          g0(:,:,iq) = g0(:,:,iq) + &
            g0_R(:,:,iR) * e_iqr(small_q(:,iq), -diffs(:,iR)) / nR
        enddo
      enddo
      !
      g0_R1 = 0._dp
      do iR = 1, size(diffs,2)
        do iq = 1, nR
          g0_R1(:,:,iR) = g0_R1(:,:,iR) + &
            g0(:,:,iq) * e_iqr(grid%xq(:,iq), diffs(:,iR))
        enddo
      enddo
      !
      print*, "Difference in g0 functions:", sum(abs(g0_R/g0_R1))/size(g0_R)
    end function
    ! open(10, file='dos_g.dat')
    ! do iw = 1, input%n_omega
    !   dos = 0._dp
    !   do iq = 1, out_grid%nqtot
    !     do ibnd = 1, S%nat3
    !       dos = dos + Tq(ibnd,ibnd,iq,iw) * wg%qw(iq)
    !     enddo
    !   enddo
    !   write(10, "(3E20.8)") wg%en(iw), dos
    ! enddo
    ! close(10)
    !
    ! open(10, file="T-real.dat")
    ! open(11, file="T-imag.dat")
    ! do iq = 1, grid%nqtot
    !   do ibnd = 1, S%nat3
    !     write(10, "(6E20.8)") REAL(T(:,ibnd,iq,2), DP)
    !     write(11, "(6E20.8)") AIMAG(T(:,ibnd,iq,2))
    !   enddo
    !   write(10, *) " "
    !   write(11, *) " "
    ! enddo
    !
  end subroutine
  !
  function sum_phases(x)
    real(dp), intent(in) :: x
    real(dp) :: sum_phases
    !
    sum_phases = 3 * intsimp(bessel_radius, 0.0_dp, x)/(x**3)
    !
  contains
    function intsimp(f, a, b)
      real(dp), external :: f
      integer, parameter :: n = 1000   ! must be even
      real(dp) :: h, intsimp
      real(dp) :: a, b
      integer :: i

      h = (b - a) / n
      intsimp = f(a) + f(b)

      do i = 1, n-1
        if (mod(i,2) == 0) then
          intsimp = intsimp + 2.0d0 * f(a + i*h)
        else
          intsimp = intsimp + 4.0d0 * f(a + i*h)
        end if
      end do

      intsimp = intsimp * h / 3.0d0

    end function
    !
    function bessel_radius(xin)
      real(dp), intent(in) :: xin
      real(dp) :: bessel_radius
      !
      bessel_radius = bessel_j0(xin) * xin**2
    end function
    !
  end function
  !
  subroutine main_defect(S, fc2, fc2_sc, grid, out_grid, input)
    use constants, only: BOHR_RADIUS_CM, RY_TO_CMM1
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(forceconst2_grid), intent(in) :: fc2
    type(q_grid), intent(in) :: grid, out_grid
    type(code_input_type), intent(in) :: input
    type(tetra_output) :: w_in, w_out
    !
    ! type(forceconst2_sc) :: TR(0:input%n_omega)
    complex(dp), dimension(S%nat3, S%nat3) :: Vqq, T_Us, Vqq2, green
    real(dp), dimension(S%nat3, grid%nqtot) :: freqs, freqs1
    ! real(dp) :: weights_flat(S%nat3*grid%nqtot)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs, out_freqs1
    complex(dp), dimension(S%nat3) :: e2

    complex(dp), dimension(S%nat3, out_grid%nqtot) :: lws_out
    complex(dp), dimension(S%nat3, out_grid%nqtot) :: lws_full

    complex(dp), dimension(S%nat3, grid%nqtot) :: lws

    ! complex(dp) :: tetra_flat(S%nat3*grid%nqtot,0:input%n_omega)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3, S%nat3, out_grid%nqtot) :: out_Us!, out_us_copy
    complex(dp), allocatable :: interp(:,:)
    complex(dp), allocatable, dimension(:,:,:,:) :: V
    complex(dp), allocatable, dimension(:,:) :: V_flat, Id_flat, R_flat, &
      GV_flat, GV2_flat, GV3_flat, GV4_flat, poly_T, inv_T
    logical, parameter :: full_born = .false.
    real(dp) :: omega, omegaq
    real(dp) :: lws_intrinsic(S%nat3, grid%nqtot)
    complex(dp) :: lw2
    integer :: iq, iqp, iw, ibnd, nR, jq, Nin, Nout, comp, compj, i
    integer ::  jbnd, ibnd1, dim_deg ! , jn1, jn2
    ! real(dp), allocatable :: freqs_sym(:,:)
    complex(dp), allocatable, dimension(:,:,:) :: Vqqs, Vqqs_out, Vms
    character(15) :: filename
    ! real(dp) :: R_grid(3,out_grid%nqtot)
    real(dp) :: V0_cm3, mult, QTmax
    real(dp) :: delta_q, const_time
    real(dp) :: concentrations(3)
    real(dp) :: spectral_function(3,input%n_omega)
    integer :: iconc, iq_star
    integer :: ierr
    character(30) :: fname
    real(dp) :: cq(3)
    complex(dp) :: self_energy(S%nat3, out_grid%nqtot, input%n_omega)
    complex(dp) :: pixel(grid%nqtot, out_grid%nqtot)
    real(dp) :: dos(input%n_omega)
    real(dp) :: vel(3,S%nat3)
    type(q_grid):: grid_sym
    complex(dp), allocatable :: Us_sym(:,:,:)
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    real(dp) :: E(input%n_omega)
    complex(dp) :: phase_factor(grid%nqtot)
    complex(dp) :: vk_plus_vm
    complex(dp), dimension(S%nat3,S%nat3) :: U_dagger, U_dagger_vm
    !
    !> full born quantities in real space
    ! complex(dp), allocatable :: VG(:,:,:,:)
    !
    concentrations = [0.0_dp, 1e-4_dp, 1e-3_dp] ! defects per cell
    V0_cm3 = S%omega * (BOHR_RADIUS_CM)**3
    print*, "V0_cm3", V0_cm3
    nR   = product(fc2%nq)
    Nin  = S%nat3*grid%nqtot
    Nout = S%nat3*out_grid%nqtot
    allocate(Vqqs(S%nat3,S%nat3,grid%nqtot))
    allocate(Vms(S%nat3,S%nat3,grid%nqtot))
    allocate(Vqqs_out(S%nat3,S%nat3,out_grid%nqtot))
    allocate(V(S%nat3,S%nat3,out_grid%nqtot,grid%nqtot))
    allocate(V_flat(Nin,Nout))
    !
    !> input files are periodic, it can be changed
    !> SC: grid type   (big_nat3,big_nat3,1)
    !> UC: grid type   (nat3,nat3,nR)
    !> RR: new SC type (nat3,nat3,nR,nR)

    ! call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)
    ! call freq_in_grid(S, fc2, grid, freqs, Us)
    !
    call freq_in_grid_degen(S, fc2, fc2_sc, out_grid, out_freqs, out_Us, out_freqs1)
    call freq_in_grid_degen(S, fc2, fc2_sc, grid, freqs, Us, freqs1)

    !
    !> max_freq is slightly larger than the maximum frequency, to be sure that
    !> maxval(out_freqs) <= max_freq
    ! max_freq = maxval(freqs) * 1.1_dp
    ! TETRA_MULTIPLIER = 1.0_dp / max_freq**2
    !
    !> tetra have already the square of freq, it can be changed with
    !> the usual delta formula after some benchmarking
    !
    ! do iw = 1, input%n_omega
    !   E(iw) = max_freq * (iw-1) / input%n_omega
    ! enddo
    !
    call set_wg(S, fc2, grid, input%n_omega, 1.5_dp, w_in)
    call print_message("end of tetra initialization")
    !
    !> serially calculate tetras for an equally spaced omega
    !> interval, after we will interpolate them (even at more than 1st order).
    !> The cycle is serial cause the tetra_weights_green is already parallelized
    ! do iq = 1, out_grid%nqtot
    !   matmul(Sd%fc, matmul(G0(input%sc_grid, omegaq**2, w_in, S, grid_sym, Us_sym). Sd%fc))
    ! !
    ! if(full_born) then
    !   call deallocate_tetra()
    !   call tetra_init_grid_sym(out_grid, S, fc2, w_out)
    !   allocate(w_out%w(S%nat3, w_out%nsym, 0:input%n_omega))
    !   do iw = 0, input%n_omega
    !     omega = max_freq*iw/input%n_omega
    !     w_out%w(:,:,iw) = tetra_weights_green(omega**2)
    !     do iq = 1, out_grid%nqtot
    !       do ibnd = 1, S%nat3
    !         if(isnan(ABS(w_out%w(ibnd,w_out%e(iq),iw)))) w_out%w(ibnd,w_out%e(iq),iw) = 0._dp
    !       enddo
    !     enddo
    !     do iq = 1, out_grid%nqtot
    !       call merge_degen(S%nat3, w_out%w(:,w_out%e(iq),iw), out_freqs(:,iq))
    !     enddo
    !   enddo
    ! endif
    !
    ! open(10, file="phdos.dat", status='replace', action='write')
    ! do iw = 1, input%n_omega
    !   omega = max_freq*(iw-1)/input%n_omega
    !   write(10, "(E17.4,100E17.4)") omega, sum(matmul(w_in%w(:,:,iw), w_in%qw))
    ! enddo
    ! call print_message("end of tetra weights calculation")
    !
    ! call full_V(S, fc2, out_grid, fc2_sc, V, out_us, grid, Us)
    ! V_flat = flatten_RR_cmplx(V)
    !
    do iq = 1, out_grid%nq
      call fc2_sc%r2q(out_grid%xq(:,iq))
      call fc2_sc%r2q(out_grid%xq(:,iq), Vqqs_out(:,:,iq))
    enddo
    !
    ! do ibnd = 6, 8
    !   do jbnd = 7,7
    !     open(1000+jbnd+10*ibnd, file="Vqqs_"//trim(int_to_char(ibnd))// &
    !       "_"//trim(int_to_char(jbnd))//".dat", status='replace', action='write')
    !   enddo
    ! enddo
    !
    ! call read_lw('Ph-pz/lw_NK.30x30x1_T300_s5.out', 91, S%nat3, lws_intrinsic)
    !
    allocate(interp(S%nat3, grid%nqtot))
    do iq = 1, out_grid%nq
      phase_factor = 0._dp
      call fc2_sc%r2q(out_grid%xq(:,iq))
      U_dagger = conjg(transpose(out_Us(:,:,iq)))
      U_dagger_vm = U_dagger
      U_dagger_vm(:,4:) = 0._dp
      do jq = 1, grid%nqtot
        ! vel = velocity(S, fc2, grid%xq(:,jq))
        ! delta_q = S%tpiba * norm2(grid%xq(:,jq) - out_grid%xq(:,iq))
        ! do jbnd = 1, S%nat3
        !   QTmax = delta_q * norm2(vel(:,jbnd)) / lws_intrinsic(jbnd,w_in%e(jq)) * RY_TO_CMM1
        !   if (QTmax < 2._dp) &
        !   print"(4E20.8)", norm2(grid%xq(:,jq)), norm2(vel(:,jbnd)), sum_phases(QTmax)
        ! enddo
        call fc2_sc%r2q(grid%xq(:,jq), Vqqs(:,:,jq))
        do i = 1, size(fc2_sc%defects,2)
          phase_factor(jq) = phase_factor(jq) + &
            e_iqr(grid%xq(:,jq)-out_grid%xq(:,iq), S%tau(:,fc2_sc%defects(1,i))-fc2_sc%taudef)
        enddo
        phase_factor(jq) = phase_factor(jq) / size(fc2_sc%defects,2)
        Vqqs(:,:,jq) = matmul(U_dagger, matmul(Vqqs(:,:,jq), Us(:,:,jq)))
        Vms(:,:,jq) = matmul(U_dagger_vm, Us(:,:,jq))
        ! if (norm2(grid%xq(:,jq) - out_grid%xq(:,iq)) < 1e-10_dp) then
        !   call fftinterp_mat2(out_grid%xq(:,iq), S, fc2, Vqq)
        !   Vqqs(:,:,jq) = Vqqs(:,:,jq) - Vqq
        ! endif
        pixel(jq,iq) = sum(Vqqs(:,:,jq))
      enddo
      !
      do iw = 1, input%n_omega
        if(mod(iw, input%n_omega/10) == 0) &
          print"(A,A,I3,A)", input%calculation, " progress ", NINT(100 * REAL(iw,DP) / input%n_omega), "%"
        omega = w_in%max_f*(iw-1)/input%n_omega
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
            ! Vqq = Vqqs(:,:,jq)
            ! if(allocated(fc2_sc%inclusion_eig)) then
            !   do concurrent(i=1:size(fc2_sc%inclusion_eig))!, ABS(omegaq**2 - fc2_sc%inclusion_eig(i)) > 1e-10)
            !     Vqq = Vqq + outer_product2(fc2_sc%inclusion_Dnx_out(:,i,iq) / &
            !       (omegaq**2 - fc2_sc%inclusion_eig(i)), fc2_sc%inclusion_Dnx_in(:,i,jq))
            !   enddo
            ! endif
            do jbnd = 1, S%nat3
              ! if(norm2(grid%xq(:,jq) - out_grid%xq(:,iq)) < 1e-10_dp .and. &
              !   jbnd == ibnd) cycle
              ! vk_plus_vm = Vms(ibnd,jbnd,jq) * fc2_sc%mass_ratios(1) * omega**2
              ! vk_plus_vm = vk_plus_vm * phase_factor(jq)
              lws_out(ibnd,iq) = lws_out(ibnd,iq) + abs(Vqqs(ibnd,jbnd,jq))**2 * interp(jbnd,w_in%e(jq))
              ! if (ibnd == 6 .and. jbnd == 7 .and. iw == 2548) write(1000+jbnd+ibnd*10, "(6E20.8)") &
              !   cryst2cart(grid%xq(:,jq), S%at, -1), real(lw2), &
              !   AIMAG(interp(jbnd,w_in%e(jq))), aimag(lw2 * interp(jbnd,w_in%e(jq)))
              ! if (ibnd == 7 .and. jbnd == 7 .and. iw == 2548) write(1000+jbnd+ibnd*10, "(6E20.8)") &
              !   cryst2cart(grid%xq(:,jq), S%at, -1), real(lw2), &
              !   AIMAG(interp(jbnd,w_in%e(jq))), aimag(lw2 * interp(jbnd,w_in%e(jq)))
              ! if (ibnd == 8 .and. jbnd == 7 .and. iw == 2560) write(1000+jbnd+ibnd*10, "(6E20.8)") &
              !   cryst2cart(grid%xq(:,jq), S%at, -1), real(lw2), &
              !   AIMAG(interp(jbnd,w_in%e(jq))), aimag(lw2 * interp(jbnd,w_in%e(jq)))
              ! sum_phases(QTmax) ** 2
              ! sum_phases(delta_q * norm2(velocity(S, fc2, out_grid%xq(:,iq))) / lws_intrinsic(ibnd,iq))
            enddo
            ! lws_out(ibnd,iq) = lws_out(ibnd,iq) + lw2
            ! if (ibnd == 6 .and. iw == 2551) write(1236, "(4E20.8)") cryst2cart(grid%xq(:,jq), S%at, -1), aimag(lw2)
            ! if (ibnd == 6 .and. iw == 2568) write(1237, "(4E20.8)") cryst2cart(grid%xq(:,jq), S%at, -1), aimag(lw2)
            ! if (ibnd == 8 .and. iw == 2560) write(1238, "(4E20.8)") cryst2cart(grid%xq(:,jq), S%at, -1), aimag(lw2)
          enddo
        enddo
        call merge_degen(S%nat3, lws_out(:,iq), out_freqs(:,iq))
        !
        if(trim(input%calculation) == 'spfdef') then
          do iconc = 1, size(concentrations)
            mult = concentrations(iconc)
            call tetra_init_sym_cmplx(out_grid, S, out_freqs**2 + mult*lws_out)
            spectral_function(iconc, iw) = &
              sum(matmul(AIMAG(tetra_weights_green_cmplx(omega**2)), out_grid%w))
            ! do iq = 1, out_grid%nqtot
            !   do ibnd = 1, S%nat3
            !     write(11, "(3e14.5,I3,3F7.2,E13.4)") out_freqs(ibnd,iq), lws_out(ibnd,iq), ibnd, out_grid%xq(:,iq), out_grid%w(iq)
            !   enddo
            ! enddo
            ! close(11)
          enddo
        elseif((trim(input%calculation) == 'self') ) then
          ! if (out_grid%nqtot /= 1) &
          ! call errore("main_defect", "you can calculate self-energy only in one q-point at once", 1)
          self_energy(:,iq,iw) = lws_out(:,iq)
        endif
        ! call mpi_bsum(S%nat3, out_grid%nqtot, lws_out)
      enddo
    enddo
    !
    ! do ibnd = 6, 8
    !   do jq = 7,7
    !     close(1000+jq+10*ibnd)
    !   enddo
    ! enddo
    ! !
    open(10, file="dos.dat", status='replace', action='write')
    do iw = 1, input%n_omega
      write(10, "(10E17.8)") w_in%max_f*(iw-1)/input%n_omega, AIMAG(matmul(w_in%w(:,:,iw), w_in%qw(:)))
    enddo
    close(10)
    !
    ! !
    ! write(filename, "(A,I1,A)") "V-real", fc2%nq(1), "p.dat"
    ! open(10, file=filename, status='replace', action='write')
    ! write(filename, "(A,I1,A)") "V-imag", fc2%nq(1), "p.dat"
    ! open(11, file=filename, status='replace', action='write')
    ! do iq = 1, out_grid%nq
    !   write(10, "(10000E15.5)") REAL(pixel(:,iq),DP)
    !   write(11, "(10000E15.5)") AIMAG(pixel(:,iq))
    ! enddo
    ! close(10)
    ! close(11)




    ! lws = 0._dp
    ! lws_out = 0._dp
    ! do iq = 1, out_grid%nq
    !   iqp = iq !+ out_grid%iq0
    !   call fc2_sc%r2q(out_grid%xq(:,iq))
    !   do jq = 1, grid%nq
    !     call fc2_sc%r2q(grid%xq(:,jq), Vqqs(:,:,jq))
    !     if (ALL(ABS(grid%xq(:,jq) - out_grid%xq(:,iq)) < 1e-10_dp)) then
    !       call fftinterp_mat2(out_grid%xq(:,jq), S, fc2, Vqq)
    !       Vqqs(:,:,jq) = Vqqs(:,:,jq) - Vqq
    !     endif
    !   enddo
    !   ibnd = 1
    !   do !ibnd
    !     do ibnd1 = ibnd+1, S%nat3
    !       if(.not. near(out_freqs(ibnd,iq) + out_freqs1(ibnd,iq), &
    !         out_freqs(ibnd1,iq) + out_freqs1(ibnd1,iq))) exit
    !     enddo
    !     dim_deg = ibnd1-ibnd
    !     ! print*, "deg?", iq, ibnd, dim_deg
    !     omegaq = out_freqs(ibnd,iqp)
    !     interp = interp1_matrix(w_in%w, omegaq*input%n_omega/max_freq)
    !     ! interp = tetra_weights(ibnd,iq,iw)
    !     Vqq2 = 0._dp
    !     do jq = 1, grid%nq
    !       green = 0._dp
    !       do jbnd = 1, S%nat3
    !         if(near(norm2(out_grid%xq(:,iq)-grid%xq(:,jq))) .and. near(omegaq, freqs(jbnd,jq))) cycle
    !         ! if(near(norm2(out_grid%xq(:,iq)+grid%xq(:,jq))) .and. near(omegaq, freqs(jbnd,jq))) cycle
    !         SELECT CASE(input%delta_approx)
    !          CASE('gauss')
    !           green = green + outer_product(Us(:,jbnd,jq)) * &
    !             f_gauss(omegaq**2-freqs(jbnd,jq)**2, 1e-10_dp)
    !          CASE('tetra')
    !           green = green + outer_product(Us(:,jbnd,jq)) * &
    !             interp(jbnd,w_in%e(jq))
    !         END SELECT
    !       enddo
    !       Vqq2 = Vqq2 + matmul(Vqqs(:,:,jq),matmul(green, conjg(transpose(Vqqs(:,:,jq)))))
    !     enddo
    !     if (dim_deg > 1) then
    !       call diag_degen_cmplx(S%nat3, dim_deg, Vqq2, &
    !         out_Us(:,ibnd:ibnd+dim_deg-1,iqp), e2(ibnd:ibnd+dim_deg-1))
    !       lws_out(ibnd:ibnd+dim_deg-1,iqp) = e2(ibnd:ibnd+dim_deg-1) / omegaq
    !     else
    !       lws_out(ibnd,iqp) = braket(out_us(:,ibnd,iqp),Vqq2) / omegaq
    !     endif
    !     ibnd = ibnd + dim_deg
    !     if(ibnd > S%nat3) exit
    !   enddo
    ! enddo
    ! call mpi_bsum(S%nat3, out_grid%nqtot, lws_out)
    ! call write_file(out_freqs, lws_out, 'degen-1B.dat', out_grid%type)
    ! call print_message("end of degen calculation")
    !
    ! out_us = us


    ! do iw = my_id, input%n_omega, num_procs
    !   lws_out = 0._dp
    !   if(mod(iw, input%n_omega/10) == 0) &
    !     print"(A,A,I3,A)", input%calculation, " progress ", 100 * iw / input%n_omega, "%"
    !   omega = max_freq*iw/input%n_omega
    !   ! write(filename, '(A5,I3.3,A4)') "self-", iw, '.dat'
    !   ! open(11, file=filename, status='unknown')
    !   do iq = 1, out_grid%nq
    !     iqp = iq + out_grid%iq0
    !     do ibnd = 1, S%nat3
    !       ! lws_out(ibnd,iq) = lws_out(ibnd,iq) + braket(out_us(:,ibnd,iqp), Vqqs_out(:,:,iq))
    !       ! omegaq = out_freqs(ibnd,iqp)
    !       do jq = 1, grid%nqtot
    !         do jbnd = 1, S%nat3
    !           lws_out(ibnd,iq) = lws_out(ibnd,iq) + &
    !             ABS(V(jbnd,ibnd,jq,iq))**2 * w_in%w(jbnd,w_in%e(jq),iw)
    !         enddo
    !       enddo
    !       ! if (ABS(aimag(lws_out(ibnd,iq))) < 1e-15_dp) lws_out(ibnd,iq) = CMPLX(REAL(lws_out(ibnd,iq), dp), 0._dp, dp)
    !     enddo
    !   enddo
    !   if(out_grid%nqtot > 1) then
    !     do iconc = 1, size(concentrations)
    !       mult = concentrations(iconc)
    !       call tetra_init_sym_cmplx(out_grid, S, TETRA_MULTIPLIER * (out_freqs**2 + mult*lws_out))
    !       spectral_function(iconc, iw) = &
    !         sum(matmul(AIMAG(tetra_weights_green_cmplx(TETRA_MULTIPLIER * omega**2)), out_grid%w)) * &
    !         TETRA_MULTIPLIER
    !       ! do iq = 1, out_grid%nqtot
    !       !   do ibnd = 1, S%nat3
    !       !     write(11, "(3e14.5,I3,3F7.2,E13.4)") out_freqs(ibnd,iq), lws_out(ibnd,iq), ibnd, out_grid%xq(:,iq), out_grid%w(iq)
    !       !   enddo
    !       ! enddo
    !       ! close(11)
    !     enddo
    !   else
    !     self_energy(:,iw) = lws_out(:,1)
    !   endif
    ! enddo
    !
    select case(trim(input%calculation))
     case("lw")
      if(trim(out_grid%type) == 'path') then
        write(filename, '(A4,I2.2,A4)') 'path', grid%n(1), '.dat'
      else
        write(filename, '(I2.2,I2.2,A4)') out_grid%n(1), grid%n(1), '.dat'
      endif
      call write_file(out_freqs, lws_out, filename, out_grid%type)
      call print_message("end of 1B calculation")
     case("spfdef")
      open(10, file="spectral_function.dat", status='replace', action='write')
      do iw = 0, input%n_omega
        omega = w_in%max_f*iw/input%n_omega
        write(10, "(100E17.4)") omega, spectral_function(:,iw)
      enddo
     case("self")
      open(10, file="self-energy-def.dat", status='replace', action='write')
      write(10, "(100E25.8)") out_freqs(:,1)
      do iw = 1, input%n_omega
        write(10, "(1200E17.4)") sum(aimag(self_energy(:,:,iw)), dim=1)
      enddo
      close(10)
    end select
    !
    if(full_born) then
      deallocate(V, V_flat)
      allocate(V(S%nat3, S%nat3, out_grid%nqtot, out_grid%nqtot))
      allocate(V_flat(Nout, Nout))
      allocate(GV_flat(Nout,Nout))
      ! allocate(GV2_flat, GV3_flat, GV4_flat, source=GV_flat)
      allocate(inv_T, R_flat, Id_flat, source=GV_flat)
      !
      call full_V(S, out_grid, fc2_sc, V_flat, out_us)
      !
      Id_flat = id_mat(S%nat3*out_grid%nqtot)
      !
      deallocate(interp)
      allocate(interp(S%nat3, out_grid%nqtot))
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          comp = (iq-1)*S%nat3 + ibnd
          omegaq = out_freqs(ibnd,iq)
          interp = interp1_matrix(w_out%w, omegaq*input%n_omega/w_in%max_f)
          do jq = 1, out_grid%nqtot
            do jbnd = 1, S%nat3
              compj = (jq-1)*S%nat3 + jbnd
              GV_flat(compj,comp) = V_flat(compj,comp) * &
                interp(jbnd,w_out%e(jq))
            enddo
          enddo
        enddo
      enddo
      ! print*, sum(GV_flat), sum(V_flat)

      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, GV_flat, &
      !   Nout, GV_flat, Nout, 0._dp, R_flat, Nout)
      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, GV2_flat, &
      !   Nout, GV_flat, Nout, 1._dp, R_flat, Nout)
      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, GV3_flat, &
      !   Nout, GV_flat, Nout, 1._dp, R_flat, Nout)
      ! ! R_flat = GV_flat + GV2_flat + GV3_flat + GV4_flat

      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, V_flat, &
      !   Nout, R_flat, Nout, 0._dp, poly_T, Nout)
      !
      ! > These are the only two lines that enable Full Born
      R_flat = Id_flat - GV_flat
      call invzmat(S%nat3*out_grid%nqtot, R_flat)
      call zgemm('N', 'N', Nout, Nout, Nout, cmplx(1._dp, 0._dp,dp), V_flat, &
        Nout, R_flat, Nout, cmplx(0._dp, 0._dp, dp), inv_T, Nout)
      !
      do iqp = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          comp = ibnd + (iqp-1)*S%nat3
          omegaq = out_freqs(ibnd,iqp)
          ! lws_out(ibnd,iqp) = poly_T(comp,comp) / omegaq
          lws_full(ibnd,iqp) = inv_T(comp,comp) / omegaq
        enddo
      enddo
      call print_message("end of born calculation")
      !
    endif
    !
    if(ionode) then

      if(full_born) then
        call write_file(out_freqs, lws_full, 'FB.dat', out_grid%type)
        call write_file(out_freqs, lws_out, '4B.dat', out_grid%type)
      endif
    endif
    !
    !
    if (full_born) deallocate(V_flat, GV_flat)
    if (full_born) deallocate(w_out%w)
    ! deallocate(V, V_flat)
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
    open(10, file=filename_, status='replace', action='write')
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
  subroutine write_file_raja(freqs_, lws_, filename_, type)
    real(dp), intent(in) :: freqs_(:,:)
    complex(dp), intent(in) :: lws_(:,:)
    character(*), intent(in) :: filename_
    character(*), intent(in) :: type
    !
    integer :: iq_, ibnd_
    !
    open(10, file=filename_, status='replace', action='write')
    do iq_ = 1, size(lws_, 2)
      write(10, "(100E15.5)") -AIMAG(lws_(:,iq_))
    enddo
    close(10)
    !
    open(10, file="freq.dat", status='replace', action='write')
    do iq_ = 1, size(freqs_, 2)
      write(10, "(100E15.5)") freqs_(:,iq_)
    enddo
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
end module
