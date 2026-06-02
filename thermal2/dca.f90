module dca
  use iso_fortran_env, only: int64
  use kinds, only: dp
  use symm_q_mat, only : apply_sym, apply_sym_q
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg, equiv_grid
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils, only: bz2simple, grid_vec_cart, &
    e_iqr, index2v, cryst2cart, v2index, freq_in_grid, diag_cmplx, diag, &
    id_mat
  use defutils, only : flatten_RR_cmplx, unflatten_RR_cmplx, &
    quter_cmplx, fftinterp_mat2_cmplx, quter_R, trace
  use defect, only : tetra_from_self_cart, write_spf_ndiag, write_spf, &
    write_self, write_dos, tetra_from_self_diag, find_where, tetra_from_self, enlarge_R
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize, setup_simple_grid
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi, pi
  use quter_defect, only : forceconst2_sc, inside_ws, full_born_center_images
  use functions, only: f_gauss
  use fc3_interpolate, only: forceconst3, sparse, d3_mixed, sum_R3
  use merge_degenerate, only: merge_degen
  use test_print, only: allclose
  use simtet, only: tetra_init_sym_cmplx, tetra_weights_green_cmplx
  use ph_velocity, only : velocity
  use constants, only : RY_TO_CMM1
  use symm_base, only : symm_matrices => s
  use EPW_utilities, only : mix_broyden_full
  use mpi_thermal, only : mpi_bsum, ionode, num_procs, my_id, mpi_broadcast
  use quter_module, only : quter
  !
  type :: support_site_type
    complex(dp), allocatable :: V_modes(:,:)
    real(dp), allocatable :: V_lambda(:)
    real(dp), allocatable :: R_cart(:,:), diff_cart(:,:)
    integer, allocatable :: iR_large(:,:), dof_R(:), dof_cart(:)
  end type
contains
  subroutine init_random_seed()
    implicit none
    integer :: n
    integer, allocatable :: seed(:)

    call random_seed(size=n)
    allocate(seed(n))

    seed = 12345
    call random_seed(put=seed)

  end subroutine
  !
  function sinc_3d(k, R)
    real(dp), intent(in) :: k(3), R(3)
    real(dp) :: sinc_3d
    real(dp) :: kr(3), sk(3)
    integer :: i
    !
    do i = 1, 3
      kr(i) = tpi * k(i) * R(i)
      if (abs(kr(i)) < 1e-8_dp) then
        sk(i) = 1.0_dp
      else
        sk(i) = sin(kr(i)) / kr(i)
      endif
    enddo
    sinc_3d = sk(1) * sk(2) * sk(3)
  end function
  !
  subroutine find_where_real(R, list_of_R, idx) !res(idx)
    real(dp), intent(in) :: R(3)
    real(dp), intent(in) :: list_of_R(:,:)
    integer, intent(out) :: idx
    integer :: iR
    !
    idx = -1
    do iR = 1, size(list_of_R,2)
      if (all(abs(R - list_of_R(:,iR)) < 1e-8_dp)) then
        idx = iR
        exit
      endif
    enddo
    if (idx == -1) print*, "Error in find_where_real: R not found", R
  end subroutine
  !
  subroutine check_hermitian(name, matq)
    character(*), intent(in) :: name
    complex(dp), intent(in) :: matq(:,:,:)
    integer :: iq
    real(dp) :: max_diff
    !
    max_diff = 0.0_dp
    do iq = 1, size(matq,3)
      max_diff = max(max_diff, maxval(abs(matq(:,:,iq) - conjg(transpose(matq(:,:,iq))))))
    enddo
    if (ionode) print"(A,A,E20.8)", name, "-> maximum deviation from Hermiticity:", max_diff
  end subroutine
  !
  subroutine symmetrize_mat_cmplx(equiv, vec)
    integer, intent(in) :: equiv(:)
    complex(dp), intent(inout) :: vec(:,:,:)
    !
    integer :: iq, Nq_sym
    complex(dp), allocatable :: vec_sym(:,:,:)
    integer, allocatable :: N_star(:)
    !
    Nq_sym = maxval(equiv)
    allocate(vec_sym(size(vec,1), size(vec,2), Nq_sym))
    allocate(N_star(Nq_sym))
    vec_sym = 0.0_dp
    N_star = 0
    !
    do iq = 1, size(vec,3)
      vec_sym(:,:,equiv(iq)) = vec_sym(:,:,equiv(iq)) + vec(:,:,iq)
      N_star(equiv(iq)) = N_star(equiv(iq)) + 1
    enddo
    !
    do iq = 1, size(vec,3)
      vec(:,:,iq) = vec_sym(:,:,equiv(iq)) / N_star(equiv(iq))
    enddo
  end subroutine
  !
  subroutine inside_ws_dca(NQ, R_list, r, nrws, rws, weights, ind_out)
    integer, intent(in) :: NQ
    real(dp), intent(in) :: R_list(:,:)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: rws(:,:)
    !
    real(dp), allocatable, intent(out) :: weights(:)
    integer, allocatable, intent(out) :: ind_out(:)
    !
    real(dp), allocatable :: weights_(:)
    integer, allocatable :: ind_out_(:)
    integer :: i, nweights, nrws
    real(dp) :: wg, wg_tot
    real(dp), external :: wsweight
    !
    allocate(weights_(size(R_list, 2)))
    allocate(ind_out_(size(R_list, 2)))
    wg_tot = 0._dp
    nweights = 0
    do i = 1, size(R_list, 2)
      wg = wsweight(r+R_list(:,i),rws,nrws)
      if (wg > 1e-8) then
        wg_tot = wg_tot + wg
        nweights = nweights + 1
        weights_(nweights) = wg
        ind_out_(nweights) = i
      endif
    enddo
    !
    wg_tot = wg_tot / NQ
    !
    if (ABS(wg_tot-1)>1e-8) then
      print*, "wg_tot is", wg_tot
      call errore("inside_ws", "sum of weights is not 1", 1)
    endif
    allocate(weights(nweights), ind_out(nweights))
    weights = weights_(:nweights)
    ind_out = ind_out_(:nweights)
    deallocate(weights_, ind_out_)
  end subroutine
  !
  integer(int64) function complex_mem_bytes(nel) result(bytes)
    integer(int64), intent(in) :: nel
    bytes = 16_int64 * nel
  end function
  !
  integer(int64) function real_mem_bytes(nel) result(bytes)
    integer(int64), intent(in) :: nel
    bytes = 8_int64 * nel
  end function
  !
  integer(int64) function int_mem_bytes(nel) result(bytes)
    integer(int64), intent(in) :: nel
    bytes = 4_int64 * nel
  end function
  !
  real(dp) function bytes_to_gib(bytes) result(gib)
    integer(int64), intent(in) :: bytes
    gib = real(bytes, dp) / 1073741824._dp
  end function
  !
  integer(int64) function mem_available_bytes() result(bytes)
    character(len=256) :: line, key
    integer :: unit_, ios
    integer(int64) :: kb
    !
    bytes = 0_int64
    open(newunit=unit_, file="/proc/meminfo", status="old", action="read", iostat=ios)
    if(ios /= 0) return
    do
      read(unit_, "(A)", iostat=ios) line
      if(ios /= 0) exit
      read(line, *, iostat=ios) key, kb
      if(ios /= 0) cycle
      if(trim(key) == "MemAvailable:") then
        bytes = kb * 1024_int64
        exit
      endif
    enddo
    close(unit_)
  end function
  !
  subroutine check_allocation_fits(where_, bytes_)
    character(*), intent(in) :: where_
    integer(int64), intent(in) :: bytes_
    integer(int64) :: available_
    character(len=256) :: msg_
    !
    if(bytes_ <= 0_int64) then
      call errore("dca memory check", "invalid or overflowing allocation size for "//trim(where_), 1)
    endif
    available_ = mem_available_bytes()
    if(ionode) then
      if(available_ > 0_int64) then
        print"(A,A,A,F10.3,A,F10.3,A)", "Memory check [", trim(where_), "]: ", &
          bytes_to_gib(bytes_), " GiB requested, ", bytes_to_gib(available_), " GiB available"
      else
        print"(A,A,A,F10.3,A)", "Memory check [", trim(where_), "]: ", &
          bytes_to_gib(bytes_), " GiB requested"
      endif
    endif
    if(available_ > 0_int64 .and. real(bytes_, dp) > 0.80_dp * real(available_, dp)) then
      write(msg_, "(A,F10.3,A,F10.3,A)") trim(where_)//" requests ", &
        bytes_to_gib(bytes_), " GiB; available memory is ", bytes_to_gib(available_), " GiB"
      call errore("dca memory check", trim(msg_), 1)
    endif
  end subroutine
  !
  subroutine dca_selfnrg(S, S_sc, input, fc2, fc2_sc, out_grid)
    type(ph_system_info), intent(in) :: S, S_sc
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    !! centered
    type(forceconst2_sc), intent(inout) :: fc2_sc
    !! not centered
    type(q_grid) :: in_grid
    !! symmetric and scattered, normal order
    type(q_grid), intent(in) :: out_grid
    !! symmetric and not scattered, normal order
    ! complex(dp), allocatable, intent(in) :: self_fb(:,:,:,:)
    ! real(dp), allocatable, intent(in) :: xR_fb(:,:)
    type(tetra_output) :: wg, wg_out
    !
    complex(dp), allocatable, dimension(:) :: delta_in, delta_out, self_diff
    complex(dp), allocatable, dimension(:,:) :: &
      Gf_conf, Gf_avg, out_den_weights, Gf0i, dv, df, Vi, V_born, self_cluster_flat, &
      support_V, support_green, support_g0i, support_self, support_self_before, &
      support_V_modes, support_green_mode, support_self_before_mode, G_clean_flat
    complex(dp), allocatable, dimension(:,:,:) :: &
      self_fine, U, UT, UT_fine, U_fine, UT_out, U_out, self_prime, D, G, &
      self_out_diag, V_conf, self_before, self_next, G0i_cluster, Gi_coarse, &
      self_outp_diag, self_sym, weights, c_out_atom, c_out_full, c_out_right, &
      self_support_out, self_support_R, support_green_R, support_mode_fine, &
      support_mode_cluster, support_mode_out
    complex(dp), allocatable, dimension(:,:,:,:) :: &
      Gi_conf, self_QQ, phase_mat, V, self_out_grid, c_in, c_out, &
      self_outp_grid, self_R, Vout_cluster
    complex(dp), allocatable, dimension(:,:,:,:,:) :: Vqqs, self_uf
    real(dp), allocatable, dimension(:,:) :: xq, R, f, out_freqs, support_diff_cart, support_R_cart, &
      defect_shift_cart
    real(dp), allocatable, dimension(:) :: support_V_lambda
    integer, allocatable :: pos(:), kq(:,:), support_iR_large(:,:), &
      support_dof_R(:), support_dof_cart(:), defect_pos(:,:,:)
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, n_eq_sites, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ, kkq, na1, na2, j1, j2, &
      cluster_mesh(3), Npos, ntot, iqp, N_ITER_TOT, iR, jR, RL_iter, Q_mesh(3), &
      idx(S%nat3), locations(3), jq_center, r1, r2, ipos_cluster, max_conf_defects, &
      max_conf_rank, support_mode_total
    real(dp) :: conc, defect_conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX, max_diff, max_diff_coarse
    real(dp) :: dos(input%n_omega), shift(3), simulated_conc, xqq(3), eta
    logical :: conv, low_concentration, use_compressed_QQ, use_factorized_configs
    complex(dp) :: A(S%nat3,S%nat3), B(S%nat3,S%nat3), A_re(S%nat3,S%nat3), A_im(S%nat3,S%nat3), &
      ialpha(S%nat3,S%nat3), eigc(S%nat3), &
      ee(S%nat3,product(input%sc_grid)), ees(S%nat3,product(input%sc_grid))
    integer, allocatable :: N_sites(:,:), c_equiv(:), iq_of(:), window_ind(:,:), &
      ind(:), equiv_full(:), window_count(:), equiv(:)
    character(len=100) :: filename
    type(q_grid) :: c_grid, in_grid_full
    type(forceconst2_sc) :: fc2_sc_centered
    type(support_site_type), allocatable :: support_site(:)
    complex(dp), allocatable :: proj(:,:), proj_fine(:,:)
    real(dp), allocatable :: eig_proj(:)
    complex(dp), allocatable :: self_uncoarsed(:,:)
    !
    real(dp), allocatable :: img_xR(:,:,:,:), img_weight(:,:,:), diff(:,:)
    integer :: img_nR(S%nat, S%nat)
    complex(dp), allocatable :: W(:,:,:)
    complex(dp), allocatable :: uncoarse(:,:,:,:)
    !
    if(input%calculation == "test") call init_random_seed()
    !
    MAXITER = 100
    ABS_TOLERANCE = 1e-17_dp
    REL_TOLERANCE = 1e-8_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 4
    conc = input%conc ! example concentration
    cluster_mesh = input%sc_grid
    Nc = product(cluster_mesh) !* size(fc2_sc%defects,2)
    n_eq_sites = size(fc2_sc%defects,2)
    defect_conc = conc * n_eq_sites
    NSAMPLES = 100
    NQ = product(input%nk_in) / Nc
    Q_mesh = input%nk_in / cluster_mesh
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true.)
    shift = xq(:,v2index([1,1,1],cluster_mesh))
    shift = cryst2cart(shift, S%at, -1)
    shift = shift * (Q_mesh-1) / 2 / Q_mesh
    shift = cryst2cart(shift, S%bg, 1)
    if (all(input%nk_in == cluster_mesh)) shift = 0._dp
    ! do iq = 1, size(xq,2)
    !   xq(:,iq) = xq(:,iq) + shift
    ! enddo
    !
    call setup_simple_grid(S%bg, input%nk_in(1), input%nk_in(2), input%nk_in(3), in_grid_full, -shift)
    call q_grid_copy(in_grid_full, in_grid)
    call in_grid%symmetrize(S)
    call equiv_grid(in_grid, S, equiv)
    if(num_procs > 1) call in_grid%scatter()
    !if(num_procs > 1) call in_grid_full%scatter()
    ! call in_grid_full%symmetrize(S)
    ! call equiv_grid(in_grid_full, S, equiv_full)
    ! call in_grid_full%destroy()
    ! call setup_simple_grid(S%bg, in_grid%n(1), in_grid%n(2), in_grid%n(3), in_grid_full)
    !
    call setup_simple_grid(S%bg, cluster_mesh(1), cluster_mesh(2), cluster_mesh(3), c_grid)
    call c_grid%symmetrize(S)
    call equiv_grid(c_grid, S, c_equiv)
    call c_grid%destroy()
    !
    if (any(mod(in_grid%n, cluster_mesh) /= 0)) &
      call errore("dca_selfnrg", "grid size not multiple of fc2 grid size", 1)
    R = grid_vec_cart(cluster_mesh, S%at)
    if(input%calculation == 'dos' .or. input%calculation == 'test') &
      call set_wg(S, fc2, out_grid, input%n_omega, wg_out)
    call set_wg(S, fc2, in_grid, input%n_omega, wg)
    !
    call check_allocation_fits("DCA harmonic/grid buffers", &
      real_mem_bytes(int(S%nat3, int64) * int(in_grid_full%nqtot, int64)) + &
      complex_mem_bytes(2_int64 * int(S%nat3, int64) * int(S%nat3, int64) * int(Nc, int64)) + &
      complex_mem_bytes(3_int64 * int(S%nat3, int64) * int(S%nat3, int64) * int(in_grid_full%nqtot, int64)))
    allocate(f(S%nat3,in_grid_full%nqtot), U(S%nat3,S%nat3,Nc), UT(S%nat3,S%nat3,Nc))
    allocate(UT_fine(S%nat3,S%nat3,in_grid_full%nqtot), U_fine(S%nat3,S%nat3,in_grid_full%nqtot))
    allocate(D(S%nat3,S%nat3,in_grid_full%nqtot))
    call freq_in_grid(S, fc2, in_grid_full, f, U_fine)
    do iq = 1, in_grid_full%nq
      iqp = iq + in_grid_full%iq0
      call fftinterp_mat2(in_grid_full%xq(:,iq), S, fc2, D(:,:,iqp))
      if(norm2(in_grid_full%xq(:,iq)) < 1e-8_dp) then
        U_fine(:,:,iq) = D(:,:,iqp)
        call mat2_diag(S%nat3, U_fine(:,:,iq), f(:,1))
      endif
    enddo
    call mpi_bsum(S%nat3, S%nat3, in_grid_full%nqtot, D)
    do iq = 1, in_grid_full%nqtot
      UT_fine(:,:,iq) = conjg(transpose(U_fine(:,:,iq)))
    enddo
    deallocate(f)
    allocate(f(S%nat3, Nc))
    do iq = 1, Nc
      call freq_phq_safe(xq(:,iq), S, fc2, f(:,iq), U(:,:,iq))
      UT(:,:,iq) = conjg(transpose(U(:,:,iq)))
    enddo
    if(ionode) print*, "DCA cluster size:", Nc
    if(ionode) print*, "number of configurations to be averaged:", NSAMPLES
    !
    call check_allocation_fits("DCA output buffers", &
      complex_mem_bytes(int(S%nat3, int64) * int(out_grid%nqtot, int64) * int(input%n_omega, int64)) + &
      complex_mem_bytes(3_int64 * int(S%nat3, int64) * int(S%nat3, int64) * int(out_grid%nqtot, int64)) + &
      real_mem_bytes(int(S%nat3, int64) * int(out_grid%nqtot, int64)) + &
      complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * &
      int(out_grid%nqtot, int64) * int(input%n_omega, int64)))
    allocate(self_out_diag(S%nat3,out_grid%nqtot,input%n_omega))
    allocate(UT_out(S%nat3,S%nat3,out_grid%nqtot), U_out(S%nat3,S%nat3,out_grid%nqtot))
    allocate(out_freqs(S%nat3,out_grid%nqtot), self_out_grid(S%nat3,S%nat3,out_grid%nqtot,input%n_omega))
    allocate(self_support_out(S%nat3,S%nat3,out_grid%nqtot))
    call freq_in_grid(S, fc2, out_grid, out_freqs, U_out)
    do iq = 1, out_grid%nqtot
      UT_out(:,:,iq) = conjg(transpose(U_out(:,:,iq)))
    enddo
    self_out_grid = 0.0_dp
    self_out_diag = 0.0_dp
    !
    call check_allocation_fits("DCA common SCF buffers", &
      complex_mem_bytes(4_int64 * int(S%nat3, int64) * int(S%nat3, int64) * int(Nc, int64)) + &
      complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(in_grid_full%nqtot, int64)) + &
      complex_mem_bytes(2_int64 * int(S%nat3, int64) * int(S%nat3, int64) * &
      int(Nc, int64) * int(MEMORY, int64)) + &
      complex_mem_bytes(2_int64 * int(S%nat3, int64) * int(S%nat3, int64) * int(Nc, int64)) + &
      int_mem_bytes(int(NQ, int64) * int(Nc, int64) + int(Nc, int64)) + &
      complex_mem_bytes(int(S%nat3, int64) * int(out_grid%nqtot, int64)))
    allocate(Gi_coarse(S%nat3, S%nat3, Nc))
    allocate(G0i_cluster(S%nat3, S%nat3, Nc))
    allocate(self_before(S%nat3, S%nat3, Nc))
    allocate(self_next(S%nat3, S%nat3, Nc))
    allocate(kq(NQ, Nc))!, big_iq(grid_scat%nqtot))
    allocate(out_den_weights(S%nat3, out_grid%nqtot))
    allocate(pos(Nc))
    allocate(G(S%nat3, S%nat3, in_grid_full%nqtot))
    !
    allocate(df(S%nat3**2*Nc, MEMORY))
    allocate(dv(S%nat3**2*Nc, MEMORY))
    allocate(delta_in(S%nat3**2*Nc))
    allocate(delta_out(S%nat3**2*Nc))
    df = 0._dp
    dv = 0._dp
    delta_in = 0._dp
    delta_out = 0._dp
    !
    fc2_sc_centered = fc2_sc
    call fc2_sc_centered%center(cluster_mesh, S)
    call prepare_low_concentration_support(fc2_sc_centered, support_V, &
      support_iR_large, support_diff_cart, support_R_cart, &
      support_dof_R, support_dof_cart)
    allocate(defect_shift_cart(3, n_eq_sites))
    do idef = 1, n_eq_sites
      defect_shift_cart(:,idef) = S%tau(:,fc2_sc%defects(1,idef)) - &
        S%tau(:,fc2_sc%defects(1,1))
    enddo
    use_compressed_QQ = input%dca_v_rank > 0
    call project_support_potential(support_V, support_iR_large, support_diff_cart, &
      support_dof_R, support_dof_cart, defect_conc, out_grid%xq, self_support_out)
    !
    allocate(N_sites(n_eq_sites, NSAMPLES))
    do idef = 1, n_eq_sites
      call assign_defect_counts(NSAMPLES, Nc, conc, N_sites(idef,:))
    enddo
    call mpi_broadcast(n_eq_sites, NSAMPLES, N_sites)
    !
    ntot = sum(N_sites)
    simulated_conc = real(ntot,dp) / (Nc * n_eq_sites * NSAMPLES)
    if(ionode) print"(A,E15.4)", "simulated concentration:", simulated_conc
    low_concentration = .false.
    if(all(N_sites < 2)) then
      low_concentration = .true.
      if(ionode) print*, "Using low concentration approximation (at most 1 defect per configuration)"
      simulated_conc = conc
    endif
    use_factorized_configs = .not. low_concentration .and. use_compressed_QQ
    if(use_factorized_configs .and. ionode) &
      print*, "Using support-factor multi-defect DCA solver with rank:", input%dca_v_rank
    if(use_factorized_configs) then
      allocate(support_site(n_eq_sites))
      do idef = 1, n_eq_sites
        call prepare_equivalent_site_support_modes(idef, input%dca_v_rank, support_site(idef))
      enddo
      support_mode_total = site_mode_count(support_site)
    elseif(use_compressed_QQ) then
      call prepare_support_V_modes(support_V, input%dca_v_rank, support_V_modes, support_V_lambda)
    endif
    if(.not. low_concentration) then
      allocate(defect_pos(Nc, n_eq_sites, NSAMPLES))
      defect_pos = 0
      do it = 1+my_id, NSAMPLES, num_procs
        do idef = 1, n_eq_sites
          Npos = N_sites(idef, it)
          call sample_canonical(Nc, Npos, pos)
          if(Npos > 0) defect_pos(1:Npos, idef, it) = pos(1:Npos)
        enddo
      enddo
      call mpi_bsum(Nc, n_eq_sites, NSAMPLES, defect_pos)
    endif
    if(low_concentration .and. use_compressed_QQ) then
      call check_allocation_fits("DCA support-mode moment buffers", &
        complex_mem_bytes(int(S%nat3, int64) * int(size(support_V_lambda), int64) * &
        (int(in_grid_full%nqtot, int64) + int(Nc, int64))) + &
        complex_mem_bytes(2_int64 * int(size(support_V_lambda), int64) * &
        int(size(support_V_lambda), int64)))
      allocate(support_mode_fine(S%nat3, size(support_V_lambda), in_grid_full%nqtot))
      allocate(support_mode_cluster(S%nat3, size(support_V_lambda), Nc))
      allocate(support_green_mode(size(support_V_lambda), size(support_V_lambda)))
      allocate(support_self_before_mode(size(support_V_lambda), size(support_V_lambda)))
      call prepare_support_mode_phase(support_V_modes, support_R_cart, &
        support_dof_R, support_dof_cart, in_grid_full%xq, support_mode_fine)
      call prepare_support_mode_phase(support_V_modes, support_R_cart, &
        support_dof_R, support_dof_cart, xq, support_mode_cluster)
    elseif(use_factorized_configs) then
      call check_allocation_fits("DCA site support-mode output projector", &
        complex_mem_bytes(int(S%nat3, int64) * int(support_mode_total, int64) * &
        (int(Nc, int64) + int(out_grid%nqtot, int64))) + &
        complex_mem_bytes(int(support_mode_total, int64) * int(support_mode_total, int64)))
      allocate(support_mode_cluster(S%nat3, support_mode_total, Nc))
      allocate(support_mode_out(S%nat3, support_mode_total, out_grid%nqtot))
      allocate(support_self_before_mode(support_mode_total, support_mode_total))
      call prepare_site_support_mode_phase(support_site, xq, support_mode_cluster)
      call prepare_site_support_mode_phase(support_site, out_grid%xq, support_mode_out)
    elseif(low_concentration .and. NQ > 1) then
      call check_allocation_fits("DCA support-moment dense buffers", &
        complex_mem_bytes(4_int64 * int(size(support_V,1), int64) * int(size(support_V,1), int64)) + &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * &
        int(size(support_diff_cart,2), int64)))
      allocate(support_green(size(support_V,1), size(support_V,1)))
      allocate(support_g0i(size(support_V,1), size(support_V,1)))
      allocate(support_self(size(support_V,1), size(support_V,1)))
      allocate(support_self_before(size(support_V,1), size(support_V,1)))
      allocate(support_green_R(S%nat3, S%nat3, size(support_diff_cart,2)))
    endif
    if(.not. low_concentration) then
      if(use_factorized_configs) then
        max_conf_defects = maxval(sum(N_sites, dim=1))
        max_conf_rank = 0
        do it = 1, NSAMPLES
          max_conf_rank = max(max_conf_rank, config_factor_rank(support_site, N_sites(:,it)))
        enddo
        call check_allocation_fits("DCA support-factor multi-defect buffers", &
          complex_mem_bytes(5_int64 * int(S%nat3, int64) * int(Nc, int64) * &
          int(S%nat3, int64) * int(Nc, int64)) + &
          complex_mem_bytes(3_int64 * int(S%nat3, int64) * int(Nc, int64) * &
          int(max(1, max_conf_rank), int64)) + &
          complex_mem_bytes(int(max(1, max_conf_rank), int64) * &
          int(max(1, max_conf_rank), int64)))
        if(ionode) then
          print*, "Maximum defects per sampled configuration:", max_conf_defects
          print*, "Maximum support-factor rank per configuration:", max_conf_rank
        endif
        allocate(Gf0i(S%nat3*Nc, S%nat3*Nc))
        allocate(Gf_conf(S%nat3*Nc, S%nat3*Nc))
        allocate(Gf_avg(S%nat3*Nc, S%nat3*Nc))
        allocate(self_cluster_flat(S%nat3*Nc, S%nat3*Nc))
        allocate(G_clean_flat(S%nat3*Nc, S%nat3*Nc))
        if(trim(input%calculation) == 'test') &
          call diagnose_factorized_single_defects(support_mode_cluster, support_site)
      else
        call center_V(xq, S, S_sc, fc2_sc, Vqqs)
        call check_allocation_fits("DCA non-low-concentration QQ' buffers", &
          complex_mem_bytes(5_int64 * int(S%nat3, int64) * int(Nc, int64) * &
          int(S%nat3, int64) * int(Nc, int64)) + &
          complex_mem_bytes(int(S%nat3, int64) * int(Nc, int64) * &
          int(S%nat3, int64) * int(Nc, int64) * &
          int(NSAMPLES, int64)) + &
          complex_mem_bytes(int(Nc, int64) * int(Nc, int64) * &
          int(n_eq_sites, int64) * int(NSAMPLES, int64)))
        allocate(Gf0i(S%nat3*Nc, S%nat3*Nc))
        allocate(Gf_conf(S%nat3*Nc, S%nat3*Nc))
        allocate(Gf_avg(S%nat3*Nc, S%nat3*Nc))
        allocate(self_cluster_flat(S%nat3*Nc, S%nat3*Nc))
        allocate(V_born(S%nat3*Nc, S%nat3*Nc))
        allocate(V_conf(S%nat3*Nc, S%nat3*Nc, NSAMPLES))
        allocate(phase_mat(Nc, Nc, n_eq_sites, NSAMPLES))
        !
        V_born = flatten_RR_cmplx(Vqqs(:,:,:,:,1))
        phase_mat = 0.0_dp
        do it = 1, NSAMPLES
          do idef = 1, n_eq_sites
            Npos = N_sites(idef, it)
            do ipos = 1, Npos
              ipos_cluster = defect_pos(ipos, idef, it)
              do k = 1, Nc
                do j = 1, Nc
                  phase_mat(j,k,idef,it) = phase_mat(j,k,idef,it) + &
                    e_iqr(xq(:,k)-xq(:,j), R(:,ipos_cluster))
                enddo
              enddo
            enddo
          enddo
        enddo
        !
        V_conf = 0.0_dp
        do it = 1, NSAMPLES
          do jq = 1, Nc
            do iq = 1, Nc
              do idef = 1, n_eq_sites
                V_conf((iq-1)*S%nat3+1:iq*S%nat3,(jq-1)*S%nat3+1:jq*S%nat3,it) = &
                  V_conf((iq-1)*S%nat3+1:iq*S%nat3,(jq-1)*S%nat3+1:jq*S%nat3,it) + &
                  Vqqs(:,:,iq,jq,idef) * phase_mat(iq,jq,idef,it) / Nc
              enddo
            enddo
          enddo
          V_conf(:,:,it) = (V_conf(:,:,it)+conjg(transpose(V_conf(:,:,it)))) / 2.0_dp
        enddo
        deallocate(phase_mat, Vqqs)
      endif
    endif
    !
    !
    !> construction of G0_coarse and G0i_coarse, which are averaged over the small
    !> patch in which self-energy is assumed constant
    do iq = 1, Nc
      do jq = 1, NQ
        kq(jq,iq) = v2index(index2v(jq, Q_mesh) + &
          index2v(iq, cluster_mesh) * Q_mesh, in_grid%n)
      enddo
    enddo
    !
    ! The diagonal-Q projector can be tested with pair_diff_projector, but for
    ! this centered support it is underdetermined: many true differences alias
    ! onto the same cluster residue.  Keep the full c_out interpolation active.
    ! allocate(W(Nc, S%nat, S%nat))
    ! allocate(uncoarse(S%nat3, S%nat3, Nc, Nc))
    ! uncoarse = 0._dp
    !
    ! W = 0._dp
    ! do na2 = 1, S%nat
    !   do na1 = 1, S%nat
    !     do iR = 1, img_nR(na1,na2)
    !       jR = v2index(bz2simple(NINT(cryst2cart(img_xR(:,iR,na1,na2), S%bg, -1)), cluster_mesh), cluster_mesh)
    !       do jq = 1, NQ
    !         W(jR,na1,na2) = W(jR,na1,na2) + &
    !           e_iqr(-in_grid_full%xq(:,kq(jq,1)), img_xR(:,iR,na1,na2)) * img_weight(iR,na1,na2) / NQ
    !       enddo
    !     enddo
    !   enddo
    ! enddo
    !
    ! do jq = 1, Nc
    !   do iq = 1, Nc
    !     do na2 = 1, S%nat
    !       do na1 = 1, S%nat
    !         do iR = 1, Nc
    !           if(abs(W(iR,na1,na2)) < 1e-10_dp) print*, W(iR,na1,na2), "is very small, check q-grid and image construction"
    !           uncoarse(3*(na1-1)+1:3*na1,3*(na2-1)+1:3*na2, iq, jq) = &
    !             uncoarse(3*(na1-1)+1:3*na1,3*(na2-1)+1:3*na2, iq, jq) + &
    !             e_iqr(-xq(:,jq) + xq(:,iq), R(:,iR)) / Nc / W(iR,na1,na2)
    !         enddo
    !       enddo
    !     enddo
    !   enddo
    ! enddo
    !
    !>
    if(.not. low_concentration .and. .not. use_factorized_configs) then
      call center_V_mixed(out_grid%xq, xq, Vout_cluster)
      call fit_potential_cout_full(c_out_full, V_born, Vout_cluster)
      deallocate(Vout_cluster)
      call center_V_mixed(xq, out_grid%xq, Vout_cluster)
      call fit_potential_cout_right(c_out_right, V_born, Vout_cluster)
      deallocate(Vout_cluster)
    endif
    !
    self_next = 1._dp
    self_before = 0._dp
    !
    dos = 0._dp
    N_ITER_TOT = 0
    if(ionode) print*, "Starting DCA self-energy calculation..."
    if(ionode .and. low_concentration .and. NQ > 1 .and. .not. use_compressed_QQ) &
      print*, "Using support-moment DCA bath for the low-concentration solver."
    if(ionode .and. low_concentration .and. use_compressed_QQ) &
      print*, "Using compressed support-moment low-concentration solver with rank:", size(support_V_lambda)
    if(ionode) print*, ""
    do iw = 3, input%n_omega
      do iq = 1, in_grid_full%nqtot
        G(:,:,iq) = matmul(U_fine(:,:,iq), matmul(diag_cmplx(1/wg%w(:,wg%e(iq),iw)), UT_fine(:,:,iq)))
        ! call invzmat(S%nat3, G(:,:,iq))
      enddo
      !
      conv = .false.
      do sc_iter = 1, MAXITER
        Gi_coarse = 0.0_dp
        if(low_concentration .and. use_compressed_QQ) support_green_mode = 0._dp
        if(low_concentration .and. NQ > 1 .and. .not. use_compressed_QQ) support_green_R = 0._dp
        do iq = 1, Nc
          do jq = 1, NQ
            kkq = kq(jq,iq)
            A = G(:,:,kkq) - self_before(:,:,iq)
            call invzmat(S%nat3, A)
            Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) + A
            if(low_concentration .and. use_compressed_QQ) &
              call add_support_mode_green(A, support_mode_fine(:,:,kkq), support_green_mode)
            if(low_concentration .and. NQ > 1 .and. .not. use_compressed_QQ) &
              call add_support_moment_green_R(A, in_grid_full%xq(:,kkq), &
              support_diff_cart, support_green_R)
          enddo
          Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) / real(NQ, dp)
          call invzmat(S%nat3, Gi_coarse(:,:,iq))
          G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_before(:,:,iq)
        enddo
        !
        if(low_concentration) then
          if(use_compressed_QQ) then
            call support_self_mode_from_cluster(self_before, support_mode_cluster, &
              support_self_before_mode)
            call low_concentration_support_self_compressed_moment(support_green_mode, &
              support_self_before_mode, support_V_modes, support_V_lambda, support_iR_large, &
              support_dof_R, support_dof_cart, self_support_R)
          elseif(NQ == 1) then
            call low_concentration_support_self_R(G0i_cluster, support_V, &
              support_iR_large, support_diff_cart, support_R_cart, &
              support_dof_R, support_dof_cart, self_support_R)
          else
            call expand_support_matrix(support_green_R, support_iR_large, &
              support_dof_R, support_dof_cart, support_green)
            call support_self_from_cluster(self_before, support_diff_cart, support_iR_large, &
              support_dof_R, support_dof_cart, support_self_before)
            support_g0i = support_green
            call invzmat(size(support_V,1), support_g0i)
            support_g0i = support_g0i + support_self_before
            support_green = support_g0i
            call invzmat(size(support_V,1), support_green)
            call low_concentration_support_self_dense(support_green, support_V, &
              support_iR_large, support_dof_R, support_dof_cart, self_support_R, support_self)
          endif
          call project_support_self(self_support_R, support_diff_cart, xq, self_next)
          call project_support_self(self_support_R, support_diff_cart, out_grid%xq, self_support_out)
        else
          Gf0i = 0._dp
          do iq = 1, Nc
            Gf0i((iq-1)*S%nat3+1:iq*S%nat3, (iq-1)*S%nat3+1:iq*S%nat3) = &
              G0i_cluster(:,:,iq)
          enddo
          if(use_factorized_configs) then
            G_clean_flat = 0._dp
            do iq = 1, Nc
              A = G0i_cluster(:,:,iq)
              call invzmat(S%nat3, A)
              G_clean_flat((iq-1)*S%nat3+1:iq*S%nat3, &
                (iq-1)*S%nat3+1:iq*S%nat3) = A
            enddo
            Gf_avg = 0.0_dp
            do it = 1+my_id, NSAMPLES, num_procs
              call factorized_config_green_Q(G_clean_flat, support_site, &
                defect_pos(:,:,it), N_sites(:,it), Gf_conf)
              Gf_avg = Gf_avg + Gf_conf
            enddo
            call mpi_bsum(S%nat3*Nc, S%nat3*Nc, Gf_avg)
            Gf_avg = Gf_avg / real(NSAMPLES, dp)
            Gf_conf = Gf_avg
            call invzmat(S%nat3*Nc, Gf_conf)
            self_cluster_flat = Gf0i - Gf_conf
          else
            Gf_avg = 0.0_dp
            do it = 1+my_id, NSAMPLES, num_procs
              Gf_conf = Gf0i - V_conf(:,:,it)
              call invzmat(S%nat3*Nc, Gf_conf)
              Gf_avg = Gf_avg + Gf_conf
            enddo
            call mpi_bsum(S%nat3*Nc, S%nat3*Nc, Gf_avg)
            Gf_avg = Gf_avg / real(NSAMPLES, dp)
            Gf_conf = Gf_avg
            call invzmat(S%nat3*Nc, Gf_conf)
            self_cluster_flat = Gf0i - Gf_conf
          endif
          do iq = 1, Nc
            self_next(:,:,iq) = self_cluster_flat((iq-1)*S%nat3+1:iq*S%nat3, &
              (iq-1)*S%nat3+1:iq*S%nat3)
          enddo
        endif
        !
        call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_next, c_equiv, xq, .true.)
        max_diff = maxval(abs(self_before - self_next))
        if(all(abs(real(self_before - self_next, dp)) < ABS_TOLERANCE + abs(real(self_next, dp)) * REL_TOLERANCE) .and. &
          all(abs(aimag(self_before - self_next)) < ABS_TOLERANCE + abs(aimag(self_next)) * REL_TOLERANCE)) conv = .true.
        if(conv) then
          self_before = self_next
          exit
        endif
        !
        delta_out = reshape(self_next, [S%nat3**2*Nc])
        call mix_broyden_full(S%nat3**2*Nc, delta_out, delta_in, &
          ALPHA_MIX, sc_iter, MEMORY, df, dv)
        self_before = reshape(delta_in, [S%nat3, S%nat3, Nc])
      enddo ! self-energy SC cycle
      !
      N_ITER_TOT = N_ITER_TOT + min(sc_iter, MAXITER)
      if(conv) then
        if(ionode) print"(A,I4,A,I4,A)", "frequency ", iw, " converged in ", sc_iter, " iterations."
      else
        if(ionode) print"(A,I4,A,I4,A,E15.3)", "frequency ", iw, &
          " NOT converged in ", MAXITER, " iterations. Max diff: ", max_diff
      end if
      if(use_factorized_configs) then
        call cluster_self_to_support_modes(self_cluster_flat, support_mode_cluster, &
          support_self_before_mode)
      endif
      !
      do iq = 1, out_grid%nqtot
        if(low_concentration) then
          A = self_support_out(:,:,iq)
        elseif(use_factorized_configs) then
          call support_mode_self_at_q(support_self_before_mode, support_mode_out(:,:,iq), A)
        else
          Gf_conf = self_cluster_flat - defect_conc * V_born
          A = self_support_out(:,:,iq) + &
            matmul(c_out_full(:,:,iq), matmul(Gf_conf, c_out_right(:,:,iq)))
        endif
        call apply_sym_q(S, out_grid%xq(:,iq), A)
        self_out_grid(:,:,iq,iw) = matmul(UT_out(:,:,iq), matmul(A, U_out(:,:,iq)))
        do i = 1, S%nat3
          self_out_diag(i,iq,iw) = self_out_grid(i,i,iq,iw)
        enddo
      enddo
    enddo ! frequency loop
    !
    if(ionode) print*, "Average number of iterations per frequency:", real(N_ITER_TOT,dp) / input%n_omega
    select case(input%calculation)
     case ('spf-def')
      call write_spf_ndiag('spf-dca-ndiag.dat', wg%en, self_out_grid, out_freqs)
      call write_spf('spf-dca.dat', wg%en, self_out_diag, out_freqs)
      call write_self('self-dca.dat', wg%en, self_out_diag)
     case('self')
      call write_self('self-dca.dat', wg%en, self_out_diag)
      ! call write_self('selfp-dca.dat', wg%en, self_outp_diag)
     case('dos')
      do iw = 1+my_id, input%n_omega, num_procs
        call tetra_from_self(S, out_grid, out_freqs, self_out_grid(:,:,:,iw), wg_out%en(iw)**2, out_den_weights)
        dos(iw) = sum(matmul(AIMAG(out_den_weights), wg_out%qw)) * product(out_grid%n)
      enddo
      call mpi_bsum(input%n_omega, dos)
      call write_dos('dos-dca.dat', wg_out%en, dos)
     case("test")
      do iw = 1+my_id, input%n_omega, num_procs
        call tetra_from_self(S, out_grid, out_freqs, self_out_grid(:,:,:,iw), wg_out%en(iw)**2, out_den_weights)
        dos(iw) = sum(matmul(AIMAG(out_den_weights), wg_out%qw)) * product(out_grid%n)
      enddo
      call mpi_bsum(input%n_omega, dos)
      call write_dos('dos-dca-test.dat', wg_out%en, dos)
      call write_spf_ndiag('spf-dca-ndiag-test.dat', wg%en, self_out_grid(:,:,:1,:), out_freqs(:,:1))
      call write_spf('spf-dca-test.dat', wg%en, self_out_diag(:,:1,:), out_freqs(:,:1))
      call write_self('self-dca-test.dat', wg%en, self_out_diag(:,:1,:))
     case default
      if(ionode) print*, "WARNING: unknown DCA calculation type"
    end select
    !
  contains
    subroutine low_concentration_cluster_self(G0i_flat_, V_cluster_flat_, c_cluster_, self_flat_)
      ! Fast cluster-DCA update:
      ! <G> = (1-c) g + c (g^-1 - V)^-1, using the same unnormalised
      ! Fourier convention as full_born_center.
      complex(dp), intent(in) :: G0i_flat_(:,:), V_cluster_flat_(:,:)
      real(dp), intent(in) :: c_cluster_
      complex(dp), intent(out) :: self_flat_(:,:)
      !
      complex(dp), allocatable :: G_clean_(:,:), G1_(:,:), Gavg_(:,:)
      integer :: nflat_
      !
      nflat_ = S%nat3 * Nc
      if(size(G0i_flat_, 1) /= nflat_ .or. size(G0i_flat_, 2) /= nflat_) &
        call errore("low_concentration_cluster_self", "bad G0 size", 1)
      if(size(V_cluster_flat_, 1) /= nflat_ .or. size(V_cluster_flat_, 2) /= nflat_) &
        call errore("low_concentration_cluster_self", "bad V size", 1)
      if(size(self_flat_, 1) /= nflat_ .or. size(self_flat_, 2) /= nflat_) &
        call errore("low_concentration_cluster_self", "bad Sigma size", 1)
      !
      allocate(G_clean_(nflat_,nflat_), G1_(nflat_,nflat_), Gavg_(nflat_,nflat_))
      G_clean_ = G0i_flat_
      call invzmat(nflat_, G_clean_)
      G1_ = G0i_flat_ - V_cluster_flat_
      call invzmat(nflat_, G1_)
      Gavg_ = (1._dp - c_cluster_) * G_clean_ + c_cluster_ * G1_
      self_flat_ = Gavg_
      call invzmat(nflat_, self_flat_)
      self_flat_ = G0i_flat_ - self_flat_
      !
      deallocate(G_clean_, G1_, Gavg_)
    end subroutine
    !
    subroutine factorized_config_green_Q(G_clean_, support_site_, config_pos_, config_counts_, G_conf_)
      ! Multi-defect configuration solve with
      !   V_conf(Q,Q') = B(Q) Lambda B(Q')^dagger.
      ! The site-specific support modes are already centered at the proper
      ! equivalent atom; here we only translate them by the sampled cluster R.
      complex(dp), intent(in) :: G_clean_(:,:)
      type(support_site_type), intent(in) :: support_site_(:)
      integer, intent(in) :: config_pos_(:,:), config_counts_(:)
      complex(dp), intent(out) :: G_conf_(:,:)
      !
      complex(dp), allocatable :: Bfac_(:,:), GB_(:,:), right_(:,:), K_(:,:)
      real(dp), allocatable :: lambda_col_(:)
      real(dp) :: center_R_(3), support_R_(3), inv_sqrt_Nc_, lambda_scale_
      complex(dp) :: phase_
      integer :: nflat_, N_, M_, Mconf_, col_, idef_, ipos_, im_, i_, iq_, idx_
      integer :: center_idx_
      !
      nflat_ = S%nat3 * Nc
      Mconf_ = config_factor_rank(support_site_, config_counts_)
      if(size(G_clean_, 1) /= nflat_ .or. size(G_clean_, 2) /= nflat_) &
        call errore("factorized_config_green_Q", "bad clean Green size", 1)
      if(size(G_conf_, 1) /= nflat_ .or. size(G_conf_, 2) /= nflat_) &
        call errore("factorized_config_green_Q", "bad configuration Green size", 1)
      if(size(support_site_) /= size(config_counts_)) &
        call errore("factorized_config_green_Q", "site/count size mismatch", 1)
      !
      G_conf_ = G_clean_
      if(Mconf_ == 0) return
      do idef_ = 1, size(support_site_)
        lambda_scale_ = maxval(abs(support_site_(idef_)%V_lambda))
        if(lambda_scale_ <= 0._dp) call errore("factorized_config_green_Q", "zero support spectrum", 1)
        if(any(abs(support_site_(idef_)%V_lambda) < 1e-14_dp * lambda_scale_)) &
          call errore("factorized_config_green_Q", "retained zero support eigenvalue", 1)
      enddo
      !
      allocate(Bfac_(nflat_, Mconf_), GB_(nflat_, Mconf_), &
        right_(Mconf_, nflat_), K_(Mconf_, Mconf_))
      allocate(lambda_col_(Mconf_))
      Bfac_ = 0._dp
      inv_sqrt_Nc_ = 1._dp / sqrt(real(Nc, dp))
      col_ = 0
      do idef_ = 1, size(config_counts_)
        N_ = size(support_site_(idef_)%dof_R)
        M_ = size(support_site_(idef_)%V_lambda)
        if(size(support_site_(idef_)%V_modes, 1) /= N_ .or. &
          size(support_site_(idef_)%V_modes, 2) /= M_) &
          call errore("factorized_config_green_Q", "bad site support mode size", 1)
        if(size(support_site_(idef_)%dof_cart) /= N_) &
          call errore("factorized_config_green_Q", "bad site active map", 1)
        if(any(support_site_(idef_)%dof_R < 1) .or. &
          any(support_site_(idef_)%dof_R > size(support_site_(idef_)%R_cart, 2))) &
          call errore("factorized_config_green_Q", "active R out of range", 1)
        if(any(support_site_(idef_)%dof_cart < 1) .or. &
          any(support_site_(idef_)%dof_cart > S%nat3)) &
          call errore("factorized_config_green_Q", "active cartesian index out of range", 1)
        do ipos_ = 1, config_counts_(idef_)
          center_idx_ = config_pos_(ipos_, idef_)
          if(center_idx_ < 1 .or. center_idx_ > Nc) &
            call errore("factorized_config_green_Q", "sampled defect position out of range", 1)
          center_R_ = R(:,center_idx_)
          do im_ = 1, M_
            col_ = col_ + 1
            lambda_col_(col_) = support_site_(idef_)%V_lambda(im_)
            do i_ = 1, N_
              support_R_ = center_R_ + support_site_(idef_)%R_cart(:,support_site_(idef_)%dof_R(i_))
              do iq_ = 1, Nc
                idx_ = support_site_(idef_)%dof_cart(i_) + (iq_ - 1) * S%nat3
                phase_ = e_iqr(xq(:,iq_), -support_R_)
                Bfac_(idx_,col_) = Bfac_(idx_,col_) + &
                  inv_sqrt_Nc_ * phase_ * support_site_(idef_)%V_modes(i_,im_)
              enddo
            enddo
          enddo
        enddo
      enddo
      if(col_ /= Mconf_) call errore("factorized_config_green_Q", "bad factor count", 1)
      !
      GB_ = matmul(G_clean_, Bfac_)
      right_ = matmul(conjg(transpose(Bfac_)), G_clean_)
      K_ = -matmul(conjg(transpose(Bfac_)), GB_)
      do col_ = 1, Mconf_
        K_(col_,col_) = K_(col_,col_) + 1._dp / lambda_col_(col_)
      enddo
      call invzmat(Mconf_, K_)
      G_conf_ = G_conf_ + matmul(GB_, matmul(K_, right_))
      !
      deallocate(Bfac_, GB_, right_, K_, lambda_col_)
    end subroutine
    !
    subroutine prepare_low_concentration_support(fc2sc_, V_support_, iR_large_, diff_cart_, R_cart_, &
      dof_R_, dof_cart_)
      ! Build the frequency-independent centered support used by full_born_center.
      type(forceconst2_sc), intent(in) :: fc2sc_
      complex(dp), allocatable, intent(out) :: V_support_(:,:)
      integer, allocatable, intent(out) :: iR_large_(:,:)
      real(dp), allocatable, intent(out) :: diff_cart_(:,:), R_cart_(:,:)
      integer, allocatable, intent(out) :: dof_R_(:), dof_cart_(:)
      !
      integer, pointer :: R_list_(:,:), diff_list_large_(:,:)
      complex(dp), allocatable :: V_full_(:,:)
      logical, allocatable :: active_(:)
      integer, allocatable :: active_idx_(:)
      integer :: iR1_, iR2_, iR_, nR_, nR_large_, i_, j_, N_, nactive_, idx_
      integer :: Rdiff_(3)
      real(dp), parameter :: support_tol_ = 1e-14_dp
      !
      allocate(R_list_(3,1))
      R_list_(:,1) = fc2sc_%yR2(:,1)
      nR_ = 1
      do iR2_ = 1, fc2sc_%n_R2
        call find_where(fc2sc_%yR2(:,iR2_), R_list_, iR_)
        if(iR_ == -1) call enlarge_R(fc2sc_%yR2(:,iR2_), R_list_, nR_)
        do iR1_ = 1, fc2sc_%n_R1(iR2_)
          call find_where(fc2sc_%yR1(:,iR1_,iR2_), R_list_, iR_)
          if(iR_ == -1) call enlarge_R(fc2sc_%yR1(:,iR1_,iR2_), R_list_, nR_)
        enddo
      enddo
      !
      nR_large_ = 1
      allocate(diff_list_large_(3,1))
      diff_list_large_(:,1) = [0,0,0]
      do i_ = 1, nR_
        do j_ = 1, nR_
          Rdiff_ = R_list_(:,i_) - R_list_(:,j_)
          call find_where(Rdiff_, diff_list_large_, iR_)
          if(iR_ == -1) call enlarge_R(Rdiff_, diff_list_large_, nR_large_)
        enddo
      enddo
      !
      allocate(iR_large_(nR_,nR_))
      allocate(diff_cart_(3,nR_large_))
      allocate(R_cart_(3,nR_))
      do iR_ = 1, nR_
        R_cart_(:,iR_) = cryst2cart(real(R_list_(:,iR_), dp), S%at, 1)
      enddo
      do iR_ = 1, nR_large_
        diff_cart_(:,iR_) = cryst2cart(real(diff_list_large_(:,iR_), dp), S%at, 1)
      enddo
      do j_ = 1, nR_
        do i_ = 1, nR_
          call find_where(R_list_(:,i_) - R_list_(:,j_), diff_list_large_, iR_)
          iR_large_(i_,j_) = iR_
        enddo
      enddo
      !
      N_ = S%nat3 * nR_
      call check_allocation_fits("DCA centered support V_full", &
        complex_mem_bytes(int(N_, int64) * int(N_, int64)))
      allocate(V_full_(N_,N_))
      V_full_ = 0._dp
      do iR2_ = 1, fc2sc_%n_R2
        call find_where(fc2sc_%yR2(:,iR2_), R_list_, j_)
        do iR1_ = 1, fc2sc_%n_R1(iR2_)
          call find_where(fc2sc_%yR1(:,iR1_,iR2_), R_list_, i_)
          V_full_((i_-1)*S%nat3+1:i_*S%nat3, (j_-1)*S%nat3+1:j_*S%nat3) = &
            cmplx(fc2sc_%fc(:,:,iR1_,iR2_), 0._dp, dp)
        enddo
      enddo
      !
      allocate(active_(N_), active_idx_(N_))
      active_ = .false.
      do i_ = 1, N_
        active_(i_) = any(abs(V_full_(i_,:)) > support_tol_) .or. &
          any(abs(V_full_(:,i_)) > support_tol_)
      enddo
      nactive_ = count(active_)
      if(nactive_ == 0) call errore("prepare_low_concentration_support", "empty active support", 1)
      !
      j_ = 0
      do i_ = 1, N_
        if(.not. active_(i_)) cycle
        j_ = j_ + 1
        active_idx_(j_) = i_
      enddo
      call check_allocation_fits("DCA active centered support V", &
        complex_mem_bytes(int(nactive_, int64) * int(nactive_, int64)) + &
        int_mem_bytes(2_int64 * int(nactive_, int64)))
      allocate(V_support_(nactive_,nactive_))
      allocate(dof_R_(nactive_), dof_cart_(nactive_))
      do i_ = 1, nactive_
        idx_ = active_idx_(i_)
        dof_R_(i_) = (idx_ - 1) / S%nat3 + 1
        dof_cart_(i_) = mod(idx_ - 1, S%nat3) + 1
        do j_ = 1, nactive_
          V_support_(i_,j_) = V_full_(idx_, active_idx_(j_))
        enddo
      enddo
      if(ionode) print"(A,I8,A,I8)", &
        "Low-concentration active support dimension: ", nactive_, " / ", N_
      !
      deallocate(V_full_, active_, active_idx_)
      deallocate(R_list_, diff_list_large_)
    end subroutine
    !
    subroutine prepare_support_V_modes(V_support_, requested_rank_, V_modes_, V_lambda_)
      ! Diagonalize the centered Hermitian support potential and keep the
      ! strongest modes by absolute eigenvalue.  This gives
      ! V_support ~= V_modes * diag(V_lambda) * V_modes^dagger.
      complex(dp), intent(in) :: V_support_(:,:)
      integer, intent(in) :: requested_rank_
      complex(dp), allocatable, intent(out) :: V_modes_(:,:)
      real(dp), allocatable, intent(out) :: V_lambda_(:)
      !
      complex(dp), allocatable :: eigvec_(:,:), work_(:)
      complex(dp) :: work_query_(1)
      real(dp), allocatable :: eigval_(:), rwork_(:)
      logical, allocatable :: used_(:)
      integer :: N_, M_, info_, lwork_, im_, i_, best_
      real(dp) :: best_abs_, total_norm_, kept_norm_
      external :: zheev
      !
      N_ = size(V_support_, 1)
      if(size(V_support_, 2) /= N_) call errore("prepare_support_V_modes", "bad support V", 1)
      M_ = min(max(requested_rank_, 1), N_)
      !
      call check_allocation_fits("DCA support-V eigensolver dense copy", &
        complex_mem_bytes(int(N_, int64) * int(N_, int64)) + &
        real_mem_bytes(int(N_, int64) + int(max(1,3*N_-2), int64)))
      allocate(eigvec_(N_,N_), eigval_(N_), rwork_(max(1,3*N_-2)))
      eigvec_ = V_support_
      call zheev('V', 'U', N_, eigvec_, N_, eigval_, work_query_, -1, rwork_, info_)
      call errore("prepare_support_V_modes", "ZHEEV workspace query failed", abs(info_))
      lwork_ = max(1, int(real(work_query_(1), dp)))
      call check_allocation_fits("DCA support-V eigensolver workspace", &
        complex_mem_bytes(int(lwork_, int64)))
      allocate(work_(lwork_))
      call zheev('V', 'U', N_, eigvec_, N_, eigval_, work_, lwork_, rwork_, info_)
      call errore("prepare_support_V_modes", "ZHEEV failed", abs(info_))
      !
      call check_allocation_fits("DCA compressed support-V modes", &
        complex_mem_bytes(int(N_, int64) * int(M_, int64)) + &
        real_mem_bytes(int(M_, int64)) + int_mem_bytes(int(N_, int64)))
      allocate(V_modes_(N_,M_), V_lambda_(M_), used_(N_))
      used_ = .false.
      total_norm_ = sum(eigval_**2)
      kept_norm_ = 0._dp
      do im_ = 1, M_
        best_ = 1
        best_abs_ = -1._dp
        do i_ = 1, N_
          if(used_(i_)) cycle
          if(abs(eigval_(i_)) > best_abs_) then
            best_abs_ = abs(eigval_(i_))
            best_ = i_
          endif
        enddo
        used_(best_) = .true.
        V_modes_(:,im_) = eigvec_(:,best_)
        V_lambda_(im_) = eigval_(best_)
        kept_norm_ = kept_norm_ + eigval_(best_)**2
      enddo
      !
      if(ionode) then
        print"(A,I8,A,I8)", "Compressed centered V modes: ", M_, " / ", N_
        if(total_norm_ > 0._dp) print"(A,F12.8)", &
          "Captured Frobenius-norm weight of V: ", kept_norm_ / total_norm_
        print"(A,ES16.8)", "Largest retained |lambda|: ", maxval(abs(V_lambda_))
        print"(A,ES16.8)", "Smallest retained |lambda|: ", minval(abs(V_lambda_))
      endif
      !
      deallocate(eigvec_, eigval_, rwork_, work_, used_)
    end subroutine
    !
    integer function site_mode_count(support_site_)
      type(support_site_type), intent(in) :: support_site_(:)
      !
      integer :: isite_
      !
      site_mode_count = 0
      do isite_ = 1, size(support_site_)
        site_mode_count = site_mode_count + size(support_site_(isite_)%V_lambda)
      enddo
    end function
    !
    integer function config_factor_rank(support_site_, config_counts_)
      type(support_site_type), intent(in) :: support_site_(:)
      integer, intent(in) :: config_counts_(:)
      !
      integer :: isite_
      !
      if(size(support_site_) /= size(config_counts_)) &
        call errore("config_factor_rank", "site/count size mismatch", 1)
      config_factor_rank = 0
      do isite_ = 1, size(support_site_)
        config_factor_rank = config_factor_rank + &
          config_counts_(isite_) * size(support_site_(isite_)%V_lambda)
      enddo
    end function
    !
    subroutine prepare_equivalent_site_support_modes(site_, requested_rank_, support_site_)
      integer, intent(in) :: site_, requested_rank_
      type(support_site_type), intent(inout) :: support_site_
      !
      type(forceconst2_sc) :: fc_site_
      complex(dp), allocatable :: V_support_(:,:)
      !
      if(ionode) print*, "Preparing symmetry-centered support modes for site:", site_
      call center_equivalent_defect_fc(site_, fc_site_)
      call prepare_low_concentration_support(fc_site_, V_support_, &
        support_site_%iR_large, support_site_%diff_cart, support_site_%R_cart, &
        support_site_%dof_R, support_site_%dof_cart)
      call prepare_support_V_modes(V_support_, requested_rank_, &
        support_site_%V_modes, support_site_%V_lambda)
      deallocate(V_support_)
      call fc_site_%deallocate()
    end subroutine
    !
    subroutine project_support_potential(V_support_, iR_large_, diff_cart_, &
      dof_R_, dof_cart_, prefactor_, xq_target_, self_cart_)
      complex(dp), intent(in) :: V_support_(:,:)
      integer, intent(in) :: iR_large_(:,:)
      real(dp), intent(in) :: diff_cart_(:,:), xq_target_(:,:)
      integer, intent(in) :: dof_R_(:), dof_cart_(:)
      real(dp), intent(in) :: prefactor_
      complex(dp), intent(out) :: self_cart_(:,:,:)
      !
      complex(dp), allocatable :: self_R_(:,:,:)
      integer :: i_, j_, N_, nR_large_
      !
      N_ = size(dof_R_)
      nR_large_ = size(diff_cart_, 2)
      if(size(dof_cart_) /= N_) call errore("project_support_potential", "bad active map", 1)
      if(size(V_support_, 1) /= N_ .or. size(V_support_, 2) /= N_) &
        call errore("project_support_potential", "bad support V", 1)
      !
      call check_allocation_fits("DCA projected support potential", &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(nR_large_, int64)))
      allocate(self_R_(S%nat3,S%nat3,nR_large_))
      self_R_ = 0._dp
      do j_ = 1, N_
        do i_ = 1, N_
          self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) = &
            self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) + &
            prefactor_ * V_support_(i_,j_)
        enddo
      enddo
      call project_support_self(self_R_, diff_cart_, xq_target_, self_cart_)
      deallocate(self_R_)
    end subroutine
    !
    subroutine add_support_moment_green_R(gq_, xq_fine_, diff_cart_, G_R_)
      ! Add exp(i q d) G(q) to every unique support difference d = R-R'.
      ! The dense support matrix is expanded only once after the fine-q loop.
      complex(dp), intent(in) :: gq_(:,:)
      real(dp), intent(in) :: xq_fine_(3), diff_cart_(:,:)
      complex(dp), intent(inout) :: G_R_(:,:,:)
      !
      integer :: iR_
      !
      if(size(G_R_, 1) /= S%nat3 .or. size(G_R_, 2) /= S%nat3 .or. &
        size(G_R_, 3) /= size(diff_cart_, 2)) &
        call errore("add_support_moment_green_R", "bad real-space support G", 1)
      do iR_ = 1, size(diff_cart_, 2)
        G_R_(:,:,iR_) = G_R_(:,:,iR_) + gq_ * e_iqr(xq_fine_, diff_cart_(:,iR_))
      enddo
    end subroutine
    !
    subroutine expand_support_matrix(mat_R_, iR_large_, dof_R_, dof_cart_, mat_support_)
      complex(dp), intent(in) :: mat_R_(:,:,:)
      integer, intent(in) :: iR_large_(:,:), dof_R_(:), dof_cart_(:)
      complex(dp), intent(out) :: mat_support_(:,:)
      !
      integer :: i_, j_, N_
      !
      N_ = size(dof_R_)
      if(size(dof_cart_) /= N_) call errore("expand_support_matrix", "bad active map", 1)
      if(size(mat_support_, 1) /= N_ .or. size(mat_support_, 2) /= N_) &
        call errore("expand_support_matrix", "bad support matrix", 1)
      !
      do j_ = 1, N_
        do i_ = 1, N_
          mat_support_(i_,j_) = mat_R_(dof_cart_(i_),dof_cart_(j_), &
            iR_large_(dof_R_(i_),dof_R_(j_)))
        enddo
      enddo
    end subroutine
    !
    subroutine support_self_from_cluster(self_cluster_, diff_cart_, iR_large_, dof_R_, dof_cart_, self_support_)
      ! Periodize the current diagonal DCA self-energy to the support.  The
      ! 1/Nc is the inverse transform paired with project_support_self.
      complex(dp), intent(in) :: self_cluster_(:,:,:)
      real(dp), intent(in) :: diff_cart_(:,:)
      integer, intent(in) :: iR_large_(:,:), dof_R_(:), dof_cart_(:)
      complex(dp), intent(out) :: self_support_(:,:)
      !
      complex(dp), allocatable :: self_R_(:,:,:)
      integer :: iq_, iR_, N_
      !
      N_ = size(dof_R_)
      if(size(dof_cart_) /= N_) call errore("support_self_from_cluster", "bad active map", 1)
      if(size(self_support_, 1) /= N_ .or. size(self_support_, 2) /= N_) &
        call errore("support_self_from_cluster", "bad support Sigma", 1)
      !
      call check_allocation_fits("DCA cluster self-energy support projection", &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(size(diff_cart_,2), int64)))
      allocate(self_R_(S%nat3,S%nat3,size(diff_cart_,2)))
      self_R_ = 0._dp
      do iR_ = 1, size(diff_cart_, 2)
        do iq_ = 1, Nc
          self_R_(:,:,iR_) = self_R_(:,:,iR_) + &
            e_iqr(xq(:,iq_), diff_cart_(:,iR_)) * self_cluster_(:,:,iq_) / real(Nc, dp)
        enddo
      enddo
      call expand_support_matrix(self_R_, iR_large_, dof_R_, dof_cart_, self_support_)
      deallocate(self_R_)
    end subroutine
    !
    subroutine prepare_support_mode_phase(V_modes_, R_cart_, dof_R_, dof_cart_, xq_grid_, Bmode_)
      ! B_alpha,m(q) = sum_i exp(-i q R_i) U_i,m for active support rows
      ! with cartesian index alpha.  This projects support moments without
      ! forming the dense support matrix.
      complex(dp), intent(in) :: V_modes_(:,:)
      real(dp), intent(in) :: R_cart_(:,:), xq_grid_(:,:)
      integer, intent(in) :: dof_R_(:), dof_cart_(:)
      complex(dp), intent(out) :: Bmode_(:,:,:)
      !
      complex(dp) :: phase_
      integer :: iq_, i_, N_, M_, nq_
      !
      N_ = size(dof_R_)
      M_ = size(V_modes_, 2)
      nq_ = size(xq_grid_, 2)
      if(size(dof_cart_) /= N_) call errore("prepare_support_mode_phase", "bad active map", 1)
      if(size(V_modes_, 1) /= N_) call errore("prepare_support_mode_phase", "bad mode matrix", 1)
      if(size(Bmode_, 1) /= S%nat3 .or. size(Bmode_, 2) /= M_ .or. size(Bmode_, 3) /= nq_) &
        call errore("prepare_support_mode_phase", "bad phase buffer", 1)
      if(any(dof_R_ < 1) .or. any(dof_R_ > size(R_cart_, 2))) &
        call errore("prepare_support_mode_phase", "active R out of range", 1)
      if(any(dof_cart_ < 1) .or. any(dof_cart_ > S%nat3)) &
        call errore("prepare_support_mode_phase", "active cartesian index out of range", 1)
      !
      Bmode_ = 0._dp
      do iq_ = 1, nq_
        do i_ = 1, N_
          phase_ = e_iqr(xq_grid_(:,iq_), -R_cart_(:,dof_R_(i_)))
          Bmode_(dof_cart_(i_),:,iq_) = Bmode_(dof_cart_(i_),:,iq_) + &
            phase_ * V_modes_(i_,:)
        enddo
      enddo
    end subroutine
    !
    subroutine prepare_site_support_mode_phase(support_site_, xq_grid_, Bmode_)
      type(support_site_type), intent(in) :: support_site_(:)
      real(dp), intent(in) :: xq_grid_(:,:)
      complex(dp), intent(out) :: Bmode_(:,:,:)
      !
      complex(dp), allocatable :: Bsite_(:,:,:)
      integer :: isite_, m0_, M_, nq_
      !
      nq_ = size(xq_grid_, 2)
      if(size(Bmode_, 1) /= S%nat3 .or. size(Bmode_, 2) /= site_mode_count(support_site_) .or. &
        size(Bmode_, 3) /= nq_) call errore("prepare_site_support_mode_phase", "bad phase buffer", 1)
      Bmode_ = 0._dp
      m0_ = 0
      do isite_ = 1, size(support_site_)
        M_ = size(support_site_(isite_)%V_lambda)
        call check_allocation_fits("DCA site support-mode phase", &
          complex_mem_bytes(int(S%nat3, int64) * int(M_, int64) * int(nq_, int64)))
        allocate(Bsite_(S%nat3, M_, nq_))
        call prepare_support_mode_phase(support_site_(isite_)%V_modes, &
          support_site_(isite_)%R_cart, support_site_(isite_)%dof_R, &
          support_site_(isite_)%dof_cart, xq_grid_, Bsite_)
        Bmode_(:,m0_+1:m0_+M_,:) = Bsite_
        m0_ = m0_ + M_
        deallocate(Bsite_)
      enddo
    end subroutine
    !
    subroutine prepare_shifted_support_mode_phase(V_modes_, R_cart_, dof_R_, dof_cart_, &
      defect_shift_cart_, xq_grid_, Bmode_)
      complex(dp), intent(in) :: V_modes_(:,:)
      real(dp), intent(in) :: R_cart_(:,:), defect_shift_cart_(:,:), xq_grid_(:,:)
      integer, intent(in) :: dof_R_(:), dof_cart_(:)
      complex(dp), intent(out) :: Bmode_(:,:,:)
      !
      complex(dp), allocatable :: Bref_(:,:,:)
      complex(dp) :: phase_
      integer :: M_, nsite_, nq_, isite_, iq_, m0_
      !
      M_ = size(V_modes_, 2)
      nsite_ = size(defect_shift_cart_, 2)
      nq_ = size(xq_grid_, 2)
      if(size(defect_shift_cart_, 1) /= 3) &
        call errore("prepare_shifted_support_mode_phase", "bad site shifts", 1)
      if(size(Bmode_, 1) /= S%nat3 .or. size(Bmode_, 2) /= M_ * nsite_ .or. &
        size(Bmode_, 3) /= nq_) call errore("prepare_shifted_support_mode_phase", "bad phase buffer", 1)
      !
      call check_allocation_fits("DCA shifted support-mode reference phase", &
        complex_mem_bytes(int(S%nat3, int64) * int(M_, int64) * int(nq_, int64)))
      allocate(Bref_(S%nat3, M_, nq_))
      call prepare_support_mode_phase(V_modes_, R_cart_, dof_R_, dof_cart_, xq_grid_, Bref_)
      do iq_ = 1, nq_
        do isite_ = 1, nsite_
          m0_ = (isite_ - 1) * M_
          phase_ = e_iqr(xq_grid_(:,iq_), -defect_shift_cart_(:,isite_))
          Bmode_(:,m0_+1:m0_+M_,iq_) = phase_ * Bref_(:,:,iq_)
        enddo
      enddo
      deallocate(Bref_)
    end subroutine
    !
    subroutine cluster_self_to_support_modes(self_flat_, Bmode_cluster_, self_mode_)
      complex(dp), intent(in) :: self_flat_(:,:), Bmode_cluster_(:,:,:)
      complex(dp), intent(out) :: self_mode_(:,:)
      !
      complex(dp), allocatable :: Phi_(:,:), Gram_(:,:), invGram_(:,:), tmp_(:,:), proj_(:,:), approx_(:,:)
      real(dp) :: inv_sqrt_Nc_, diag_max_, reg_, res_norm_, target_norm_
      integer :: M_, nflat_, iq_, imode_
      !
      nflat_ = S%nat3 * Nc
      M_ = size(Bmode_cluster_, 2)
      if(size(self_flat_, 1) /= nflat_ .or. size(self_flat_, 2) /= nflat_) &
        call errore("cluster_self_to_support_modes", "bad cluster self-energy", 1)
      if(size(Bmode_cluster_, 1) /= S%nat3 .or. size(Bmode_cluster_, 3) /= Nc) &
        call errore("cluster_self_to_support_modes", "bad support-mode phase", 1)
      if(size(self_mode_, 1) /= M_ .or. size(self_mode_, 2) /= M_) &
        call errore("cluster_self_to_support_modes", "bad mode self-energy", 1)
      !
      call check_allocation_fits("DCA support-mode output fit", &
        complex_mem_bytes(int(nflat_, int64) * int(M_, int64)) + &
        complex_mem_bytes(4_int64 * int(M_, int64) * int(M_, int64)) + &
        complex_mem_bytes(int(nflat_, int64) * int(nflat_, int64)))
      allocate(Phi_(nflat_,M_), Gram_(M_,M_), invGram_(M_,M_), &
        tmp_(nflat_,M_), proj_(M_,M_), approx_(nflat_,nflat_))
      inv_sqrt_Nc_ = 1._dp / sqrt(real(Nc, dp))
      do iq_ = 1, Nc
        Phi_((iq_-1)*S%nat3+1:iq_*S%nat3,:) = inv_sqrt_Nc_ * Bmode_cluster_(:,:,iq_)
      enddo
      Gram_ = matmul(conjg(transpose(Phi_)), Phi_)
      diag_max_ = 0._dp
      do imode_ = 1, M_
        diag_max_ = max(diag_max_, abs(Gram_(imode_,imode_)))
      enddo
      reg_ = max(1e-18_dp, 1e-12_dp * diag_max_)
      invGram_ = Gram_
      do imode_ = 1, M_
        invGram_(imode_,imode_) = invGram_(imode_,imode_) + cmplx(reg_, 0._dp, dp)
      enddo
      call invzmat(M_, invGram_)
      tmp_ = matmul(self_flat_, Phi_)
      proj_ = matmul(conjg(transpose(Phi_)), tmp_)
      self_mode_ = matmul(invGram_, matmul(proj_, invGram_))
      if(ionode) then
        approx_ = matmul(Phi_, matmul(self_mode_, conjg(transpose(Phi_))))
        res_norm_ = sum(abs(approx_ - self_flat_)**2)
        target_norm_ = sum(abs(self_flat_)**2)
        print"(A,E12.4,A,E12.4)", "support-mode Sigma fit residual: ", &
          sqrt(res_norm_ / max(target_norm_, 1e-300_dp)), " regularization: ", reg_
      endif
      deallocate(Phi_, Gram_, invGram_, tmp_, proj_, approx_)
    end subroutine
    !
    subroutine support_mode_self_at_q(self_mode_, Bmode_q_, self_q_)
      complex(dp), intent(in) :: self_mode_(:,:), Bmode_q_(:,:)
      complex(dp), intent(out) :: self_q_(:,:)
      !
      if(size(Bmode_q_, 2) /= size(self_mode_, 1) .or. size(self_mode_, 1) /= size(self_mode_, 2)) &
        call errore("support_mode_self_at_q", "bad mode dimensions", 1)
      if(size(Bmode_q_, 1) /= size(self_q_, 1) .or. size(Bmode_q_, 1) /= size(self_q_, 2)) &
        call errore("support_mode_self_at_q", "bad q self dimensions", 1)
      self_q_ = matmul(Bmode_q_, matmul(self_mode_, conjg(transpose(Bmode_q_))))
    end subroutine
    !
    subroutine diagnose_factorized_single_defects(Bmode_cluster_, support_site_)
      complex(dp), intent(in) :: Bmode_cluster_(:,:,:)
      type(support_site_type), intent(in) :: support_site_(:)
      !
      complex(dp), allocatable :: Vqqs_ref_(:,:,:,:,:), Vblock_(:,:), tmp_(:,:)
      integer :: M_, nsite_, idef_, iq_, jq_, m0_, imode_
      real(dp) :: err2_, ref2_, maxerr_, rel_
      !
      nsite_ = size(support_site_)
      if(size(Bmode_cluster_, 1) /= S%nat3 .or. &
        size(Bmode_cluster_, 2) /= site_mode_count(support_site_) .or. &
        size(Bmode_cluster_, 3) /= Nc) call errore("diagnose_factorized_single_defects", &
        "bad shifted support modes", 1)
      if(ionode) print*, "Diagnosing factorized single-defect V(Q,Q') against dense center_V..."
      call center_V(xq, S, S_sc, fc2_sc, Vqqs_ref_)
      allocate(Vblock_(S%nat3,S%nat3))
      m0_ = 0
      do idef_ = 1, nsite_
        M_ = size(support_site_(idef_)%V_lambda)
        allocate(tmp_(S%nat3,M_))
        err2_ = 0._dp
        ref2_ = 0._dp
        maxerr_ = 0._dp
        do jq_ = 1, Nc
          do iq_ = 1, Nc
            tmp_ = Bmode_cluster_(:,m0_+1:m0_+M_,iq_)
            do imode_ = 1, M_
              tmp_(:,imode_) = tmp_(:,imode_) * support_site_(idef_)%V_lambda(imode_)
            enddo
            Vblock_ = matmul(tmp_, conjg(transpose(Bmode_cluster_(:, &
              m0_ + 1:m0_ + M_,jq_)))) / real(Nc, dp)
            maxerr_ = max(maxerr_, maxval(abs(Vblock_ - &
              Vqqs_ref_(:,:,iq_,jq_,idef_) / real(Nc, dp))))
            err2_ = err2_ + sum(abs(Vblock_ - &
              Vqqs_ref_(:,:,iq_,jq_,idef_) / real(Nc, dp))**2)
            ref2_ = ref2_ + sum(abs(Vqqs_ref_(:,:,iq_,jq_,idef_) / real(Nc, dp))**2)
          enddo
        enddo
        rel_ = sqrt(err2_ / max(ref2_, 1e-300_dp))
        if(ionode) print"(A,I4,A,E12.4,A,E12.4)", &
          "factorized V(Q,Q') residual site ", idef_, ": fro=", rel_, " max=", maxerr_
        m0_ = m0_ + M_
        deallocate(tmp_)
      enddo
      deallocate(Vqqs_ref_, Vblock_)
    end subroutine
    !
    subroutine add_support_mode_green(gq_, Bmode_q_, Gmode_)
      complex(dp), intent(in) :: gq_(:,:), Bmode_q_(:,:)
      complex(dp), intent(inout) :: Gmode_(:,:)
      !
      complex(dp) :: tmp_(size(Bmode_q_,1), size(Bmode_q_,2))
      integer :: M_
      !
      M_ = size(Bmode_q_, 2)
      if(size(Bmode_q_, 1) /= S%nat3) call errore("add_support_mode_green", "bad mode phase", 1)
      if(size(gq_, 1) /= S%nat3 .or. size(gq_, 2) /= S%nat3) &
        call errore("add_support_mode_green", "bad Green block", 1)
      if(size(Gmode_, 1) /= M_ .or. size(Gmode_, 2) /= M_) &
        call errore("add_support_mode_green", "bad mode Green", 1)
      !
      tmp_ = matmul(gq_, Bmode_q_)
      Gmode_ = Gmode_ + matmul(conjg(transpose(Bmode_q_)), tmp_)
    end subroutine
    !
    subroutine support_self_mode_from_cluster(self_cluster_, Bmode_cluster_, self_mode_)
      complex(dp), intent(in) :: self_cluster_(:,:,:), Bmode_cluster_(:,:,:)
      complex(dp), intent(out) :: self_mode_(:,:)
      !
      complex(dp) :: tmp_(size(Bmode_cluster_,1), size(Bmode_cluster_,2))
      integer :: iq_, M_
      !
      M_ = size(Bmode_cluster_, 2)
      if(size(Bmode_cluster_, 1) /= S%nat3 .or. size(Bmode_cluster_, 3) /= Nc) &
        call errore("support_self_mode_from_cluster", "bad mode phase", 1)
      if(size(self_cluster_, 1) /= S%nat3 .or. size(self_cluster_, 2) /= S%nat3 .or. &
        size(self_cluster_, 3) /= Nc) call errore("support_self_mode_from_cluster", "bad Sigma", 1)
      if(size(self_mode_, 1) /= M_ .or. size(self_mode_, 2) /= M_) &
        call errore("support_self_mode_from_cluster", "bad mode Sigma", 1)
      !
      self_mode_ = 0._dp
      do iq_ = 1, Nc
        tmp_ = matmul(self_cluster_(:,:,iq_), Bmode_cluster_(:,:,iq_))
        self_mode_ = self_mode_ + &
          matmul(conjg(transpose(Bmode_cluster_(:,:,iq_))), tmp_) / real(Nc, dp)
      enddo
    end subroutine
    !
    subroutine low_concentration_support_self_compressed_moment(G_int_mode_, self_before_mode_, &
      V_modes_, V_lambda_, iR_large_, dof_R_, dof_cart_, self_R_)
      complex(dp), intent(in) :: G_int_mode_(:,:), self_before_mode_(:,:), V_modes_(:,:)
      real(dp), intent(in) :: V_lambda_(:)
      integer, intent(in) :: iR_large_(:,:), dof_R_(:), dof_cart_(:)
      complex(dp), allocatable, intent(out) :: self_R_(:,:,:)
      !
      complex(dp), allocatable :: G_bath_mode_(:,:), Kmode_(:,:), Taumode_(:,:), row_tau_(:)
      complex(dp) :: sigma_ij_
      integer :: i_, j_, im_, N_, M_, nR_large_
      !
      N_ = size(dof_R_)
      M_ = size(V_lambda_)
      nR_large_ = maxval(iR_large_)
      if(size(dof_cart_) /= N_) &
        call errore("low_concentration_support_self_compressed_moment", "bad active map", 1)
      if(size(V_modes_, 1) /= N_ .or. size(V_modes_, 2) /= M_) &
        call errore("low_concentration_support_self_compressed_moment", "bad mode matrix", 1)
      if(size(G_int_mode_, 1) /= M_ .or. size(G_int_mode_, 2) /= M_) &
        call errore("low_concentration_support_self_compressed_moment", "bad mode Green", 1)
      if(size(self_before_mode_, 1) /= M_ .or. size(self_before_mode_, 2) /= M_) &
        call errore("low_concentration_support_self_compressed_moment", "bad mode Sigma", 1)
      if(any(dof_cart_ < 1) .or. any(dof_cart_ > S%nat3)) &
        call errore("low_concentration_support_self_compressed_moment", &
        "active cartesian index out of range", 1)
      !
      call check_allocation_fits("DCA compressed support-moment solver", &
        complex_mem_bytes(3_int64 * int(M_, int64) * int(M_, int64)) + &
        complex_mem_bytes(int(M_, int64)) + &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(nR_large_, int64)))
      allocate(G_bath_mode_(M_,M_), Kmode_(M_,M_), Taumode_(M_,M_), row_tau_(M_))
      !
      G_bath_mode_ = G_int_mode_
      call invzmat(M_, G_bath_mode_)
      G_bath_mode_ = G_bath_mode_ + self_before_mode_
      call invzmat(M_, G_bath_mode_)
      !
      Kmode_ = id_mat(M_)
      do j_ = 1, M_
        Kmode_(:,j_) = Kmode_(:,j_) - (1._dp - defect_conc) * G_bath_mode_(:,j_) * V_lambda_(j_)
      enddo
      call invzmat(M_, Kmode_)
      Taumode_ = Kmode_
      do i_ = 1, M_
        Taumode_(i_,:) = defect_conc * V_lambda_(i_) * Taumode_(i_,:)
      enddo
      !
      allocate(self_R_(S%nat3,S%nat3,nR_large_))
      self_R_ = 0._dp
      do i_ = 1, N_
        row_tau_ = 0._dp
        do im_ = 1, M_
          row_tau_(:) = row_tau_(:) + V_modes_(i_,im_) * Taumode_(im_,:)
        enddo
        do j_ = 1, N_
          sigma_ij_ = sum(row_tau_(:) * conjg(V_modes_(j_,:)))
          self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) = &
            self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) + &
            sigma_ij_
        enddo
      enddo
      !
      deallocate(G_bath_mode_, Kmode_, Taumode_, row_tau_)
    end subroutine
    !
    subroutine low_concentration_support_self_dense(G_bath_, V_support_, iR_large_, &
      dof_R_, dof_cart_, self_R_, self_support_)
      ! Support-space low-concentration solver:
      !   Sigma = c V [1 - (1-c) G_bath V]^{-1}
      ! where G_bath contains the DCA moment-preserving Weiss bath.
      complex(dp), intent(in) :: G_bath_(:,:), V_support_(:,:)
      integer, intent(in) :: iR_large_(:,:)
      integer, intent(in) :: dof_R_(:), dof_cart_(:)
      complex(dp), allocatable, intent(out) :: self_R_(:,:,:)
      complex(dp), intent(out) :: self_support_(:,:)
      !
      complex(dp), allocatable :: work_(:,:)
      integer :: i_, j_, N_, nR_large_
      !
      N_ = size(dof_R_)
      nR_large_ = maxval(iR_large_)
      if(size(dof_cart_) /= N_) call errore("low_concentration_support_self_dense", "bad active map", 1)
      if(size(G_bath_, 1) /= N_ .or. size(G_bath_, 2) /= N_) &
        call errore("low_concentration_support_self_dense", "bad support G", 1)
      if(size(V_support_, 1) /= N_ .or. size(V_support_, 2) /= N_) &
        call errore("low_concentration_support_self_dense", "bad support V", 1)
      if(size(self_support_, 1) /= N_ .or. size(self_support_, 2) /= N_) &
        call errore("low_concentration_support_self_dense", "bad support Sigma", 1)
      !
      call check_allocation_fits("DCA dense low-concentration solver work", &
        complex_mem_bytes(int(N_, int64) * int(N_, int64)) + &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(nR_large_, int64)))
      allocate(work_(N_,N_))
      work_ = id_mat(N_) - (1._dp - defect_conc) * matmul(G_bath_, V_support_)
      call invzmat(N_, work_)
      self_support_ = matmul(defect_conc * V_support_, work_)
      !
      allocate(self_R_(S%nat3,S%nat3,nR_large_))
      self_R_ = 0._dp
      do j_ = 1, N_
        do i_ = 1, N_
          self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) = &
            self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) + &
            self_support_(i_,j_)
        enddo
      enddo
      !
      deallocate(work_)
    end subroutine
    !
    subroutine low_concentration_support_self_R(G0i_cluster_, V_support_, iR_large_, diff_cart_, R_cart_, &
      dof_R_, dof_cart_, self_R_)
      ! Stable full-support form of the low-concentration average.  It is
      ! algebraically equivalent to Sigma = g^-1 - [(1-c)g+cG1]^-1, but avoids
      ! inverting the nearly singular support-space <G> directly.  The clean
      ! support Green function has rank <= Nc*nat3, so use Woodbury instead of
      ! an Nsupport x Nsupport inverse.
      complex(dp), intent(in) :: G0i_cluster_(:,:,:)
      complex(dp), intent(in) :: V_support_(:,:)
      integer, intent(in) :: iR_large_(:,:)
      real(dp), intent(in) :: diff_cart_(:,:), R_cart_(:,:)
      integer, intent(in) :: dof_R_(:), dof_cart_(:)
      complex(dp), allocatable, intent(out) :: self_R_(:,:,:)
      !
      complex(dp), allocatable :: A_low_(:,:), FV_(:,:), GV_(:,:), BA_(:,:), &
        K_(:,:), VA_(:,:), Sigma_(:,:)
      complex(dp) :: gQ_(S%nat3,S%nat3,Nc), block_(S%nat3,S%nat3)
      integer :: iR_, nR_, nR_large_, i_, j_, iq_, N_, nflat_, idx_, row0_
      real(dp), parameter :: support_tol_ = 1e-14_dp
      !
      do iq_ = 1, Nc
        block_ = G0i_cluster_(:,:,iq_)
        call invzmat(S%nat3, block_)
        gQ_(:,:,iq_) = block_
      enddo
      !
      nR_ = size(iR_large_, 1)
      nR_large_ = size(diff_cart_, 2)
      N_ = size(dof_R_)
      nflat_ = S%nat3 * Nc
      if(size(iR_large_, 2) /= nR_) call errore("low_concentration_support_self_R", "bad support map", 1)
      if(size(R_cart_, 2) /= nR_) call errore("low_concentration_support_self_R", "bad R list", 1)
      if(size(dof_cart_) /= N_) call errore("low_concentration_support_self_R", "bad active support map", 1)
      if(any(dof_R_ < 1) .or. any(dof_R_ > nR_)) &
        call errore("low_concentration_support_self_R", "active R out of range", 1)
      if(any(dof_cart_ < 1) .or. any(dof_cart_ > S%nat3)) &
        call errore("low_concentration_support_self_R", "active cartesian index out of range", 1)
      if(size(V_support_, 1) /= N_ .or. size(V_support_, 2) /= N_) &
        call errore("low_concentration_support_self_R", "bad support V", 1)
      !
      call check_allocation_fits("DCA full-support low-concentration solver", &
        complex_mem_bytes(3_int64 * int(N_, int64) * int(nflat_, int64)) + &
        complex_mem_bytes(2_int64 * int(nflat_, int64) * int(nflat_, int64)) + &
        complex_mem_bytes(int(N_, int64) * int(N_, int64)) + &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(nR_large_, int64)))
      allocate(A_low_(N_,nflat_), FV_(nflat_,N_), GV_(nflat_,N_), &
        BA_(nflat_,nflat_), K_(nflat_,nflat_), VA_(N_,nflat_), Sigma_(N_,N_))
      A_low_ = 0._dp
      FV_ = 0._dp
      VA_ = 0._dp
      !
      ! g0_support = A_low * G(Q) * A_low^dagger, using the same
      ! unnormalised Fourier convention as full_born_center.
      do i_ = 1, N_
        do iq_ = 1, Nc
          idx_ = dof_cart_(i_) + (iq_ - 1) * S%nat3
          A_low_(i_,idx_) = e_iqr(xq(:,iq_), R_cart_(:,dof_R_(i_)))
        enddo
      enddo
      !
      ! Sparse products with V_support: FV = A_low^dagger V and VA = V A_low.
      do j_ = 1, N_
        do i_ = 1, N_
          if(abs(V_support_(i_,j_)) <= support_tol_) cycle
          do iq_ = 1, Nc
            FV_(dof_cart_(i_) + (iq_ - 1) * S%nat3,j_) = &
              FV_(dof_cart_(i_) + (iq_ - 1) * S%nat3,j_) + &
              e_iqr(xq(:,iq_), -R_cart_(:,dof_R_(i_))) * V_support_(i_,j_)
            VA_(i_,dof_cart_(j_) + (iq_ - 1) * S%nat3) = &
              VA_(i_,dof_cart_(j_) + (iq_ - 1) * S%nat3) + &
              V_support_(i_,j_) * e_iqr(xq(:,iq_), R_cart_(:,dof_R_(j_)))
          enddo
        enddo
      enddo
      !
      GV_ = 0._dp
      do iq_ = 1, Nc
        row0_ = (iq_ - 1) * S%nat3
        GV_(row0_+1:row0_+S%nat3,:) = &
          matmul(gQ_(:,:,iq_), FV_(row0_+1:row0_+S%nat3,:))
      enddo
      !
      BA_ = matmul(GV_, A_low_)
      K_ = id_mat(nflat_) - (1._dp - defect_conc) * BA_
      call invzmat(nflat_, K_)
      Sigma_ = defect_conc * (V_support_ + &
        (1._dp - defect_conc) * matmul(VA_, matmul(K_, GV_)))
      !
      allocate(self_R_(S%nat3,S%nat3,nR_large_))
      self_R_ = 0._dp
      do j_ = 1, N_
        do i_ = 1, N_
          self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) = &
            self_R_(dof_cart_(i_),dof_cart_(j_),iR_large_(dof_R_(i_),dof_R_(j_))) + &
            Sigma_(i_,j_)
        enddo
      enddo
      !
      deallocate(A_low_, FV_, GV_, BA_, K_, VA_, Sigma_)
    end subroutine
    !
    subroutine project_support_self(self_R_, diff_cart_, xq_target_, self_cart_)
      complex(dp), intent(in) :: self_R_(:,:,:)
      real(dp), intent(in) :: diff_cart_(:,:)
      real(dp), intent(in) :: xq_target_(:,:)
      complex(dp), intent(out) :: self_cart_(:,:,:)
      !
      integer :: iq_, iR_, ntarget_, nR_
      !
      ntarget_ = size(xq_target_, 2)
      nR_ = size(diff_cart_, 2)
      if(size(self_R_, 3) /= nR_) &
        call errore("project_support_self", "inconsistent support sizes", 1)
      if(size(self_cart_, 1) /= S%nat3 .or. size(self_cart_, 2) /= S%nat3 .or. &
        size(self_cart_, 3) /= ntarget_) &
        call errore("project_support_self", "wrong output size", 1)
      !
      self_cart_ = 0._dp
      do iq_ = 1, ntarget_
        do iR_ = 1, nR_
          self_cart_(:,:,iq_) = self_cart_(:,:,iq_) + &
            self_R_(:,:,iR_) * e_iqr(xq_target_(:,iq_), -diff_cart_(:,iR_))
        enddo
      enddo
    end subroutine
    !
    subroutine pair_diff_projector(c_pair, fc2sc, xq_cluster, grid)
      ! Build the diagonal-Q interpolation on the same finite support used by
      ! full_born_center, but separately for each cartesian matrix element.
      complex(dp), allocatable, intent(out) :: c_pair(:,:,:,:)
      type(forceconst2_sc), intent(in) :: fc2sc
      real(dp), intent(in) :: xq_cluster(:,:)
      type(q_grid), intent(in) :: grid
      !
      integer, allocatable :: R_list(:,:), diff_pair(:,:,:,:), n_diff(:,:)
      integer, allocatable :: res_count(:)
      logical, allocatable :: row_has(:,:), col_has(:,:)
      complex(dp), allocatable :: gram(:,:), inv_gram(:,:), h(:), res_phase(:,:)
      real(dp), allocatable :: res_cart(:,:)
      integer :: max_R, max_diff_pair, nR, iR1_, iR2_, iR_, jR_, idx_
      integer :: i_, j_, kQ_, lQ_, iq_out_, iD_, max_n_diff
      integer :: ires_, Rmod_(3), max_n_res, max_alias
      integer :: Rdiff_(3)
      real(dp) :: Rdiff_cart_(3), support_tol_, diag_max_, reg_
      !
      support_tol_ = 1e-14_dp
      max_R = fc2sc%n_R2 + sum(fc2sc%n_R1)
      allocate(R_list(3,max_R))
      nR = 0
      do iR2_ = 1, fc2sc%n_R2
        call add_unique_R_int(R_list, nR, fc2sc%yR2(:,iR2_), idx_)
        do iR1_ = 1, fc2sc%n_R1(iR2_)
          call add_unique_R_int(R_list, nR, fc2sc%yR1(:,iR1_,iR2_), idx_)
        enddo
      enddo
      !
      allocate(row_has(S%nat3,nR), col_has(S%nat3,nR))
      row_has = .false.
      col_has = .false.
      do iR2_ = 1, fc2sc%n_R2
        call find_where(fc2sc%yR2(:,iR2_), R_list(:,:nR), jR_)
        do iR1_ = 1, fc2sc%n_R1(iR2_)
          call find_where(fc2sc%yR1(:,iR1_,iR2_), R_list(:,:nR), iR_)
          do i_ = 1, S%nat3
            if(any(abs(fc2sc%FC(i_,:,iR1_,iR2_)) > support_tol_)) row_has(i_,iR_) = .true.
          enddo
          do j_ = 1, S%nat3
            if(any(abs(fc2sc%FC(:,j_,iR1_,iR2_)) > support_tol_)) col_has(j_,jR_) = .true.
          enddo
        enddo
      enddo
      !
      max_diff_pair = nR * nR
      allocate(diff_pair(3,max_diff_pair,S%nat3,S%nat3))
      allocate(n_diff(S%nat3,S%nat3))
      n_diff = 0
      do j_ = 1, S%nat3
        do i_ = 1, S%nat3
          do jR_ = 1, nR
            if(.not. col_has(j_,jR_)) cycle
            do iR_ = 1, nR
              if(.not. row_has(i_,iR_)) cycle
              Rdiff_ = R_list(:,iR_) - R_list(:,jR_)
              call add_unique_R_int(diff_pair(:,:,i_,j_), n_diff(i_,j_), Rdiff_, idx_)
            enddo
          enddo
        enddo
      enddo
      max_n_diff = maxval(n_diff)
      if(ionode) print*, "Max pair-specific R differences:", max_n_diff
      !
      allocate(c_pair(S%nat3,S%nat3,Nc,grid%nqtot))
      c_pair = 0._dp
      allocate(gram(Nc,Nc), inv_gram(Nc,Nc), h(Nc))
      allocate(res_count(Nc), res_phase(Nc,grid%nqtot), res_cart(3,Nc))
      do ires_ = 1, Nc
        res_cart(:,ires_) = cryst2cart(real(index2v(ires_, cluster_mesh), dp), S%at, 1)
      enddo
      max_n_res = 0
      max_alias = 0
      do j_ = 1, S%nat3
        do i_ = 1, S%nat3
          if(n_diff(i_,j_) == 0) cycle
          res_count = 0
          res_phase = 0._dp
          do iD_ = 1, n_diff(i_,j_)
            Rmod_ = bz2simple(diff_pair(:,iD_,i_,j_), cluster_mesh)
            ires_ = v2index(Rmod_, cluster_mesh)
            res_count(ires_) = res_count(ires_) + 1
            Rdiff_cart_ = cryst2cart(real(diff_pair(:,iD_,i_,j_),dp), S%at, 1)
            do iq_out_ = 1, grid%nqtot
              res_phase(ires_,iq_out_) = res_phase(ires_,iq_out_) + &
                e_iqr(-grid%xq(:,iq_out_), Rdiff_cart_)
            enddo
          enddo
          max_n_res = max(max_n_res, count(res_count > 0))
          max_alias = max(max_alias, maxval(res_count))
          !
          gram = 0._dp
          do lQ_ = 1, Nc
            do kQ_ = 1, Nc
              do ires_ = 1, Nc
                if(res_count(ires_) == 0) cycle
                gram(kQ_,lQ_) = gram(kQ_,lQ_) + &
                  res_count(ires_) * e_iqr(xq_cluster(:,lQ_) - xq_cluster(:,kQ_), res_cart(:,ires_))
              enddo
            enddo
          enddo
          diag_max_ = 0._dp
          do kQ_ = 1, Nc
            diag_max_ = max(diag_max_, abs(gram(kQ_,kQ_)))
          enddo
          reg_ = max(1e-18_dp, 1e-12_dp * diag_max_)
          inv_gram = gram
          do kQ_ = 1, Nc
            inv_gram(kQ_,kQ_) = inv_gram(kQ_,kQ_) + cmplx(reg_, 0._dp, dp)
          enddo
          call invzmat(Nc, inv_gram)
          do iq_out_ = 1, grid%nqtot
            h = 0._dp
            do lQ_ = 1, Nc
              do ires_ = 1, Nc
                if(res_count(ires_) == 0) cycle
                h(lQ_) = h(lQ_) + res_phase(ires_,iq_out_) * &
                  e_iqr(xq_cluster(:,lQ_), res_cart(:,ires_))
              enddo
            enddo
            c_pair(i_,j_,:,iq_out_) = matmul(h, inv_gram)
          enddo
        enddo
      enddo
      if(ionode) print*, "Max pair-specific cluster residues:", max_n_res
      if(ionode) print*, "Max alias count per residue:", max_alias
      deallocate(R_list, row_has, col_has, diff_pair, n_diff, gram, inv_gram, h, &
        res_count, res_phase, res_cart)
    end subroutine
    !
    subroutine add_unique_R_int(R_list, nR, R_in, idx_out)
      integer, intent(inout) :: R_list(:,:)
      integer, intent(inout) :: nR
      integer, intent(in) :: R_in(3)
      integer, intent(out) :: idx_out
      !
      integer :: i_
      !
      do i_ = 1, nR
        if(all(R_list(:,i_) == R_in)) then
          idx_out = i_
          return
        endif
      enddo
      nR = nR + 1
      if(nR > size(R_list,2)) call errore("add_unique_R_int", "R list is too short", 1)
      R_list(:,nR) = R_in
      idx_out = nR
    end subroutine
    !
    subroutine diff_list(fc2sc, at, diff_cart)
      type(forceconst2_sc), intent(in) :: fc2sc
      real(dp), intent(in) :: at(3, 3)
      integer, pointer :: R_list(:,:)
      integer, pointer :: diff_cryst(:,:)!, ir_large(:,:)
      integer :: nR, iR1, iR2, RR(3), nr_large
      real(dp), allocatable :: diff_cart(:,:)
      !
      allocate(R_list(3,1))
      R_list(:,1) = fc2sc%yR2(:,1)
      nR = 1
      !
      do iR2 = 1, fc2sc%n_R2
        call find_where(fc2sc%yR2(:,iR2), R_list, iR)
        if(iR == -1) call enlarge_R(fc2sc%yR2(:,iR2), R_list, nR)
        do iR1 = 1, fc2sc%n_R1(iR2)
          call find_where(fc2sc%yR1(:,iR1,iR2), R_list, iR)
          if(iR == -1) call enlarge_R(fc2sc%yR1(:,iR1,iR2), R_list, nR)
        enddo
      enddo
      if(ionode) print*, "Number of unique R vectors: ", nR
      !
      nR_large = 1
      allocate(diff_cryst(3,1))
      diff_cryst(:,1) = [0,0,0]
      do i = 1, nR
        do j = 1, nR
          RR = R_list(:,i) - R_list(:,j)
          call find_where(RR, diff_cryst, iR)
          ! iR_large(i,j) = iR
          if (iR == -1) then
            call enlarge_R(RR, diff_cryst, nR_large)
            ! ir_large(i,j) = nR_large
          endif
        enddo
      enddo
      print*, "Number of unique R differences: ", nR_large
      !
      allocate(diff_cart(3,nR_large))
      diff_cart = cryst2cart(real(diff_cryst,dp), at, 1)
    end subroutine
    !
    subroutine center_equivalent_defect_fc(site_, fc_centered_)
      use symm_base, only : irt, nsym, ft, symm_mat => s, invs
      !
      integer, intent(in) :: site_
      type(forceconst2_sc), intent(inout) :: fc_centered_
      !
      real(dp), allocatable :: tens4_(:,:,:,:,:,:), work_(:,:,:,:,:,:)
      real(dp) :: tau_crys_(3), tau_rot_crys_(3), G_crys_(3)
      integer, allocatable :: G_atom_(:,:)
      integer, dimension(3,3) :: SM_, SMT_, SM1_, SMT1_
      integer :: N_, count_, isym_map_, isym_, inv_map_
      integer :: nai_, naj_, naii_, najj_, iR_, jR_, iiR_, jjR_
      integer :: Ri_(3), Rj_(3)
      !
      if(site_ < 1 .or. site_ > size(fc2_sc%defects, 2)) &
        call errore("center_equivalent_defect_fc", "bad equivalent site", 1)
      N_ = product(fc2_sc%nq)
      call fc_centered_%allocate(S, S_sc, fc2_sc%nq)
      call check_allocation_fits("DCA equivalent-site centered support", &
        real_mem_bytes(2_int64 * 9_int64 * int(S%nat, int64) * int(S%nat, int64) * &
        int(N_, int64) * int(N_, int64)))
      allocate(tens4_(3, S%nat, 3, S%nat, N_, N_))
      allocate(work_(3, S%nat, 3, S%nat, N_, N_))
      allocate(G_atom_(3, S%nat))
      !
      tens4_ = reshape(fc2_sc%fc, [3, S%nat, 3, S%nat, N_, N_])
      call transform_tns4_cart(tens4_, -1)
      count_ = 0
      do isym_ = 1, nsym
        if(irt(isym_, fc2_sc%defects(1,1)) == fc2_sc%defects(1,site_)) then
          isym_map_ = isym_
          inv_map_ = invs(isym_map_)
          SM1_ = symm_mat(:,:,inv_map_)
          SMT1_ = transpose(SM1_)
          SM_ = symm_mat(:,:,isym_map_)
          SMT_ = transpose(SM_)
          count_ = count_ + 1
        endif
      enddo
      if(count_ == 0) call errore("center_equivalent_defect_fc", "could not map defect site", 1)
      !
      do nai_ = 1, S%nat
        naii_ = irt(isym_map_, nai_)
        tau_crys_ = cryst2cart(S%tau(:,nai_), S%bg, -1)
        tau_rot_crys_ = cryst2cart(S%tau(:,naii_), S%bg, -1)
        G_crys_ = matmul(SMT_, tau_crys_) - ft(:,isym_map_) - tau_rot_crys_
        G_atom_(:,nai_) = nint(G_crys_)
        if(any(abs(G_crys_ - real(G_atom_(:,nai_), dp)) > 1e-6_dp)) then
          print*, "G_crys:", G_crys_, " rounded:", G_atom_(:,nai_)
          call errore("center_equivalent_defect_fc", "non-integer atom-dependent symmetry shift", 1)
        endif
      enddo
      !
      work_ = 0._dp
      do nai_ = 1, S%nat
        naii_ = irt(isym_map_, nai_)
        do naj_ = 1, S%nat
          najj_ = irt(isym_map_, naj_)
          do iR_ = 1, N_
            do jR_ = 1, N_
              Ri_ = matmul(SMT_, index2v(iR_, fc2_sc%nq)) + G_atom_(:,nai_)
              Rj_ = matmul(SMT_, index2v(jR_, fc2_sc%nq)) + G_atom_(:,naj_)
              iiR_ = v2index(bz2simple(Ri_, fc2_sc%nq), fc2_sc%nq)
              jjR_ = v2index(bz2simple(Rj_, fc2_sc%nq), fc2_sc%nq)
              work_(:,naii_,:,najj_,iiR_,jjR_) = work_(:,naii_,:,najj_,iiR_,jjR_) + &
                matmul(matmul(SM1_, tens4_(:,nai_,:,naj_,iR_,jR_)), SMT1_)
            enddo
          enddo
        enddo
      enddo
      !
      call transform_tns4_cart(work_, 1)
      fc_centered_%fc = reshape(work_, [3*S%nat, 3*S%nat, N_, N_])
      fc_centered_%taudef = fc2_sc%taudef + S%tau(:,fc2_sc%defects(1,site_)) - &
        S%tau(:,fc2_sc%defects(1,1))
      call fc_centered_%center(fc2_sc%nq, S)
      deallocate(tens4_, work_, G_atom_)
    end subroutine
    !
    subroutine center_V_mixed(xq_left, xq_right, Vmix)
      ! Build V(q_left, q_right) after applying the same centering used for
      ! the cluster V(Q,Q') entering Gi_avg.
      use symm_base, only : irt, nsym, ft, symm_mat => s, invs
      !
      real(dp), intent(in) :: xq_left(:,:), xq_right(:,:)
      complex(dp), allocatable, intent(out) :: Vmix(:,:,:,:)
      !
      type(forceconst2_sc) :: fc_temp
      real(dp), allocatable :: tens4(:,:,:,:,:,:), work(:,:,:,:,:,:)
      real(dp) :: trans(3), tau_crys(3), tau_rot_crys(3), G_crys(3)
      integer, allocatable :: G_atom(:,:)
      integer, dimension(3,3) :: SM, SMT, SM1, SMT1
      integer :: N_, nql_, nqr_, site_, count_
      integer :: isym_map_, isym_, inv_map_
      integer :: nai_, naj_, naii_, najj_, iR_, jR_, iiR_, jjR_, iq_, jq_
      integer :: Ri_(3), Rj_(3)
      !
      nql_ = size(xq_left, 2)
      nqr_ = size(xq_right, 2)
      N_ = product(fc2_sc%nq)
      call check_allocation_fits("DCA mixed centered V interpolation", &
        complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(nql_, int64) * int(nqr_, int64)) + &
        real_mem_bytes(2_int64 * 9_int64 * int(S%nat, int64) * int(S%nat, int64) * &
        int(N_, int64) * int(N_, int64)))
      allocate(Vmix(S%nat3, S%nat3, nql_, nqr_))
      allocate(tens4(3, S%nat, 3, S%nat, N_, N_))
      allocate(work(3, S%nat, 3, S%nat, N_, N_))
      allocate(G_atom(3, S%nat))
      !
      tens4 = reshape(fc2_sc%fc, [3, S%nat, 3, S%nat, N_, N_])
      call transform_tns4_cart(tens4, -1)
      !
      site_ = 1
      call fc_temp%allocate(S, S_sc, fc2_sc%nq)
      count_ = 0
      do isym_ = 1, nsym
        if (irt(isym_, fc2_sc%defects(1,1)) == fc2_sc%defects(1,site_)) then
          isym_map_ = isym_
          inv_map_ = invs(isym_map_)
          SM1 = symm_mat(:,:,inv_map_)
          SMT1 = transpose(SM1)
          SM = symm_mat(:,:,isym_map_)
          SMT = transpose(SM)
          trans = cryst2cart(ft(:,isym_map_), S%at, 1)
          count_ = count_ + 1
        endif
      enddo
      if (count_ == 0) call errore("center_V_mixed", "could not map defect site", 1)
      !
      do nai_ = 1, S%nat
        naii_ = irt(isym_map_, nai_)
        tau_crys = cryst2cart(S%tau(:,nai_), S%bg, -1)
        tau_rot_crys = cryst2cart(S%tau(:,naii_), S%bg, -1)
        G_crys = matmul(SMT, tau_crys) - ft(:,isym_map_) - tau_rot_crys
        G_atom(:,nai_) = nint(G_crys)
        if(any(abs(G_crys - real(G_atom(:,nai_), dp)) > 1e-6_dp)) then
          print*, "G_crys:", G_crys, " rounded:", G_atom(:,nai_)
          call errore("center_V_mixed", "non-integer atom-dependent symmetry shift", 1)
        endif
      enddo
      !
      work = 0.0_dp
      do nai_ = 1, S%nat
        naii_ = irt(isym_map_, nai_)
        do naj_ = 1, S%nat
          najj_ = irt(isym_map_, naj_)
          do iR_ = 1, N_
            do jR_ = 1, N_
              Ri_ = matmul(SMT, index2v(iR_, fc2_sc%nq)) + G_atom(:,nai_)
              Rj_ = matmul(SMT, index2v(jR_, fc2_sc%nq)) + G_atom(:,naj_)
              iiR_ = v2index(bz2simple(Ri_, fc2_sc%nq), fc2_sc%nq)
              jjR_ = v2index(bz2simple(Rj_, fc2_sc%nq), fc2_sc%nq)
              work(:,naii_,:,najj_,iiR_,jjR_) = work(:,naii_,:,najj_,iiR_,jjR_) + matmul( &
                matmul(SM1, tens4(:,nai_,:,naj_,iR_,jR_)), SMT1)
            enddo
          enddo
        enddo
      enddo
      !
      call transform_tns4_cart(work, 1)
      fc_temp%fc = reshape(work, [3*S%nat, 3*S%nat, N_, N_])
      fc_temp%taudef = fc2_sc%taudef + S%tau(:,fc2_sc%defects(1,site_)) - &
        S%tau(:, fc2_sc%defects(1,1))
      call fc_temp%center(fc2_sc%nq, S)
      !
      do iq_ = 1, nql_
        call fc_temp%r2q(xq_left(:,iq_))
        do jq_ = 1, nqr_
          call fc_temp%r2q(xq_right(:,jq_), Vmix(:,:,iq_,jq_))
        enddo
      enddo
      !
      call fc_temp%deallocate()
      deallocate(tens4, work, G_atom)
    end subroutine
    !
    subroutine fit_potential_cout_full(c_out_, V_cluster_flat, V_target)
      ! Full row-space interpolation for the centered potential:
      ! V_target(q,P) ~= C(q) V_cluster(Q,P), with C acting on
      ! the combined cartesian/cluster row index.
      complex(dp), allocatable, intent(out) :: c_out_(:,:,:)
      complex(dp), intent(in) :: V_cluster_flat(:,:)
      complex(dp), intent(in) :: V_target(:,:,:,:)
      !
      complex(dp), allocatable :: metric(:,:), inv_metric(:,:), target_flat(:,:), approx(:,:)
      integer :: nrow_, n_target_, itarget_, iP_, j_
      real(dp) :: diag_max_, reg_, res_norm_, target_norm_, rel_, max_rel_
      !
      nrow_ = size(V_cluster_flat, 1)
      n_target_ = size(V_target, 3)
      if(size(V_cluster_flat, 2) /= nrow_) &
        call errore("fit_potential_cout_full", "cluster V is not square", 1)
      if(nrow_ /= S%nat3 * Nc) &
        call errore("fit_potential_cout_full", "cluster V has wrong size", 1)
      if(size(V_target, 4) /= Nc) &
        call errore("fit_potential_cout_full", "V_target right leg is not on the cluster", 1)
      !
      call check_allocation_fits("DCA full left interpolation fit", &
        complex_mem_bytes(int(S%nat3, int64) * int(nrow_, int64) * int(n_target_, int64)) + &
        complex_mem_bytes(2_int64 * int(nrow_, int64) * int(nrow_, int64)) + &
        complex_mem_bytes(2_int64 * int(S%nat3, int64) * int(nrow_, int64)))
      allocate(c_out_(S%nat3, nrow_, n_target_))
      allocate(metric(nrow_, nrow_), inv_metric(nrow_, nrow_))
      allocate(target_flat(S%nat3, nrow_), approx(S%nat3, nrow_))
      !
      metric = matmul(V_cluster_flat, conjg(transpose(V_cluster_flat)))
      diag_max_ = 0._dp
      do iP_ = 1, nrow_
        diag_max_ = max(diag_max_, abs(metric(iP_,iP_)))
      enddo
      reg_ = max(1e-18_dp, 1e-12_dp * diag_max_)
      inv_metric = metric
      do iP_ = 1, nrow_
        inv_metric(iP_,iP_) = inv_metric(iP_,iP_) + cmplx(reg_, 0._dp, dp)
      enddo
      call invzmat(nrow_, inv_metric)
      !
      max_rel_ = 0._dp
      do itarget_ = 1, n_target_
        do iP_ = 1, Nc
          do j_ = 1, S%nat3
            target_flat(:, j_ + (iP_-1)*S%nat3) = V_target(:,j_,itarget_,iP_)
          enddo
        enddo
        c_out_(:,:,itarget_) = matmul( &
          matmul(target_flat, conjg(transpose(V_cluster_flat))), inv_metric)
        approx = matmul(c_out_(:,:,itarget_), V_cluster_flat)
        !
        res_norm_ = sum(abs(approx - target_flat)**2)
        target_norm_ = sum(abs(target_flat)**2)
        rel_ = sqrt(res_norm_ / max(target_norm_, 1e-300_dp))
        max_rel_ = max(max_rel_, rel_)
      enddo
      if(ionode) print"(A,E12.4,A,E12.4)", &
        "full c_out fit max relative residual: ", max_rel_, " regularization: ", reg_
      !
      deallocate(metric, inv_metric, target_flat, approx)
    end subroutine
    !
    subroutine fit_potential_cout_right(c_right_, V_cluster_flat, V_target)
      ! Full column-space interpolation for the centered potential:
      ! V_target(P,q) ~= V_cluster(P,Q) C_right(q).
      complex(dp), allocatable, intent(out) :: c_right_(:,:,:)
      complex(dp), intent(in) :: V_cluster_flat(:,:)
      complex(dp), intent(in) :: V_target(:,:,:,:)
      !
      complex(dp), allocatable :: metric(:,:), inv_metric(:,:), target_flat(:,:), approx(:,:)
      integer :: nrow_, n_target_, itarget_, iP_, j_
      real(dp) :: diag_max_, reg_, res_norm_, target_norm_, rel_, max_rel_
      !
      nrow_ = size(V_cluster_flat, 1)
      n_target_ = size(V_target, 4)
      if(size(V_cluster_flat, 2) /= nrow_) &
        call errore("fit_potential_cout_right", "cluster V is not square", 1)
      if(nrow_ /= S%nat3 * Nc) &
        call errore("fit_potential_cout_right", "cluster V has wrong size", 1)
      if(size(V_target, 3) /= Nc) &
        call errore("fit_potential_cout_right", "V_target left leg is not on the cluster", 1)
      !
      call check_allocation_fits("DCA full right interpolation fit", &
        complex_mem_bytes(int(nrow_, int64) * int(S%nat3, int64) * int(n_target_, int64)) + &
        complex_mem_bytes(2_int64 * int(nrow_, int64) * int(nrow_, int64)) + &
        complex_mem_bytes(2_int64 * int(nrow_, int64) * int(S%nat3, int64)))
      allocate(c_right_(nrow_, S%nat3, n_target_))
      allocate(metric(nrow_, nrow_), inv_metric(nrow_, nrow_))
      allocate(target_flat(nrow_, S%nat3), approx(nrow_, S%nat3))
      !
      metric = matmul(conjg(transpose(V_cluster_flat)), V_cluster_flat)
      diag_max_ = 0._dp
      do iP_ = 1, nrow_
        diag_max_ = max(diag_max_, abs(metric(iP_,iP_)))
      enddo
      reg_ = max(1e-18_dp, 1e-12_dp * diag_max_)
      inv_metric = metric
      do iP_ = 1, nrow_
        inv_metric(iP_,iP_) = inv_metric(iP_,iP_) + cmplx(reg_, 0._dp, dp)
      enddo
      call invzmat(nrow_, inv_metric)
      !
      max_rel_ = 0._dp
      do itarget_ = 1, n_target_
        do iP_ = 1, Nc
          do j_ = 1, S%nat3
            target_flat(j_ + (iP_-1)*S%nat3, :) = V_target(j_, :, iP_, itarget_)
          enddo
        enddo
        c_right_(:,:,itarget_) = matmul(inv_metric, &
          matmul(conjg(transpose(V_cluster_flat)), target_flat))
        approx = matmul(V_cluster_flat, c_right_(:,:,itarget_))
        !
        res_norm_ = sum(abs(approx - target_flat)**2)
        target_norm_ = sum(abs(target_flat)**2)
        rel_ = sqrt(res_norm_ / max(target_norm_, 1e-300_dp))
        max_rel_ = max(max_rel_, rel_)
      enddo
      if(ionode) print"(A,E12.4,A,E12.4)", &
        "right c_out fit max relative residual: ", max_rel_, " regularization: ", reg_
      !
      deallocate(metric, inv_metric, target_flat, approx)
    end subroutine
    !
    subroutine transform_tns4_cart(tns4, flag)
      real(dp), intent(inout) :: tns4(:,:,:,:,:,:)
      integer, intent(in) :: flag
      !
      integer :: nai_, naj_, iR_, jR_
      !
      do nai_ = 1, S%nat
        do naj_ = 1, S%nat
          do iR_ = 1, size(tns4, 5)
            do jR_ = 1, size(tns4, 6)
              call transform_mat_cart(tns4(:,nai_,:,naj_,iR_,jR_), flag)
            enddo
          enddo
        enddo
      enddo
    end subroutine
    !
    subroutine transform_mat_cart(mat, flag)
      real(dp), intent(inout) :: mat(:,:)
      integer, intent(in) :: flag
      real(dp) :: U_(3,3)
      !
      if (flag == 1) then
        U_ = S%bg
      else
        U_ = transpose(S%at)
      endif
      mat = matmul(U_, matmul(mat, transpose(U_)))
    end subroutine
    !
    subroutine scalar_fourier_basis(c_Qq, grid)
      complex(dp), allocatable, intent(out) :: c_Qq(:,:)
      type(q_grid), intent(in) :: grid
      !
      integer :: iq_, jq_, iR_, n_target
      !
      n_target = size(grid%xq, 2)
      allocate(c_Qq(Nc, n_target))
      c_Qq = 0._dp
      do jq_ = 1, n_target
        do iq_ = 1, Nc
          do iR_ = 1, Nc
            c_Qq(iq_,jq_) = c_Qq(iq_,jq_) + &
              e_iqr(xq(:,iq_) - grid%xq(:,jq_), R(:,iR_)) / Nc
          enddo
        enddo
      enddo
    end subroutine
    !
    subroutine atom_fourier_basis(c_Qq, grid)
      complex(dp), allocatable, intent(out) :: c_Qq(:,:,:)
      type(q_grid), intent(in) :: grid
      !
      integer :: iq_, jq_, na_, iR_, n_target
      integer :: l1_, l2_, l3_, nRbig_, nRout_, far_
      real(dp) :: atws_(3,3), dist_(3), wg_, totalweight_
      real(dp), allocatable :: Rbig_(:,:)
      integer, parameter :: nrwsx_ = 5000
      integer :: nrws_
      real(dp) :: rws_(0:3,nrwsx_)
      real(dp), external :: wsweight
      !
      far_ = 2
      n_target = size(grid%xq, 2)
      allocate(c_Qq(S%nat, Nc, n_target))
      c_Qq = 0._dp
      !
      atws_(:,1) = cluster_mesh(1) * S%at(:,1)
      atws_(:,2) = cluster_mesh(2) * S%at(:,2)
      atws_(:,3) = cluster_mesh(3) * S%at(:,3)
      call wsinit(rws_, nrwsx_, nrws_, atws_)
      !
      nRbig_ = (2*far_*cluster_mesh(1)+1) * &
        (2*far_*cluster_mesh(2)+1) * (2*far_*cluster_mesh(3)+1)
      allocate(Rbig_(3,nRbig_))
      nRbig_ = 0
      do l1_ = -far_*cluster_mesh(1), far_*cluster_mesh(1)
        do l2_ = -far_*cluster_mesh(2), far_*cluster_mesh(2)
          do l3_ = -far_*cluster_mesh(3), far_*cluster_mesh(3)
            nRbig_ = nRbig_ + 1
            Rbig_(:,nRbig_) = S%at(:,1)*l1_ + S%at(:,2)*l2_ + S%at(:,3)*l3_
          enddo
        enddo
      enddo
      !
      do na_ = 1, S%nat
        totalweight_ = 0._dp
        nRout_ = 0
        do iR_ = 1, nRbig_
          dist_ = Rbig_(:,iR_) + S%tau(:,na_) - fc2_sc%taudef
          wg_ = wsweight(dist_, rws_, nrws_)
          if(wg_ > 0._dp) then
            nRout_ = nRout_ + 1
            totalweight_ = totalweight_ + wg_ / real(Nc, dp)
            do jq_ = 1, n_target
              do iq_ = 1, Nc
                c_Qq(na_,iq_,jq_) = c_Qq(na_,iq_,jq_) + &
                  e_iqr(xq(:,iq_) - grid%xq(:,jq_), Rbig_(:,iR_)) * wg_ / Nc
              enddo
            enddo
          endif
        enddo
        if(abs(totalweight_ - 1._dp) > 1e-8_dp) then
          print*, totalweight_, na_, nRout_
          call errore("atom_fourier_basis", "wrong totalweight", 1)
        endif
      enddo
      !
      deallocate(Rbig_)
    end subroutine
    !
    subroutine fourier_basis(c_Qq, grid)
      complex(dp), allocatable, intent(out) :: c_Qq(:,:,:,:)
      type(q_grid), intent(in) :: grid
      !
      integer :: iq_, jq_, na1_, na2_, iR_, n_target
      !
      n_target = size(grid%xq, 2)
      allocate(c_Qq(S%nat3, S%nat3, Nc, n_target))
      c_Qq = 0._dp
      do jq_ = 1, n_target
        do iq_ = 1, Nc
          do na2_ = 1, S%nat
            do na1_ = 1, S%nat
              do iR_ = 1, img_nR(na1_,na2_)
                c_Qq(3*(na1_-1)+1:3*na1_,3*(na2_-1)+1:3*na2_, iq_, jq_) = &
                  c_Qq(3*(na1_-1)+1:3*na1_,3*(na2_-1)+1:3*na2_, iq_, jq_) + &
                  e_iqr(xq(:,iq_) - grid%xq(:,jq_), img_xR(:,iR_,na1_,na2_)) * &
                  img_weight(iR_,na1_,na2_) / Nc
              enddo
            enddo
          enddo
        enddo
      enddo
    end subroutine
    !
  end subroutine
  !
  subroutine assign_defect_counts(NSAMPLES, Nc, conc, N_sites)
    implicit none
    integer, intent(in) :: NSAMPLES, Nc
    real(dp), intent(in) :: conc
    integer, intent(out) :: N_sites(NSAMPLES)
    real(dp) :: probs(0:Nc), cumsum(0:Nc+1), r
    integer :: i, j

    ! Compute log probs, then cumulative
    probs = 0.0_dp
    do i = 0, Nc
      probs(i) = log_binom(Nc, i, conc)  ! ln[ binom * c^i * (1-c)^{Nc-i} ]
    end do
    cumsum(0) = 0.0_dp
    do i = 1, Nc+1
      cumsum(i) = cumsum(i-1) + exp(probs(i-1))
    end do
    cumsum(Nc+1) = 1.0_dp + 1.0e-12_dp  ! numerical safety

    ! Assign via inverse CDF
    do j = 1, NSAMPLES
      call random_number(r)
      ! Find i such that cumsum(i) <= r < cumsum(i+1)
      do i = 0, Nc
        if (r < cumsum(i+1)) then
          N_sites(j) = i
          exit
        end if
      end do
    end do
  end subroutine
!
  function real_to_randint(rr)
    real(dp), intent(in) :: rr
    real(dp) :: x
    integer :: real_to_randint
    !
    call random_number(x)
    if(x > rr - floor(rr)) then
      real_to_randint = floor(rr)
    else
      real_to_randint = ceiling(rr)
    end if
    !
  end function
  !
  function log_binom(n, k, c)
    implicit none
    integer, intent(in) :: n, k
    real(dp), intent(in) :: c
    integer :: i
    real(dp) :: lb, log_binom
    !
    log_binom = log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - log_gamma(real(n-k+1,dp)) + k*log(c) + (n-k)*log(1-c)
    ! log_binom = exp(lb)
  end function
!
  subroutine sample_canonical(N, Npos, pos)
    integer, intent(in) :: N, Npos
    integer, intent(out) :: pos(N)
    !
    real(dp) :: r
    integer :: i, j, tmp
    !
    if (Npos > N) call errore("sample_canonical", "Npos cannot be larger than N", 1)
    do i = 1, N
      pos(i) = i
    enddo
    do i = 1, Npos
      call random_number(r)
      j = i + floor(r * (N - i + 1))
      tmp = pos(i)
      pos(i) = pos(j)
      pos(j) = tmp
    enddo
  end subroutine
  !
  subroutine sample_binomial(N, c, pos, npos)
    integer, intent(in) :: N
    real(dp), intent(in) :: c
    integer, intent(out) :: pos(N)
    integer, intent(out) :: npos
    !
    integer :: i
    real(dp) :: r
    !
    !
    npos = 0
    do i = 1, N
      call random_number(r)
      if (r < c) then
        npos = npos + 1
        pos(npos) = i
      end if
    enddo
    !
  end subroutine
!
  SUBROUTINE center_V(xq, S, S_sc, fc2_sc, Vqqs)
    !-----------------------------------------------------------------------
    ! Apply ONE symmetry that maps atom1 -> atom2 to center potential there
    !
    USE kinds,     ONLY : DP
    USE symm_base, ONLY : irt, nsym, ft, symm_mat => s, invs
    use symme, only : cart_to_crys, crys_to_cart
    !
    real(dp), intent(in) :: xq(:,:)
    type(ph_system_info), INTENT(IN) :: S, S_sc
    type(forceconst2_sc), INTENT(IN) :: fc2_sc
    complex(dp), allocatable, intent(out) :: Vqqs(:,:,:,:,:)
    !
    real(dp):: trans(3)
    integer, dimension(3,3) :: SM, SMT, SM1, SMT1
    integer, allocatable :: G_atom(:,:)
    type(forceconst2_sc) :: fc_temp
    INTEGER :: isym_map, isym, nai, naj, naii, najj, iR1, iR2
    integer :: site, N, iR, jR, iiR, jjR, iq, jq, inv_map, Nq, count
    integer :: Ri(3), Rj(3)
    real(dp) :: tau_crys(3), tau_rot_crys(3), G_crys(3)
    real(dp), allocatable :: tens4(:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: work(:,:,:,:,:,:)
    !
    !
    Nq = size(xq,2)
    N = product(fc2_sc%nq)
    call check_allocation_fits("DCA centered cluster V(Q,Q')", &
      complex_mem_bytes(int(S%nat3, int64) * int(S%nat3, int64) * int(Nq, int64) * &
      int(Nq, int64) * int(size(fc2_sc%defects,2), int64)) + &
      real_mem_bytes(2_int64 * 9_int64 * int(S%nat, int64) * int(S%nat, int64) * &
      int(N, int64) * int(N, int64)))
    allocate(Vqqs(S%nat3,S%nat3,Nq,Nq,size(fc2_sc%defects,2)))
    allocate(tens4(3,S%nat,3,S%nat, N, N))
    allocate(G_atom(3,S%nat))
    ALLOCATE( work, source=tens4 )
    !
    tens4 = reshape(fc2_sc%fc, [3,S%nat,3,S%nat,N,N])
    !
    call transform_tns4(tens4, -1)
    !
    do site = 1, size(fc2_sc%defects,2)
      call fc_temp%allocate(S, S_sc, fc2_sc%nq)
      count = 0

      DO isym = 1, nsym
        IF ( irt(isym, fc2_sc%defects(1,1)) == fc2_sc%defects(1,site) ) THEN
          isym_map = isym
          inv_map = invs(isym_map)
          SM1 = symm_mat(:,:,inv_map)
          SMT1 = transpose(SM1)
          SM = symm_mat(:,:,isym_map)
          SMT = transpose(SM)
          trans = cryst2cart(ft(:,isym_map), S%at, 1)
          count = count + 1
        END IF
      END DO
      !
      ! In crystal coordinates, QE's direct-space operation is
      !   tau' = transpose(SM) tau - ft.
      ! When atom nai is mapped to naii, the basis may also cross a cell
      ! boundary:
      !   transpose(SM) tau_nai - ft = tau_naii + G_atom(:,nai).
      ! That integer vector must be included when rotating the two lattice
      ! vector indices of V.
      do nai = 1, S%nat
        naii = irt(isym_map, nai)
        tau_crys = cryst2cart(S%tau(:,nai), S%bg, -1)
        tau_rot_crys = cryst2cart(S%tau(:,naii), S%bg, -1)
        G_crys = matmul(SMT, tau_crys) - ft(:,isym_map) - tau_rot_crys
        G_atom(:,nai) = NINT(G_crys)
        if(any(abs(G_crys - real(G_atom(:,nai), dp)) > 1e-6_dp)) then
          print*, "G_crys:", G_crys, " rounded:", G_atom(:,nai)
          call errore("center_V", "non-integer atom-dependent symmetry shift", 1)
        endif
      enddo
      !
      work = 0.0_DP
      DO nai = 1, S%nat
        naii = irt(isym_map, nai)
        DO naj = 1, S%nat
          najj = irt(isym_map, naj)
          do iR = 1, N
            do jR = 1, N
              Ri = matmul(SMT, index2v(iR, fc2_sc%nq)) + G_atom(:,nai)
              Rj = matmul(SMT, index2v(jR, fc2_sc%nq)) + G_atom(:,naj)
              iiR = v2index(bz2simple(Ri, fc2_sc%nq), fc2_sc%nq)
              jjR = v2index(bz2simple(Rj, fc2_sc%nq), fc2_sc%nq)
              work(:,naii,:,najj,iiR,jjR) = work(:,naii,:,najj,iiR,jjR) + matmul( &
                matmul( SM1, tens4(:,nai,:,naj,iR,jR) ), SMT1 )
            END DO
          END DO
        END DO
      END DO
      !
      call transform_tns4(work, 1)
      !
      fc_temp%fc = reshape( work, [3*S%nat, 3*S%nat, N, N] )
      fc_temp%taudef = fc2_sc%taudef + S%tau(:,fc2_sc%defects(1,site)) - S%tau(:, fc2_sc%defects(1,1))
      call fc_temp%center(fc2_sc%nq, S)
      !
      ! call translate_xR(fc_temp, trans)
      do iq = 1, Nq
        call fc_temp%r2q(xq(:,iq))
        do jq = 1, Nq
          call fc_temp%r2q(xq(:,jq), Vqqs(:,:,iq,jq,site))
        enddo
      enddo
      ! call translate_xR(fc_temp, -trans)
      !
      call fc_temp%deallocate()
    enddo
    !
    DEALLOCATE( work )
    DEALLOCATE( G_atom )
    !
  contains
    subroutine translate_xR(fc, t)
      real(dp), intent(in) :: t(3)
      type(forceconst2_sc), intent(inout) :: fc
      !
      do iR2 = 1, N
        fc%xR2(:,iR2) = fc%xR2(:,iR2) + t
        do iR1 = 1, N
          fc%xR1(:,iR1,iR2) = fc%xR1(:,iR1,iR2) + t
        enddo
      enddo
    end subroutine
    !
    subroutine transform_tns4(tns4, flag)
      real(dp), intent(inout) :: tns4(:,:,:,:,:,:)
      integer, intent(in) :: flag
      !
      do nai=1,S%nat; do naj=1,S%nat; do iR=1,N; do jR=1,N
              CALL transform_mat( tns4(:,nai,:,naj,iR,jR), flag )
            enddo; enddo; enddo; enddo
    end subroutine
    !
    subroutine transform_mat(mat, flag)
      real(dp), intent(inout) :: mat(:,:)
      integer, intent(in) :: flag
      real(dp) :: U(3,3)
      !
      if (flag == 1) then
        U = S%bg
      else
        U = transpose(S%at)
      end if
      mat = matmul( U, matmul( mat, transpose(U) ) )
    end subroutine
  END SUBROUTINE
  !
end module
