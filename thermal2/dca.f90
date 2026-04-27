module dca
  use kinds, only: dp
  use symm_q_mat, only : apply_sym
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg, equiv_grid
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils, only: bz2simple, grid_vec_cart, &
    e_iqr, index2v, cryst2cart, v2index, freq_in_grid, diag_cmplx, diag, &
    id_mat
  use defutils, only : flatten_RR_cmplx, unflatten_RR_cmplx, &
    quter_cmplx, fftinterp_mat2_cmplx
  use defect, only : tetra_from_self_cart, write_spf_ndiag, write_spf, &
    write_self, write_dos, tetra_from_self_diag, find_where, tetra_from_self
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize, setup_simple_grid
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi, pi
  use quter_defect, only : forceconst2_sc, inside_ws
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
    complex(dp), allocatable, dimension(:,:) :: den_weights, &
      Gf_conf, Gf_avg, overlap, &
      out_den_weights, self_diag, &
      w2_plus_self, Gf0i, dv, df, self_full_mat, self_full_mat_out
    complex(dp), allocatable, dimension(:,:,:) :: den_UL, den_UR, &
      self_fine, U, UT, UT_fine, U_fine, UT_out, U_out, &
      self_out_diag, V_conf, self_in, self_out, G0i_cluster, Gi_coarse, &
      self_R, self_coarse, rest, self_RL, self_den, D, self_full, weights
    complex(dp), allocatable, dimension(:,:,:,:) :: G_avg, Gi_conf, &
      Gi_avg, phase_mat, V, self_out_grid
    complex(dp), allocatable, dimension(:,:,:,:,:) :: Vqqs, &
      self_uf
    real(dp), allocatable, dimension(:,:) :: self_xR, xq, R, f, out_freqs, &
      window_weights, xq_patch
    integer, allocatable :: pos(:), kq(:,:)!, big_iq(:)
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, n_eq_sites, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ, kkq, &
      cluster_mesh(3), Npos, ntot, iqp, N_ITER_TOT, iR, RL_iter, Q_mesh(3), idx(S%nat3), locations(3)
    real(dp) :: conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX, max_diff, max_diff_coarse
    real(dp) :: dos(input%n_omega), shift(3), simulated_conc, xqq(3)
    logical :: conv, low_concentration
    complex(dp) :: A(S%nat3,S%nat3), eta, ialpha(S%nat3,S%nat3), eigc(S%nat3)
    real(dp) :: eig(S%nat3), bg_patch(3,3)
    integer, allocatable :: N_sites(:,:), c_equiv(:), iq_of(:), window_ind(:,:), &
      ind(:), equiv_full(:), window_count(:), equiv(:)
    character(len=100) :: filename
    type(q_grid) :: c_grid, in_grid_full
    complex(dp), allocatable :: proj(:,:), proj_fine(:,:)
    real(dp), allocatable :: eig_proj(:)
    complex(dp), allocatable :: self_uncoarsed(:,:)
    !
    if(input%calculation == "test") call init_random_seed()
    !
    MAXITER = 100
    ABS_TOLERANCE = 1e-12_dp
    REL_TOLERANCE = 1e-5_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 4
    conc = input%conc ! example concentration
    cluster_mesh = input%sc_grid
    Nc = product(cluster_mesh) !* size(fc2_sc%defects,2)
    n_eq_sites = size(fc2_sc%defects,2)
    NSAMPLES = 100
    NQ = product(input%nk_in) / Nc
    Q_mesh = input%nk_in / cluster_mesh
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true.)
    shift = xq(:,v2index([1,1,1],cluster_mesh))
    shift = cryst2cart(shift, S%at, -1)
    shift = shift * (Q_mesh-1) / 2 / Q_mesh
    shift = cryst2cart(shift, S%bg, 1)
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
    ! do i = 1, 3
    !   bg_patch(:,i) = S%bg(:,i) / cluster_mesh(i)
    ! enddo
    ! xq_patch = grid_vec_cart(in_grid%n*3, S%bg*3, divide=.true., center=.true.)
    !
    if (any(mod(in_grid%n, cluster_mesh) /= 0)) &
      call errore("dca_selfnrg", "grid size not multiple of fc2 grid size", 1)
    R = grid_vec_cart(cluster_mesh, S%at)
    if(input%calculation == 'dos' .or. input%calculation == 'test') &
      call set_wg(S, fc2, out_grid, input%n_omega, wg_out)
    call set_wg(S, fc2, in_grid, input%n_omega, wg)
    !
    allocate(UT_fine(S%nat3,S%nat3,in_grid%nqtot), U_fine(S%nat3,S%nat3,in_grid%nqtot))
    allocate(f(S%nat3,in_grid%nqtot), U(S%nat3,S%nat3,Nc), UT(S%nat3,S%nat3,Nc))
    allocate(D(S%nat3,S%nat3,in_grid%nqtot))
    call freq_in_grid(S, fc2, in_grid, f, U_fine)
    do iq = 1, in_grid%nq
      iqp = iq + in_grid%iq0
      call fftinterp_mat2(in_grid%xq(:,iq), S, fc2, D(:,:,iqp))
    enddo
    call mpi_bsum(S%nat3, S%nat3, in_grid%nqtot, D)
    do iq = 1, in_grid%nqtot
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
    allocate(self_out_diag(S%nat3,out_grid%nqtot,input%n_omega))
    allocate(UT_out(S%nat3,S%nat3,out_grid%nqtot), U_out(S%nat3,S%nat3,out_grid%nqtot))
    allocate(out_freqs(S%nat3,out_grid%nqtot), self_out_grid(S%nat3,S%nat3,out_grid%nqtot,input%n_omega))
    call freq_in_grid(S, fc2, out_grid, out_freqs, U_out)
    do iq = 1, out_grid%nqtot
      UT_out(:,:,iq) = conjg(transpose(U_out(:,:,iq)))
    enddo
    self_out_grid = 0.0_dp
    self_out_diag = 0.0_dp
    !
    allocate(proj(Nc, Nc))
    allocate(Gf0i(S%nat3*Nc, S%nat3*Nc))
    allocate(V_conf(S%nat3*Nc, S%nat3*Nc, NSAMPLES))
    allocate(self_diag(S%nat3, in_grid%nqtot))
    allocate(G_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(Gi_conf(S%nat3, S%nat3, Nc, Nc))
    allocate(Gf_conf(S%nat3*Nc, S%nat3*Nc))
    allocate(Gf_avg(S%nat3*Nc, S%nat3*Nc))
    allocate(Gi_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(Gi_coarse(S%nat3, S%nat3, Nc))
    allocate(G0i_cluster(S%nat3, S%nat3, Nc))
    allocate(self_in(S%nat3, S%nat3, Nc))
    allocate(self_uf(3, 3, S%nat, S%nat, Nc))
    allocate(self_fine(S%nat3, S%nat3, in_grid%nqtot))
    allocate(self_out(S%nat3, S%nat3, Nc))
    allocate(V(S%nat3, S%nat3, Nc, Nc))
    allocate(phase_mat(Nc, Nc, n_eq_sites, NSAMPLES))
    allocate(kq(NQ, Nc))!, big_iq(grid_scat%nqtot))
    allocate(den_weights(S%nat3, in_grid%nqtot))
    allocate(out_den_weights(S%nat3, out_grid%nqtot))
    allocate(den_UL(S%nat3, S%nat3, in_grid%nqtot))
    allocate(den_UR(S%nat3, S%nat3, in_grid%nqtot))
    allocate(overlap(S%nat3, in_grid%nqtot))
    ! allocate(den_eig(S%nat3, in_grid_full%nqtot))
    allocate(pos(Nc))
    allocate(w2_plus_self(S%nat3, in_grid%nqtot))
    allocate(iq_of(in_grid%nqtot))
    allocate(self_coarse(S%nat3, S%nat3, Nc))
    allocate(rest(S%nat3, S%nat3, in_grid_full%nqtot))
    allocate(self_RL(S%nat3, S%nat3, in_grid_full%nqtot))
    allocate(self_den(S%nat3, S%nat3, in_grid_full%nqtot))
    allocate(self_full_mat(S%nat3**2, in_grid_full%nqtot))
    allocate(self_full_mat_out(S%nat3**2, in_grid_full%nqtot))
    !
    allocate(df(S%nat3**2*Nc, MEMORY))
    allocate(dv(S%nat3**2*Nc, MEMORY))
    allocate(delta_in(S%nat3**2*Nc))
    allocate(delta_out(S%nat3**2*Nc))
    allocate(self_diff(S%nat3**2*Nc))
    allocate(self_full(S%nat3, S%nat3, in_grid_full%nqtot))
    allocate(weights(S%nat3, S%nat3, in_grid_full%nqtot))
    allocate(self_uncoarsed(S%nat3, input%n_omega))
    df = 0._dp
    dv = 0._dp
    delta_in = 0._dp
    delta_out = 0._dp
    self_uncoarsed = 0._dp
    !
    call center_V(xq, S, S_sc, fc2_sc, Vqqs)
    ! do idef = 1, n_eq_sites
    !   do iq = 1, Nc
    !     do jq = 1, Nc
    !       Vqqs(:,:,iq,jq,idef) = matmul(UT(:,:,iq), matmul(Vqqs(:,:,iq,jq,idef),U(:,:,jq))) / Nc
    !     enddo
    !   enddo
    ! enddo
    ! !
    ! allocate(window_count(in_grid_full%nqtot))
    ! window_count = 0
    ! allocate(window_weights(NQ*3, in_grid_full%nqtot))
    ! allocate(window_ind(NQ*3, in_grid_full%nqtot))
    ! call wsinit(rws, nrwsx, nrws, bg_patch)
    ! do iq = 1, in_grid_full%nqtot
    !   call inside_ws_dca(NQ, xq_patch, -in_grid_full%xq(:,iq), nrws, rws, wgt, ind)
    !   do i = 1, size(ind)
    !     xqq = cryst2cart(xq_patch(:,ind(i))+in_grid_full%xq(:,iq), S%at, -1) * in_grid%n
    !     kkq = v2index(bz2simple(NINT(xqq), in_grid%n), in_grid%n)
    !     window_count(iq) = window_count(iq) + 1
    !     window_weights(window_count(iq), iq) = wgt(i)
    !     window_ind(window_count(iq), iq) = kkq
    !   enddo
    ! enddo
    ! !> construction of of the phase, used to translate the potential to
    ! !> random configurations inside the cluster
    ntot = 0
    phase_mat = 0.0_dp
    !
    allocate(N_sites(n_eq_sites, NSAMPLES))
    do idef = 1, n_eq_sites
      call assign_defect_counts(NSAMPLES, Nc, conc, N_sites(idef,:))
    enddo
    call mpi_broadcast(n_eq_sites, NSAMPLES, N_sites)
    !
    do it = 1+my_id, NSAMPLES, num_procs
      do idef = 1, n_eq_sites
        ! call sample_binomial(Nc, conc, pos, Npos)
        Npos = N_sites(idef, it)
        call sample_canonical(Nc, Npos, pos)
        ntot = ntot + Npos
        do ipos = 1, Npos
          do k = 1, Nc
            do j = 1, Nc
              phase_mat(j,k,idef,it) = phase_mat(j,k,idef,it) + &
                e_iqr(xq(:,k)-xq(:,j), R(:,pos(ipos)))
              !       V(:,:,k,j,it) = V(:,:,k,j,it) + Vqqs(:,:,k,j,idef) * &
              !         phase_mat(k,j,pos(ipos)) / Nc
              !       ! e_iqr(xq(:,j)-xq(:,k), R(:,pos(ipos))) / Nc
            enddo
          enddo
        enddo
      enddo
    enddo
    call mpi_bsum(Nc, Nc, n_eq_sites, NSAMPLES, phase_mat)
    call mpi_bsum(ntot)
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
    !
    simulated_conc = real(ntot,dp) / (Nc * n_eq_sites * NSAMPLES)
    if(ionode) print"(A,E15.4)", "simulated concentration:", simulated_conc
    low_concentration = .false.
    if(all(N_sites < 2)) then
      low_concentration = .true.
      if(ionode) print*, "Using low concentration approximation (at most 1 defect per configuration)"
      simulated_conc = conc
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
    ! proj = 0._dp
    ! do jq = 1, Nc
    !   do iq = 1, Nc
    !     do i = 1, NQ
    !       do iR = 1, size(fc2%xR,2)
    !         proj(iq,jq) = proj(iq,jq) + e_iqr(xq(:,jq) - in_grid_full%xq(:,kq(i,iq)), fc2%xR(:,iR))
    !       enddo
    !     enddo
    !   enddo
    ! enddo
    ! !
    ! print*, maxval(abs(aimag(proj)))
    ! allocate(eig_proj(Nc))
    ! call mat2_diag(Nc, proj, eig_proj)
    ! do iq = 1, Nc
    !   print*, eig_proj(iq)
    ! enddo

    allocate(proj_fine(in_grid_full%nqtot, in_grid_full%nqtot))
    proj_fine = 0._dp
    do jq = 1, in_grid_full%nqtot
      do iq = 1, in_grid_full%nqtot
        xqq = cryst2cart(in_grid_full%xq(:,jq) - in_grid_full%xq(:,iq), S%at, -1) * &
          in_grid_full%n
        if(all(2*NINT(ABS(xqq)) + 1 <= Q_mesh)) proj_fine(iq,jq) = 1._dp / NQ
      enddo
    enddo
    !
    call invzmat(in_grid_full%nqtot, proj_fine)
    proj_fine = transpose(proj_fine)
    ! allocate(eig_proj(in_grid_full%nqtot))
    ! call mat2_diag(in_grid_full%nqtot, proj_fine, eig_proj)
    ! open(10, file='eig_proj.dat')
    ! do iq = 1, in_grid_full%nqtot
    !   write(10, *) eig_proj(iq)
    ! enddo
    ! close(10)
    ! stop 1
    ! do jq = 1, NQ
    !   print*, cryst2cart(in_grid_full%xq(:,kq(jq,1)), S%at, -1)
    ! enddo
    !>
    self_out = 1._dp
    self_in = 0._dp
    ! allocate(df(S%nat3**2*Nc, MEMORY))
    ! allocate(dv(S%nat3**2*Nc, MEMORY))
    ! allocate(delta_in(S%nat3**2*Nc))
    ! allocate(delta_out(S%nat3**2*Nc))
    ! delta_out = 0._dp
    ! delta_in = 0._dp
    ! df = 0._dp
    ! dv = 0._dp
    !
    dos = 0._dp
    N_ITER_TOT = 0
    self_fine = 0._dp
    self_coarse = 0._dp
    !
    if(ionode) print*, "Starting DCA self-energy calculation..."
    if(ionode) print*, ""
    do iw = 1, input%n_omega
      max_diff = huge(1.0_dp)
      do sc_iter = 1, MAXITER
        !
        if (.not. allocated(weights)) allocate(weights(S%nat3, S%nat3, in_grid%nqtot))
        call tetra_from_self_cart(S, in_grid, D + self_fine, wg%en(iw)**2, weights, mpi=.true.)
        call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, weights, equiv, in_grid_full%xq, .false.)
        !
        Gi_coarse = 0.0_dp
        do iq = 1, Nc
          do jq = 1, NQ
            kkq = kq(jq,iq)
            Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) + weights(:,:,kkq) * Nc
          enddo
          ! A = (Gi_coarse(:,:,iq) - conjg(transpose(Gi_coarse(:,:,iq)))) / cmplx(0._dp, 2.0_dp)
          ! call mat2_diag(S%nat3, A, eig)
          ! print"(6E15.4)", eig
          !
          call invzmat(S%nat3, Gi_coarse(:,:,iq))
          G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_coarse(:,:,iq)
        enddo
        !
        ! call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, G0i_cluster, c_equiv, xq, .true.)
        deallocate(weights)
        !
        Gf0i = 0.0_dp
        do iq = 1, Nc
          Gf0i((iq-1)*S%nat3+1:iq*S%nat3,(iq-1)*S%nat3+1:iq*S%nat3) = G0i_cluster(:,:,iq)
        enddo
        if(low_concentration) then
          Gi_conf = 0.0_dp
          do iq = 1, Nc
            A = G0i_cluster(:,:,iq)
            call invzmat(S%nat3, A)
            Gi_conf(:,:,iq,iq) = A
          enddo
          Gf_conf = flatten_RR_cmplx(Gi_conf)
          Gf_avg = Gf_conf * (1 - conc * n_eq_sites * Nc)
          do idef = 1, n_eq_sites
            Gi_conf = 0.0_dp
            do iq = 1, Nc
              Gi_conf(:,:,iq,iq) = G0i_cluster(:,:,iq)
              do jq = 1, Nc
                Gi_conf(:,:,iq,jq) = Gi_conf(:,:,iq,jq) - Vqqs(:,:,iq,jq,idef) / Nc
              enddo
            enddo
            Gf_conf = flatten_RR_cmplx(Gi_conf)
            call invzmat(S%nat3*Nc, Gf_conf)
            Gf_avg = Gf_avg + Gf_conf * (conc * Nc)
          enddo
        else
          Gf_avg = 0.0_dp
          do it = 1+my_id, NSAMPLES, num_procs
            ! > construction of G_conf for a given configuration
            ! Gi_conf = 0.0_dp
            !
            Gf_conf = Gf0i - V_conf(:,:,it)
            call invzmat(S%nat3*Nc, Gf_conf)
            !
            Gf_avg = Gf_avg + Gf_conf
            !
          enddo
          call mpi_bsum(S%nat3*Nc, S%nat3*Nc, Gf_avg)
          Gf_avg = Gf_avg / real(NSAMPLES, dp)
        endif
        !
        !> self energy is G0_cluster^-1 - <G>^-1
        Gi_avg = unflatten_RR_cmplx(Gf_avg, Nc, Nc)
        do iq = 1, Nc
          call invzmat(S%nat3, Gi_avg(:,:,iq,iq))
          self_out(:,:,iq) = G0i_cluster(:,:,iq) - Gi_avg(:,:,iq,iq)
        enddo
        !
        ! do iq = 1, Nc
        !   A = self_out(:,:,iq)
        !   call mat2_diag(S%nat3, A, self_in(:,1,iq))
        ! enddo
        ! !
        ! call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_out, c_equiv, xq)

        ! do iq = 1, Nc
        !   call mat2_diag(S%nat3, self_out(:,:,iq), self_in(:,2,iq))
        !   print*, sum(abs(self_in(:,2,iq))), sum(abs(self_in(:,1,iq))), c_equiv(iq)
        ! enddo
        ! print*, "-----------------------------------"
        ! call symmetrize_mat_cmplx(c_equiv, self_out)
        !>
        !
        ! do iq = 1, Nc
        !   call mat2_diag(S%nat3, self_out(:,:,iq), eigc)
        !   delta_out((iq-1)*S%nat3+1:iq*S%nat3) = eigc
        ! enddo
        delta_out = reshape(self_out, [S%nat3**2*Nc])
        call mix_broyden_full(S%nat3**2*Nc, delta_out, delta_in, &
          ALPHA_MIX, sc_iter, MEMORY, df, dv)
        self_in = reshape(delta_in, [S%nat3, S%nat3, Nc])
        !
        call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_in, c_equiv, xq, .true.)

        if(all(abs(real(self_in - self_coarse, dp)) < ABS_TOLERANCE + abs(real(self_coarse, dp)) * REL_TOLERANCE) .and. &
          all(abs(aimag(self_in - self_coarse)) < ABS_TOLERANCE + abs(aimag(self_coarse)) * REL_TOLERANCE)) exit
        ! do iq = 1, Nc
        !   self_in(:,:,iq) = matmul(self_out(:,:,iq), matmul( &
        !     diag_cmplx(delta_in((iq-1)*S%nat3+1:iq*S%nat3)), conjg(transpose(self_out(:,:,iq)))))
        ! enddo

        ! do iq = 1, Nc
        !   A = self_in(:,:,iq)
        !   call mat2_diag(S%nat3, A, eig)
        !   print*, sum(abs(eig)), c_equiv(iq)
        ! enddo
        ! do iq = 1, Nc
        !   call mat2_diag(S%nat3, self_out(:,:,iq), self_in(:,2,iq))
        !   print*, sum(abs(self_in(:,2,iq))), sum(abs(self_in(:,1,iq))), c_equiv(iq)
        ! enddo
        ! print*, "-----------------------------------"
        ! self_in = self_out
        ! max_diff = 0.0_dp
        ! do iq = 1, Nc
        !   A = (self_in(:,:,iq) - conjg(transpose(self_in(:,:,iq)))) / cmplx(0._dp, 2.0_dp)
        !   call mat2_diag(S%nat3, A, eig)
        !   max_diff = max(max_diff, maxval(eig))
        ! enddo
        ! if(ionode) print"(A,E15.4)", "maximum eig of self_in", max_diff
        !
        !
        ! ialpha = id_mat(S%nat3) * cmplx(0._dp, 1e-10_dp, dp)
        do iq = 1, Nc
          ! A = matmul(U(:,:,iq), matmul(self_in(:,:,iq), UT(:,:,iq)))
          ! self_in(:,:,iq) = A
          ! A = id_mat(S%nat3) * wg%en(iw)**2 - A + ialpha
          ! A = A - ialpha
          ! call invzmat(S%nat3, A)
          self_uf(:,:,:,:,iq) = unflatten_RR_cmplx(self_in(:,:,iq), S%nat, S%nat)
        enddo
        !
        CALL quter_cmplx(input%sc_grid(1), input%sc_grid(2), input%sc_grid(3), &
          S%nat, S%tau, S%at, S%bg, self_uf, xq, self_R, self_xR, 2)
        !
        ! if(ionode) print"(A,I4.4,A,I4.4,A,E15.4,A,E15.4)", "iw =", iw, ", sc_iter =", sc_iter, ", &
        !   sum(abs(self_in)) =", sum(abs(self_in)) /S%nat3**2 / Nc, &
        !   ", sum(abs(self_out)) =", sum(abs(self_out)) /S%nat3**2 / Nc
        ! do iR = 1, size(self_xR,2)
        !   self_R(:,:,iR) = self_R(:,:,iR) / sinc_3d(shift, self_xR(:,iR))
        ! enddo
        !
        self_fine = 0.0_dp
        do iq = 1, in_grid%nq
          iqp = iq + in_grid%iq0
          call fftinterp_mat2_cmplx(in_grid%xq(:,iq), S, self_R, self_xR, self_fine(:,:,iqp))
          ! call invzmat(S%nat3, A)
          ! ! self_fine(:,:,iqp) = ialpha + id_mat(S%nat3) * wg%en(iw)**2 - A
          ! self_fine(:,:,iqp) = A + ialpha
          ! self_fine(:,:,iqp) = matmul(UT_fine(:,:,iqp), matmul(A, U_fine(:,:,iqp)))
        enddo
        call mpi_bsum(S%nat3, S%nat3, in_grid%nqtot, self_fine)
        !
        ! write(filename, "(A,I4.4,A)") "self_R_", sc_iter, ".dat"
        ! open(10, file=filename)
        ! do iR = 1, size(self_xR,2)
        !   do i = 1, S%nat
        !     do j = 1, S%nat
        !       write(10, "(3E20.8)") norm2(self_xR(:,iR) + S%tau(:,i)- S%tau(:,j)), sum(abs(self_R((i-1)*3+1:i*3,(j-1)*3+1:j*3,iR)))
        !     enddo
        !   enddo
        ! enddo
        ! close(10)
        !
        ! max_diff = 0.0_dp
        ! do iq = 1, in_grid%nqtot
        !   A = (self_fine(:,:,iq) - conjg(transpose(self_fine(:,:,iq)))) / cmplx(0._dp, 2.0_dp)
        !   call mat2_diag(S%nat3, A, eig)
        !   max_diff = max(max_diff, maxval(eig))
        ! enddo
        ! if(ionode) print*, "maximum eig of self_fine", max_diff
        !
        ! RL_iter = 0
        ! do iq = 1, in_grid_full%nqtot
        !   self_RL(:,:,iq) = self_fine(:,:,wg%e(iq))
        ! enddo
        !
        ! !> RL algorithm
        ! do
        !   self_den = 0._dp
        !   do iq = 1, in_grid_full%nqtot
        !     do jq = 1, window_count(iq)
        !       self_den(:,:,iq) = self_den(:,:,iq) + self_RL(:,:,window_ind(jq,iq)) * window_weights(jq,iq) / NQ
        !     enddo
        !     ! Pack Re-ratio → real(rest), Im-ratio → aimag(rest)
        !     rest(:,:,iq) = cmplx( &
        !       merge(real(self_fine(:,:,wg%e(iq)),dp)  / real(self_den(:,:,iq),dp),  1._dp, &
        !       abs(real(self_den(:,:,iq),dp))  >= 1e-3_dp*maxval(real(self_den(:,:,iq),dp))), &
        !       merge(aimag(self_fine(:,:,wg%e(iq))) / aimag(self_den(:,:,iq)), 1._dp, &
        !       abs(aimag(self_den(:,:,iq))) >= 1e-3_dp*maxval(aimag(self_den(:,:,iq)))), dp)
        !   enddo
        !   !
        !   self_den = 0._dp
        !   do iq = 1, in_grid_full%nqtot
        !     do jq = 1, window_count(iq)
        !       self_den(:,:,iq) = self_den(:,:,iq) + rest(:,:,window_ind(jq,iq)) * window_weights(jq,iq) / NQ
        !     enddo
        !   enddo
        !   print*, sum(self_den) / S%nat3**2 / in_grid_full%nqtot
        !   !
        !   do iq = 1, in_grid_full%nqtot
        !     self_RL(:,:,iq) = cmplx( &
        !       real(self_RL(:,:,iq),dp)  * real(self_den(:,:,iq),dp), &
        !       aimag(self_RL(:,:,iq)) * aimag(self_den(:,:,iq)), dp)
        !   enddo
        !   print*, "-", sum(self_RL)
        !   if(all(abs(real(self_den,dp)  - 1._dp) < 1e-5_dp) .and. &
        !     all(abs(aimag(self_den) - 1._dp) < 1e-5_dp)) exit
        !   RL_iter = RL_iter + 1
        !   if(RL_iter > 100) then
        !     if(ionode) print*, "Warning: inner loop not converging", &
        !       maxval(abs(real(self_den,dp)-1._dp)), maxval(abs(aimag(self_den)-1._dp))
        !     exit
        !   endif
        ! enddo !> RL loop
        !
        ! jq = v2index((Q_mesh - 1) / 2, Q_mesh)
        ! iq = 10
        if(allocated(self_full)) deallocate(self_full)
        allocate(self_full, source=self_fine)
        self_full = self_fine
        call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_full, equiv, in_grid_full%xq, .false.)
        ! print*, in_grid_full%xq(:,kq(jq,iq))
        ! print*, sum(self_full(:,:,kq(jq,iq)))
        ! print*, sum(self_in(:,:,iq))
        self_full_mat = reshape(self_full, [S%nat3**2, in_grid_full%nqtot])
        self_full_mat_out = 0._dp
        call zgemm('N','N', S%nat3**2, in_grid_full%nqtot, in_grid_full%nqtot, (1._dp, 0._dp), &
          self_full_mat, S%nat3**2, proj_fine, in_grid_full%nqtot, (1._dp, 0._dp), self_full_mat_out, S%nat3**2)
        self_full = reshape(self_full_mat_out, [S%nat3, S%nat3, in_grid_full%nqtot])
        ! call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_full, equiv, in_grid_full%xq, .true.)
        do iq = 1, in_grid%nqtot
          self_fine(:,:,iq) = self_full(:,:,wg%e(iq))
        enddo
        !
        self_coarse = 0.0_dp
        do iq = 1, Nc
          do jq = 1, NQ
            self_coarse(:,:,iq) = self_coarse(:,:,iq) + self_full(:,:,kq(jq,iq)) / NQ
          enddo
          ! self_coarse(:,:,iq) = matmul(UT(:,:,iq), matmul(self_coarse(:,:,iq), U(:,:,iq)))
        enddo
        !
        ! print*, sum(self_coarse(:,:,10))
        ! call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_coarse, c_equiv, xq, .true.)
        ! do iq = 1, Nc
        !   A = (self_out(:,:,iq) - conjg(transpose(self_out(:,:,iq)))) / cmplx(0._dp, 2.0_dp)
        !   call mat2_diag(S%nat3, A, eig)
        !   if(ionode) print*, eig
        ! enddo
        ! do iq = 1, in_grid%nqtot
        !   self_fine(:,:,iq) = matmul(UT_fine(:,:,iq), matmul(self_fine(:,:,iq), U_fine(:,:,iq)))
        ! enddo
        ! A = self_in(:,:,10)
        ! call mat2_diag(S%nat3, A, eigc)
        ! eig = real(eigc, dp)
        ! idx = 0
        ! CALL hpsort(S%nat3, eig, idx)
        ! eigc = eigc(idx)
        ! if(ionode) print"(A,12E20.8)", "in    ",  eigc
        ! A = self_coarse(:,:,10)
        ! call mat2_diag(S%nat3, A, eigc)
        ! eig = real(eigc, dp)
        ! idx = 0
        ! CALL hpsort(S%nat3, eig, idx)
        ! eigc = eigc(idx)
        ! if(ionode) print"(A,12E20.8)", "coarse",  eigc
        !
      enddo ! self-energy SC cycle
      N_ITER_TOT = N_ITER_TOT + sc_iter - 1
      if(sc_iter - 1 /= MAXITER) then
        if(ionode) print"(A,I4,A,I4,A)", "frequency ", iw, " converged in ", sc_iter-1, " iterations."
      else
        if(ionode) print"(A,I4,A,I4,A,E15.3)", "frequency ", iw, &
          " NOT converged in ", sc_iter-1, " iterations. Max diff: ", max_diff
      end if
      !
      jq = v2index((Q_mesh - 1) / 2, Q_mesh)
      kkq = kq(jq,1)
      call freq_phq_safe(in_grid_full%xq(:,kkq), S, fc2, f(:,1), U(:,:,1))
      A = self_full(:,:,kkq)
      A = matmul(conjg(transpose(U(:,:,1))), matmul(A, U(:,:,1)))
      do i = 1, S%nat3
        self_uncoarsed(i,iw) = A(i,i)
      enddo
      !
      do iq = 1, out_grid%nqtot
        call fftinterp_mat2_cmplx(out_grid%xq(:,iq), S, self_R, self_xR, A)
        ! call invzmat(S%nat3, A)
        ! A = A + ialpha
        self_out_grid(:,:,iq,iw) = matmul(UT_out(:,:,iq), matmul(A, U_out(:,:,iq))) * &
          conc / simulated_conc
        do i = 1, S%nat3
          self_out_diag(i,iq,iw) = self_out_grid(i,i,iq,iw)
        enddo
      enddo
    enddo ! frequency loop

    !
    open(115, file='self_uncoarsed.dat')
    do iw = 1, input%n_omega
      if(ionode) write(115, "(I4, 12E20.8)") wg%en(iw), self_uncoarsed(:,iw)
    enddo
    close(115)
    !

    if(ionode) print*, "Average number of iterations per frequency:", real(N_ITER_TOT,dp) / input%n_omega
    select case(input%calculation)
     case ('spf-def')
      call write_spf_ndiag('spf-dca-ndiag.dat', wg%en, self_out_grid, out_freqs)
      call write_spf('spf-dca.dat', wg%en, self_out_diag, out_freqs)
      call write_self('self-dca.dat', wg%en, self_out_diag)
     case('self')
      call write_self('self-dca.dat', wg%en, self_out_diag)
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
    type(forceconst2_sc) :: fc_temp
    INTEGER :: isym_map, isym, nai, naj, naii, najj, iR1, iR2
    integer :: site, N, iR, jR, iiR, jjR, iq, jq, inv_map, Nq
    real(dp), allocatable :: tens4(:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: work(:,:,:,:,:,:)
    !
    !
    Nq = size(xq,2)
    N = product(fc2_sc%nq)
    allocate(Vqqs(S%nat3,S%nat3,Nq,Nq,size(fc2_sc%defects,2)))
    allocate(tens4(3,S%nat,3,S%nat, N, N))
    ALLOCATE( work, source=tens4 )
    !
    tens4 = reshape(fc2_sc%fc, [3,S%nat,3,S%nat,N,N])
    !
    call transform_tns4(tens4, -1)
    !
    do site = 1, size(fc2_sc%defects,2)
      call fc_temp%allocate(S, S_sc, fc2_sc%nq)
      DO isym = 1, nsym
        IF ( irt(isym, fc2_sc%defects(1,1)) == fc2_sc%defects(1,site) ) THEN
          isym_map = isym
          inv_map = invs(isym_map)
          SM1 = symm_mat(:,:,inv_map)
          SMT1 = transpose(SM1)
          SM = symm_mat(:,:,isym_map)
          SMT = transpose(SM)
          trans = cryst2cart(ft(:,isym_map), S%at, 1)
          EXIT
        END IF
      END DO
      !
      work = 0.0_DP
      DO nai = 1, S%nat
        naii = irt(isym_map, nai)
        DO naj = 1, S%nat
          najj = irt(isym_map, naj)
          do iR = 1, N
            do jR = 1, N
              iiR = v2index(bz2simple(matmul(SM, index2v(iR, fc2_sc%nq)), fc2_sc%nq), fc2_sc%nq)
              jjR = v2index(bz2simple(matmul(SM, index2v(jR, fc2_sc%nq)), fc2_sc%nq), fc2_sc%nq)
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
      fc_temp%taudef = S%tau(:,fc2_sc%defects(1,site)) - S%tau(:, fc2_sc%defects(1,1))
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
