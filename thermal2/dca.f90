module dca
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg, equiv_grid
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils, only: bz2simple, grid_vec_cart, &
    e_iqr, index2v, cryst2cart, v2index, freq_in_grid, diag_cmplx, diag
  use defutils, only : flatten_RR_cmplx, unflatten_RR_cmplx, &
    quter_cmplx, fftinterp_mat2_cmplx
  use defect, only : tetra_from_self, write_spf_ndiag, write_spf, &
    write_self, write_dos, tetra_from_self_diag
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize, setup_simple_grid
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi, pi
  use quter_defect, only : forceconst2_sc
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
  subroutine dca_selfnrg(S, S_sc, input, fc2, fc2_sc, in_grid, out_grid)
    type(ph_system_info), intent(in) :: S, S_sc
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    !! centered
    type(forceconst2_sc), intent(inout) :: fc2_sc
    !! not centered
    type(q_grid), intent(in) :: in_grid
    !! symmetric and scattered, normal order
    type(q_grid), intent(in) :: out_grid
    !! symmetric and not scattered, normal order
    ! complex(dp), allocatable, intent(in) :: self_fb(:,:,:,:)
    ! real(dp), allocatable, intent(in) :: xR_fb(:,:)
    type(tetra_output) :: wg, wg_out
    !
    type(forceconst2_grid) :: self_R
    complex(dp), allocatable, dimension(:) :: delta_in, delta_out, self_diff
    complex(dp), allocatable, dimension(:,:) :: den_weights, &
      den_eig, Gf_conf, Gf_avg, overlap, &
      out_den_weights, self_diag, &
      w2_plus_self, VL, VR, Gf0i, dv, df
    complex(dp), allocatable, dimension(:,:,:) :: den_UL, den_UR, &
      self_fine, U, UT, UT_fine, U_fine, UT_out, U_out, &
      self_out_diag, V_conf, self_in, self_out, G0i_cluster, Gi_coarse
    complex(dp), allocatable, dimension(:,:,:,:) :: G_avg, Gi_conf, &
      Gi_avg, phase_mat, V, self_out_grid
    complex(dp), allocatable, dimension(:,:,:,:,:) :: Vqqs, &
      self_uf
    real(dp), allocatable, dimension(:,:) :: self_xR, xq, R, f, out_freqs
    integer, allocatable :: pos(:), kq(:,:)!, big_iq(:)
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, n_eq_sites, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ, kkq, &
      cluster_mesh(3), Npos, ntot, iqp, N_ITER_TOT, iR
    real(dp) :: conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX, max_diff
    real(dp) :: dos(input%n_omega), shift(3), simulated_conc
    logical :: conv, low_concentration
    complex(dp) :: A(S%nat3,S%nat3), eta
    real(dp) :: eig(S%nat3)
    integer, allocatable :: N_sites(:,:), c_equiv(:)
    character(len=100) :: filename
    type(q_grid) :: c_grid
    real(dp), allocatable :: w2(:)
    !
    if(input%calculation == "test") call init_random_seed()
    !
    MAXITER = 1000
    ABS_TOLERANCE = 1e-10_dp * input%conc
    REL_TOLERANCE = 1e-6_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 4
    conc = input%conc ! example concentration
    cluster_mesh = [4,4,1]
    Nc = product(cluster_mesh) !* size(fc2_sc%defects,2)
    n_eq_sites = size(fc2_sc%defects,2)
    NSAMPLES = 100
    NQ = product(in_grid%n) / Nc
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true.)
    shift = xq(:,v2index([1,1,1],cluster_mesh)) / 2
    if (NQ == 1) shift = 0.0_dp
    do iq = 1, size(xq,2)
      xq(:,iq) = xq(:,iq) + shift
    enddo
    !
    call setup_simple_grid(S%bg, cluster_mesh(1), cluster_mesh(2), cluster_mesh(3), c_grid, shift)
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
    allocate(UT_fine(S%nat3,S%nat3,in_grid%nqtot), U_fine(S%nat3,S%nat3,in_grid%nqtot))
    allocate(f(S%nat3,in_grid%nqtot), U(S%nat3,S%nat3,Nc), UT(S%nat3,S%nat3,Nc))
    call freq_in_grid(S, fc2, in_grid, f, U_fine)
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
    allocate(den_eig(S%nat3, in_grid%nqtot))
    allocate(pos(Nc))
    allocate(w2_plus_self(S%nat3, in_grid%nqtot))

    allocate(VL, source=Gf_avg)
    allocate(VR, source=Gf_avg)
    allocate(w2(S%nat3*Nc))
    !
    allocate(df(S%nat3**2*Nc, MEMORY))
    allocate(dv(S%nat3**2*Nc, MEMORY))
    allocate(delta_in(S%nat3**2*Nc))
    allocate(delta_out(S%nat3**2*Nc))
    allocate(self_diff(S%nat3**2*Nc))
    df = 0._dp
    dv = 0._dp
    delta_in = 0._dp
    delta_out = 0._dp
    !
    call center_V(xq, S, S_sc, fc2_sc, Vqqs)
    do idef = 1, n_eq_sites
      do iq = 1, Nc
        do jq = 1, Nc
          Vqqs(:,:,iq,jq,idef) = matmul(UT(:,:,iq), matmul(Vqqs(:,:,iq,jq,idef),U(:,:,jq))) / Nc
        enddo
      enddo
    enddo
    !
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
              Vqqs(:,:,iq,jq,idef) * phase_mat(iq,jq,idef,it)
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
        kq(jq,iq) = v2index(index2v(jq, in_grid%n / cluster_mesh) + &
          index2v(iq, cluster_mesh) * in_grid%n / cluster_mesh, in_grid%n)
      enddo
    enddo
    !
    !>
    !
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
    !
    print*, in_grid%xq
    if(ionode) print*, "Starting DCA self-energy calculation..."
    if(ionode) print*, ""
    do iw = 10, input%n_omega
      do sc_iter = 1, MAXITER
        Gi_coarse = 0.0_dp
        do iq = 1, Nc
          do jq = 1, NQ
            kkq = wg%e(kq(jq,iq))
            ! print*, kkq, kq(jq,iq), NINT(cryst2cart(in_grid%xq(:,kkq), S%at, -1) * in_grid%n)
            A = diag(wg%f(:,kkq)**2 + wg%en(iw)**2) + self_in(:,:,iq)
            call mat2_diag(S%nat3, A, eig)
            Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) - matmul(A, matmul(diag(1._dp/eig), conjg(transpose(A))))
            ! Gi_coarse(:,iq) = Gi_coarse(:,iq) + den_weights(:,kkq)
          enddo
          A = Gi_coarse(:,:,iq) / NQ
          call mat2_diag(S%nat3, A, eig)
          Gi_coarse(:,:,iq) = matmul(A, matmul(diag(1._dp/eig), conjg(transpose(A))))
          ! Gi_coarse(:,:,iq) = (Gi_coarse(:,:,iq) + conjg(transpose(Gi_coarse(:,:,iq)))) / 2.0_dp
          ! Gi_coarse(:,iq) = 1._dp / Gi_coarse(:,iq) / Nc
          G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_in(:,:,iq)
        enddo
        Gf0i = 0.0_dp
        do iq = 1, Nc
          Gf0i((iq-1)*S%nat3+1:iq*S%nat3,(iq-1)*S%nat3+1:iq*S%nat3) = G0i_cluster(:,:,iq)
        enddo
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
        !
        !> self energy is G0_cluster^-1 - <G>^-1
        Gi_avg = unflatten_RR_cmplx(Gf_avg, Nc, Nc)
        do iq = 1, Nc
          A = Gi_avg(:,:,iq,iq)
          call mat2_diag(S%nat3, A, eig)
          Gi_avg(:,:,iq,iq) = matmul(A, matmul(diag(1._dp/eig), conjg(transpose(A))))
          self_out(:,:,iq) = G0i_cluster(:,:,iq) - Gi_avg(:,:,iq,iq)
        enddo
        call symmetrize_mat_cmplx(c_equiv, self_out)
        !>
        !
        do iq = 1, Nc
          self_out(:,:,iq) = (self_out(:,:,iq) + conjg(transpose(self_out(:,:,iq)))) / 2.0_dp
        enddo
        delta_out = reshape(self_out, [S%nat3**2*Nc])
        !
        call mix_broyden_full(S%nat3**2*Nc, delta_out, delta_in, &
          ALPHA_MIX, sc_iter, MEMORY, df, dv)
        self_in = reshape(delta_in, [S%nat3, S%nat3, Nc])
        !
        conv = .true.
        max_diff = 0._dp
        do iq = 1, S%nat3**2*Nc
            max_diff = max(max_diff, ABS(delta_out(iq)))
            if( ABS(delta_out(iq)) > REL_TOLERANCE * &
              abs(delta_in(iq)) + ABS_TOLERANCE) conv = .false.
        enddo
        if(ionode) print*, iw, sc_iter, max_diff, sum(abs(self_in)) /S%nat3**2 / Nc, sum(abs(self_out)) /S%nat3**2 / Nc
        if (conv) exit
      enddo ! self-energy SC cycle
      N_ITER_TOT = N_ITER_TOT + sc_iter - 1
      if(conv) then
        if(ionode) print"(A,I4,A,I4,A)", "frequency ", iw, " converged in ", sc_iter-1, " iterations."
      else
        if(ionode) print"(A,I4,A,I4,A,E15.3)", "frequency ", iw, &
          " NOT converged in ", sc_iter-1, " iterations. Max diff: ", max_diff
      end if
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
    integer :: i
    !
    do i = 1, Npos
      call random_number(r)
      pos(i) = 1 + floor(r * N)
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
