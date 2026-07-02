module dca
  use kinds, only: dp
  use symm_q_mat, only : apply_sym_q, apply_sym
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg, equiv_grid
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils, only: bz2simple, grid_vec_cart, &
    e_iqr, index2v, cryst2cart, v2index, freq_in_grid, diag_cmplx, diag, &
    id_mat, zgemm_N
  use defutils, only : flatten_RR_cmplx, unflatten_RR_cmplx, &
    quter_cmplx, fftinterp_mat2_cmplx, quter_R
  use defect, only : tetra_from_self_cart, write_spf_ndiag, write_spf, &
    write_self, write_dos, tetra_from_self_diag, find_where, tetra_from_self
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize, setup_simple_grid
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi, pi
  use quter_defect, only : forceconst2_sc, inside_ws, write_fc2_sc_RR
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
  use nist_isotopes_db, only : get_natural_isotopes
  use more_constants, only : MASS_DALTON_TO_RY
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
  subroutine flip_positive_imag_eigs(n, mat)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: mat(n,n)
    !
    complex(dp) :: eig_l(n,n), eig_r(n,n), eig(n), overlap(n)
    integer :: i, j, k
    logical :: changed
    !
    eig_l = mat
    call mat2_diag(n, eig_l, eig_r, eig)
    changed = .false.
    do i = 1, n
      if(aimag(eig(i)) > 0._dp) then
        eig(i) = cmplx(real(eig(i), dp), -aimag(eig(i)), kind=dp)
        changed = .true.
      endif
    enddo
    if(changed) then
      do i = 1, n
        overlap(i) = dot_product(eig_l(:,i), eig_r(:,i))
      enddo
      mat = 0._dp
      do k = 1, n
        do j = 1, n
          do i = 1, n
            mat(i,j) = mat(i,j) + eig(k) * conjg(eig_l(j,k)) * eig_r(i,k) / overlap(k)
          enddo
        enddo
      enddo
    endif
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
    complex(dp), allocatable, dimension(:) :: delta_in, delta_out
    complex(dp), allocatable, dimension(:,:) :: &
      Gf_conf, Gf_avg, out_den_weights, Gf0i, dv, df, Gf0
    complex(dp), allocatable, dimension(:,:,:) :: &
      self_fine, U, UT, UT_fine, U_fine, UT_out, U_out, &
      self_out_diag, self_before, self_next, G0i_cluster, Gi_coarse, &
      Gi, VM, VK, V, phases
    complex(dp), allocatable, dimension(:,:,:,:) :: &
      Gi_avg, phase_mat, mass_phase_mat, self_out_grid, c_out
    complex(dp), allocatable, dimension(:,:,:,:,:) :: Vqqs
    real(dp), allocatable, dimension(:,:) :: xq, R, f, out_freqs, site_mass_eps
    real(dp), allocatable, dimension(:) :: host_isotope_mass, host_isotope_conc, &
      host_isotope_cdf, impurity_isotope_mass, impurity_isotope_conc, impurity_isotope_cdf
    integer, allocatable, dimension(:) :: pos, c_equiv, iq_of, ind, equiv_full, equiv
    integer, allocatable, dimension(:,:) :: kq
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, n_eq_sites, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ, kkq, na1, na2, j1, j2, &
      cluster_mesh(3), iqp, N_ITER_TOT, iR, jR, Q_mesh(3), ntot
    real(dp) :: conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX, max_diff, max_diff_coarse
    real(dp) :: dos(input%n_omega), shift(3), simulated_conc
    logical :: conv, single_sample
    complex(dp) :: A(S%nat3,S%nat3), eigc(S%nat3)
    logical, allocatable :: impurity_site(:,:,:)
    integer :: host_atom, impurity_atom, host_type, na_def
    real(dp) :: host_reference_mass, host_natural_mass, impurity_natural_mass, output_scale
    character(len=2) :: host_element, impurity_element
    type(q_grid) :: c_grid, in_grid_full
    !
    integer :: defect_idef(S%nat)
    real(dp), allocatable :: img_xR(:,:,:,:), img_weight(:,:,:), xR(:,:)
    integer :: img_nR(S%nat, S%nat), n_tot_sites
    !
    if(input%calculation == "test") call init_random_seed()
    !
    MAXITER = 100
    ABS_TOLERANCE = 1e-12_dp
    REL_TOLERANCE = 1e-5_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 4
    conc = input%conc ! example concentration
    cluster_mesh = input%dca_grid
    Nc = product(cluster_mesh) !* size(fc2_sc%defects,2)
    n_eq_sites = size(fc2_sc%defects,2)
    NSAMPLES = 100
    NQ = product(input%nk_in) / Nc
    Q_mesh = input%nk_in / cluster_mesh
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true.)
    xR = grid_vec_cart(cluster_mesh, S%at)
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
    if(num_procs > 1) call in_grid_full%scatter()
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
      call set_wg(S, fc2, out_grid, input, wg_out)
    call set_wg(S, fc2, in_grid, input, wg, &
      skip_w0 = (fc2_sc%def_type == "inclusion"))
    !
    allocate(UT_fine(S%nat3,S%nat3,in_grid_full%nqtot), U_fine(S%nat3,S%nat3,in_grid_full%nqtot))
    allocate(f(S%nat3,in_grid_full%nqtot), U(S%nat3,S%nat3,Nc), UT(S%nat3,S%nat3,Nc))
    call freq_in_grid(S, fc2, in_grid_full, f, U_fine)
    do iq = 1, in_grid_full%nq
      iqp = iq + in_grid_full%iq0
      if(norm2(in_grid_full%xq(:,iq)) < 1e-10_dp) then
        call fftinterp_mat2(in_grid_full%xq(:,iq), S, fc2, U_fine(:,:,iqp))
        call mat2_diag(S%nat3, U_fine(:,:,iqp), f(:,iqp))
      endif
      UT_fine(:,:,iqp) = conjg(transpose(U_fine(:,:,iqp)))
    enddo
    deallocate(f)
    allocate(f(S%nat3, Nc))
    do iq = 1, Nc
      call freq_phq_safe(xq(:,iq), S, fc2, f(:,iq), U(:,:,iq))
      UT(:,:,iq) = conjg(transpose(U(:,:,iq)))
    enddo
    !
    if(input%isotope_scattering) then
      if(fc2_sc%def_type /= "substitution") &
        call errore("dca_selfnrg", "isotope_scattering currently requires a substitutional defect", 1)
      ! All entries in defects are symmetry-equivalent copies of the same
      ! substitution, so they share one host and one impurity distribution.
      host_atom = fc2_sc%defects(1,1)
      impurity_atom = fc2_sc%defects(3,1)
      host_type = S%ityp(host_atom)
      host_element = S%atm(host_type)
      if(len_trim(input%impurity_element) > 0) then
        impurity_element = input%impurity_element
      else
        impurity_element = S_sc%atm(S_sc%ityp(impurity_atom))
      endif
      if(len_trim(impurity_element) == 0) &
        call errore("dca_selfnrg", &
        "missing impurity symbol in defect FC file; set impurity_element in definput", 1)
      if(any(S%ityp(fc2_sc%defects(1,:)) /= host_type)) &
        call errore("dca_selfnrg", "equivalent defect sites have different host species", 1)
      ! NIST masses are in daltons.  The reference mass remains exactly the
      ! average host mass supplied by the user and stored internally by QE.
      call prepare_isotope_distribution(host_element, host_isotope_mass, &
        host_isotope_conc, host_isotope_cdf)
      call prepare_isotope_distribution(impurity_element, impurity_isotope_mass, &
        impurity_isotope_conc, impurity_isotope_cdf)
      host_reference_mass = S%amass(host_type) / MASS_DALTON_TO_RY
      ! host_natural_mass = sum(host_isotope_mass * host_isotope_conc)
      ! impurity_natural_mass = sum(impurity_isotope_mass * impurity_isotope_conc)
      if(ionode) then
        print*, "Explicit isotope scattering enabled"
        ! print"(A,A,A,F12.6)", "Host ", trim(host_element), &
        !   " reference mass (user average): ", host_reference_mass
        ! call print_isotope_distribution("Host", host_element, host_isotope_mass, host_isotope_conc)
        ! call print_isotope_distribution("Impurity", impurity_element, &
        ! impurity_isotope_mass, impurity_isotope_conc)
        ! if(abs(host_natural_mass-host_reference_mass) > 1.e-3_dp) &
        !   print"(A,F12.6,A,F12.6)", "WARNING: NIST natural host average ", &
        !     host_natural_mass, " differs from the user reference mass ", host_reference_mass
        ! print"(A,F12.6)", "Impurity natural average mass: ", impurity_natural_mass
      endif
    endif
    !
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
    allocate(Gf0(S%nat3*Nc, S%nat3*Nc))
    allocate(Gf_conf(S%nat3*Nc, S%nat3*Nc))
    allocate(Gf_avg(S%nat3*Nc, S%nat3*Nc))
    allocate(Gi_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(Gi_coarse(S%nat3, S%nat3, Nc))
    allocate(G0i_cluster(S%nat3, S%nat3, Nc))
    allocate(self_fine(S%nat3, S%nat3, in_grid%nqtot))
    allocate(self_before(S%nat3, S%nat3, Nc))
    allocate(self_next(S%nat3, S%nat3, Nc))
    allocate(Gi(S%nat3, S%nat3, in_grid_full%nqtot))

    allocate(kq(NQ, Nc))!, big_iq(grid_scat%nqtot))
    allocate(out_den_weights(S%nat3, out_grid%nqtot))
    allocate(pos(Nc))
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
    call center_V(xq, S, S_sc, fc2_sc, Vqqs)
    ! Explicit isotope disorder is present on every site, so it must always be
    ! configuration-averaged even when the impurity concentration is dilute.
    single_sample = Nc * n_eq_sites * conc < 1._dp
    !
    if(single_sample) then
      if(input%isotope_scattering) &
        call errore("dca_selfnrg", "single_sample approximation is incompatible with isotope_scattering", 1)
      if(ionode) print*, "Using low concentration approximation (at most 1 defect per configuration)"
      simulated_conc = conc
      allocate(VK(S%nat3*Nc, S%nat3*Nc,1))
      allocate(VM(S%nat3*Nc, S%nat3*Nc,1))
      allocate(V(S%nat3*Nc, S%nat3*Nc,1))
      VK(:,:,1) = flatten_RR_cmplx(Vqqs(:,:,:,:,1))
      VM = 0._dp
      na_def = fc2_sc%defects(1,1)
      do jq = 1, Nc
        do iq = 1, Nc
          do j1 = 1, 3
            VM((iq-1)*S%nat3+(na_def-1)*3+j1,(jq-1)*S%nat3+(na_def-1)*3+j1,1) = &
              fc2_sc%eps / Nc
          enddo
        enddo
      enddo
    else
      allocate(phases(Nc, Nc, Nc))
      allocate(phase_mat(Nc, Nc, S%nat, NSAMPLES))
      allocate(mass_phase_mat(Nc, Nc, S%nat, NSAMPLES))
      allocate(impurity_site(S%nat,Nc, NSAMPLES))
      allocate(site_mass_eps(S%nat,Nc))
      phase_mat = 0.0_dp
      mass_phase_mat = 0.0_dp
      impurity_site = .false.
      site_mass_eps = 0
      !
      defect_idef = 0
      do idef = 1, n_eq_sites
        defect_idef(fc2_sc%defects(1,idef)) = idef
      enddo
      !
      do iR = 1, Nc
        do k = 1, Nc
          do j = 1, Nc
            phases(j,k,iR) = e_iqr(xq(:,k)-xq(:,j), R(:,iR))
          enddo
        enddo
      enddo
      !
      do
        impurity_site = .false.
        do it = 1 + my_id, NSAMPLES, num_procs
          do iR = 1, Nc
            do idef = 1, n_eq_sites
              na_def = fc2_sc%defects(1,idef)
              call sample_defect(conc, impurity_site(na_def,iR,it))
            enddo
          enddo
        enddo

        ntot = count(impurity_site)
        call mpi_bsum(ntot)
        if (ntot == nint(real(NSAMPLES*Nc*n_eq_sites, dp) * conc)) exit
      enddo
      !
      if(ionode) print*, "found sampling with concentration:", &
        ntot / real(NSAMPLES*Nc*n_eq_sites, dp)

      do it = 1+my_id, NSAMPLES, num_procs
        do ipos = 1, Nc
          do na1 = 1, S%nat
            if(input%isotope_scattering) then
              if(impurity_site(na1,ipos,it)) then
                site_mass_eps(na1,ipos) = 1._dp - sample_isotope_mass( &
                  impurity_isotope_mass, impurity_isotope_cdf) / host_reference_mass
              else
                site_mass_eps(na1,ipos) = 1._dp - sample_isotope_mass( &
                  host_isotope_mass, host_isotope_cdf) / host_reference_mass
              endif
            else
              if(impurity_site(na1,ipos,it)) site_mass_eps(na1,ipos) = fc2_sc%eps
            endif
          enddo
        enddo
        !
        do ipos = 1, Nc
          do na1 = 1, S%nat
            do k = 1, Nc
              do j = 1, Nc
                if(impurity_site(na1,ipos,it)) then
                  phase_mat(j,k,na1,it) = phase_mat(j,k,na1,it) + &
                    phases(j,k,ipos)
                endif
                !
                mass_phase_mat(j,k,na1,it) = mass_phase_mat(j,k,na1,it) + &
                  site_mass_eps(na1,ipos) * phases(j,k,ipos)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(phases)
      !
      allocate(VK(S%nat3*Nc, S%nat3*Nc, NSAMPLES))
      allocate(VM(S%nat3*Nc, S%nat3*Nc, NSAMPLES))
      allocate(V(S%nat3*Nc, S%nat3*Nc, NSAMPLES))
      VK = 0._dp
      VM = 0._dp
      do it = 1+my_id, NSAMPLES, num_procs
        do jq = 1, Nc
          do iq = 1, Nc
            do na1 = 1, S%nat
              do j1 = 1, 3
                VM((iq-1)*S%nat3+(na1-1)*3+j1,(jq-1)*S%nat3+(na1-1)*3+j1,it) = &
                  VM((iq-1)*S%nat3+(na1-1)*3+j1,(jq-1)*S%nat3+(na1-1)*3+j1,it) + &
                  mass_phase_mat(iq,jq,na1,it) / Nc
              enddo
              idef = defect_idef(na1)
              if(idef == 0 .or. (idef /= 1 .and. single_sample)) cycle
              VK((iq-1)*S%nat3+1:iq*S%nat3,(jq-1)*S%nat3+1:jq*S%nat3,it) = &
                VK((iq-1)*S%nat3+1:iq*S%nat3,(jq-1)*S%nat3+1:jq*S%nat3,it) + &
                Vqqs(:,:,iq,jq,idef) * phase_mat(iq,jq,na1,it)
            enddo
          enddo
        enddo
      enddo
      !
    endif
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
    call quter_R(cluster_mesh, S%nat, S%tau, S%at, S%bg, img_xR, img_nR, img_weight)
    !
    call fourier_basis(c_out, out_grid)
    !
    !>
    self_next = 0._dp
    self_before = 0._dp
    !
    dos = 0._dp
    N_ITER_TOT = 0
    self_fine = 0._dp
    !
    if(ionode) print*, "Starting DCA self-energy calculation..."
    if(ionode) print*, ""
    do iw = 1, input%n_omega
      if(fc2_sc%def_type == "inclusion" .and. iw == 1) cycle
      V = VK + VM * wg%en(iw)**2
      Gi = 0._dp
      do iq = 1, in_grid_full%nq
        iqp = iq + in_grid_full%iq0
        Gi(:,:,iqp) = matmul(U_fine(:,:,iqp), &
          matmul(diag_cmplx(1/(wg%w(:,wg%e(iqp),iw)* Nc * NQ)), UT_fine(:,:,iqp)))
      enddo
      call mpi_bsum(S%nat3, S%nat3, in_grid_full%nqtot, Gi)
      do sc_iter = 1, MAXITER
        !
        if(ionode) print*, sc_iter
        Gf0 = 0._dp
        Gf0i = 0._dp
        Gi_coarse = 0._dp
        G0i_cluster = 0._dp
        do iq = 1+my_id, Nc, num_procs
          do jq = 1, NQ
            kkq = kq(jq,iq)
            A = Gi(:,:,kkq) - self_before(:,:,iq)
            call invzmat(S%nat3, A)
            Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) + A / NQ
          enddo
          call invzmat(S%nat3, Gi_coarse(:,:,iq))
          G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_before(:,:,iq)
          Gf0i((iq-1)*S%nat3+1:iq*S%nat3,(iq-1)*S%nat3+1:iq*S%nat3) = G0i_cluster(:,:,iq)
          if(single_sample) then
            A = G0i_cluster(:,:,iq)
            call invzmat(S%nat3, A)
            Gf0((iq-1)*S%nat3+1:iq*S%nat3,(iq-1)*S%nat3+1:iq*S%nat3) = A
          endif
        enddo
        call mpi_bsum(S%nat3, S%nat3, Nc, G0i_cluster)
        call mpi_bsum(S%nat3*Nc, S%nat3*Nc, Gf0i)
        if(single_sample) call mpi_bsum(S%nat3*Nc, S%nat3*Nc, Gf0)
        !
        if(single_sample) then
          Gf_conf = Gf0i - V(:,:,1)
          call invzmat(S%nat3*Nc, Gf_conf)
          Gf_avg = Gf0 * (1-Nc*conc*n_eq_sites) + Gf_conf * (Nc * conc * n_eq_sites)
          !> self energy is G0_cluster^-1 - <G>^-1
        else
          Gf_avg = 0.0_dp
          do it = 1+my_id, NSAMPLES, num_procs
            ! > construction of G_conf for a given configuration
            ! Gi_conf = 0.0_dp
            !
            Gf_conf = Gf0i - V(:,:,it)
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
          self_next(:,:,iq) = G0i_cluster(:,:,iq) - Gi_avg(:,:,iq,iq)
          call apply_sym_q(S, xq(:,iq), self_next(:,:,iq))
          call flip_positive_imag_eigs(S%nat3, self_next(:,:,iq))
        enddo
        !
        ! call apply_sym(S%at, S%bg, S%nat, S%ityp, S%tau, self_next, c_equiv, xq, .true.)
        ! self_before = self_next
        ! exit
        delta_out = reshape(self_next, [S%nat3**2*Nc])
        call mix_broyden_full(S%nat3**2*Nc, delta_out, delta_in, &
          ALPHA_MIX, sc_iter, MEMORY, df, dv)
        self_before = reshape(delta_in, [S%nat3, S%nat3, Nc])
        !
        !
        if(all(abs(real(self_before - self_next, dp)) < ABS_TOLERANCE + abs(real(self_next, dp)) * REL_TOLERANCE) .and. &
          all(abs(aimag(self_before - self_next)) < ABS_TOLERANCE + abs(aimag(self_next)) * REL_TOLERANCE)) exit
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
      do iq = 1, out_grid%nqtot
        A = 0._dp
        do jq = 1+my_id, Nc, num_procs
          A = A + self_before(:,:,jq) * c_out(:,:,jq,iq)
        enddo
        call mpi_bsum(S%nat3, S%nat3, A)
        call apply_sym_q(S, out_grid%xq(:,iq), A)
        ! call invzmat(S%nat3, A)
        ! A = A + ialpha
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
     case('dos')
      do iw = 1+my_id, input%n_omega, num_procs
        call tetra_from_self(S, out_grid, out_freqs, self_out_grid(:,:,:,iw), wg%en(iw)**2, out_den_weights)
        dos(iw) = sum(matmul(AIMAG(out_den_weights), wg_out%qw)) * product(out_grid%n)
      enddo
      call mpi_bsum(input%n_omega, dos)
      call write_dos('dos-dca.dat', wg%en, dos)
     case("test")
      do iw = 1+my_id, input%n_omega, num_procs
        call tetra_from_self(S, out_grid, out_freqs, self_out_grid(:,:,:,iw), wg%en(iw)**2, out_den_weights)
        dos(iw) = sum(matmul(AIMAG(out_den_weights), wg_out%qw)) * product(out_grid%n)
      enddo
      call mpi_bsum(input%n_omega, dos)
      call write_dos('dos-dca-test.dat', wg%en, dos)
      call write_spf_ndiag('spf-dca-ndiag-test.dat', wg%en, self_out_grid(:,:,:1,:), out_freqs(:,:1))
      call write_spf('spf-dca-test.dat', wg%en, self_out_diag(:,:1,:), out_freqs(:,:1))
      call write_self('self-dca-test.dat', wg%en, self_out_diag(:,:1,:))
     case default
      if(ionode) print*, "WARNING: unknown DCA calculation type"
    end select
    !
  contains
    !
    subroutine prepare_isotope_distribution(element_name, isotope_mass, isotope_conc, isotope_cdf)
      character(len=*), intent(in) :: element_name
      real(dp), allocatable, intent(out) :: isotope_mass(:), isotope_conc(:), isotope_cdf(:)
      integer :: ii
      !
      call get_natural_isotopes(element_name, isotope_mass, isotope_conc)
      allocate(isotope_cdf(size(isotope_conc)))
      isotope_cdf(1) = isotope_conc(1)
      do ii = 2, size(isotope_conc)
        isotope_cdf(ii) = isotope_cdf(ii-1) + isotope_conc(ii)
      enddo
      isotope_cdf(size(isotope_cdf)) = 1._dp
    end subroutine prepare_isotope_distribution
    !
    function sample_isotope_mass(isotope_mass, isotope_cdf) result(sampled_mass)
      real(dp), intent(in) :: isotope_mass(:), isotope_cdf(:)
      real(dp) :: sampled_mass, random_value
      integer :: ii
      !
      call random_number(random_value)
      sampled_mass = isotope_mass(size(isotope_mass))
      do ii = 1, size(isotope_mass)
        if(random_value < isotope_cdf(ii)) then
          sampled_mass = isotope_mass(ii)
          exit
        endif
      enddo
    end function sample_isotope_mass
    !
    subroutine print_isotope_distribution(label, element_name, isotope_mass, isotope_conc)
      character(len=*), intent(in) :: label, element_name
      real(dp), intent(in) :: isotope_mass(:), isotope_conc(:)
      integer :: ii
      !
      print"(A,A,A)", trim(label), " isotope distribution for ", trim(element_name)
      do ii = 1, size(isotope_mass)
        if(isotope_conc(ii) > 0._dp) &
          print"(2X,F12.6,2X,F10.6)", isotope_mass(ii), isotope_conc(ii)
      enddo
    end subroutine print_isotope_distribution
    !
    subroutine fourier_basis(c_Qq, grid)
      complex(dp), allocatable, intent(out) :: c_Qq(:,:,:,:)
      type(q_grid), intent(in) :: grid
      !
      allocate(c_Qq(S%nat3, S%nat3, Nc, grid%nqtot))
      c_Qq = 0._dp
      do jq = 1, grid%nqtot
        do iq = 1, Nc
          do na2 = 1, S%nat
            do na1 = 1, S%nat
              do iR = 1, img_nR(na1,na2)
                c_Qq(3*(na1-1)+1:3*na1,3*(na2-1)+1:3*na2, iq, jq) = &
                  c_Qq(3*(na1-1)+1:3*na1,3*(na2-1)+1:3*na2, iq, jq) + &
                  e_iqr(xq(:,iq)- grid%xq(:,jq), img_xR(:,iR,na1,na2)) * img_weight(iR,na1,na2) / Nc
              enddo
            enddo
          enddo
        enddo
      enddo
    end subroutine
    !
    subroutine sample_defect(c, defect_here)
      real(dp), intent(in) :: c
      logical, intent(out) :: defect_here
      !
      real(dp) :: random_value
      !
      call random_number(random_value)
      if(random_value < c) then
        defect_here = .true.
      else
        defect_here = .false.
      endif
    end subroutine sample_defect
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
    real(dp):: trans(3), shift_real(3)
    integer, dimension(3,3) :: SM, SMT, SM1, SMT1
    integer, dimension(3,S%nat) :: atom_shift
    type(forceconst2_sc) :: fc_temp
    INTEGER :: isym_map, isym, nai, naj, naii, najj, iR1, iR2
    integer :: site, N, iR, jR, iiR, jjR, iq, jq, inv_map, Nq
    real(dp), allocatable :: tens4(:,:,:,:,:,:)
    real(dp) :: tau_cryst(3,S%nat)
    REAL(DP), ALLOCATABLE :: work(:,:,:,:,:,:)
    character(100) :: filename
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
    do nai = 1, S%nat
      tau_cryst(:,nai) = cryst2cart(S%tau(:,nai), S%bg, -1)
    enddo
    !
    do site = 1, size(fc2_sc%defects,2)
      call fc_temp%allocate(S, S_sc, fc2_sc%nq)
      if (site == 1) then
        fc_temp%fc = fc2_sc%fc
        fc_temp%taudef = fc2_sc%taudef
      else
        isym_map = -1
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
        if (isym_map == -1) call errore("center_V", "could not find symmetry mapping equivalent defect site", site)
        do nai = 1, S%nat
          naii = irt(isym_map, nai)
          shift_real = matmul(transpose(real(SM, dp)), tau_cryst(:,nai)) - &
            ft(:,isym_map) - tau_cryst(:,naii)
          atom_shift(:,nai) = nint(shift_real)
          if (any(abs(shift_real - real(atom_shift(:,nai), dp)) > 1e-8_dp)) &
            call errore("center_V", "symmetry atom shift is not a lattice vector", nai)
        enddo
        !
        work = 0.0_DP
        DO nai = 1, S%nat
          naii = irt(isym_map, nai)
          DO naj = 1, S%nat
            najj = irt(isym_map, naj)
            do iR = 1, N
              do jR = 1, N
                iiR = v2index(bz2simple(matmul(SMT, index2v(iR, fc2_sc%nq)) + &
                  atom_shift(:,nai), fc2_sc%nq), fc2_sc%nq)
                jjR = v2index(bz2simple(matmul(SMT, index2v(jR, fc2_sc%nq)) + &
                  atom_shift(:,naj), fc2_sc%nq), fc2_sc%nq)
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
        fc_temp%taudef = S%tau(:,fc2_sc%defects(1,site))
      endif
      call fc_temp%center(fc2_sc%nq, S)
      write(filename, "(A,I3.3,A)") "fc2_sc", site, ".dat"
      call write_fc2_sc_RR(S, fc_temp, trim(filename))
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
    Vqqs = Vqqs / real(Nq, dp)
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
