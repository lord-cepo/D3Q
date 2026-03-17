module dca
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils, only: v2index_n, index2v_n, bz2simple, grid_vec_cart, &
    e_iqr, index2v, cryst2cart, v2index, freq_in_grid
  use defutils, only : flatten_RR_cmplx, unflatten_RR_cmplx, &
    quter_cmplx, fftinterp_mat2_cmplx
  use defect, only : tetra_from_self, write_spf_ndiag, write_spf, write_self, write_dos
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
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
  use mpi_thermal, only : mpi_bsum, ionode, num_procs, my_id
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
  subroutine dca_selfnrg(S, S_sc, input, fc2, fc2_sc, grid, out_grid)
    type(ph_system_info), intent(in) :: S, S_sc
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid, out_grid
    type(tetra_output) :: wg, wg_out
    !
    complex(dp), allocatable, dimension(:) :: delta_in, delta_out, G__
    complex(dp), allocatable, dimension(:,:) :: den_weights, &
      den_eig, Gf_conf, Gf_avg, df, dv, overlap, V__, I_gV__, G_avg__
    complex(dp), allocatable, dimension(:,:,:) :: den_UL, den_UR, &
      G0i_cluster, G_coarse, Gi_coarse, self_fine, self_R, U, UT, &
      UT_fine, U_fine, UT_out, U_out, self_out_diag
    complex(dp), allocatable, dimension(:,:,:,:) :: G_avg, Gi_conf, &
      Gi_avg, self_in, self_out, phase_mat, V, self_out_grid
    complex(dp), allocatable, dimension(:,:,:,:,:) :: Vqqs, &
      self_uf
    real(dp), allocatable, dimension(:,:) :: self_xR, xq, R, f, out_freqs
    integer, allocatable :: pos(:), kq(:,:)!, big_iq(:)
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, n_eq_sites, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ, kkq, &
      cluster_mesh(3), Npos, ntot
    real(dp) :: conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX, max_diff
    real(dp) :: dos(input%n_omega), shift(3)
    logical :: conv
    complex(dp) :: A(S%nat3,S%nat3)
    !
    if(input%calculation == "test") call init_random_seed()
    !
    MAXITER = 50
    ABS_TOLERANCE = 1e-13_dp
    REL_TOLERANCE = 1e-3_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 4
    conc = input%conc ! example concentration
    cluster_mesh = input%sc_grid
    Nc = product(cluster_mesh) !* size(fc2_sc%defects,2)
    n_eq_sites = size(fc2_sc%defects,2)
    NSAMPLES = max(nint(1e2_dp / Nc / conc / n_eq_sites), 1)
    NQ = product(grid%n) / Nc
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true., natural=.true.)
    shift = xq(:,v2index_n([1,1,1],cluster_mesh)) / 2
    if (NQ == 1) shift = 0.0_dp
    do iq = 1, size(xq,2)
      xq(:,iq) = xq(:,iq) + shift
    enddo
    if (any(mod(grid%n, cluster_mesh) /= 0)) &
      call errore("dca_selfnrg", "grid size not multiple of fc2 grid size", 1)
    R = grid_vec_cart(cluster_mesh, S%at, natural=.true.)
    call set_wg(S, fc2, grid, input%n_omega, wg)
    call set_wg(S, fc2, out_grid, input%n_omega, wg_out)
    !
    allocate(UT_fine(S%nat3,S%nat3,grid%nqtot), U_fine(S%nat3,S%nat3,grid%nqtot))
    allocate(f(S%nat3,Nc), U(S%nat3,S%nat3,Nc), UT(S%nat3,S%nat3,Nc))
    do iq = 1, grid%nqtot
      call freq_phq_safe(grid%xq(:,iq), S, fc2, f(:,1), U_fine(:,:,iq))
      UT_fine(:,:,iq) = conjg(transpose(U_fine(:,:,iq)))
    enddo
    do iq = 1, Nc
      call freq_phq_safe(xq(:,iq), S, fc2, f(:,iq), U(:,:,iq))
      UT(:,:,iq) = conjg(transpose(U(:,:,iq)))
    enddo
    if(ionode) print*, "DCA cluster size:", Nc, "number of configurations to be averaged:", NSAMPLES
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
    allocate(G_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(Gi_conf(S%nat3, S%nat3, Nc, Nc))
    allocate(Gf_conf(S%nat3*Nc, S%nat3*Nc))
    allocate(Gf_avg(S%nat3*Nc, S%nat3*Nc))
    allocate(Gi_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(G_coarse(S%nat3, S%nat3, Nc))
    allocate(Gi_coarse(S%nat3, S%nat3, Nc))
    allocate(G0i_cluster(S%nat3, S%nat3, Nc))
    allocate(self_in(S%nat3, S%nat3, Nc, input%n_omega))
    allocate(self_uf(3, 3, S%nat, S%nat, Nc))
    allocate(self_fine(S%nat3, S%nat3, grid%nqtot))
    allocate(self_out(S%nat3, S%nat3, Nc, input%n_omega))
    allocate(V(S%nat3, S%nat3, Nc, Nc))
    allocate(phase_mat(Nc, Nc, n_eq_sites, NSAMPLES))
    allocate(kq(NQ, Nc))!, big_iq(grid%nqtot))
    allocate(den_weights(S%nat3, grid%nqtot))
    allocate(den_UL(S%nat3, S%nat3, grid%nqtot))
    allocate(den_UR(S%nat3, S%nat3, grid%nqtot))
    allocate(overlap(S%nat3, grid%nqtot))
    allocate(den_eig(S%nat3, grid%nqtot))
    allocate(pos(Nc))
    allocate(G__(S%nat3*Nc))
    allocate(V__(S%nat3*Nc, S%nat3*Nc))
    allocate(I_gV__(S%nat3*Nc, S%nat3*Nc))
    allocate(G_avg__(S%nat3*Nc, S%nat3*Nc))

    ! allocate(Npos(NSAMPLES))

    ! deallocate(fc2_sc%defects)
    ! allocate(fc2_sc%defects(3, 1))
    ! fc2_sc%defects(:,1) = [1,1,1]

    if(ionode) print*, NSAMPLES, "DCA samples to be used"
    call center_V(xq, S, S_sc, fc2_sc, cluster_mesh, Vqqs)
    do idef = 1, n_eq_sites
      do iq = 1, Nc
        do jq = 1, Nc
          Vqqs(:,:,iq,jq,idef) = matmul(UT(:,:,iq), matmul(Vqqs(:,:,iq,jq,idef),U(:,:,jq))) / Nc
        enddo
      enddo
    enddo
    V__ = flatten_RR_cmplx(Vqqs(:,:,:,:,1))

    !
    ! !> construction of of the phase, used to translate the potential to
    ! !> random configurations inside the cluster

    ntot = 0
    phase_mat = 0.0_dp
    !
    do it = 1, NSAMPLES
      do idef = 1, n_eq_sites
        call sample_binomial(Nc, conc, pos, Npos)
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
    if(ionode) print"(A,E15.4)", "simulated concentration:", real(ntot,dp) / (Nc * n_eq_sites * NSAMPLES)
    !
    !
    !> construction of G0_coarse and G0i_coarse, which are averaged over the small
    !> patch in which self-energy is assumed constant
    do iq = 1, Nc
      do jq = 1, NQ
        kq(jq,iq) = v2index_n(index2v_n(jq, grid%n / cluster_mesh) + &
          index2v_n(iq, cluster_mesh) * grid%n / cluster_mesh, grid%n)
        ! big_iq(wg%e(kq(jq,iq))) = iq
      enddo
    enddo
    !
    !>
    !
    self_out = 1._dp
    self_in = 0._dp
    allocate(df(S%nat3**2*Nc, MEMORY))
    allocate(dv(S%nat3**2*Nc, MEMORY))
    allocate(delta_in(S%nat3**2*Nc))
    allocate(delta_out(S%nat3**2*Nc))
    delta_out = 0._dp
    !
    dos = 0._dp
    !
    if(ionode) print*, "Starting DCA self-energy calculation..."
    do iw = 1+my_id, input%n_omega, num_procs
      ! do iq = 1, Nc
      !   do i = 1, S%nat3
      !     g__(i+S%nat3*(iq-1)) = wg%w(i, wg%e(kq(1,iq)), iw) * Nc
      !   enddo
      !   ! print*, sum(abs(wg%en(iw)**2 - f(:,iq)**2 - 1/g__(S%nat3*(iq-1)+1:S%nat3*iq))) / sum(abs(wg%en(iw)**2 - f(:,iq)**2))
      ! enddo
      df = 0._dp
      dv = 0._dp
      delta_in = 0._dp
      do sc_iter = 1, MAXITER
        conv = .true.
        max_diff = 0._dp
        do iq = 1, Nc
          do i = 1, S%nat3
            max_diff = max(max_diff, ABS(self_out(i,i,iq,iw)-self_in(i,i,iq,iw)))
            if( ABS(self_out(i,i,iq,iw)-self_in(i,i,iq,iw)) > REL_TOLERANCE * &
              abs(self_in(i,i,iq,iw)) + ABS_TOLERANCE) conv = .false.
          enddo
        enddo
        print*, iw, sc_iter, max_diff
        if (conv) exit
        !> construction of G0_cluster from self-energy
        !
        do iq = 1, Nc
          A = matmul(U(:,:,iq), matmul(self_in(:,:,iq,iw), UT(:,:,iq)))
          self_uf(:,:,:,:,iq) = unflatten_RR_cmplx(A, S%nat, S%nat)
        enddo
        CALL quter_cmplx(input%sc_grid(1), input%sc_grid(2), input%sc_grid(3), &
          S%nat, S%tau, S%at, S%bg, self_uf, xq, self_R, self_xR, 2)
        !
        ! write(filename, "(A,I2.2,A)") "self_R_", iw, ".dat"
        ! open(10, file=filename)
        ! do iR = 1, size(self_xR,2)
        !   do i = 1, S%nat
        !     do j = 1, S%nat
        !       write(10, "(3E20.8)") norm2(self_xR(:,iR) + S%tau(:,i) - S%tau(:,j)), &
        !         sum(abs(self_R((i-1)*3+1:i*3,(j-1)*3+1:j*3,iR)))
        !     enddo
        !   enddo
        ! enddo
        ! close(10)
        do iq = 1, grid%nqtot
          call fftinterp_mat2_cmplx(grid%xq(:,iq), S, self_R, self_xR, self_fine(:,:,iq))
          self_fine(:,:,iq) = matmul(UT_fine(:,:,iq), matmul(self_fine(:,:,iq), U_fine(:,:,iq)))
          where(aimag(self_fine(:,:,iq)) > 0._dp) self_fine(:,:,iq) = conjg(self_fine(:,:,iq))
        enddo
        !
        call tetra_from_self(S, grid, wg%f, self_fine, wg%en(iw)**2, &
          den_weights, den_UL, den_UR, overlap)
        !}
        do iq = 1, grid%nqtot
          call merge_degen(S%nat3, den_weights(:,iq), wg%f(:,iq))
        enddo
        Gi_coarse = 0.0_dp
        do iq = 1, Nc
          do jq = 1, NQ
            kkq = wg%e(kq(jq,iq))
            do k = 1, S%nat3
              do j = 1, S%nat3
                do i = 1, S%nat3
                  Gi_coarse(i,j,iq) = Gi_coarse(i,j,iq) + &
                    den_UR(i,k,kkq) * conjg(den_UL(j,k,kkq)) * &
                    den_weights(k,kkq) / overlap(k,kkq) * Nc
                enddo
              enddo
            enddo
          enddo
          call invzmat(S%nat3, Gi_coarse(:,:,iq))
          G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_in(:,:,iq,iw)
        enddo
        !>
        !
        !> construction of the flatten average Gf_avg over niter configurations
        !----------------------------------------------------------------------
        Gf_avg = 0.0_dp
        do it = 1, NSAMPLES
          ! > construction of G_conf for a given configuration
          Gi_conf = 0.0_dp
          V = 0.0_dp
          do idef = 1, n_eq_sites
            do jq = 1, Nc
              do iq = 1, Nc
                V(:,:,iq,jq) = Vqqs(:,:,iq,jq,idef) * phase_mat(iq,jq,idef,it)
              enddo
            enddo
          enddo
          !
          do jq = 1, Nc
            Gi_conf(:,:,jq,jq) = G0i_cluster(:,:,jq)
            do iq = 1, Nc
              Gi_conf(:,:,iq,jq) = Gi_conf(:,:,iq,jq) - V(:,:,iq,jq)
            enddo
          enddo
          !>
          !
          Gf_conf = flatten_RR_cmplx(Gi_conf)
          call invzmat(S%nat3*Nc, Gf_conf)
          !
          Gf_avg = Gf_avg + Gf_conf
          !
        enddo
        Gf_avg = Gf_avg / real(NSAMPLES, dp)
        !----------------------------------------------------------------------

        ! Gi_conf = 0.0_dp
        ! do iq = 1, Nc
        !   Gi_conf(:,:,iq,iq) = G0i_cluster(:,:,iq)
        ! enddo
        ! Gf_conf = flatten_RR_cmplx(Gi_conf)
        ! call invzmat(S%nat3*Nc, Gf_conf)
        ! Gf_avg = Gf_conf * (1 - conc * n_eq_sites * Nc)
        ! ! Gf_avg = diag_cmplx(g__ * (1 - n_eq_sites * Nc * conc))
        ! do idef = 1, n_eq_sites
        !   Gi_conf = 0.0_dp
        !   do iq = 1, Nc
        !     Gi_conf(:,:,iq,iq) = G0i_cluster(:,:,iq)
        !     do jq = 1, Nc
        !       Gi_conf(:,:,iq,jq) = Gi_conf(:,:,iq,jq) - Vqqs(:,:,iq,jq,idef)
        !     enddo
        !   enddo
        !   Gf_conf = flatten_RR_cmplx(Gi_conf)
        !   call invzmat(S%nat3*Nc, Gf_conf)
        !   Gf_avg = Gf_avg + Gf_conf * (conc * Nc)
        ! enddo
        !>
        !
        !> self energy is G0_cluster^-1 - <G>^-1
        Gi_avg = unflatten_RR_cmplx(Gf_avg, Nc, Nc)
        do iq = 1, Nc
          call invzmat(S%nat3, Gi_avg(:,:,iq,iq))
          self_out(:,:,iq,iw) = G0i_cluster(:,:,iq) - Gi_avg(:,:,iq,iq)
          !diag(wg%en(iw)**2 - f(:,iq)**2) - Gi_avg(:,:,iq,iq)
        enddo
        !>
        !
        delta_out = reshape(self_out(:,:,:,iw), [S%nat3**2*Nc])
        !
        call mix_broyden_full(S%nat3**2*Nc, delta_out, delta_in, &
          ALPHA_MIX, sc_iter, MEMORY, df, dv)
        self_in(:,:,:,iw) = reshape(delta_in, [S%nat3, S%nat3, Nc])
      enddo ! self-energy SC cycle
      ! print*, self_out(4,4,1,iw), self_out(5,5,2,iw)
      do iq = 1, out_grid%nqtot
        call fftinterp_mat2_cmplx(out_grid%xq(:,iq), S, self_R, self_xR, self_out_grid(:,:,iq,iw))
        self_out_grid(:,:,iq,iw) = matmul(UT_out(:,:,iq), matmul(self_out_grid(:,:,iq,iw), U_out(:,:,iq)))
        where(aimag(self_out_grid(:,:,iq,iw)) > 0._dp) self_out_grid(:,:,iq,iw) = conjg(self_out_grid(:,:,iq,iw))
        do i = 1, S%nat3
          self_out_diag(i, iq, iw) = self_out_grid(i,i,iq,iw)
        enddo
      enddo
    enddo ! frequency loop
    call mpi_bsum(S%nat3, S%nat3, out_grid%nqtot, input%n_omega, self_out_grid)
    call mpi_bsum(S%nat3, out_grid%nqtot, input%n_omega, self_out_diag)
    !
    select case(input%calculation)
     case ('spf-def')
      call write_spf_ndiag('spf-dca-ndiag.dat', wg%en, self_out_grid, out_freqs)
      call write_spf('spf-dca.dat', wg%en, self_out_diag, out_freqs)
      call write_self('self-dca.dat', wg%en, self_out_diag)
     case('self')
      call write_self('self-dca.dat', wg%en, self_out_diag)
     case('dos')
      do iw = 1+my_id, input%n_omega, num_procs
        call tetra_from_self(S, out_grid, out_freqs, self_out_grid(:,:,:,iw), wg_out%en(iw)**2, den_weights)
        dos(iw) = sum(matmul(AIMAG(den_weights), wg_out%qw)) * product(out_grid%n)
      enddo
      call mpi_bsum(input%n_omega, dos)
      call write_dos('dos-dca.dat', wg_out%en, dos)
     case("test")
      do iw = 1+my_id, input%n_omega, num_procs
        call tetra_from_self(S, out_grid, out_freqs, self_out_grid(:,:,:,iw), wg_out%en(iw)**2, den_weights)
        dos(iw) = sum(matmul(AIMAG(den_weights), wg_out%qw)) * product(out_grid%n)
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
  SUBROUTINE center_V(xq, S, S_sc, fc2_sc, cluster_mesh, Vqqs)
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
    integer, dimension(3), intent(in) :: cluster_mesh
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
    call fc_temp%allocate(S, S_sc, fc2_sc%nq)
    !
    Nq = product(cluster_mesh)
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
      call fc_temp%center(fc2_sc%nq, S)
      !
      ! call translate_xR(fc_temp, trans)
      do iq = 1, size(xq,2)
        call fc_temp%r2q(xq(:,iq))
        do jq = 1, size(xq,2)
          call fc_temp%r2q(xq(:,jq), Vqqs(:,:,iq,jq,site))
        enddo
      enddo
      ! call translate_xR(fc_temp, -trans)
      !
    enddo
    !
    call fc_temp%deallocate()
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
