module dca
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  use defect, only : tetra_from_self
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
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
  use symm_base, only : symm_matrices => s
  use EPW_utilities, only : mix_broyden_full
  use mpi_thermal, only : mpi_bsum, ionode, num_procs, my_id
contains
  subroutine dca_selfnrg(S, S_sc, input, fc2, fc2_sc, grid)
    type(ph_system_info), intent(in) :: S, S_sc
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: grid
    type(tetra_output) :: wg
    !
    complex(dp), allocatable, dimension(:) :: delta_in, delta_out
    complex(dp), allocatable, dimension(:,:) :: den_weights, &
      den_eig, Gf_conf, Gf_avg, df, dv, overlap
    complex(dp), allocatable, dimension(:,:,:) :: den_UL, den_UR, &
      G0i_cluster, G_coarse, Gi_coarse, self_fine, self_R
    complex(dp), allocatable, dimension(:,:,:,:) :: G_avg, Gi_conf, &
      Gi_avg, self_in, self_out
    complex(dp), allocatable, dimension(:,:,:,:,:) :: Vqqs, &
      self_uf, V, phase_mat
    real(dp), allocatable, dimension(:,:) :: self_xR, xq, R
    integer, allocatable :: pos(:,:), kq(:,:), big_iq(:), Npos(:)
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, n_eq_sites, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ, kkq, cluster_mesh(3)
    real(dp) :: conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX, max_diff
    real(dp) :: dos(input%n_omega)
    logical :: conv
    !
    NSAMPLES = 300
    MAXITER = 50
    ABS_TOLERANCE = 1e-13_dp
    REL_TOLERANCE = 1e-4_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 4
    conc = 1e-3_dp  ! example concentration
    cluster_mesh = fc2%nq
    Nc = product(cluster_mesh) !* size(fc2_sc%defects,2)
    NQ = product(grid%n) / Nc
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true., natural=.true.)
    if (any(mod(grid%n, cluster_mesh) /= 0)) &
      call errore("dca_selfnrg", "grid size not multiple of fc2 grid size", 1)
    R = grid_vec_cart(cluster_mesh, S%at, natural=.true.)
    call set_wg(S, fc2, grid, input%n_omega, wg)

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
    allocate(V(S%nat3, S%nat3, Nc, Nc, NSAMPLES))
    allocate(phase_mat(S%nat3, S%nat3, Nc, Nc, Nc))
    allocate(kq(NQ, Nc), big_iq(grid%nqtot))
    allocate(den_weights(S%nat3, grid%nqtot))
    allocate(den_UL(S%nat3, S%nat3, grid%nqtot))
    allocate(den_UR(S%nat3, S%nat3, grid%nqtot))
    allocate(overlap(S%nat3, grid%nqtot))
    allocate(den_eig(S%nat3, grid%nqtot))
    allocate(pos(Nc, NSAMPLES))
    allocate(Npos(NSAMPLES))

    ! deallocate(fc2_sc%defects)
    ! allocate(fc2_sc%defects(3, 1))
    ! fc2_sc%defects(:,1) = [1,1,1]
    n_eq_sites = size(fc2_sc%defects,2)
    call center_V(S, S_sc, fc2_sc, cluster_mesh, Vqqs)
    !
    !> construction of of the phase, used to translate the potential to
    !> random configurations inside the cluster
    do i = 1, Nc
      do j = 1, Nc
        do k = 1, Nc
          phase_mat(:,:,k,j,i) = e_iqr(xq(:,j)-xq(:,k), R(:,i))
        enddo
      enddo
    enddo

    V = 0.0_dp
    do it = 1, NSAMPLES
      do idef = 1, n_eq_sites
        call sample_binomial(Nc, conc/n_eq_sites, pos(:,it), Npos(it))
        do ipos = 1, Npos(it)
          V(:,:,:,:,it) = V(:,:,:,:,it) + &
            Vqqs(:,:,:,:,idef) * phase_mat(:,:,:,:,pos(ipos,it)) / Nc
        enddo
      enddo
    enddo
    !
    !
    !> construction of G0_coarse and G0i_coarse, which are averaged over the small
    !> patch in which self-energy is assumed constant
    do iq = 1, Nc
      do jq = 1, NQ
        kq(jq,iq) = v2index_n(index2v_n(jq, grid%n / cluster_mesh) + &
          index2v_n(iq, cluster_mesh) * grid%n / cluster_mesh, grid%n)
        big_iq(wg%e(kq(jq,iq))) = iq
      enddo
    enddo
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
    print*, "Starting DCA self-energy calculation..."
    do iw = 70+my_id, input%n_omega, num_procs
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
          self_uf(:,:,:,:,iq) = unflatten_RR_cmplx(self_in(:,:,iq,iw), S%nat, S%nat)
        enddo
        CALL quter_cmplx(input%sc_grid(1), input%sc_grid(2), input%sc_grid(3), &
          S%nat, S%tau, S%at, S%bg, self_uf, xq, self_R, self_xR, 2)
        do iq = 1, grid%nqtot
          call fftinterp_mat2_cmplx(grid%xq(:,iq), S, self_R, self_xR, self_fine(:,:,iq))
        enddo
        !
        call tetra_from_self(S, grid, wg%f, self_fine, wg%en(iw)**2, &
          den_weights, den_UL, den_UR, overlap)
        !}
        Gi_coarse = 0.0_dp
        do iq = 1, Nc
          do jq = 1, NQ
            kkq = wg%e(kq(jq,iq))
            do i = 1, S%nat3
              do j = 1, S%nat3
                do k = 1, S%nat3
                  Gi_coarse(i,j,iq) = Gi_coarse(i,j,iq) + &
                    den_UR(i,k,kkq) * conjg(den_UL(j,k,kkq)) * &
                    den_weights(k,kkq) / overlap(k,kkq)
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
        ! Gf_avg = 0.0_dp
        ! do it = 1, NSAMPLES
        !   !> construction of V for a given configuration
        !   !>
        !   !
        !   !> construction of G_conf for a given configuration
        !   G_conf = 0.0_dp
        !   do iq = 1, Nc
        !     G_conf(:,:,iq,iq) = G0i_cluster(:,:,iq)
        !     do jq = 1, Nc
        !       G_conf(:,:,iq,jq) = G_conf(:,:,iq,jq) - V(:,:,iq,jq,it)
        !     enddo
        !   enddo
        !   !>
        !   !
        !   Gf_conf = flatten_RR_cmplx(G_conf)
        !   call invzmat(S%nat3*Nc, Gf_conf)
        !   !
        !   Gf_avg = Gf_avg + Gf_conf
        !   !
        ! enddo
        ! Gf_avg = Gf_avg / real(NSAMPLES, dp)
        Gf_avg = 0.0_dp
        Gi_conf = 0.0_dp
        do iq = 1, Nc
          Gi_conf(:,:,iq,iq) = G0i_cluster(:,:,iq)
        enddo
        Gf_conf = flatten_RR_cmplx(Gi_conf)
        call invzmat(S%nat3*Nc, Gf_conf)
        Gf_avg = Gf_avg + Gf_conf * (1 - n_eq_sites * Nc * conc)
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
        !>
        !
        !> self energy is G0_cluster^-1 - <G>^-1
        call invzmat(S%nat3*Nc, Gf_avg)
        Gi_avg = unflatten_RR_cmplx(Gf_avg, Nc, Nc)
        do iq = 1, Nc
          self_out(:,:,iq,iw) = G0i_cluster(:,:,iq) - Gi_avg(:,:,iq,iq)
        enddo
        !>
        !
        delta_out = reshape(self_out(:,:,:,iw), [S%nat3**2*Nc])
        !
        call mix_broyden_full(S%nat3**2*Nc, delta_out, delta_in, &
          ALPHA_MIX, sc_iter, MEMORY, df, dv)
        self_in(:,:,:,iw) = reshape(delta_in, [S%nat3, S%nat3, Nc])
      enddo ! self-energy SC cycle
      do iq = 1, grid%nqtot
        dos(iw) = dos(iw) + aimag(sum(den_weights(:,iq))) * wg%qw(iq)
      enddo
    enddo ! frequency loop
    call mpi_bsum(input%n_omega, dos)
    !
    open(110, file="dos_dca.dat")
    do iw = 1, input%n_omega
      if(ionode) WRITE(110, "(3E20.8)") wg%en(iw) * RY_TO_CMM1, &
        - dos(iw) / pi * product(grid%n) * 2 * wg%en(iw) / RY_TO_CMM1
    enddo
    close(110)
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
  SUBROUTINE center_V( S, S_sc, fc2_sc, cluster_mesh, Vqqs)
    !-----------------------------------------------------------------------
    ! Apply ONE symmetry that maps atom1 -> atom2 to center potential there
    !
    USE kinds,     ONLY : DP
    USE symm_base, ONLY : irt, nsym, ft, symm_mat => s, invs
    use symme, only : cart_to_crys, crys_to_cart
    !
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
    real(dp), allocatable :: xq(:,:)
    !
    call fc_temp%allocate(S, S_sc, fc2_sc%nq)
    !
    Nq = product(cluster_mesh)
    N = product(fc2_sc%nq)
    xq = grid_vec_cart(cluster_mesh, S%bg, divide=.true., natural=.true.)
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
