module dca
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
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
  use constants, only : RY_TO_CMM1
  use symm_base, only : symm_matrices => s
  use EPW_utilities, only : mix_broyden_full
contains
  subroutine dca_selfnrg(S, S_sc, input, fc2, fc2_sc, grid)
    type(ph_system_info), intent(in) :: S, S_sc
    type(code_input_type), intent(in) :: input
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc), intent(in) :: fc2_sc
    type(q_grid), intent(in) :: grid
    type(tetra_output) :: wg
    !
    complex(dp), allocatable :: V(:,:,:,:), phase_mat(:,:,:,:,:)
    ! complex(dp), allocatable :: G0_coarse(:,:,:)
    ! complex(dp), allocatable :: G0i_coarse(:,:,:)
    complex(dp), allocatable :: den_weights(:,:), den_U(:,:,:)
    real(dp), allocatable :: den_eig(:,:)
    complex(dp), allocatable :: G0i_cluster(:,:,:)
    complex(dp), allocatable :: G_coarse(:,:,:)
    complex(dp), allocatable :: G_avg(:,:,:,:)
    complex(dp), allocatable :: Gi_coarse(:,:,:)
    complex(dp), allocatable :: self_in(:,:,:,:)
    complex(dp), allocatable :: Vqqs(:,:,:,:,:)
    complex(dp), allocatable :: G_conf(:,:,:,:), Gf_conf(:,:)
    complex(dp), allocatable :: Gf_avg(:,:), Gi_avg(:,:,:,:)
    complex(dp), allocatable :: self_out(:,:,:,:)
    complex(dp), allocatable :: delta_in(:)
    complex(dp), allocatable :: delta_out(:)
    complex(dp), allocatable :: df(:,:), dv(:,:)
    complex(dp), allocatable :: A(:,:)
    integer :: Nc, idef, ipos, i, j, k, iq, jq, iw, ndef, &
      it, NSAMPLES, MAXITER, sc_iter, MEMORY, NQ
    real(dp) :: conc, ABS_TOLERANCE, REL_TOLERANCE, ALPHA_MIX
    integer, allocatable :: pos(:), kq(:,:)
    real(dp), allocatable :: xq(:,:), R(:,:), small_xq(:,:)
    !
    NSAMPLES = 50
    MAXITER = 100
    ABS_TOLERANCE = 1e-8_dp
    REL_TOLERANCE = 1e-3_dp
    ALPHA_MIX = 0.3_dp
    MEMORY = 8
    conc = 0.05_dp  ! example concentration
    ndef = size(fc2_sc%defects,2)
    Nc = product(fc2%nq) !* size(fc2_sc%defects,2)
    NQ = grid%nqtot / Nc
    xq = grid_vec_cart(fc2%nq, S%at, divide=.true.)
    if (any(mod(grid%nq, fc2%nq) /= 0)) &
      call errore("dca_selfnrg", "grid size not multiple of fc2 grid size", 1)
    small_xq = grid_vec_cart(grid%n / fc2%nq, S%at)
    do iq = 1, size(small_xq,2)
      small_xq(:,iq) = small_xq(:,iq) / real(grid%n, dp)
    enddo
    R = grid_vec_cart(fc2%nq, S%at)
    call set_wg(S, fc2, grid, input%n_omega, wg)

    ! call tetra_init_sym_cmplx(grid, S, cmplx(wg%f, 0._dp, dp))
    !
    ! allocate(G0_coarse(S%nat3, Nc, input%n_omega))
    ! allocate(G0i_coarse(S%nat3, Nc, input%n_omega))
    allocate(G_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(G_conf(S%nat3, S%nat3, Nc, Nc))
    allocate(Gf_conf(S%nat3*Nc, S%nat3*Nc))
    allocate(Gf_avg(S%nat3*Nc, S%nat3*Nc))
    allocate(Gi_avg(S%nat3, S%nat3, Nc, Nc))
    allocate(G_coarse(S%nat3, S%nat3, Nc))
    allocate(Gi_coarse(S%nat3, S%nat3, Nc))
    allocate(G0i_cluster(S%nat3, S%nat3, Nc))
    allocate(self_in(S%nat3, S%nat3, Nc, input%n_omega))
    allocate(self_out(S%nat3, S%nat3, Nc, input%n_omega))
    allocate(V(S%nat3, S%nat3, Nc, Nc))
    allocate(A(S%nat3, S%nat3))
    allocate(phase_mat(S%nat3, S%nat3, Nc, Nc, Nc))
    allocate(kq(NQ, Nc))
    allocate(den_weights(S%nat3, NQ))
    allocate(den_U(S%nat3, S%nat3, NQ))
    allocate(den_eig(S%nat3, NQ))
    call center_V(S, S_sc, fc2_sc, Vqqs)
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
    !
    ! G0_coarse = 0.0_dp
    !
    !> construction of G0_coarse and G0i_coarse, which are averaged over the small
    !> patch in which self-energy is assumed constant
    ! do iw = 1, input%n_omega
    do iq = 1, Nc
      do jq = 1, NQ
        kq(jq,iq) = v2index_n(index2v_n(jq, grid%n / fc2%nq) + &
          index2v_n(iq, fc2%nq) * grid%n / fc2%nq, grid%n)
        ! do kq = 1, grid%nq
        !   if(norm2(small_xq(:,jq) + xq(:,iq) - grid%xq(:,kq)) < 1e-8_dp) then
        ! G0_coarse(:,iq,iw) = G0_coarse(:,iq,iw) + wg%w(:,wg%e(kq(jq,iq)),iw)
        ! endif
        ! enddo
        !     enddo
        !     G0i_coarse(:,iq,iw) = 1.0_dp / G0_coarse(:,iq,iw)
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
    delta_in = 0._dp
    !
    print*, "Starting DCA self-energy calculation..."
    do iw = 2, input%n_omega
      df = 0._dp
      dv = 0._dp
      do sc_iter = 1, MAXITER
        print*, iw, sc_iter, SUM(ABS(self_out(:,:,:,iw)-self_in(:,:,:,iw))) / &
          sum(ABS(self_in(:,:,:,iw))+1e-20_dp)
        if (all(abs(self_out(:,:,:,iw)-self_in(:,:,:,iw)) < REL_TOLERANCE * &
          abs(self_in(:,:,:,iw)) + ABS_TOLERANCE)) exit
        !> construction of G0_cluster from self-energy
        Gi_coarse = 0.0_dp
        do iq = 1, Nc
          do jq = 1, NQ
            den_U(:,:,jq) = diag(wg%f(:,wg%e(kq(jq,iq)))**2) + self_in(:,:,iq,iw)
            call mat2_diag(S%nat3, den_U(:,:,jq), den_eig(:,jq))
            ! A = diag(wg%en(iw)**2 - wg%f(:,wg%e(kq(jq,iq)))**2) - self_in(:,:,iq,iw)
            ! ! A = diag_cmplx(1._dp/wg%w(:,wg%e(kq(jq,iq)),iw)) - self_in(:,:,iq,iw)
            ! call invzmat(S%nat3, A)
            ! Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) + A
            ! G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_in(:,:,iq,iw)
          enddo
          call tetra_init(grid%n / fc2%nq, S%bg, den_eig, opt=.false.)
          den_weights = tetra_weights_green(wg%en(iw)**2)
          do jq = 1, NQ
            do i = 1, S%nat3
              Gi_coarse(:,:,iq) = Gi_coarse(:,:,iq) + &
                outer_product(den_U(:,i,jq)) * den_weights(i,jq)
            enddo
          enddo
          call invzmat(S%nat3, Gi_coarse(:,:,iq))
          G0i_cluster(:,:,iq) = Gi_coarse(:,:,iq) + self_in(:,:,iq,iw)
        enddo
        !>
        !
        !> construction of the flatten average Gf_avg over niter configurations
        Gf_avg = 0.0_dp
        do it = 1, NSAMPLES
          !> construction of V for a given configuration
          V = 0.0_dp
          do idef = 1, ndef
            call sample_binomial(Nc, conc/ndef, pos)
            do ipos = 1, size(pos)
              V = V + Vqqs(:,:,:,:,idef) * phase_mat(:,:,:,:,pos(ipos))
            enddo
          enddo
          !>
          !
          !> construction of G_conf for a given configuration
          G_conf = 0.0_dp
          do iq = 1, Nc
            G_conf(:,:,iq,iq) = G0i_cluster(:,:,iq)
            do jq = 1, Nc
              G_conf(:,:,iq,jq) = G_conf(:,:,iq,jq) - V(:,:,iq,jq)
            enddo
          enddo
          !>
          !
          Gf_conf = flatten_RR_cmplx(G_conf)
          call invzmat(S%nat3*Nc, Gf_conf)
          !
          Gf_avg = Gf_avg + Gf_conf
          !
        enddo
        Gf_avg = Gf_avg / real(NSAMPLES, dp)
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
      if (sc_iter == MAXITER + 1) &
        print*, "DCA not converged for frequency ", iw
    enddo ! frequency loop
    !
  end subroutine
!
  subroutine sample_binomial(N, c, pos)
    integer, intent(in) :: N
    real(dp), intent(in) :: c
    integer, allocatable, intent(out) :: pos(:)
    !
    integer :: i, npos
    real(dp) :: r
    integer :: pos_(N)
    !
    npos = 0
    do i = 1, N
      call random_number(r)
      if (r < c) then
        npos = npos + 1
        pos_(npos) = i
      end if
    enddo
    !
    allocate(pos(npos))
    pos = pos_(:npos)
    !
  end subroutine
!
  SUBROUTINE center_V( S, S_sc, fc2_sc, Vqqs)
    !-----------------------------------------------------------------------
    ! Apply ONE symmetry that maps atom1 -> atom2 to center potential there
    !
    USE kinds,     ONLY : DP
    USE symm_base, ONLY : irt, nsym, ft, symm_mat => s
    use symme, only : cart_to_crys, crys_to_cart
    !
    type(ph_system_info), INTENT(IN) :: S, S_sc
    type(forceconst2_sc), INTENT(IN) :: fc2_sc
    complex(dp), allocatable, intent(out) :: Vqqs(:,:,:,:,:)
    !
    real(dp):: trans(3)
    integer :: SM(3,3), SMT(3,3)
    type(forceconst2_sc) :: fc_temp
    INTEGER :: isym_map, isym, nai, naj, nbi, nbj, iR1, iR2
    integer :: site, N, iR, jR, iiR, jjR, iq, jq
    real(dp), allocatable :: tens4(:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: work(:,:,:,:,:,:)
    real(dp), allocatable :: xq(:,:)
    !
    call fc_temp%allocate(S, S_sc, fc2_sc%nq)
    !
    N = product(fc2_sc%nq)
    xq = grid_vec_cart(fc2_sc%nq, S%at, divide=.true.)
    allocate(Vqqs(S%nat3,S%nat3,N,N,size(fc2_sc%defects,2)))
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
          SM = symm_mat(:,:,isym_map)
          SMT = transpose(SM)
          trans = cryst2cart(ft(:,isym_map), S%at, 1)
          EXIT
        END IF
      END DO
      !
      work = 0.0_DP
      DO nai = 1, S%nat
        nbi = irt(isym_map, nai)
        DO naj = 1, S%nat
          nbj = irt(isym_map, naj)
          do iR = 1, N
            do jR = 1, N
              iiR = v2index(bz2simple(matmul(SM, index2v(iR, fc2_sc%nq)), fc2_sc%nq), fc2_sc%nq)
              jjR = v2index(bz2simple(matmul(SM, index2v(jR, fc2_sc%nq)), fc2_sc%nq), fc2_sc%nq)
              work(:,nai,:,naj,iiR,jjR) = work(:,nai,:,naj,iiR,jjR) + matmul( &
                matmul( SM, tens4(:,nbi,:,nbj,iR,jR) ), SMT )
            END DO
          END DO
        END DO
      END DO
      !
      call transform_tns4(work, 1)
      !
      fc_temp%fc = reshape( work, [3*S%nat, 3*S%nat, N, N] )
      !
      call translate_xR(fc_temp, trans)
      do iq = 1, size(xq,2)
        call fc_temp%r2q(xq(:,iq))
        do jq = 1, size(xq,2)
          call fc_temp%r2q(xq(:,jq), Vqqs(:,:,iq,jq,site))
        enddo
      enddo
      call translate_xR(fc_temp, -trans)
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
