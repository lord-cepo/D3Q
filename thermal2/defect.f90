module defect
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  ! only: outer_product, freq_in_grid, interp1_matrix, &
  ! id_mat, braket, index2v, v2index, e_iqr, interp1_tns4, grid_vec
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, setup_grid, q_grid_copy
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi
  use quter_defect
  use functions, only: f_gauss
  use fc3_interpolate, only: forceconst3, sparse, d3_mixed, sum_R3
  use merge_degenerate, only: merge_degen
  ! use tetra_raja
  !
  implicit none
  !
contains
  !
  subroutine main_defect(S, Sd, fc2, fc2d, grid, out_grid, input)
    type(ph_system_info), intent(in):: S, Sd
    type(forceconst2_grid), intent(in) :: fc2, fc2d
    type(q_grid), intent(in) :: grid, out_grid
    type(code_input_type), intent(in) :: input
    !
    type(q_grid) :: grid_sym
    type(sparse) :: fc3
    type(forceconst2_grid) :: fc2_centered
    type(d3_mixed) :: Dqr
    type(ph_system_info) :: S_sc, S3
    type(forceconst2_sc) :: fc2_sc
    complex(dp), dimension(S%nat3, S%nat3) :: Vqq, T_Us, temp
    real(dp) :: freqs(S%nat3, grid%nqtot)
    ! real(dp) :: weights_flat(S%nat3*grid%nqtot)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs
    complex(dp) :: lws_full_iw(S%nat3, out_grid%nqtot, 0:input%n_omega)
    complex(dp), dimension(S%nat3, out_grid%nqtot) :: lws, lws_full
    complex(dp) :: tetra_flat(S%nat3*grid%nqtot,0:input%n_omega)
    complex(dp), dimension(S%nat3*grid%nqtot) :: interp_flat, weights
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    complex(dp) :: tetra_weights(S%nat3, grid%nqtot, 0:input%n_omega)
    complex(dp) :: interp(S%nat3, grid%nqtot)
    complex(dp), allocatable, dimension(:,:,:,:) :: V
    complex(dp), allocatable, dimension(:,:) :: V_flat, Id_flat, R_flat, &
      VG_flat, VG2_flat, VG3_flat, VG4_flat, PP, PP2
    logical, parameter :: full_born = .false.
    real(dp) :: max_freq, omega, omegaq, R_def(3), mass_def
    real(dp) :: mass_matrix(S%nat3, S%nat3), mean, max_diff_fc2d, max_diff_rel
    integer :: iq, iqp, iw, ibnd, nR, jq, Nin, Nout, comp, R1, R2
    integer :: iR_def, na_def, j, compj, jbnd, ir, jr, jqp ! , jn1, jn2
    complex(dp) :: phase_def
    real(dp), allocatable :: freqs_sym(:,:)
    real(dp), dimension(3) :: xq1, xq2
    real(dp) :: fc2_transp(S%nat3, S%nat3, product(fc2%nq), product(fc2%nq))
    real(dp) :: mean_fc2d(Sd%nat3, Sd%nat3)
    integer, allocatable :: yR2_out(:,:), yR3_out(:,:,:), index2(:), index3(:,:)
    ! !> elphbolt tetra
    ! INTEGER, ALLOCATABLE :: tetra(:,:), tetracount(:), tetramap(:,:,:)
    ! REAL(DP), ALLOCATABLE :: evals(:,:), tetraevals(:,:,:)
    !
    nR   = product(fc2%nq)
    Nin  = S%nat3*grid%nqtot
    Nout = S%nat3*out_grid%nqtot
    ! allocate(V(S%nat3,S%nat3,grid%nqtot,out_grid%nqtot))
    ! allocate(V_flat(Nin,Nout))
    allocate(PP(S%nat3,grid%nqtot), PP2(S%nat3,grid%nqtot))
    ! allocate(V_flat, source=V_flat)
    if(full_born) allocate(R_flat, VG_flat, VG2_flat, VG3_flat, VG4_flat, source=V_flat)
    !
    !> construct mass_matrix, it will be put in the correct R position
    !> in real space, in reciprocal space the translation becomes a phase
    ! CALL build_mass_ratios(S, Sd, fc2%nq, mass_def, iR_def, na_def)
    ! mass_matrix = 0.0_dp
    ! do j = 1, 3
    !   mass_matrix(j+3*(na_def),j+3*(na_def)) = mass_def / nR
    ! enddo
    !> needed to construct the phase in reciprocal space
    ! R_def = REAL(index2v(iR_def, fc2%nq), DP)
    CALL cryst_to_cart(1, R_def, S%at, 1)
    !
    !> input files are periodic, it can be changed
    CALL fc2_recenter(S, fc2, fc2_centered, 2)
    call fc3%read(input%file_mat3, S3, .true.)
    !
    max_diff_fc2d = 0._dp
    do ibnd = 1, Sd%nat3
      do jbnd = ibnd+1, Sd%nat3
        mean = (fc2d%fc(ibnd,jbnd,1) + fc2d%fc(jbnd,ibnd,1)) / 2.0_dp
        if (ABS(fc2d%fc(ibnd,jbnd,1) - fc2d%fc(jbnd,ibnd,1)) > max_diff_fc2d) then
          max_diff_fc2d = ABS(fc2d%fc(ibnd,jbnd,1) - fc2d%fc(jbnd,ibnd,1))
          max_diff_rel = max_diff_fc2d / ABS(mean)
        endif
        mean_fc2d(ibnd,jbnd) = mean
        mean_fc2d(jbnd,ibnd) = mean
      enddo
    enddo
    print*, "max_diff_fc2d", max_diff_fc2d
    print*, "max_diff_rel", max_diff_rel
    !> SC: grid type   (big_nat3,big_nat3,1)
    !> UC: grid type   (nat3,nat3,nR)
    !> RR: new SC type (nat3,nat3,nR,nR)
    ! print*, fc2_sc%fc(1,3,2,2), fc2%fc(1,3,1), fc2_sc%xR1(:,2)

    ! do R1 = 1, 2
    !   do R2 = 1, nR
    !     do ibnd = 1, 6
    !       do jbnd = 1, 6
    !         if (fc2_sc%fc(ibnd, jbnd, R2, R1) - fc2_sc%fc(jbnd, ibnd, R1, R2) > 1e-10_dp) then
    !           print*, "R1", R1, "R2", R2, "ibnd", ibnd, "jbnd", jbnd
    !           print*, fc2_sc%fc(ibnd, jbnd, R2, R1), fc2_sc%fc(jbnd, ibnd, R1, R2)
    !         endif
    !       enddo
    !     enddo
    !   enddo
    ! enddo
    ! print*, ABS(fc2_transp - fc2_sc%fc) < 1e-10_dp
    ! do ir = 1, nR
    !   do jr = 1, nR
    !     if (ir > 2 .or. jr > 2) fc2_sc%fc(:,:,ir,jr) = 0.0_dp
    !     ! do ibnd = 1, S%nat3
    !     !   do jbnd = 1, S%nat3
    !     !     if(ibnd /= jbnd) fc2_sc%fc(ibnd,jbnd,ir,jr) = 0.0_dp
    !     !   enddo
    !     ! enddo
    !   enddo
    ! enddo

    !> VKR should have the following symmetry:
    !> VKR(na1,na2,i,j) == VKR(na2,na1,j,i) (CHECKED)

    !> centering procedure gives different output
    call S_uc2sc(S, fc2%nq, S_sc)
    S_sc%ityp(1) = 2
    CALL fc2_sc%allocate(S, fc2%nq)
    fc2_sc%fc = fc_sc2RR(fc2%nq, S, Sd, fc2d%fc) - &
      fc_uc2RR(fc2%nq, S, fc2%fc)
    call center2(fc2_sc, fc2%nq, S, Sd)
    ! call matd2RR(fc3, S, fc2_sc)
    ! call div_mass_fcsc(S, Sd, fc2_sc)
    !
    call freq_in_grid_degen(S, fc2_centered, fc2_Sc, out_grid, out_freqs, out_Us)
    do iq = 1, out_grid%nq
      call merge_degen(S%nat3, out_freqs(:,iq), out_freqs(:,iq))
    enddo
    call freq_in_grid_degen(S, fc2_centered, fc2_sc, grid, freqs, Us)
    do iq = 1, grid%nq
      call merge_degen(S%nat3, freqs(:,iq), freqs(:,iq))
    enddo
    !
    !> max_freq is slightly larger than the maximum frequency, to be sure that
    !> maxval(out_freqs) <= max_freq
    max_freq = maxval(freqs) * 1.1_dp
    !
    !> tetra have already the square of freq, it can be changed with
    !> the usual delta formula after some benchmarking
    call q_grid_copy(grid, grid_sym)
    ! if (grid%symmetrized) then
    call grid_sym%symmetrize(S)
    allocate(freqs_sym(S%nat3, grid_sym%nqtot))
    call freq_in_grid(S, fc2_centered, grid_sym, freqs_sym)
    call tetra_init_sym(grid_sym, S, freqs_sym**2, .true.)
    ! else
    !   call tetra_init(grid%n, S%bg, freqs**2, .true.)
    ! endif
    !
    call print_message("end of tetra initialization")
    !
    ! ALLOCATE(tetra(6*grid%nqtot, 4), tetracount(grid%nqtot), tetramap(2, grid%nqtot, 24))
    ! ALLOCATE(evals(grid%nqtot,S%nat3), tetraevals(6*grid%nqtot,S%nat3, 4))
    ! evals = TRANSPOSE(freqs)
    ! CALL form_tetrahedra_3d(grid%nqtot, grid%n, tetra, tetracount, tetramap)
    ! CALL fill_tetrahedra_3d(grid%nqtot, S%nat3, tetra, freqs**2, tetraevals)
    ! !
    !> serially calculate tetras for an equally spaced omega
    !> interval, after we will interpolate them (even at more than 1st order).
    !> The cycle is serial cause the tetra_weights_green is already parallelized
    do iw = 0, input%n_omega
      omega = max_freq*iw/input%n_omega
      tetra_weights(:,:,iw) = tetra_weights_green(omega**2)
      tetra_flat(:,iw) = reshape(tetra_weights(:,:,iw), [S%nat3*grid%nqtot])
    enddo
    !
    call print_message("end of tetra weights calculation")
    !
    ! V = 0.0_dp
    ! do iq = 1, grid%nq
    !   iqp = iq + grid%iq0
    !   call fc2_sc%interpolate( grid%xq(:,iq), S)
    !   T_Us = CONJG(TRANSPOSE(Us(:,:,iqp)))
    !   do jq = 1, out_grid%nqtot
    !     call fc2_sc%interpolate( out_grid%xq(:,jq), S, Vqq)
    !     ! call zgemm('N', 'N', S%nat3, S%nat3, S%nat3, 1._dp, Vqq, &
    !     !   S%nat3, out_Us(:,:,jq), S%nat3, 0._dp, temp, S%nat3)
    !     ! call zgemm('C', 'N', S%nat3, S%nat3, S%nat3, 1._dp, Us(:,:,iqp), &
    !     !   S%nat3, temp, S%nat3, 0._dp, V(:,:,iqp,jq), S%nat3)
    !     V(:,:,iqp,jq) = matmul(T_Us, matmul(Vqq, out_Us(:,:,jq)))
    !   enddo
    ! enddo
    ! call mpi_bsum(S%nat3, S%nat3, grid%nqtot, out_grid%nqtot, V)
    ! !
    ! V_flat = reshape_RR_cmplx(V)
    !
    ! call print_message("end of braket calculation")
    ! !
    ! lws_full = 0._dp
    ! lws_full_iw = 0._dp
    ! if(full_born) then
    !   if (out_grid%nqtot /= grid%nqtot) &
    !     call errore("defect", "out_grid%nqtot /= grid%nqtot", 1)
    !   !
    !   Id_flat = id_mat(S%nat3*out_grid%nqtot)
    !   do iw = my_id, input%n_omega, num_procs
    !     omega = max_freq*iw/input%n_omega
    !     !
    !     ! V_flat = V_flat
    !     do comp = 1, out_grid%nqtot * S%nat3
    !       VG_flat(:,comp) = V_flat(:,comp) * tetra_flat(:,iw)
    !     enddo

    !     !
    !     !> These are the only two lines that enable Full Born
    !     ! M_flat = Id_flat - M_flat
    !     ! call invzmat(S%nat3*out_grid%nqtot, M_flat)
    !     ! M_flat = M_flat + matmul(M_flat, M_flat)
    !     ! M_flat = matmul(V_flat, M_flat)
    !     !
    !     call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, VG_flat, &
    !       Nin, VG_flat, Nin, 0._dp, VG2_flat, Nin)
    !     call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, VG2_flat, &
    !       Nin, V_flat, Nin, 0._dp, VG3_flat, Nin)
    !     call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, VG3_flat, &
    !       Nin, V_flat, Nin, 0._dp, VG4_flat, Nin)
    !     R_flat = VG_flat + VG2_flat + VG3_flat + VG4_flat
    !     call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, V_flat, &
    !       Nin, R_flat, Nin, 0._dp, VG_flat, Nin)
    !     !
    !     do iqp = 1, out_grid%nqtot
    !       do ibnd = 1, S%nat3
    !         comp = ibnd + (iqp-1)*S%nat3
    !         lws_full_iw(ibnd,iqp,iw) = VG_flat(comp,comp)
    !       enddo
    !     enddo
    !   enddo
    !   call mpi_bsum(S%nat3, out_grid%nqtot, input%n_omega+1, lws_full_iw)
    !   !
    !   call print_message("end of born calculation")
    !   !
    ! endif
    ! !
    lws = 0._dp
    lws_full = 0._dp

    ibnd = 1
    iq = 10
    call fc2_sc%interpolate(out_grid%xq(:,iq), S)
    call fc2_sc%interpolate(grid%xq(:,iq), S, Vqq)
    print*, braket(out_Us(:,ibnd,iq), Vqq, Us(:,ibnd,iq))
    call fc2_sc%interpolate(grid%xq(:,iq), S)
    call fc2_sc%interpolate(out_grid%xq(:,iq), S, Vqq)
    print*, braket(Us(:,ibnd,iq), Vqq, out_Us(:,ibnd,iq))



    ! weights = reshape(spread(grid%w, 1, S%nat3), [S%nat3*grid%nqtot])
    do iq = 1, out_grid%nq
      iqp = iq !+ out_grid%iq0
      ! call fc3%sum_R2(out_grid%xq(:,iq), S%nat3, Dqr, .true.)
      call fc2_sc%interpolate(out_grid%xq(:,iq), S)
      do ibnd = 1, S%nat3
        ! comp = ibnd + (iqp-1)*S%nat3
        omegaq = out_freqs(ibnd,iqp)
        ! lws_full(ibnd,iqp) = interp1_scl(lws_full_iw(ibnd,iqp,:), omegaq*input%n_omega/max_freq) / omegaq
        PP = 0._dp
        do jq = 1, grid%nq
          jqp = jq + grid%iq0
          call fc2_sc%interpolate(grid%xq(:,jq), S, vqq)
          ! call sum_R3(S, grid%xq(:,jq), Dqr, Vqq)
          ! call interp_at_once(fc2_sc, out_grid%xq(:,iq), grid%xq(:,jq), S, Vqq)
          do jbnd = 1, S%nat3
            if(near(norm2(out_grid%xq(:,iq)-grid%xq(:,jq))) .and. near(omegaq, freqs(jbnd,jqp))) cycle
            PP(jbnd,jqp) = braket(out_Us(:,ibnd,iq), Vqq, Us(:,jbnd,jqp)) !* f_gauss(omegaq-freqs(ibnd,jqp), 1e-4_dp)
            ! PP2(jbnd,jqp) = braket(out_Us(:,ibnd,iq), temp, Us(:,jbnd,jqp))
          enddo
        enddo
        ! do jq = 1, grid%nq
        !   do jbnd = 1, S%nat3
        !     interp_flat(jbnd + (jq-1)*S%nat3) = &
        !       delta_fn_tetra(out_freqs(ibnd,iqp)**2, jq, jbnd, grid%n, tetramap, tetracount, tetraevals)
        !   enddo
        ! enddo
        ! lws_full(ibnd,iqp) = SUM(ABS(V_flat(:,comp))**2 * interp_flat) / omegaq
        interp = interp1_matrix(tetra_weights, omegaq*input%n_omega/max_freq)
        ! interp_flat = interp1_vector(tetra_flat, omegaq*input%n_omega/max_freq)
        lws(ibnd,iqp) = SUM( ABS(PP)**2 * interp) / omegaq /grid%nqtot
        ! lws_full(ibnd,iqp) = SUM(ABS(PP2)**2 * interp) / omegaq /grid%nqtot
      enddo
    enddo
    call mpi_bsum(S%nat3, out_grid%nqtot, lws)
    where (abs(lws) < 1e-20_dp)
      lws = 0._dp
    endwhere
    ! call mpi_bsum(S%nat3, out_grid%nqtot, lws_full)
    !
    call print_message("writing defect linewidths to file")
    !
    ! write freqs and lws to file
    if(ionode) then
      open(10, file='defect.dat', status='unknown')
      do iq = 1, out_grid%nqtot
        ! do ibnd = 1, S%nat3
          write(10, "(6e14.5)") AIMAG(lws(:,iq))
        ! enddo
      enddo
      close(10)
    endif
    !
    ! if(ionode) then
    !   open(10, file='defect_symm.dat', status='unknown')
    !   do iq = 1, out_grid%nqtot
    !     do ibnd = 1, S%nat3
    !       write(10, "(3e14.5,I3)") out_freqs(ibnd,iq), lws_full(ibnd,iq), ibnd
    !     enddo
    !   enddo
    !   close(10)
    ! endif
    !
    if (full_born) deallocate(R_flat, VG_flat, VG2_flat, VG3_flat, VG4_flat)
    ! deallocate(V, V_flat)
  end subroutine
  !
end module
