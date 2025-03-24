module defect
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  ! only: outer_product, freq_in_grid, interp1_matrix, &
  ! id_mat, braket, index2v, v2index, e_iqr, interp1_tns4, grid_vec
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, setup_grid
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi
  use quter_defect
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
    type(forceconst2_grid) :: fc2_centered
    type(ph_system_info) :: S_sc
    type(forceconst2_sc) :: fc2_sc
    complex(dp), dimension(S%nat3, S%nat3) :: Vqq, T_Us
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
    complex(dp), allocatable, dimension(:,:,:,:) :: V
    complex(dp), allocatable, dimension(:,:) :: V_flat, Id_flat, R_flat, &
      VG_flat, VG2_flat, VG3_flat, VG4_flat
    logical, parameter :: full_born = .false.
    real(dp) :: max_freq, omega, omegaq, R_def(3), mass_def
    real(dp) :: mass_matrix(S%nat3, S%nat3)
    integer :: iq, iqp, iw, ibnd, nR, jq, Nin, Nout, comp
    integer :: iR_def, na_def, j, compj, jbnd, ir, jr
    complex(dp) :: phase_def
    !
    nR   = product(fc2%nq)
    Nin  = S%nat3*grid%nqtot
    Nout = S%nat3*out_grid%nqtot
    allocate(V(S%nat3,S%nat3,grid%nqtot,out_grid%nqtot))
    allocate(V_flat(Nin,Nout))
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
    !
    !> SC: grid type   (big_nat3,big_nat3,1)
    !> UC: grid type   (nat3,nat3,nR)
    !> RR: new SC type (nat3,nat3,nR,nR)
    CALL fc2_sc%allocate(S, fc2%nq)
    fc2_sc%fc = fc_sc2RR(fc2%nq, S, Sd, fc2d%fc(:,:,1)) - &
      fc_uc2RR(fc2%nq, S, fc2%fc)

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
    call fc2_sc%center(fc2%nq, S, S_sc)
    !
    call freq_in_grid(S, fc2_centered, grid, freqs, Us)
    call freq_in_grid(S, fc2_centered, out_grid, out_freqs, out_Us)
    !
    !> max_freq is slightly larger than the maximum frequency, to be sure that
    !> maxval(out_freqs) <= max_freq
    max_freq = maxval(freqs) * 1.1_dp
    !
    !> tetra have already the square of freq, it can be changed with
    !> the usual delta formula after some benchmarking
    if (grid%symmetrized) then
      call tetra_init_sym(grid, S, freqs**2)
    else
      call tetra_init(grid%n, S%bg, freqs**2)
    endif
    !
    call print_message("end of tetra initialization")
    !
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
    V = 0.0_dp
    do iq = 1, grid%nq
      iqp = iq + grid%iq0
      call fc2_sc%interpolate( grid%xq(:,iq), S)
      T_Us = CONJG(TRANSPOSE(Us(:,:,iqp)))
      do jq = 1, out_grid%nqtot
        call fc2_sc%interpolate( out_grid%xq(:,jq), S, Vqq)
        V(:,:,iqp,jq) = matmul(T_Us, &
          matmul(Vqq, out_Us(:,:,jq)))
      enddo
    enddo
    call mpi_bsum(S%nat3, S%nat3, grid%nqtot, out_grid%nqtot, V)
    WRITE(0,*) sum(V)
    !
    V_flat = reshape_RR_cmplx(V)
    !
    call print_message("end of braket calculation")
    !
    lws_full = 0._dp
    lws_full_iw = 0._dp
    if(full_born) then
      if (out_grid%nqtot /= grid%nqtot) &
        call errore("defect", "out_grid%nqtot /= grid%nqtot", 1)
      !
      Id_flat = id_mat(S%nat3*out_grid%nqtot)
      do iw = my_id, input%n_omega, num_procs
        omega = max_freq*iw/input%n_omega
        !
        ! V_flat = V_flat
        do comp = 1, out_grid%nqtot * S%nat3
          VG_flat(:,comp) = V_flat(:,comp) * tetra_flat(:,iw)
        enddo

        !
        !> These are the only two lines that enable Full Born
        ! M_flat = Id_flat - M_flat
        ! call invzmat(S%nat3*out_grid%nqtot, M_flat)
        ! M_flat = M_flat + matmul(M_flat, M_flat)
        ! M_flat = matmul(V_flat, M_flat)
        !
        call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, VG_flat, &
          Nin, VG_flat, Nin, 0._dp, VG2_flat, Nin)
        call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, VG2_flat, &
          Nin, V_flat, Nin, 0._dp, VG3_flat, Nin)
        call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, VG3_flat, &
          Nin, V_flat, Nin, 0._dp, VG4_flat, Nin)
        R_flat = VG_flat + VG2_flat + VG3_flat + VG4_flat
        call zgemm('N', 'N', Nin, Nin, Nin, 1._dp, V_flat, &
          Nin, R_flat, Nin, 0._dp, VG_flat, Nin)
        !
        do iqp = 1, out_grid%nqtot
          do ibnd = 1, S%nat3
            comp = ibnd + (iqp-1)*S%nat3
            lws_full_iw(ibnd,iqp,iw) = VG_flat(comp,comp)
          enddo
        enddo
      enddo
      call mpi_bsum(S%nat3, out_grid%nqtot, input%n_omega+1, lws_full_iw)
      !
      call print_message("end of born calculation")
      !
    endif
    !
    lws = 0._dp
    lws_full = 0._dp

    ! weights = reshape(spread(grid%w, 1, S%nat3), [S%nat3*grid%nqtot])
    do iq = 1, out_grid%nq
      iqp = iq ! + out_grid%iq0
      do ibnd = 1, S%nat3
        comp = ibnd + (iqp-1)*S%nat3
        omegaq = out_freqs(ibnd,iqp)

        lws_full(ibnd,iqp) = interp1_scl(lws_full_iw(ibnd,iqp,:), omegaq*input%n_omega/max_freq) / omegaq

        interp_flat = interp1_vector(tetra_flat, omegaq*input%n_omega/max_freq)
        lws(ibnd,iqp) = SUM(ABS(V_flat(:,comp))**2 * interp_flat ) / omegaq /grid%nqtot
      enddo
    enddo
    ! call mpi_bsum(S%nat3, out_grid%nqtot, lws)
    ! call mpi_bsum(S%nat3, out_grid%nqtot, lws_full)
    !
    call print_message("writing defect linewidths to file")
    !
    ! write freqs and lws to file
    if(ionode) then
      open(10, file='defect.dat', status='unknown')
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          write(10, "(3e14.5,I3)") out_freqs(ibnd,iq), lws(ibnd,iq), ibnd
        enddo
      enddo
      close(10)
    endif
    !
    if(ionode .and. full_born) then
      open(10, file='defect_FB.dat', status='unknown')
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          write(10, "(3e14.5,I3)") out_freqs(ibnd,iq), lws_full(ibnd,iq), ibnd
        enddo
      enddo
      close(10)
    endif
    !
    if (full_born) deallocate(R_flat, VG_flat, VG2_flat, VG3_flat, VG4_flat)
    deallocate(V, V_flat)
  end subroutine
  !
end module
