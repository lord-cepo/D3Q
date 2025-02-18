module defect
  use kinds, only: dp
  use thtetra, only: tetra_init, tetra_weights_green
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
  function green_function(S, grid, tetra_weights, Us, R)
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    complex(dp), intent(in) :: tetra_weights(S%nat3, grid%nqtot)
    complex(dp), intent(in) :: Us(S%nat3, S%nat3, grid%nqtot)
    real(dp), intent(in), optional :: R(3)
    !
    complex(dp) :: green_function(S%nat3, S%nat3), phase
    integer :: iq, ibnd, iqp
    !
    green_function = 0.0_dp
    do iq = 1, grid%nq
      iqp = iq + grid%iq0
      if (present(R)) then
        phase = e_iqr(grid%xq(:,iq), R)
      else
        phase = 1.0_dp
      endif
      do ibnd = 1, S%nat3
        green_function = green_function + &
          outer_product(Us(:,ibnd,iqp)) * phase * tetra_weights(ibnd,iqp)
      enddo
    enddo
    !
    if (grid%scattered) call mpi_bsum(S%nat3, S%nat3, green_function)
  end function
  !
  subroutine main_defect(S, Sd, fc2, fc2d, grid, out_grid, input)
    type(ph_system_info), intent(in):: S, Sd
    type(forceconst2_grid), intent(in) :: fc2, fc2d
    type(q_grid), intent(in) :: grid, out_grid
    type(code_input_type), intent(in) :: input
    !
    type(forceconst2_grid) :: fc2_centered
    type(forceconst2_sc) :: fc2_sc
    complex(dp), dimension(S%nat3, S%nat3) :: VK, T_Us
    real(dp) :: freqs(S%nat3, grid%nqtot)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs, lws
    real(dp) :: lws_full_iw(S%nat3, out_grid%nqtot, 0:input%n_omega)
    real(dp) :: lws_full(S%nat3, out_grid%nqtot)
    complex(dp) :: tetra_flat(S%nat3*grid%nqtot,0:input%n_omega)
    real(dp) :: interp_flat(S%nat3*grid%nqtot)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    complex(dp) :: tetra_weights(S%nat3, grid%nqtot, 0:input%n_omega)
    complex(dp), dimension(S%nat3,S%nat3,grid%nqtot,out_grid%nqtot) :: M_K, M_M
    complex(dp), dimension(S%nat3*grid%nqtot,S%nat3*out_grid%nqtot) :: M_flat, Id_flat, MK_flat, MM_flat, V_flat
    logical, parameter :: full_born = .true.
    real(dp) :: max_freq, omega, omegaq, R_def(3), mass_def
    real(dp) :: mass_matrix(S%nat3, S%nat3)
    integer :: iq, iqp, iw, ibnd, nR, jq
    integer :: iR_def, na_def, j, comp
    ! real(dp) :: eigenvalues(S%nat3*out_grid%nqtot)
    complex(dp) :: phase_def
    !
    nR = product(fc2%nq)
    !
    !> construct mass_matrix, it will be put in the correct R position
    !> in real space, in reciprocal space the translation becomes a phase
    CALL build_mass_ratios(S, Sd, fc2%nq, mass_def, iR_def, na_def)
    mass_matrix = 0.0_dp
    do j = 1, 3
      mass_matrix(j+3*(na_def),j+3*(na_def)) = mass_def
    enddo
    !> needed to construct the phase in reciprocal space
    R_def = REAL(index2v(iR_def, fc2%nq), DP)
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
    !> VKR should have the following symmetry:
    !> VKR(na1,na2,i,j) == VKR(na2,na1,j,i) (CHECKED)

    !> centering procedure gives different output
    ! call fc2_sc%center(fc2%nq, S, Sd)
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
    call tetra_init(grid%n, S%bg, freqs**2)
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
    do iq = 1, grid%nq
      iqp = iq + out_grid%iq0
      call fc2_sc%interpolate( grid%xq(:,iq), S)
      T_Us = CONJG(TRANSPOSE(Us(:,:,iqp)))
      do jq = 1, out_grid%nqtot
        call fc2_sc%interpolate( - out_grid%sxq(:,jq), S, VK)
        phase_def = e_iqr(out_grid%sxq(:,jq) - grid%xq(:,iq), R_def)
        M_K(:,:,iqp,jq) = matmul(T_Us, &
          matmul(VK, out_Us(:,:,jq)))
        !
        M_M(:,:,iqp,jq) = phase_def * matmul(T_Us, &
          matmul(mass_matrix, out_Us(:,:,jq)))
      enddo
    enddo
    MK_flat = reshape_RR_cmplx(M_K)
    MM_flat = reshape_RR_cmplx(M_M)

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
        V_flat = MK_flat + omega**2 * MM_flat
        !
        do comp = 1, out_grid%nqtot * S%nat3
          M_flat(:,comp) = V_flat(:,comp) * tetra_flat(:,iw)
        enddo
        !
        ! if (iw /= 0) then
        !   call mat2_diag(S%nat3*out_grid%nqtot, M_flat, eigenvalues)
        !   do comp = 1, size(eigenvalues)
        !     if (abs(eigenvalues(comp)) > 1._dp) then
        !       print *, "outside", iw, eigenvalues(comp)
        !     endif
        !   enddo
        ! endif
        ! print*, "siamo a ", iw

        !> These are the only two lines that enable Full Born
        M_flat = Id_flat - M_flat
        call invzmat(S%nat3*out_grid%nqtot, M_flat)
        ! M_flat = M_flat + matmul(M_flat, M_flat)
        !
        M_flat = matmul(V_flat, M_flat)
        !
        do iqp = 1, out_grid%nqtot
          do ibnd = 1, S%nat3
            comp = ibnd + (iqp-1)*S%nat3
            lws_full_iw(ibnd,iqp,iw) = -AIMAG(M_flat(comp,comp))
          enddo
        enddo
      enddo
      ! call mpi_bsum(input%n_omega+1, S%nat3, out_grid%nqtot, lws_full_iw)
      !
      call print_message("end of born calculation")
      !
    endif
    !
    lws = 0._dp
    lws_full = 0._dp
    do iq = 1, out_grid%nq
      iqp = iq + out_grid%iq0
      do ibnd = 1, S%nat3
        comp = ibnd + (iqp-1)*S%nat3
        omegaq = out_freqs(ibnd,iqp)

        lws_full(ibnd,iqp) = interp1_scl(lws_full_iw(ibnd,iqp,:), omegaq*input%n_omega/max_freq)

        interp_flat = -AIMAG(interp1_vector(tetra_flat, omegaq*input%n_omega/max_freq))
        lws(ibnd,iqp) = SUM(ABS(MK_flat(:,comp) + MM_flat(:,comp) * omegaq**2)**2 * interp_flat)
      enddo
    enddo
    call mpi_bsum(S%nat3, out_grid%nqtot, lws)
    call mpi_bsum(S%nat3, out_grid%nqtot, lws_full)
    !
    call print_message("writing defect linewidths to file")
    !
    ! write freqs and lws to file
    if(ionode) then
      open(10, file='defect.dat', status='unknown')
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          write(10, *) out_freqs(ibnd,iq), lws(ibnd,iq), ibnd
        enddo
      enddo
      close(10)
    endif
    !
    if(ionode .and. full_born) then
      open(10, file='defect_FB.dat', status='unknown')
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          write(10, *) out_freqs(ibnd,iq), lws_full(ibnd,iq), ibnd
        enddo
      enddo
      close(10)
    endif
    !
  end subroutine
  !
end module
