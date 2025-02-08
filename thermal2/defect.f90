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
    complex(dp), dimension(S%nat3, S%nat3) :: VK
    real(dp) :: freqs(S%nat3, grid%nqtot)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs, lws, lwsR
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    complex(dp) :: tetra_weights(S%nat3, grid%nqtot, 0:input%n_omega)
    complex(dp) :: tetra_interp(S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3,S%nat3,product(fc2%nq),product(fc2%nq),0:input%n_omega) :: TsR, GsR
    complex(dp), dimension(S%nat3,S%nat3,product(fc2%nq),product(fc2%nq)) :: TR, GR, temp
    real(dp), dimension(S%nat3,S%nat3,product(fc2%nq),product(fc2%nq)) :: VR, VKR
    complex(dp), dimension(S%nat3, product(fc2%nq)) :: UR
    logical, parameter :: full = .false.
    real(dp) :: max_freq, omega, omegaq, R(3), R_def(3), mass_def
    real(dp) :: mass_matrix(S%nat3, S%nat3)
    integer :: iq, iw, ibnd, nR, iR, jR, jbnd, jq, iRin
    integer :: iR_def, na_def, j
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
    VKR = fc_sc2RR(fc2%nq, S, Sd, fc2d%fc) - fc_uc2RR(S, fc2%nq, fc2%fc)
    CALL fc2_sc%allocate(S, fc2%nq)
    fc2_sc%fc = VKR
    !> VKR should have the following symmetry:
    !> VKR(na1,na2,i,j) == VKR(na2,na1,j,i) (CHECKED)

    !> centering procedure gives different output
    ! call fc2_sc%center(fc2%nq, S, Sd)
    !
    call freq_in_grid(S, fc2_centered, grid, freqs, Us)
    CALL freq_in_grid(S, fc2_centered, out_grid, out_freqs, out_Us)
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
    do iw = 1, input%n_omega
      omega = max_freq*iw/input%n_omega
      tetra_weights(:,:,iw) = tetra_weights_green(omega**2)

      !> builds big green function matrix in RR format
      if (full) then
        do iR = 1, nR
          do jR = 1, nR
            R = REAL(index2v(iR, fc2%nq) - index2v(jR, fc2%nq), DP)
            CALL cryst_to_cart(1, R, S%at, 1)
            GsR(:,:,iR,jR,iw) = &
              green_function(S, grid, tetra_weights(:,:,iw), Us, R)
          enddo
        enddo
      endif
    enddo

    call print_message("end of green function calculation")

    TsR = 0.0_dp
    if (full) then
      !> cycle is now parallel, for this reason I spliced it
      do iw = 1+my_id, input%n_omega, num_procs
        omega = max_freq*iw/input%n_omega
        VR = VKR ! VK
        VR(:,:,iR_def, iR_def) = VR(:,:,iR_def, iR_def) + &
          mass_matrix * omega**2 ! VM

        GR = GsR(:,:,:,:,iw)
        !> first born, temp is needed to if we want to precompute
        temp = 0.0_dp
        do iR = 1, nR
          do jR = 1, nR
            do iRin = 1, nR
              temp(:,:,iR,jR) = temp(:,:,iR,jR) + matmul( &
                VR(:,:,iR,iRin), GR(:,:,iRin,jR))
            enddo
          enddo
        enddo

        do iR = 1, nR
          do jR = 1, nR
            do iRin = 1, nR
              TsR(:,:,iR,jR,iw) = TsR(:,:,iR,jR,iw) + &
                matmul(temp(:,:,iR,iRin), VR(:,:,iRin,jR))
            enddo
          enddo
        enddo
        !
      enddo
    endif
    !
    call print_message("end of born calculation")
    !
    lws = 0.0_dp
    lwsR = 0.0_dp
    !> inner grid is parallel
    do iq = 1, out_grid%nqtot
      !> the last parameter says that we're interpolating
      !> on the ith component (1 or 2)
      call fc2_sc%interpolate( out_grid%xq(:,iq), S,  1)
      do jq = 1, grid%nq
        !> translation phase for VM
        phase_def = e_iqr(grid%xq(:,jq) - out_grid%xq(:,iq), R_def)
        !> second step of the interpolation (it's slower in the second
        !> component for centered grids, it's not that wise to keep it like this)
        call fc2_sc%interpolate( - grid%xq(:,jq), S, VK)
        !> outer band loop
        do ibnd = 4,4
          omegaq = out_freqs(ibnd,iq)
          if(omegaq < 1e-12) cycle
          !> 1D interpolation of the tetra weights
          tetra_interp = interp1_matrix(tetra_weights, omegaq*input%n_omega/max_freq)
          !> inner band loop
          do jbnd = 1, S%nat3
            lws(ibnd, iq) = lws(ibnd, iq) - &
            ABS(braket(out_Us(:,ibnd,iq), &
            VK + mass_matrix * phase_def * omegaq**2, & ! omegaq can give slightly different results than interp(VM(omega))
            Us(:,jbnd,jq + grid%iq0)))**2 * &
            AIMAG(tetra_interp(jbnd,jq + grid%iq0)) !> only IMG needed
          enddo
        enddo
      enddo

      if (full) then
        !> only outer cylce, inner is inside green function definition
        do ibnd = 1, S%nat3
          omegaq = out_freqs(ibnd,iq)
          if(omegaq < 1e-12) cycle
          !> builds long outer vector
          UR = 0.0_dp
          do iR = 1, nR
            R = REAL(index2v(iR, fc2%nq), DP)
            CALL cryst_to_cart(1, R, S%at, 1)
            UR(:, iR) = out_Us(:,ibnd,iq) * e_iqr(out_grid%xq(:,iq), R)
          enddo

          !> interpolation in RR form (4 indices)
          TR = interp1_tns4(TsR, omegaq*input%n_omega/max_freq)

          !> usual matrix multiplication with 4 indices
          do iR = 1, nR
            do jR = 1, nR
              lwsR(ibnd, iq) = lwsR(ibnd, iq) - &
                AIMAG(braket(UR(:,iR), TR(:,:,iR,jR), UR(:,jR)))
            enddo
          enddo
        enddo
      endif
    enddo
    call mpi_bsum(S%nat3, out_grid%nqtot, lws)
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
      !
      if (full) then
        open(12, file='defectR.dat', status='unknown')
        do iq = 1, out_grid%nqtot
          do ibnd = 1, S%nat3
            write(12, *) out_freqs(ibnd,iq), lwsR(ibnd,iq), ibnd
          enddo
        enddo
        close(12)
      endif
    endif
    !
  end subroutine
  !
end module
