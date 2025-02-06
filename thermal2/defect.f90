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
    type(ph_system_info) :: S, Sd
    type(forceconst2_grid) :: fc2, fc2d
    type(q_grid), intent(in) :: grid, out_grid
    type(code_input_type), intent(in) :: input
    !
    type(q_grid) :: grid_serial
    type(forceconst2_grid) :: fc2_centered
    type(forceconst2_sc) :: fc2_sc
    complex(dp), dimension(S%nat3, S%nat3) :: VM, Id, VK, VK1
    real(dp) :: freqs(S%nat3, grid%nqtot), omega2(S%nat3), omega2R(Sd%nat3)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs, lws, lws2, lws2R, lwsR, incoherent
    ! real(dp) :: dos(0:input%n_omega)
    real(dp), dimension(3) :: Ri, Rj, q
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    integer :: map_R(Sd%nat), kR, mjq
    complex(dp), dimension(S%nat3, S%nat3) :: T, G, V
    complex(dp) :: tetra_weights(S%nat3, grid%nqtot, 0:input%n_omega)
    complex(dp) :: tetra_interp(S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3, S%nat3, 0:input%n_omega) :: Ts, Ts2, Gs
    complex(dp), dimension(S%nat3,S%nat3,product(fc2%nq),product(fc2%nq),0:input%n_omega) :: TsR, Ts2R, GsR
    complex(dp), dimension(S%nat3,S%nat3,product(fc2%nq),product(fc2%nq)) :: TR, GR, VR, IdR, VKR, temp, VKK
    complex(dp), dimension(S%nat3, product(fc2%nq)) :: UR, UR1
    real(dp) :: SC(Sd%nat3, Sd%nat3), freqs1(S%nat3, grid%nqtot)
    logical, parameter :: full = .true.
    !
    ! complex(dp) :: tetra_interp(S%nat3, grid%nqtot)
    real(dp) :: max_freq, omega, omegaq, R(3) !, delta_q(3), q6(6), q3(3)
    !
    integer :: iq, iw, ibnd, iqp, nR, iR, jR, jbnd, jq !, atm_sc(S%nat, product(fc2%nq))
    integer :: iRin !, jRin
    ! integer :: na1, na2, j1, j2, na1_sc, na2_sc, jn1, jn2, R1, R2, Rint(3)
    ! integer :: R1, R2, i, j
    !
    nR = product(fc2%nq)
    Id = id_mat(S%nat3)
    ! IdR = id_mat(S%nat3*nR)
    !
    CALL fc2_recenter(S, fc2, fc2_centered, 2)
    !
    VKR = fc_gamma2RR(fc2%nq, S, Sd, fc2d%fc)
    CALL fc2_sc%allocate(S, fc2%nq)
    fc2_sc%fc = REAL(VKR, DP)
    !> VKR should have the following symmetry:
    !> VKR(na1,na2,i,j) == VKR(na2,na1,j,i) (CHECKED)
    call fc2_sc%center(fc2%nq, S, Sd)
    !
    ! call center_grid_sc(fc2d, fc2%nq, S, Sd)
    !
    call freq_in_grid(S, fc2_centered, grid, freqs, Us)
    !
    print*, "setting up serial inner grid"
    CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
      input%nk_in(2), input%nk_in(3),&
      grid_serial, scatter=.false., xq0=input%xk0_in)
    !
    !> max_freq is slightly larger than the maximum frequency, to be sure that
    !> maxval(out_freqs) <= max_freq
    max_freq = maxval(freqs) * 1.1_dp
    ! !
    call tetra_init(grid%n, S%bg, freqs**2)
    !
    do iw = 1, input%n_omega
      omega = max_freq*iw/input%n_omega
      tetra_weights(:,:,iw) = tetra_weights_green(omega**2)
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

    if(ionode) print*, "end of green function calculation"

    Ts = 0.0_dp
    Ts2 = 0.0_dp
    TsR = 0.0_dp
    Ts2R = 0.0_dp
    if (full) then
      do iw = 1+my_id, input%n_omega, num_procs
        ! omega = max_freq*iw/input%n_omega
        !> unit cell
        ! V = Id * 1e-6_dp * omega**2 + VK
        ! G = Gs(:,:,iw)
        ! !> first born
        ! Ts(:,:,iw) = matmul(matmul(V, G), V)
        ! !> full born
        ! T = Id - matmul(V, G)
        ! CALL invzmat(S%nat3, T)
        ! Ts2(:,:,iw) = matmul(T, V)

        !> supercell
        VR = VKR ! + IdR * 1e-6_dp * omega**2
        GR = GsR(:,:,:,:,iw)
        !> first bornz
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

        !> full born
        ! TR = IdR - matmul(VR, GR)
        ! CALL invzmat(S%nat3*nR, TR)
        ! Ts2R(:,:,iw) = matmul(TR, VR)
      enddo
    endif

    if(ionode) print*, "end of born calculation"
    ! !
    CALL freq_in_grid(S, fc2_centered, out_grid, out_freqs, out_Us)

    lws = 0.0_dp
    ! dos = 0.0_dp
    lws2 = 0.0_dp
    lws2R = 0.0_dp
    lwsR = 0.0_dp
    incoherent = 0.0_dp
    ! fc2%fc = fc2_sc%fc(:,:,:,1)
    ! S%ldrigid = .false.
    ! print*, SUM(ABS(fc2_sc%xR(:,:,1) - fc2%xR(:,:)))

    do iq = 1, out_grid%nq
      iqp = iq + out_grid%iq0
      call fc2_sc%interpolate( grid%xq(:,iqp), S,  1)
      do jq = 1, grid_serial%nq
        call fc2_sc%interpolate( - grid_serial%xq(:,jq), S, VK)
        do ibnd = 1, S%nat3
          omegaq = out_freqs(ibnd,iqp)
          if(omegaq < 1e-12) cycle
          tetra_interp = interp1_matrix(tetra_weights, omegaq*input%n_omega/max_freq)
          do jbnd = 1, S%nat3
            lws(ibnd, iq) = lws(ibnd, iq) - &
              ABS(braket(out_Us(:,ibnd,iqp), VK, Us(:,jbnd,jq)))**2 * &
              AIMAG(tetra_interp(jbnd,jq))
          enddo
        enddo
      enddo

      if (full) then
        do ibnd = 1, S%nat3
          UR = 0.0_dp
          do iR = 1, nR
            R = REAL(index2v(iR, fc2%nq), DP)
            CALL cryst_to_cart(1, R, S%at, 1)
            UR(:, iR) = out_Us(:,ibnd,iqp) * e_iqr(out_grid%xq(:,iq), R)
          enddo

          TR = interp1_tns4(TsR, omegaq*input%n_omega/max_freq)

          do iR = 1, nR
            do jR = 1, nR
              lwsR(ibnd, iq) = lwsR(ibnd, iq) - &
                AIMAG(braket(UR(:,iR), TR(:,:,iR,jR), UR(:,jR)))
            enddo
          enddo
        enddo
      endif
      ! call mpi_bsum(lws2(ibnd,iqp))
    enddo
    ! if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, lws)
    ! if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, lws2)
    ! if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, lwsR)
    ! if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, incoherent)


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
      ! open(12, file='defect2.dat', status='unknown')
      ! do iq = 1, out_grid%nqtot
      !   do ibnd = 1, S%nat3
      !     write(12, *) out_freqs(ibnd,iq), SQRT(ABS(lws2(ibnd,iq)))/8
      !   enddo
      ! enddo
      ! close(12)
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

      ! open(12, file='defect2R.dat', status='unknown')
      ! do iq = 1, out_grid%nqtot
      !   do ibnd = 1, S%nat3
      !     write(12, *) out_freqs(ibnd,iq), lws2R(ibnd,iq)
      !   enddo
      ! enddo
      ! close(12)
      !
    endif
    !
  end subroutine
  !
end module
