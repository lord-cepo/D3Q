module defect
  use kinds, only: dp
  use thtetra, only: tetra_init, tetra_weights_green
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils, only: outer_product, freq_in_grid, interp1_matrix, &
    id_mat, braket, index2v, v2index, e_iqr
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, setup_grid
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi
  !
  implicit none
  !
contains
  function fc2_sc(S, S_sc, sc_grid, fc2)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    type(forceconst2_grid), intent(in) :: fc2
    integer :: R1, R2, nR, j1, j2, jn1, jn2, na1, na2, na1_sc, na2_sc, R(3)
    real(dp) :: fc2_sc(S_sc%nat3, S_sc%nat3)
    integer :: atoms_sc(S%nat, product(sc_grid))
    !
    nR = product(sc_grid)
    atoms_sc = map_atm_sc(S, S_sc, sc_grid)
    do R1 = 1, nR
      do R2 = 1, nR
        R = index2v(R2, sc_grid) - index2v(R1, sc_grid)
        do j1 = 1, 3
          if (R(j1) < 0) R(j1) = R(j1) + sc_grid(j1)
        enddo
        do na1 = 1, S%nat
          do na2 = 1, S%nat
            do j1 = 1, 3
              jn1 = j1 + 3*(na1-1)
              do j2 = 1, 3
                jn2 = j2 + 3*(na2-1)
                na1_sc = atoms_sc(na1, R1)
                na2_sc = atoms_sc(na2, R2)

                fc2_sc(j1 + 3*(na1_sc-1), j2 + 3*(na2_sc-1)) = &
                  fc2%FC(jn1, jn2, v2index(R, sc_grid))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end function
  !
  !> map_atm_sc(inat, iR) = inat_sc
  function map_atm_sc(S, S_sc, sc_grid)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    !
    integer :: map_atm_sc(S%nat, PRODUCT(sc_grid))
    integer :: i, isc, iR
    real(dp) :: r_cryst(3)
    map_atm_sc = -1
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        r_cryst = S_sc%tau(:,isc)*sc_grid-S%tau(:,i)
        call cryst_to_cart(1, r_cryst, S_sc%bg, -1)
        if (NORM2(r_cryst - NINT(r_cryst))<1e-2) then
          iR = v2index(NINT(r_cryst), sc_grid)
          if (map_atm_sc(i, iR) > -1) CALL errore("map_atm_sc", "found same atom in same R", -1)
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("map_atm_sc", "R is out of bound", -ABS(iR))
          map_atm_sc(i, iR) = isc
        endif
        !
      enddo
    enddo
    if (ANY(map_atm_sc == -1)) CALL errore("map_atm_sc", "some atoms are not mapped", -1)
  end function
  !
  ! function fc00(S,  fc2d, atoms)
  !   type(ph_system_info), intent(in) :: S
  !   type(forceconst2_grid), intent(in) :: fc2d
  !   integer :: atoms(:,:)
  !   !
  !   integer :: j1, j2, na1, na2, jn1, jn2, na1_sc, na2_sc
  !   real(dp) :: fc00(S%nat3, S%nat3)

  !   do na1 = 1, S%nat
  !     do na2 = 1, S%nat
  !       do j1 = 1, 3
  !         jn1 = j1 + 3*(na1-1)
  !         do j2 = 1, 3
  !           jn2 = j2 + 3*(na2-1)
  !           na1_sc = atoms(na1,1)
  !           na2_sc = atoms(na2,1)
  !           fc00(jn1, jn2) = fc2d%FC(j1 + 3*(na1_sc-1), j2 + 3*(na2_sc-1), 1)
  !           ! S%sqrtmm1(jn1) * S%sqrtmm1(jn2)
  !           ! print"(2i9,3x,2i9,2x,1pe25.15)", j1, j2, na1, na2, current_fc
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! end function
  !
  function V_mass_diag(S)
    type(ph_system_info), intent(in) :: S
    !
    real(dp) :: V_mass_diag(S%nat3, S%nat3)
    integer :: i
    !
    V_mass_diag = 0.0_dp
    do i = 1, S%nat3
      V_mass_diag(i,i) = 1.0_dp
    enddo
  end function
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
          outer_product(Us(:,ibnd,iqp)) * tetra_weights(ibnd,iqp) * phase
      enddo
    enddo
    !
    if (grid%scattered) call mpi_bsum(S%nat3, S%nat3, green_function)
  end function
  !
  subroutine main_defect(S, Sd, fc2, fc2d, grid, out_grid, input)
    use quter_defect, only: center_sc
    type(ph_system_info), intent(in) :: S, Sd
    type(forceconst2_grid) :: fc2
    type(forceconst2_grid) :: fc2d
    type(q_grid), intent(in) :: grid, out_grid
    type(code_input_type), intent(in) :: input
    !
    type(q_grid) :: grid_serial, G_grid
    type(forceconst2_grid) :: fc2_centered
    complex(dp), dimension(S%nat3, S%nat3) :: VM, Id, VK
    real(dp) :: freqs(S%nat3, grid%nqtot), omega2(S%nat3), omega2R(Sd%nat3)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs, lws, lws2, lws2R, lwsR, incoherent
    ! real(dp) :: dos(0:input%n_omega)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp) :: out_Us(S%nat3, S%nat3, out_grid%nqtot)
    integer :: map_R(Sd%nat)
    complex(dp), dimension(S%nat3, S%nat3) :: T, G, V
    complex(dp) :: tetra_weights(S%nat3, grid%nqtot, 0:input%n_omega)
    complex(dp) :: tetra_interp(S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3, S%nat3, 0:input%n_omega) :: Ts, Ts2, Gs
    complex(dp), dimension(Sd%nat3, Sd%nat3, 0:input%n_omega) :: TsR, Ts2R, GsR
    complex(dp), dimension(Sd%nat3, Sd%nat3) :: TR, GR, VR, IdR, VKR
    complex(dp) :: UR(Sd%nat3)

    ! complex(dp) :: tetra_interp(S%nat3, grid%nqtot)
    real(dp) :: max_freq, omega, omegaq, R(3), delta_q(3), q6(6), q3(3)
    !
    integer :: iq, iw, ibnd, iqp, nR, iR, jR, jbnd, jq, atm_sc(S%nat, product(fc2%nq))
    ! integer :: na1, na2, j1, j2, na1_sc, na2_sc, jn1, jn2, R1, R2, Rint(3)
    ! integer :: R1, R2, i, j
    !
    nR = product(fc2%nq)
    Id = id_mat(S%nat3)
    IdR = id_mat(S%nat3*nR)
    !
    ! call allocate_fc2_grid(1, Sd%nat, fc2_big)
    ! iR = 0
    ! do R1 = 1, nR
    !   do R2 = 1, nR
    !     iR = iR + 1
    !     fc2_big%yR(:,iR) = [index2v(R1, fc2%nq), index2v(R2, fc2%nq)]
    !   enddo
    ! enddo
    ! call cryst_to_cart(nR**2, fc2%yR, S%at, 1)
    ! fc2_big%xR = 0.0_dp
    ! fc2_big%yR = 0.0_dp

    ! fc2%FC = 0._dp
    ! fc2%FC(1,1,3) = 1.0_dp

    ! fc2_big%FC(:,:,1) = fc2_sc(S, Sd, fc2%nq, fc2)
    CALL fc2_recenter(S, fc2, fc2_centered, 2)
    VKR = fc2d%FC(:,:,1) ! - fc2_sc(S, Sd, fc2%nq, fc2)
    call center_sc(fc2d, fc2%nq, S, Sd)
    print*, fc2d%n_R
    !
    call freq_in_grid(S, fc2_centered, grid, freqs, Us)
    ! call setup_grid("simple", S%bg, fc2%nq(1), fc2%nq(2), fc2%nq(3), G_grid, scatter=.false.)
    print*, "setting up serial inner grid"
    CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
      input%nk_in(2), input%nk_in(3),&
      grid_serial, scatter=.false., xq0=input%xk0_in)
    ! !
    ! !> max_freq is slightly larger than the maximum frequency, to be sure that
    ! !> maxval(out_freqs) <= max_freq
    max_freq = maxval(freqs) * 1.1_dp
    ! map_R = (map_sc_atm(S, Sd, fc2%nq) -1 ) / S%nat + 1
    ! !
    call tetra_init(grid%n, S%bg, freqs**2)
    !
    do iw = 1, input%n_omega
      omega = max_freq*iw/input%n_omega
      tetra_weights(:,:,iw) = tetra_weights_green(omega**2)
      ! Gs(:,:,iw) = green_function(S, grid, tetra_weights(:,:,iw), Us)
      ! do iR = 1, nR
      !   do jR = 1, nR
      !     R = REAL(index2v(iR, fc2%nq) - index2v(jR, fc2%nq), DP)
      !     CALL cryst_to_cart(1, R, S%at, 1)
      !     GsR((iR-1)*S%nat3+1:iR*S%nat3,(jR-1)*S%nat3+1:jR*S%nat3,iw) = &
      !       green_function(S, grid, tetra_weights(:,:,iw), Us, R)
      !   enddo
      !   ! R = index2v(iR, fc2%nq)
      !   ! CALL cryst_to_cart(1, R, S%at, 1)
      !   ! GsR((iR-1)*S%nat3+1:iR*S%nat3, 1:S%nat3,iw) = green_function(S, grid, tetra_weights(:,:,iw), Us,  R)
      !   ! GsR(1:S%nat3, (iR-1)*S%nat3+1:iR*S%nat3,iw) = green_function(S, grid, tetra_weights(:,:,iw), Us, -R)
      !   ! if (iR == 1) cycle
      !   ! do jR = 2, nR
      !   !   GsR((iR-1)*S%nat3+1:iR*S%nat3,(jR-1)*S%nat3+1:jR*S%nat3,iw) = &
      !   !     GsR((iR-2)*S%nat3+1:(iR-1)*S%nat3,(jR-2)*S%nat3+1:(jR-1)*S%nat3,iw)
      !   ! enddo
      ! enddo
    enddo

    if(ionode) print*, "end of green function calculation"

    ! atm_sc = map_atm_sc(S, Sd, fc2%nq)
    ! VK = fc00(S, fc2d, atm_sc) - fc2%FC(:,:,fc2%i_0)
    ! call mat2_diag(Sd%nat3, VKR, omega2R)
    ! call center_sc(VKR, fc2%nq, S, Sd)

    ! VKR = fc2d%FC(:,:,1) - fc2_big%FC(:,:,1)
    Ts = 0.0_dp
    Ts2 = 0.0_dp
    TsR = 0.0_dp
    Ts2R = 0.0_dp
    ! do iw = 1+my_id, input%n_omega, num_procs
    !   omega = max_freq*iw/input%n_omega
    !   !> unit cell
    !   ! V = Id * 1e-6_dp * omega**2 + VK
    !   ! G = Gs(:,:,iw)
    !   ! !> first born
    !   ! Ts(:,:,iw) = matmul(matmul(V, G), V)
    !   ! !> full born
    !   ! T = Id - matmul(V, G)
    !   ! CALL invzmat(S%nat3, T)
    !   ! Ts2(:,:,iw) = matmul(T, V)

    !   !> supercell
    !   VR = VKR ! + IdR * 1e-6_dp * omega**2
    !   GR = GsR(:,:,iw)
    !   !> first born
    !   TsR(:,:,iw) = matmul(matmul(VR, GR), VR)
    !   !> full born
    !   ! TR = IdR - matmul(VR, GR)
    !   ! CALL invzmat(S%nat3*nR, TR)
    !   ! Ts2R(:,:,iw) = matmul(TR, VR)
    ! enddo
    ! CALL mpi_bsum(S%nat3, S%nat3, input%n_omega, Ts)
    ! CALL mpi_bsum(S%nat3, S%nat3, input%n_omega, Ts2)
    ! CALL mpi_bsum(S%nat3*nR, S%nat3*nR, input%n_omega, TsR)
    ! CALL mpi_bsum(S%nat3*nR, S%nat3*nR, input%n_omega, Ts2R)

    if(ionode) print*, "end of born calculation"
    ! !
    CALL freq_in_grid(S, fc2_centered, out_grid, out_freqs, out_Us)

    lws = 0.0_dp
    ! dos = 0.0_dp
    lws2 = 0.0_dp
    lws2R = 0.0_dp
    lwsR = 0.0_dp
    incoherent = 0.0_dp
    do iq = 1, out_grid%nq
      iqp = iq + out_grid%iq0
      ! ! do iw = 0, input%n_omega
      ! ! omega = max_freq*iw/input%n_omega
      ! ! q6 = [-out_grid%xq(:,iq), out_grid%xq(:,iq)]
      ! ! call fftinterp_mat2(q6, S, fc2d, VK)
      ! ! call mat2_diag(S%nat3, VK, omega2)
      ! ! lws2(:,iqp) =  omega2
      do ibnd = 1, S%nat3
        omegaq = out_freqs(ibnd,iqp)
        if(omegaq < 1e-12) cycle
        tetra_interp = interp1_matrix(tetra_weights, omegaq*input%n_omega/max_freq)
        !
        do jq = 1, grid_serial%nq
          if (ibnd == 1) then
            q6 = [-out_grid%xq(:,iq), grid_serial%xq(:,jq)]
            CALL fftinterp_mat2(q6, S, fc2d, VK)
          endif
          do jbnd = 1, S%nat3
            lws(ibnd, iq) = lws(ibnd, iq) - &
              ABS(braket(out_Us(:,ibnd,iqp), VK, Us(:,jbnd,jq)))**2 * &
              AIMAG(tetra_interp(jbnd,jq))
          enddo
        enddo

        ! UR = 0.0_dp
        ! do iR = 1, nR
        !   R = REAL(index2v(iR, fc2%nq), DP)
        !   CALL cryst_to_cart(1, R, S%at, 1)
        !   UR((iR-1)*S%nat3+1:iR*S%nat3) = out_Us(:,ibnd,iqp) * e_iqr(out_grid%xq(:,iq), R)
        ! enddo

        ! TR = interp1_matrix(TsR, omegaq*input%n_omega/max_freq)
        ! lwsR(ibnd,iqp) =  - AIMAG(braket(UR,TR))

      enddo
      ! call mpi_bsum(lws2(ibnd,iqp))
    enddo
    if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, lws)
    if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, lws2)
    if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, lwsR)
    if (out_grid%scattered) CALL mpi_bsum(S%nat3, out_grid%nqtot, incoherent)


    ! write freqs and lws to file
    if(ionode) then
      open(10, file='defect.dat', status='unknown')
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          write(10, *) out_freqs(ibnd,iq), lws(ibnd,iq)
        enddo
      enddo
      close(10)
      !
      ! open(12, file='defect2.dat', status='unknown')
      ! do iq = 1, out_grid%nqtot
      !   do ibnd = 1, S%nat3
      !     write(12, *) out_freqs(ibnd,iq), SQRT(lws2(ibnd,iq)/64)
      !   enddo
      ! enddo
      ! close(12)

      ! open(12, file='defectR.dat', status='unknown')
      ! do iq = 1, out_grid%nqtot
      !   do ibnd = 1, S%nat3
      !     write(12, *) out_freqs(ibnd,iq), lwsR(ibnd,iq)
      !   enddo
      ! enddo
      ! close(12)

      ! open(12, file='defect2R.dat', status='unknown')
      ! do iq = 1, out_grid%nqtot
      !   do ibnd = 1, S%nat3
      !     write(12, *) out_freqs(ibnd,iq), lws2R(ibnd,iq)
      !   enddo
      ! enddo
      ! close(12)
      !
      ! open(12, file='incoherent.dat', status='unknown')
      ! do iq = 1, out_grid%nqtot
      !   do ibnd = 1, S%nat3
      !     write(12, *) out_freqs(ibnd,iq), incoherent(ibnd,iq)
      !   enddo
      ! enddo
      ! close(12)
    endif
    !

  end subroutine
  !
end module
