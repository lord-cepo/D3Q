module defect
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  ! only: outer_product, freq_in_grid, interp1_matrix, &
  ! id_mat, braket, index2v, v2index, e_iqr, interp1_tns4, grid_vec
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi
  use quter_defect
  use functions, only: f_gauss
  use fc3_interpolate, only: forceconst3, sparse, d3_mixed, sum_R3
  use merge_degenerate, only: merge_degen
  use test_print, only: allclose
  use simtet, only: tetra_init_sym_cmplx, tetra_weights_green_cmplx
  ! use tetra_raja
  !
  implicit none
  !
contains
  !
  subroutine tetra_init_grid_sym(grid, S, fc2, wg, grid_sym_, U_sym)
    type(q_grid), intent(in) :: grid
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(tetra_output), intent(out) :: wg
    type(q_grid), intent(out), optional :: grid_sym_
    complex(dp), allocatable, intent(out), optional :: U_sym(:,:,:)
    !
    type(q_grid) :: grid_sym
    real(dp), allocatable :: freqs_sym(:,:)
    !
    call q_grid_copy(grid, grid_sym)
    if( .not. grid_sym%symmetrized) &
      call grid_sym%symmetrize(S)
    !
    allocate(U_sym(S%nat3, S%nat3, grid_sym%nqtot))
    allocate(freqs_sym(S%nat3, grid_sym%nqtot))
    call freq_in_grid(S, fc2, grid_sym, freqs_sym, U_sym)
    call tetra_init_sym(grid_sym, S, freqs_sym**2, .false., wg)
    if (present(grid_sym_)) &
      call q_grid_copy(grid_sym, grid_sym_)
  end subroutine
  !
  ! subroutine frobenius_triangle(fc2_sc, S, filename)
  !   type(forceconst2_sc), intent(in) :: fc2_sc
  !   type(ph_system_info), intent(in) :: S
  !   character(len=*), intent(in) :: filename
  !   !
  !   integer :: ir1, ir2, na1, na2
  !   real(dp), dimension(3) :: R1, R2
  !   real(dp) :: perimeter
  !   !
  !   open(unit=10, file=filename, status="replace")
  !   do ir1 = 1, fc2_sc%n_R1
  !     do na1 = 1, S%nat
  !       R1 = fc2_sc%xR1(:,ir1) + S%tau(:,na1)
  !       do ir2 = 1, fc2_sc%n_R2(ir1)
  !         do na2 = 1, S%nat
  !           R2 = fc2_sc%xR2(:,ir2,ir1) + S%tau(:,na2)
  !           perimeter = norm2(R2 - R1) + &
  !             norm2(R1 - fc2_sc%taudef) + &
  !             norm2(R2 - fc2_sc%taudef)
  !           write(10,*) perimeter, fc2_sc%fc((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,ir2,ir1)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   close(10)

  ! end subroutine
  !
  subroutine full_V(S, fc2, out_grid, fc2_sc, V, out_Us, grid, Us)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(q_grid), target, intent(in) :: out_grid
    type(forceconst2_sc), intent(inout) :: fc2_sc
    complex(dp), allocatable, intent(out) :: V(:,:,:,:)
    complex(dp), intent(out), target, optional :: out_us(:,:,:)
    type(q_grid), intent(in), target, optional :: grid
    complex(dp), intent(in), target, optional :: Us(:,:,:)
    !
    integer :: iq, jq
    complex(dp), dimension(S%nat3, S%nat3) :: D, D1
    type(q_grid), pointer :: grid_
    complex(dp), pointer :: us_(:,:,:)
    !
    if (present(grid)) then
      grid_ => grid
      if(present(out_Us)) us_ => us
    else
      grid_ => out_grid
      if (present(out_Us)) us_ => out_Us
    endif
    allocate(V(S%nat3, S%nat3, grid_%nqtot, out_grid%nqtot))
    !
    do iq = 1, grid_%nqtot
      call fc2_sc%r2q(grid_%xq(:,iq))
      do jq = 1, out_grid%nq
        call fc2_sc%r2q(out_grid%xq(:,jq), D)
        ! if (ALL(ABS(out_grid%xq(:,jq) - grid_%xq(:,iq)) < 1e-10_dp)) then
        !   call fftinterp_mat2(out_grid%xq(:,jq), S, fc2, D1)
        !   D = D - D1
        ! endif
        if (present(out_Us)) &
          D = matmul(conjg(transpose(Us_(:,:,iq))), matmul(D, out_Us(:,:,jq)))
        V(:,:,iq,jq) = D
      enddo
    enddo
    !
    call mpi_bsum(S%nat3, S%nat3, out_grid%nqtot, grid_%nqtot, V)
  end subroutine
  !
  subroutine green_inversion(S, fc2, fc2d, fc2_sc, out_grid)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2, fc2d
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(q_grid), intent(in) :: out_grid
    ! type(code_input_type), intent(in) :: input
    !
    type(q_grid) :: grid_sym
    integer :: nR, iq, jq, ibnd, k, iR, ieta, nq, i, j, jbnd, compi, compj
    integer, parameter :: SMEARINGS = 5
    real(dp), allocatable :: Eoo(:), E0(:,:), R(:,:), lws(:,:,:), E0_sym(:,:)
    complex(dp), allocatable :: U0(:,:,:), Uoo(:,:), braks(:,:,:), &
      G_mat(:,:), S_mat(:,:), V(:,:,:,:), V_sum(:,:)
    complex(dp) :: G(SMEARINGS), brak, D(S%nat3, S%nat3)
    character(20) :: filename
    character(20) :: file_energies
    logical :: file_exists
    real(dp) :: eta(SMEARINGS)
    !
    nR = product(fc2%nq)
    nq = out_grid%nqtot
    allocate(Eoo(S%nat3*nq), E0(S%nat3, nq))
    allocate(R(3, nR))
    allocate(lws(SMEARINGS, S%nat3, nq))
    allocate(U0(S%nat3, S%nat3, nq))
    allocate(Uoo(S%nat3*nq, S%nat3*nq))
    allocate(braks(S%nat3*nR, S%nat3, nq))
    allocate(G_mat(S%nat3*nq, S%nat3*nq))
    allocate(S_mat(S%nat3*nq, S%nat3*nq))
    allocate(V(S%nat3, S%nat3, nq, nq))
    allocate(V_sum(nq, nq))
    !
    call freq_in_grid(S, fc2, out_grid, E0, U0)
    E0 = E0**2
    print*, "mean E0",  sum(E0)  / real(S%nat3*nq, dp)
    do iq = 1, nq
      print"(6E15.3)", E0(:,iq)
    enddo
    write(file_energies, '(A2,I2.2,A4)') "E-", out_grid%n(1), '.dat'
    inquire(file=file_energies, exist=file_exists)

    file_exists = .false.
    print*, file_exists
    if (file_exists) then
      open(10, file=file_energies, status='old')
      ! do j = 1, S%nat3*nq
      !   do i = 1, S%nat3*nq
      !     read(10, "(2F15.6)") Uoo(i,j)
      !   enddo
      ! enddo
      !
      do j = 1, S%nat3*nq
        read(10, "(2E15.5)") Eoo(j)
      enddo
      close(10)
      call print_message("read energies from file")
    else
      ! do iq = 1, out_grid%nq
      !   call fc2_sc%r2q(out_grid%xq(:,iq))
      !   call fc2_sc%r2q(out_grid%xq(:,iq), D)
      !   call mat2_diag(S%nat3, D, E0(:,iq))
      !   print"(6E15.3)", E0(:,iq)
      ! enddo
      !
      call full_V(S, fc2, out_grid, fc2_sc, V, U0)
      do iq = 1, nq
        ! ! call mat2_diag(S%nat3, V(:,:,iq,iq), E0(:,iq))
        call fc2_sc%r2q(out_grid%xq(:,iq))
        ! do ibnd = 1, S%nat3
        !   E0(ibnd,iq) = braket(U0(:,ibnd,iq), D)
        ! enddo
        ! print"(6E15.3)", E0(:,iq)
        do jq = 1, nq
          call fc2_sc%r2q(out_grid%xq(:,jq), D)
          V_sum(jq,iq) = sum(D)
        enddo
      enddo
      Uoo = flatten_RR_cmplx(V)
      !
      call print_message("end of interpolation")
      ! Uoo = fc2d%fc(:,:,1)
      call mat2_diag(S%nat3*nq, Uoo, Eoo)
      call print_message("end of frequency calculation")
      !
      open(10, file=file_energies, status='replace', action='write')
      ! do j = 1, S%nat3*nq
      !   do i = 1, S%nat3*nq
      !     write(10, "(2F15.6)") Uoo(i,j)
      !   enddo
      ! enddo
      !
      do j = 1, S%nat3*nq
        write(10, "(E15.5)") Eoo(j)
      enddo
      close(10)
    endif

    print*, "mean Eoo", sum(Eoo) / real(S%nat3*nq, dp)
    ! call q_grid_copy(out_grid, grid_sym)
    ! call grid_sym%symmetrize(S)
    ! allocate(E0_sym(S%nat3, grid_sym%nqtot))
    ! call freq_in_grid(S, fc2, grid_sym, E0_sym)
    ! E0_sym = E0_sym**2
    ! call tetra_init_sym(grid_sym, S, E0_sym, .false.)

    eta = [1e-4_dp, 1e-6_dp, 1e-8_dp, 1e-10_dp, 1e-12_dp]
    ! do iR = 1, nR
    !   R(:,iR) = REAL(index2v(iR, fc2%nq), dp)
    ! enddo
    ! call cryst_to_cart(nR, R, S%at, 1)
    !
    ! do iq = 1, nq
    !   do ibnd = 1, S%nat3
    !     compi = ibnd + (iq-1)*S%nat3
    !     do k = 1, S%nat3 * nq
    !       S_mat(k,compi) = dot_product(Uoo((iq-1)*S%nat3+1:iq*S%nat3,k), U0(:,ibnd,iq))
    !     enddo
    !   enddo
    ! enddo
    ! !
    ! call print_message("end of S matrix calculation")
    ! G_mat = 0._dp
    ! do iq = 1, nq
    !   do ibnd = 1, S%nat3
    !     compi = ibnd + (iq-1)*S%nat3
    !     do jq = 1, nq
    !       do jbnd = 1, S%nat3
    !         compj = jbnd + (jq-1)*S%nat3
    !         G = 0._dp
    !         do k = 1, S%nat3 * nq
    !           ! brak = 0._dp
    !           ! do iR = 1, nR
    !           !   brak = brak + &
    !           !     dot_product(U0(:,ibnd,iq), Uoo((iR-1)*S%nat3+1:iR*S%nat3, k)) * &
    !           !     e_iqr(out_grid%xq(:,iq), R(:,iR))
    !           ! enddo
    !           ! brak = dot_product(U0(:,ibnd,iq), Uoo((iq-1)*S%nat3+1:iq*S%nat3,k))
    !           ! do ieta = 1, SMEARINGS
    !           !   G(ieta) = G(ieta) + abs(brak)**2 / cmplx(1.14e-5_dp-Eoo(k), eta(ieta), dp)
    !           ! enddo
    !           G_mat(compi, compj) = G_mat(compi, compj) + &
    !             conjg(S_mat(k, compi)) * S_mat(k, compj) / cmplx(1.14e-5_dp-Eoo(k), eta(3), dp)
    !         enddo
    !       enddo
    !     enddo
    !   enddo
    ! enddo

    ! ! write(filename, '(A2,I2.2,A4)') "G-", out_grid%n(1), '.dat'
    ! ! open(10, file=filename, status='unknown')
    ! ! do iq = 1, nq
    ! !   do ibnd = 1, S%nat3
    ! !     write(10, "(e14.5,5e14.5,I3,3F7.2)") E0(ibnd,iq), lws(:,ibnd,iq), ibnd, out_grid%xq(:,iq)
    ! !   enddo
    ! ! enddo
    ! ! close(10)
    ! write(filename, '(A2,I2.2,A4)') "M-", out_grid%n(1), '.dat'
    ! open(10, file=filename, status='unknown')
    ! do i = 1, S%nat3*nq
    !   write(10, "(10000E14.5)") G_mat(:,i)
    ! enddo
    ! close(10)

    ! write(filename, '(A3,I2.2,A4)') "E0-", out_grid%n(1), '.dat'
    ! open(10, file=filename, status='unknown')
    ! do iq = 1, nq
    !   do ibnd = 1, S%nat3
    !     write(10, "(E15.5,I2,3F7.2)") E0(ibnd,iq), ibnd, out_grid%xq(:,iq)
    !   enddo
    ! enddo
    ! close(10)

    write(filename, '(A2,I2.2,A4)') "V-", out_grid%n(1), '.dat'
    open(10, file=filename, status='unknown')
    do i = 1, nq
      write(10, "(10000E14.5)") V_sum(:,i)
    enddo
    close(10)
  end subroutine
  !
  ! subroutine check_derivative_swap(fc, nat3)
  !   type(forceconst2_sc), intent(in) :: fc
  !   integer, intent(in) :: nat3
  !   !
  !   real(dp) :: mean, max_diff_fc2d, max_diff_rel
  !   integer :: ibnd, jbnd, ir1, ir2, ir1p, ir2p
  !   integer(dp) :: R(3)
  !   logical :: ok
  !   !
  !   max_diff_fc2d = 0._dp
  !   max_diff_rel = 0._dp
  !   do ir1 = 1, fc%n_R1
  !     do ir2 = 1, fc%n_R2(ir1)
  !       R = fc%yR2(:,ir2,ir1)
  !       ok = .false.
  !       do ir1p = 1, fc%n_R1
  !         if (ALL(fc%yR1(:,ir1p) == R)) then
  !           ok = .true.
  !           exit
  !         endif
  !       enddo
  !       if (.not. ok) then
  !         do ir1p = 1, fc%n_R1
  !           print"(3I2)", fc%yR1(:,ir1p)
  !         enddo
  !         print*, R
  !         call errore("check_derivative_swap", "R1 not found in R2", 1)
  !       endif
  !       R = fc%yR1(:,ir1)
  !       ok = .false.
  !       do ir2p = 1, fc%n_R2(ir1p)
  !         if (ALL(fc%yR2(:,ir2p,ir1p) == R)) then
  !           ok = .true.
  !           exit
  !         endif
  !       enddo
  !       if (.not. ok) &
  !         call errore("check_derivative_swap", "R2 not found in R1", 1)
  !       do ibnd = 1, nat3
  !         do jbnd = ibnd+1, nat3
  !           mean = (fc%fc(ibnd,jbnd,ir2,ir1) + fc%fc(jbnd,ibnd,ir2p,ir1p)) / 2.0_dp
  !           if (ABS(fc%fc(ibnd,jbnd,ir2,ir1) - fc%fc(jbnd,ibnd,ir2p,ir1p)) > max_diff_fc2d) then
  !             max_diff_fc2d = ABS(fc%fc(ibnd,jbnd,ir2,ir1) - fc%fc(jbnd,ibnd,ir2p,ir1p))
  !             max_diff_rel = max_diff_fc2d / ABS(mean)
  !           endif
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   if (ionode) print"(A,E14.4)", "max_diff_fc2d", max_diff_fc2d
  !   if (ionode) print"(A,E14.4)", "max_diff_rel", max_diff_rel
  ! end subroutine
  !
  function G0(sc_grid, omega2, wg, S,  grid, U)
    integer, intent(in) :: sc_grid(3)
    real(dp) , intent(in) :: omega2
    type(tetra_output), intent(in) :: wg
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    complex(dp), intent(in) :: U(S%nat3,S%nat3,wg%nsym)
    !<grid%S%nat3>
    complex(dp) :: G0(S%nat3*product(sc_grid), S%nat3*product(sc_grid))
    complex(dp) :: den(S%nat3, wg%nsym)
    real(dp) :: dummy(S%nat3)
    integer :: Ri, Rj, iq, ibnd, nR
    real(dp) :: R_grid(3, product(sc_grid))
    !
    den = tetra_weights_green(omega2)
    G0 = 0._dp
    nR = product(sc_grid)
    R_grid = grid_vec_cart(sc_grid, S%at)
    !
    do Ri = 1, nR
      do Rj = 1, nR
        do iq = 1, wg%nsym
          do ibnd = 1, S%nat3
            G0(S%nat3*(Ri-1)+1:S%nat3*Ri, S%nat3*(Rj-1)+1:S%nat3*Rj) =  &
              G0(S%nat3*(Ri-1)+1:S%nat3*Ri, S%nat3*(Rj-1)+1:S%nat3*Rj) + &
              outer_product(U(:,ibnd,iq)) * &
              e_iqr(grid%xq(:,iq), R_grid(:,Ri) - R_grid(:,Rj)) * &
              den(ibnd,iq) * wg%qw(iq)
          enddo
        enddo
      enddo
    enddo
    !
  end function
  !
  subroutine main_defect(S, fc2, fc2_sc, grid, out_grid, input)
    use constants, only: BOHR_RADIUS_CM, RY_TO_CMM1
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(inout) :: fc2_sc
    type(forceconst2_grid), intent(in) :: fc2
    type(q_grid), intent(in) :: grid, out_grid
    type(code_input_type), intent(in) :: input
    type(tetra_output) :: w_in, w_out
    !
    ! type(forceconst2_sc) :: TR(0:input%n_omega)
    complex(dp), dimension(S%nat3, S%nat3) :: Vqq, T_Us, Vqq2, green
    real(dp), dimension(S%nat3, grid%nqtot) :: freqs, freqs1
    ! real(dp) :: weights_flat(S%nat3*grid%nqtot)
    real(dp), dimension(S%nat3, out_grid%nqtot) :: out_freqs, out_freqs1
    complex(dp), dimension(S%nat3) :: e2

    complex(dp), dimension(S%nat3, out_grid%nqtot) :: lws_out
    complex(dp), dimension(S%nat3, out_grid%nqtot) :: lws_full

    complex(dp), dimension(S%nat3, grid%nqtot) :: lws

    ! complex(dp) :: tetra_flat(S%nat3*grid%nqtot,0:input%n_omega)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    complex(dp), dimension(S%nat3, S%nat3, out_grid%nqtot) :: out_Us!, out_us_copy
    complex(dp), allocatable :: interp(:,:)
    complex(dp), allocatable, dimension(:,:,:,:) :: V
    complex(dp), allocatable, dimension(:,:) :: V_flat, Id_flat, R_flat, &
      GV_flat, GV2_flat, GV3_flat, GV4_flat, poly_T, inv_T
    logical, parameter :: full_born = .false.
    real(dp) :: max_freq, omega, omegaq
    complex(dp) :: lw2
    integer :: iq, iqp, iw, ibnd, nR, jq, Nin, Nout, comp, compj, i
    integer ::  jbnd, ibnd1, dim_deg ! , jn1, jn2
    ! real(dp), allocatable :: freqs_sym(:,:)
    complex(dp), allocatable, dimension(:,:,:) :: Vqqs, Vqqs_out
    character(15) :: filename
    ! real(dp) :: R_grid(3,out_grid%nqtot)
    real(dp) :: V0_cm3, mult
    real(dp) :: concentrations(3)
    real(dp) :: spectral_function(3,0:input%n_omega)
    integer :: iconc, iq_star
    real(dp) :: TETRA_MULTIPLIER
    real(dp) :: cq(3)
    complex(dp) :: self_energy(S%nat3, 0:input%n_omega)
    complex(dp) :: pixel(grid%nqtot, out_grid%nqtot)

    type(q_grid):: grid_sym
    complex(dp), allocatable :: Us_sym(:,:,:)
    !
    !> full born quantities in real space
    ! complex(dp), allocatable :: VG(:,:,:,:)
    !
    concentrations = [0.0_dp, 1e-4_dp, 1e-3_dp] ! defects per cell
    V0_cm3 = S%omega * (BOHR_RADIUS_CM)**3
    print*, "V0_cm3", V0_cm3
    nR   = product(fc2%nq)
    Nin  = S%nat3*grid%nqtot
    Nout = S%nat3*out_grid%nqtot
    allocate(Vqqs(S%nat3,S%nat3,grid%nqtot))
    allocate(Vqqs_out(S%nat3,S%nat3,out_grid%nqtot))
    allocate(V(S%nat3,S%nat3,out_grid%nqtot,grid%nqtot))
    allocate(V_flat(Nin,Nout))
    !
    !> input files are periodic, it can be changed
    !> SC: grid type   (big_nat3,big_nat3,1)
    !> UC: grid type   (nat3,nat3,nR)
    !> RR: new SC type (nat3,nat3,nR,nR)

    ! call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)
    ! call freq_in_grid(S, fc2, grid, freqs, Us)
    !
    call freq_in_grid_degen(S, fc2, fc2_sc, out_grid, out_freqs, out_Us, out_freqs1)
    call freq_in_grid_degen(S, fc2, fc2_sc, grid, freqs, Us, freqs1)

    !
    !> max_freq is slightly larger than the maximum frequency, to be sure that
    !> maxval(out_freqs) <= max_freq
    max_freq = maxval(freqs) * 1.1_dp
    TETRA_MULTIPLIER = 1.0_dp / max_freq**2
    !
    !> tetra have already the square of freq, it can be changed with
    !> the usual delta formula after some benchmarking
    call tetra_init_grid_sym(grid, S, fc2, w_in, grid_sym, Us_sym)
    !
    call print_message("end of tetra initialization")
    !
    !> serially calculate tetras for an equally spaced omega
    !> interval, after we will interpolate them (even at more than 1st order).
    !> The cycle is serial cause the tetra_weights_green is already parallelized
    allocate(w_in%w(S%nat3, w_in%nsym, 0:input%n_omega))
    do iw = 0, input%n_omega
      omega = max_freq*iw/input%n_omega
      w_in%w(:,:,iw) = tetra_weights_green(omega**2)
      do iq = 1, grid%nqtot
        do ibnd = 1, S%nat3
          if(isnan(ABS(w_in%w(ibnd,w_in%e(iq),iw)))) w_in%w(ibnd,w_in%e(iq),iw) = 0._dp
        enddo
      enddo
      ! tetra_flat(:,iw) = reshape(tetra_weights(:,:,iw), [S%nat3*grid%nqtot])
      do iq = 1, grid%nqtot
        call merge_degen(S%nat3, w_in%w(:,w_in%e(iq),iw), freqs(:,iq))
      enddo
    enddo
    ! do iq = 1, out_grid%nqtot
    !   matmul(Sd%fc, matmul(G0(input%sc_grid, omegaq**2, w_in, S, grid_sym, Us_sym). Sd%fc))
    ! !
    ! if(full_born) then
    !   call deallocate_tetra()
    !   call tetra_init_grid_sym(out_grid, S, fc2, w_out)
    !   allocate(w_out%w(S%nat3, w_out%nsym, 0:input%n_omega))
    !   do iw = 0, input%n_omega
    !     omega = max_freq*iw/input%n_omega
    !     w_out%w(:,:,iw) = tetra_weights_green(omega**2)
    !     do iq = 1, out_grid%nqtot
    !       do ibnd = 1, S%nat3
    !         if(isnan(ABS(w_out%w(ibnd,w_out%e(iq),iw)))) w_out%w(ibnd,w_out%e(iq),iw) = 0._dp
    !       enddo
    !     enddo
    !     do iq = 1, out_grid%nqtot
    !       call merge_degen(S%nat3, w_out%w(:,w_out%e(iq),iw), out_freqs(:,iq))
    !     enddo
    !   enddo
    ! endif
    !
    ! open(10, file="phdos.dat", status='replace', action='write')
    ! do iw = 0, input%n_omega
    !   omega = max_freq*iw/input%n_omega
    !   write(10, "(E17.4,100E17.4)") omega, sum(matmul(w_in%w(:,:,iw), w_in%qw))
    ! enddo
    ! call print_message("end of tetra weights calculation")
    !
    call full_V(S, fc2, out_grid, fc2_sc, V, out_us, grid, Us)
    V_flat = flatten_RR_cmplx(V)
    !
    do iq = 1, out_grid%nq
      call fc2_sc%r2q(out_grid%xq(:,iq))
      call fc2_sc%r2q(out_grid%xq(:,iq), Vqqs_out(:,:,iq))
    enddo
    !
    open(1236, file="self-q6.dat", status='replace', action='write')
    open(1237, file="self-q7.dat", status='replace', action='write')
    open(1238, file="self-q8.dat", status='replace', action='write')

    allocate(interp(S%nat3, grid%nqtot))
    do iw = 0, input%n_omega
      lws_out = 0._dp
      omega = max_freq*iw/input%n_omega
      if (iw > 1 .and. trim(input%calculation) == 'lw') exit
      if(mod(iw, input%n_omega/10) == 0) &
        print"(A,A,I3,A)", input%calculation, " progress ", NINT(100 * REAL(iw,DP) / input%n_omega), "%"
      do iq = 1, out_grid%nq
        call fc2_sc%r2q(out_grid%xq(:,iq))
        do jq = 1, grid%nqtot
          call fc2_sc%r2q(grid%xq(:,jq), Vqqs(:,:,jq))
          ! if (norm2(grid%xq(:,jq) - out_grid%xq(:,iq)) < 1e-10_dp) then
          !   call fftinterp_mat2(out_grid%xq(:,iq), S, fc2, Vqq)
          !   Vqqs(:,:,jq) = Vqqs(:,:,jq) - Vqq
          ! endif
          pixel(jq,iq) = sum(Vqqs(:,:,jq))
        enddo
        !
        do ibnd = 1, S%nat3
          lws_out(ibnd,iq) = lws_out(ibnd,iq) + REAL(braket(out_us(:,ibnd,iq), Vqqs_out(:,:,iq)), dp)
          ! comp = (iq-1)*S%nat3 + ibnd
          omegaq = out_freqs(ibnd,iq)
          ! interp = tetra_weights_green(omegaq**2)
          if (trim(input%calculation) == 'lw') then
            interp = interp1_matrix(w_in%w, omegaq*input%n_omega/max_freq) / omegaq
          else
            interp = w_in%w(:,:,iw)
          endif

          do jq = 1, grid%nqtot
            Vqq = Vqqs(:,:,jq)
            ! if(allocated(fc2_sc%inclusion_eig)) then
            !   do concurrent(i=1:size(fc2_sc%inclusion_eig))!, ABS(omegaq**2 - fc2_sc%inclusion_eig(i)) > 1e-10)
            !     Vqq = Vqq + outer_product2(fc2_sc%inclusion_Dnx_out(:,i,iq) / &
            !       (omegaq**2 - fc2_sc%inclusion_eig(i)), fc2_sc%inclusion_Dnx_in(:,i,jq))
            !   enddo
            ! endif
            lw2 = 0._dp
            do jbnd = 1, S%nat3
              if(norm2(grid%xq(:,jq) - out_grid%xq(:,iq)) < 1e-10_dp .and. &
                jbnd == ibnd) cycle
              lw2 = lw2 + ABS(braket(out_us(:,ibnd,iq), Vqq, us(:,jbnd,jq)))**2 * &
                interp(jbnd,w_in%e(jq))
            enddo
            lws_out(ibnd,iq) = lws_out(ibnd,iq) + lw2
            if (ibnd == 6 .and. iw == 2481) write(1236, "(4E20.8)") cryst2cart(grid%xq(:,jq), S%at, -1), aimag(lw2)
            if (ibnd == 6 .and. iw == 2568) write(1237, "(4E20.8)") cryst2cart(grid%xq(:,jq), S%at, -1), aimag(lw2)
            if (ibnd == 8 .and. iw == 2526) write(1238, "(4E20.8)") cryst2cart(grid%xq(:,jq), S%at, -1), aimag(lw2)

          enddo
        enddo
        call merge_degen(S%nat3, lws_out(:,iq), out_freqs(:,iq))
      enddo
      close(1236)
      close(1237)
      close(1238)
      !
      if(trim(input%calculation) == 'spfdef') then
        do iconc = 1, size(concentrations)
          mult = concentrations(iconc)
          call tetra_init_sym_cmplx(out_grid, S, TETRA_MULTIPLIER * (out_freqs**2 + mult*lws_out))
          spectral_function(iconc, iw) = &
            sum(matmul(AIMAG(tetra_weights_green_cmplx(TETRA_MULTIPLIER * omega**2)), out_grid%w)) * &
            TETRA_MULTIPLIER
          ! do iq = 1, out_grid%nqtot
          !   do ibnd = 1, S%nat3
          !     write(11, "(3e14.5,I3,3F7.2,E13.4)") out_freqs(ibnd,iq), lws_out(ibnd,iq), ibnd, out_grid%xq(:,iq), out_grid%w(iq)
          !   enddo
          ! enddo
          ! close(11)
        enddo
      elseif((trim(input%calculation) == 'self') ) then
        if (out_grid%nqtot /= 1) &
          call errore("main_defect", "you can calculate self-energy only in one q-point at once", 1)
        self_energy(:,iw) = lws_out(:,1)
      endif
      call mpi_bsum(S%nat3, out_grid%nqtot, lws_out)
    enddo
    ! !
    ! write(filename, "(A,I1,A)") "V-real", fc2%nq(1), "p.dat"
    ! open(10, file=filename, status='replace', action='write')
    ! write(filename, "(A,I1,A)") "V-imag", fc2%nq(1), "p.dat"
    ! open(11, file=filename, status='replace', action='write')
    ! do iq = 1, out_grid%nq
    !   write(10, "(10000E15.5)") REAL(pixel(:,iq),DP)
    !   write(11, "(10000E15.5)") AIMAG(pixel(:,iq))
    ! enddo
    ! close(10)
    ! close(11)




    ! lws = 0._dp
    ! lws_out = 0._dp
    ! do iq = 1, out_grid%nq
    !   iqp = iq !+ out_grid%iq0
    !   call fc2_sc%r2q(out_grid%xq(:,iq))
    !   do jq = 1, grid%nq
    !     call fc2_sc%r2q(grid%xq(:,jq), Vqqs(:,:,jq))
    !     if (ALL(ABS(grid%xq(:,jq) - out_grid%xq(:,iq)) < 1e-10_dp)) then
    !       call fftinterp_mat2(out_grid%xq(:,jq), S, fc2, Vqq)
    !       Vqqs(:,:,jq) = Vqqs(:,:,jq) - Vqq
    !     endif
    !   enddo
    !   ibnd = 1
    !   do !ibnd
    !     do ibnd1 = ibnd+1, S%nat3
    !       if(.not. near(out_freqs(ibnd,iq) + out_freqs1(ibnd,iq), &
    !         out_freqs(ibnd1,iq) + out_freqs1(ibnd1,iq))) exit
    !     enddo
    !     dim_deg = ibnd1-ibnd
    !     ! print*, "deg?", iq, ibnd, dim_deg
    !     omegaq = out_freqs(ibnd,iqp)
    !     interp = interp1_matrix(w_in%w, omegaq*input%n_omega/max_freq)
    !     ! interp = tetra_weights(ibnd,iq,iw)
    !     Vqq2 = 0._dp
    !     do jq = 1, grid%nq
    !       green = 0._dp
    !       do jbnd = 1, S%nat3
    !         if(near(norm2(out_grid%xq(:,iq)-grid%xq(:,jq))) .and. near(omegaq, freqs(jbnd,jq))) cycle
    !         ! if(near(norm2(out_grid%xq(:,iq)+grid%xq(:,jq))) .and. near(omegaq, freqs(jbnd,jq))) cycle
    !         SELECT CASE(input%delta_approx)
    !          CASE('gauss')
    !           green = green + outer_product(Us(:,jbnd,jq)) * &
    !             f_gauss(omegaq**2-freqs(jbnd,jq)**2, 1e-10_dp)
    !          CASE('tetra')
    !           green = green + outer_product(Us(:,jbnd,jq)) * &
    !             interp(jbnd,w_in%e(jq))
    !         END SELECT
    !       enddo
    !       Vqq2 = Vqq2 + matmul(Vqqs(:,:,jq),matmul(green, conjg(transpose(Vqqs(:,:,jq)))))
    !     enddo
    !     if (dim_deg > 1) then
    !       call diag_degen_cmplx(S%nat3, dim_deg, Vqq2, &
    !         out_Us(:,ibnd:ibnd+dim_deg-1,iqp), e2(ibnd:ibnd+dim_deg-1))
    !       lws_out(ibnd:ibnd+dim_deg-1,iqp) = e2(ibnd:ibnd+dim_deg-1) / omegaq
    !     else
    !       lws_out(ibnd,iqp) = braket(out_us(:,ibnd,iqp),Vqq2) / omegaq
    !     endif
    !     ibnd = ibnd + dim_deg
    !     if(ibnd > S%nat3) exit
    !   enddo
    ! enddo
    ! call mpi_bsum(S%nat3, out_grid%nqtot, lws_out)
    ! call write_file(out_freqs, lws_out, 'degen-1B.dat', out_grid%type)
    ! call print_message("end of degen calculation")
    !
    ! out_us = us


    ! do iw = my_id, input%n_omega, num_procs
    !   lws_out = 0._dp
    !   if(mod(iw, input%n_omega/10) == 0) &
    !     print"(A,A,I3,A)", input%calculation, " progress ", 100 * iw / input%n_omega, "%"
    !   omega = max_freq*iw/input%n_omega
    !   ! write(filename, '(A5,I3.3,A4)') "self-", iw, '.dat'
    !   ! open(11, file=filename, status='unknown')
    !   do iq = 1, out_grid%nq
    !     iqp = iq + out_grid%iq0
    !     do ibnd = 1, S%nat3
    !       ! lws_out(ibnd,iq) = lws_out(ibnd,iq) + braket(out_us(:,ibnd,iqp), Vqqs_out(:,:,iq))
    !       ! omegaq = out_freqs(ibnd,iqp)
    !       do jq = 1, grid%nqtot
    !         do jbnd = 1, S%nat3
    !           lws_out(ibnd,iq) = lws_out(ibnd,iq) + &
    !             ABS(V(jbnd,ibnd,jq,iq))**2 * w_in%w(jbnd,w_in%e(jq),iw)
    !         enddo
    !       enddo
    !       ! if (ABS(aimag(lws_out(ibnd,iq))) < 1e-15_dp) lws_out(ibnd,iq) = CMPLX(REAL(lws_out(ibnd,iq), dp), 0._dp, dp)
    !     enddo
    !   enddo
    !   if(out_grid%nqtot > 1) then
    !     do iconc = 1, size(concentrations)
    !       mult = concentrations(iconc)
    !       call tetra_init_sym_cmplx(out_grid, S, TETRA_MULTIPLIER * (out_freqs**2 + mult*lws_out))
    !       spectral_function(iconc, iw) = &
    !         sum(matmul(AIMAG(tetra_weights_green_cmplx(TETRA_MULTIPLIER * omega**2)), out_grid%w)) * &
    !         TETRA_MULTIPLIER
    !       ! do iq = 1, out_grid%nqtot
    !       !   do ibnd = 1, S%nat3
    !       !     write(11, "(3e14.5,I3,3F7.2,E13.4)") out_freqs(ibnd,iq), lws_out(ibnd,iq), ibnd, out_grid%xq(:,iq), out_grid%w(iq)
    !       !   enddo
    !       ! enddo
    !       ! close(11)
    !     enddo
    !   else
    !     self_energy(:,iw) = lws_out(:,1)
    !   endif
    ! enddo
    !
    select case(trim(input%calculation))
     case("lw")
      if(trim(out_grid%type) == 'path') then
        write(filename, '(A4,I2.2,A4)') 'path', grid%n(1), '.dat'
      else
        write(filename, '(I2.2,I2.2,A4)') out_grid%n(1), grid%n(1), '.dat'
      endif
      call write_file(out_freqs, lws_out, filename, out_grid%type)
      call print_message("end of 1B calculation")
     case("spfdef")
      open(10, file="spectral_function.dat", status='replace', action='write')
      do iw = 0, input%n_omega
        omega = max_freq*iw/input%n_omega
        write(10, "(100E17.4)") omega, spectral_function(:,iw)
      enddo
     case("self")
      open(10, file="self-energy-def.dat", status='replace', action='write')
      write(10, "(100E25.8)") out_freqs(:,1)
      do iw = 0, input%n_omega
        omega = max_freq*iw/input%n_omega
        write(10, "(E17.4,100E17.4)") omega, self_energy(:,iw)
      enddo
      close(10)
    end select
    !
    if(full_born) then
      deallocate(V, V_flat)
      allocate(V(S%nat3, S%nat3, out_grid%nqtot, out_grid%nqtot))
      allocate(V_flat(Nout, Nout))
      allocate(GV_flat(Nout,Nout))
      ! allocate(GV2_flat, GV3_flat, GV4_flat, source=GV_flat)
      allocate(inv_T, R_flat, Id_flat, source=GV_flat)
      !
      call full_V(S, fc2, out_grid, fc2_sc, V, out_us)
      V_flat = flatten_RR_cmplx(V)
      !
      Id_flat = id_mat(S%nat3*out_grid%nqtot)
      !
      deallocate(interp)
      allocate(interp(S%nat3, out_grid%nqtot))
      do iq = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          comp = (iq-1)*S%nat3 + ibnd
          omegaq = out_freqs(ibnd,iq)
          interp = interp1_matrix(w_out%w, omegaq*input%n_omega/max_freq)
          do jq = 1, out_grid%nqtot
            do jbnd = 1, S%nat3
              compj = (jq-1)*S%nat3 + jbnd
              GV_flat(compj,comp) = V_flat(compj,comp) * &
                interp(jbnd,w_out%e(jq))
            enddo
          enddo
        enddo
      enddo
      ! print*, sum(GV_flat), sum(V_flat)

      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, GV_flat, &
      !   Nout, GV_flat, Nout, 0._dp, R_flat, Nout)
      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, GV2_flat, &
      !   Nout, GV_flat, Nout, 1._dp, R_flat, Nout)
      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, GV3_flat, &
      !   Nout, GV_flat, Nout, 1._dp, R_flat, Nout)
      ! ! R_flat = GV_flat + GV2_flat + GV3_flat + GV4_flat

      ! call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, V_flat, &
      !   Nout, R_flat, Nout, 0._dp, poly_T, Nout)
      !
      ! > These are the only two lines that enable Full Born
      R_flat = Id_flat - GV_flat
      call invzmat(S%nat3*out_grid%nqtot, R_flat)
      call zgemm('N', 'N', Nout, Nout, Nout, 1._dp, V_flat, &
        Nout, R_flat, Nout, 0._dp, inv_T, Nout)
      !
      do iqp = 1, out_grid%nqtot
        do ibnd = 1, S%nat3
          comp = ibnd + (iqp-1)*S%nat3
          omegaq = out_freqs(ibnd,iqp)
          ! lws_out(ibnd,iqp) = poly_T(comp,comp) / omegaq
          lws_full(ibnd,iqp) = inv_T(comp,comp) / omegaq
        enddo
      enddo
      call print_message("end of born calculation")
      !
    endif
    !
    if(ionode) then

      if(full_born) then
        call write_file(out_freqs, lws_full, 'FB.dat', out_grid%type)
        call write_file(out_freqs, lws_out, '4B.dat', out_grid%type)
      endif
    endif
    !
    !
    if (full_born) deallocate(V_flat, GV_flat)
    if (full_born) deallocate(w_out%w)
    ! deallocate(V, V_flat)
  contains
    subroutine write_file(freqs_, lws_, filename_, type)
      real(dp), intent(in) :: freqs_(:,:)
      complex(dp), intent(in) :: lws_(:,:)
      character(*), intent(in) :: filename_
      character(*), intent(in) :: type
      !
      integer :: iq_, ibnd_
      !
      open(10, file=filename_, status='replace', action='write')
      if (trim(type) == 'path') then
        do iq_ = 1, size(freqs_, 2)
          write(10, "(1000e14.5)") freqs_(:,iq_), -AIMAG(lws_(:,iq_))
        enddo
      else
        do iq_ = 1, size(freqs_, 2)
          do ibnd_ = 1, size(freqs_, 1)
            write(10, "(3E15.5,I2)") freqs_(ibnd_,iq_), lws_(ibnd_,iq_), ibnd_
          enddo
        enddo
      endif
      close(10)
    end subroutine
  end subroutine
  !
end module
