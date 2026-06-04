program defectp
  use defect
  USE fc3_interpolate,  ONLY : forceconst3
  use q_grids, only: setup_grid, revert_grid
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi, num_procs
  use input_fc, only: read_fc2, aux_system, div_mass_fc2, write_fc2, multiply_mass_fc2
  use asr2_module, only: impose_asr2
  use thutils, only: v2index, cryst2cart, index2v_cart
  use quter_defect
  use quter_module, only : quter
  use test_print
  use defect_eig, only : full_diag, norm_gV, norm_V
  use constants, only : RY_TO_CMM1
  use full_born, only : full_born_p, full_born_analytical
  USE parameters, ONLY : ntypx
  use defect_proj, only: project
  use dca, only : dca_selfnrg
  ! use symm_base, only : nofrac
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd, S_, S_sc
  type(forceconst2_grid) :: fc2d, fc2d_centered, fc2_centered, fc2_periodic, fc2_treated, fc2_treated_centered
  type(forceconst2_sc) :: fc2_sc, fc2_sc_centered
  type(q_grid) :: in_grid, out_grid, sym_grid, out_grid_sym_scat, in_grid_sym_scat
  type(code_input_type) :: input, input_
  class(forceconst3), pointer :: fc3, fc3_
  integer :: n_add, new_it, isc
  integer, allocatable :: map(:,:), new_ityp(:), new_ityp0(:)
  real(dp), allocatable :: new_tau(:,:), DRR(:,:,:,:), new_tau0(:,:), new_fc0(:,:,:)
  real(dp), allocatable :: D_nx_real(:,:,:)
  integer, allocatable :: inclusion_isc(:)
  complex(dp), allocatable :: inclusion_U(:,:)
  real(dp) :: c
  real(dp), dimension(3) :: d1, d2
  integer :: sc_grid(3), i, R1, R2, na1, na2, j1, j2, nR, j, iq, isc1, isc2
  real(dp), allocatable :: p(:)
  real(dp), allocatable :: freqs_sc(:)
  real(dp) :: v(3)
  real(dp), allocatable :: interp_grid(:,:)
  real(dp), allocatable :: diff_large(:,:)
  complex(dp), allocatable :: Ds(:,:,:), matq(:,:,:,:,:)
  character(len=100) :: filename
  real(dp) :: max_norm, max_freq
  real(dp), allocatable :: E(:)
  real(dp) :: eta(12)
  integer :: iw, jR, jq, iR
  complex(dp) :: mat3(3,3), eig3(3)
  complex(dp), allocatable :: V_sc(:,:)
  real(dp) :: tau(3)
  real(dp), allocatable :: R(:,:), q(:,:)
  real(dp), allocatable :: D0(:,:)
  ! integer :: wait_for_debugger
  ! integer, allocatable :: atoms(:,:)
  ! integer :: na1, na2, j1, j2, na1_sc, na2_sc, jn1, jn2, R1, R2, nR
  ! integer :: iq, i, j
  ! real(dp), allocatable :: dyn(:,:,:,:), zeu(:,:,:)
  !
  CALL start_mpi()
  !
  CALL READ_INPUT("DEF", input, out_grid, S, fc2_periodic)
  S%lrigid = .false.
  CALL fc2_recenter(S, fc2_periodic, fc2_centered, 2)
  S%lrigid = .true.
  !
  !
  if(all(input%sc_grid == -1)) then
    sc_grid = fc2_periodic%nq
  else
    sc_grid = input%sc_grid
    if(ionode) print*, "Using sc_grid = ", input%sc_grid
  end if
  nR = product(sc_grid)
  call allocate_fc2_grid(nR, S%nat, fc2_treated)
  allocate(interp_grid(3,nR))
  allocate(Ds(S%nat3, S%nat3, nR))
  allocate(matq(3,3,S%nat,S%nat,nR))
  interp_grid = grid_vec_cart(sc_grid, S%bg)
  do iq = 1, nR
    interp_grid(:,iq) = interp_grid(:,iq) / sc_grid
    call fftinterp_mat2(interp_grid(:,iq), S, fc2_centered, Ds(:,:,iq))
    do concurrent( j1=1:3, j2=1:3, na1=1:S%nat, na2=1:S%nat)
      matq(j1,j2,na1,na2,iq) = Ds(j1 + 3*(na1-1), j2 + 3*(na2-1), iq)
    enddo
  enddo
  CALL quter(sc_grid(1), sc_grid(2), sc_grid(3), S%nat, S%tau, S%at, S%bg, matq, interp_grid, fc2_treated, 0)

  CALL read_fc2(input%file_mat3, Sd, fc2d)
  ! call S_uc2sc(S, Sd, sc_grid, S_sc)
  CALL aux_system(Sd)
  ! fc2d%fc(:,:,1) = (fc2d%fc(:,:,1) + transpose(fc2d%fc(:,:,1))) / 2
  ! call impose_asr2(input%asr3, Sd%nat, fc2d)

  !> ASR imposed on RR representation of fc2d
  !>------------------------------------------------
  allocate(DRR(S%nat3, S%nat3, nR, nR))
  DRR = fc_sc2RR(sc_grid, S, Sd, fc2d%fc(:,:,1))
  if(input%asr3 /= 'no') call asr3(DRR)
  fc2d%fc(:,:,1) = fc_RR2sc(sc_grid, S, Sd, DRR)
  deallocate(DRR)
  call div_mass_fc2(Sd, fc2d)
  if(input%asr3 /= 'no') call print_message("ASR applied to fc2d")
  ! call fc2_recenter(Sd, fc2d, fc2d_centered, 2)
  !
  !--------------------------------------------------
  ! call fc2_recenter(Sd, fc2d, fc2d_centered, 2)
  ! open(138, file="band.dat")
  ! allocate(freq0(Sd%nat3))
  ! do iq = 1, out_grid%nqtot
  !   call freq_phq(out_grid%xq(:,iq), Sd, fc2d_centered, freq0)
  !   write(138,"(1000E20.8)") freq0
  ! enddo
  ! close(138)
  !
  CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
    input%nk_in(2), input%nk_in(3),&
    in_grid, scatter=.false., xq0=input%xk0_in)
  !
  call q_grid_copy(in_grid, sym_grid)
  call sym_grid%symmetrize(S)
  call q_grid_copy(sym_grid, in_grid_sym_scat)
  if(num_procs > 1) call in_grid_sym_scat%scatter()
  ! call q_grid_copy(out_grid, out_grid_sym_scat)
  ! if(.not. out_grid_sym_scat%symmetrized .and. &
  ! (out_grid_sym_scat%type == 'simple' .and. out_grid_sym_scat%type == 'grid')) &
  ! call out_grid_sym_scat%symmetrize(S)
  ! if(num_procs > 1) call out_grid_sym_scat%scatter()
  ! call q_grid_copy(sym_grid, out_grid)
  CALL fc2_sc%allocate(S, Sd, sc_grid)
  !
  ! Keep the grid in natural Fortran order.  The full-Born Green function
  ! remaps tetrahedron weights explicitly when filling the FFT buffer.
  ! call revert_grid(in_grid)

  ! call project(S, Sd, fc2_centered, fc2d_centered)
  ! call dca_selfnrg(S, input, fc2_treated, fc2_sc, in_grid, out_grid)
  !
  n_add = Sd%nat - nR*S%nat
  if (n_add < 0) then !> VACANCY ------------------------------------------
    allocate(map(S%nat, nR))
    map = map_uc2sc(S, Sd, sc_grid)
    allocate(new_tau(3,nR*S%nat))
    allocate(new_ityp(nR*S%nat))
    allocate(DRR(S%nat3, S%nat3, nR, nR))
    new_tau(:,:Sd%nat) = Sd%tau
    new_ityp(:Sd%nat) = Sd%ityp
    DRR = 0._dp
    new_it = Sd%nat + 1
    do concurrent(j1=1:3, j2=1:3, na1=1:S%nat, na2=1:S%nat, R1=1:nR, R2=1:nR, map(na1,R1) /= -1 .and. map(na2,R2) /= -1)
      DRR(j1+3*(na1-1), j2+3*(na2-1), R1, R2) = fc2d%fc(j1 + 3*(map(na1,R1)-1), j2 + 3*(map(na2,R2)-1), 1)
    enddo
    !
    do concurrent(i=1:S%nat, iR=1:nR, map(i,iR)==-1)
      new_tau(:,new_it) = (index2v_cart(iR, sc_grid, S%at) + S%tau(:,i)) / sc_grid
      new_ityp(new_it) = S%ityp(i)
      new_it = new_it + 1
    enddo
    call move_alloc(new_tau, Sd%tau)
    call move_alloc(new_ityp, Sd%ityp)
    ! call move_alloc(new_fc, fc2d%fc)
    Sd%nat = nR*S%nat
    deallocate(map)
  elseif (n_add > 0) then !> INCLUSION ------------------------------------
    ! allocate(R_map(Sd%nat))
    ! allocate(map(S%nat, nR))
    ! map = map_uc2sc(S, Sd, sc_grid)
    ! allocate(i_map(Sd%nat))
    !
    !
    allocate(new_tau0 (3,S%nat + n_add))
    allocate(new_ityp0(S%nat + n_add))
    new_tau0(:,:S%nat) = S%tau
    new_ityp0(:S%nat) = S%ityp
    allocate(new_tau (3,(S%nat + n_add)*nR))
    allocate(new_ityp((S%nat + n_add)*nR))
    new_tau(:,:Sd%nat) = Sd%tau
    new_ityp(:Sd%nat) = Sd%ityp
    S%atm(S%ntyp+1) = "void"
    Sd%atm(Sd%ntyp+1) = "void"
    !
    i = 0
    ! taudef = 0._dp
    do concurrent(isc=1:Sd%nat, fc2_sc%R_map(isc) == -1)
      i = i + 1
      tau = cryst2cart(Sd%tau(:,isc), S%bg, -1) * sc_grid
      !
      new_ityp0(S%nat + i) = S%ntyp+1 !Sd%ityp(isc)
      iR = v2index(INT(tau), sc_grid)
      fc2_sc%R_map(isc) = iR
      fc2_sc%at_map(isc) = S%nat + i
      new_tau0(:,S%nat + i) = cryst2cart( &
        (tau - index2v(iR, sc_grid)), S%at, 1)
      ! taudef = taudef + new_tau0(:,S%nat + i) / n_add
      j = 0
      do jR = 1, nR
        if (iR == jR) cycle
        j = j + 1
        new_ityp(Sd%nat + (i-1)*nR + j) = Sd%ntyp+1 !Sd%ityp(isc)
        new_tau(:,Sd%nat + (i-1)*nR + j) = cryst2cart( &
          (index2v(jR, sc_grid) + tau - index2v(iR, sc_grid)) / sc_grid, S%at, 1)
      enddo
    enddo
    !
    call move_alloc(new_tau0, S%tau)
    call move_alloc(new_tau, Sd%tau)
    call move_alloc(new_ityp0, S%ityp)
    call move_alloc(new_ityp, Sd%ityp)
    !
    ! map = map_uc2sc(S, Sd, sc_grid)
    allocate(DRR(S%nat3+n_add*3, S%nat3+n_add*3, nR, nR))
    allocate(new_fc0(S%nat3+n_add*3, S%nat3+n_add*3, nR))
    !
    DRR = 0._dp
    do concurrent(j1=1:3, j2=1:3, isc1=1:Sd%nat, isc2=1:Sd%nat)
      DRR(j1 + 3*(fc2_sc%at_map(isc1)-1), j2 + 3*(fc2_sc%at_map(isc2)-1), &
        fc2_sc%R_map(isc1), fc2_sc%R_map(isc2)) = &
        fc2d%fc(j1 + 3*(isc1-1), j2 + 3*(isc2-1), 1)
    enddo
    !
    new_fc0 = 0._dp
    new_fc0(:S%nat3, :S%nat3, :) = fc2_treated%fc
    call move_alloc(new_fc0, fc2_treated%fc)
    !
    allocate(new_fc0(S%nat3 + n_add*3, S%nat3 + n_add*3, size(fc2_centered%fc,3)))
    new_fc0 = 0._dp
    new_fc0(:S%nat3, :S%nat3, :) = fc2_centered%fc
    call move_alloc(new_fc0, fc2_centered%fc)
    !
    S%nat = S%nat + n_add
    S%nat3 = S%nat * 3
    Sd%nat = Sd%nat + (nR-1) * n_add
    Sd%nat3 = Sd%nat * 3
    S%ntyp = S%ntyp +1
    Sd%ntyp = Sd%ntyp +1
    !
  else
    allocate(DRR(S%nat3, S%nat3, nR, nR))
    DRR = fc_sc2RR(sc_grid, S, Sd, fc2d%fc)! - fc_uc2RR(fc2_treated)
  endif ! ---------------------------------------------------------------------------------------
  !
  fc2_sc%fc = DRR - fc_uc2RR(fc2_treated)
  deallocate(DRR)
  !
  call fc2_sc_centered%allocate(S, Sd, sc_grid)
  fc2_sc_centered%fc = fc2_sc%fc
  CALL fc2_sc_centered%center(sc_grid, S)
  !
  if(contain(input%mode, '1b')) &
    call main_defect(S, fc2_centered, fc2_sc_centered, in_grid, sym_grid, out_grid, input)
  if(contain(input%mode, 'fb')) &
    call full_born_center(S, input, fc2_centered, fc2_sc_centered, in_grid, sym_grid, out_grid)
  if(contain(input%mode, 'dca')) &
    call dca_selfnrg(S, Sd, input, fc2_centered, fc2_sc, out_grid)

  CALL stop_mpi()
contains
  subroutine distance_uc(fc2_centered, S, filename)
    type(forceconst2_grid), intent(in) :: fc2_centered
    type(ph_system_info), intent(in) :: S
    character(len=*), intent(in) :: filename
    integer :: na1, na2, iR, j1, j2
    real(dp) :: dist, dist2

    open(10, file=filename, status='replace')
    do iR = 1, fc2_centered%n_r
      do na1 = 1, S%nat
        do na2 = 1, S%nat
          dist = norm2(S%tau(:,na1) - S%tau(:,na2) + fc2_centered%xR(:,iR))
          dist2 = norm2(S%tau(:,na1) + fc2_centered%xR(:,iR)) + norm2(S%tau(:,na2))
          do j1 = 1, 3
            do j2 = 1, 3
              write(10, "(3E20.8, 3I5)") dist, dist2, &
                fc2_centered%fc(j1+3*(na1-1), j2+3*(na2-1), iR), &
                iR, na1, na2
            enddo
          enddo
        enddo
      enddo
    enddo
    close(10)
  end subroutine
  !
  subroutine write_fc2_sc_pixels(fc2_sc, filename)
    type(forceconst2_sc),intent(in) :: fc2_sc
    character(*), intent(in) :: filename
    !
    real(dp), allocatable :: p(:)
    !
    open(10, file=filename, status='replace')
    do R2 = 1, fc2_sc%n_r2
      allocate(p(fc2_sc%n_R1(R2)))
      do R1 = 1, fc2_sc%n_r1(R2)
        p(r1) = sum(abs(fc2_sc%fc(:,:,R1,R2)))
      enddo
      write(10, "(1000E15.5)") p
      deallocate(p)
    enddo
    close(10)
  end subroutine
  !
end program
