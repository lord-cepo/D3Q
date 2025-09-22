program defectp
  use defect
  USE fc3_interpolate,  ONLY : forceconst3
  use q_grids, only: setup_grid
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2, write_fc2
  use asr2_module, only: impose_asr2
  use thutils, only: v2index, cryst2cart, index2v_cart
  use quter_defect
  use quter_module, only : quter
  use test_print
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd, S_, S_sc
  type(forceconst2_grid) :: fc2d, fc2d_centered, fc2_centered, fc2_periodic, fc2_treated
  type(forceconst2_sc) :: fc2_sc
  type(q_grid) :: in_grid, out_grid, grid_
  type(code_input_type) :: input, input_
  class(forceconst3), pointer :: fc3, fc3_
  real(dp) :: alpha
  integer :: n_add, new_it, isc
  integer, allocatable :: map(:,:), R_map(:), i_map(:), new_ityp(:)
  real(dp), allocatable :: new_tau(:,:), new_fc(:,:,:,:)
  real(dp), allocatable :: D_nx_real(:,:,:)
  integer, allocatable :: inclusion_isc(:)
  complex(dp), allocatable :: inclusion_U(:,:)

  real(dp), dimension(3) :: d1, d2
  integer :: sc_grid(3), i, R, R1, R2, na1, na2, j1, j2, nR, j, iq, isc1, isc2
  real(dp), allocatable :: p(:)
  real(dp), allocatable :: freqs_sc(:)
  real(dp) :: v(3)
  real(dp), allocatable :: interp_grid(:,:)
  complex(dp), allocatable :: Ds(:,:,:), matq(:,:,:,:,:)
  character(len=100) :: filename
  real(dp) :: max_norm
  real(dp) :: freq0(9)
  ! integer :: wait_for_debugger
  ! integer, allocatable :: atoms(:,:)
  ! integer :: na1, na2, j1, j2, na1_sc, na2_sc, jn1, jn2, R1, R2, nR
  ! integer :: iq, i, j
  ! real(dp), allocatable :: dyn(:,:,:,:), zeu(:,:,:)
  !
  CALL start_mpi()
  !
  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  !
  ! if (ionode) then
  !   print *, "Rank 0 is waiting for debugger. Attach to PID:", getpid()
  !   wait_for_debugger = 1
  !   do while (wait_for_debugger == 1)
  !     ! Pause in the loop until debugger sets wait_for_debugger to 0
  !     call sleep(1)
  !   end do
  ! end if
  !

  ! CALL READ_INPUT("LW", input_, out_grid, S_, fc2_, fc3_)
  ! out_grid%nqtot = out_grid%nq
  ! call out_grid%destroy()
  ! call read_fc2('reference/mat2R_4periodic', S, fc2_periodic)

  CALL READ_INPUT("DEF", input, out_grid, S, fc2_periodic)
  CALL fc2_recenter(S, fc2_periodic, fc2_centered, 2)

  !
  if(all(input%sc_grid == -1)) then
    sc_grid = fc2_periodic%nq
  else
    sc_grid = input%sc_grid
    print*, "Using sc_grid = ", input%sc_grid
  end if
  nR = product(sc_grid)
  ! call allocate_fc2_grid(nR, S%nat, fc2_treated)
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

  ! do iq = 1, 10
  !   call freq_phq_safe([REAL(iq, dp)/10,0._dp,0._dp], S, fc2_centered, freq0)
  !   print"(9E20.8)", freq0
  !   call freq_phq_safe([REAL(iq, dp)/10,0._dp,0._dp], S, fc2_treated, freq0)
  !   print"(9E20.8)", freq0
  !   print*, "-----------"
  ! enddo
  !
  CALL read_fc2(input%file_mat3, Sd, fc2d)
  ! call S_uc2sc(S, Sd, sc_grid, S_sc)
  CALL aux_system(Sd)
  call impose_asr2(input%asr3, Sd%nat, fc2d)
  call div_mass_fc2(Sd, fc2d)
  ! call fc2_recenter(S_sc, fc2d, fc2d_centered, 2)
  call print_message("ASR applied to fc2d")
  !
  CALL fc2_sc%allocate(S, sc_grid)
  CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
    input%nk_in(2), input%nk_in(3),&
    in_grid, scatter=.true., xq0=input%xk0_in)
  !
  n_add = Sd%nat - nR*S%nat
  if (n_add < 0) then !> VACANCY ------------------------------------------
    allocate(map(S%nat, nR))
    map = map_uc2sc(S, Sd, sc_grid)
    allocate(new_tau(3,nR*S%nat))
    allocate(new_ityp(nR*S%nat))
    allocate(new_fc(S%nat3, S%nat3, nR, nR))
    new_tau(:,:Sd%nat) = Sd%tau
    new_ityp(:Sd%nat) = Sd%ityp
    new_fc = 0._dp
    new_it = Sd%nat + 1
    do concurrent(j1=1:3, j2=1:3, na1=1:S%nat, na2=1:S%nat, R1=1:nR, R2=1:nR, map(na1,R1) /= -1 .and. map(na2,R2) /= -1)
      new_fc(j1+3*(na1-1), j2+3*(na2-1), R1, R2) = fc2d%fc(j1 + 3*(map(na1,R1)-1), j2 + 3*(map(na2,R2)-1), 1)
    enddo
    !
    do concurrent(i=1:S%nat, R=1:nR, map(i,R)==-1)
      new_tau(:,new_it) = (index2v_cart(R, sc_grid, S%at) + S%tau(:,i)) / sc_grid
      new_ityp(new_it) = S%ityp(i)
      new_it = new_it + 1
    enddo
    fc2_sc%fc = new_fc - fc_uc2RR(fc2_periodic)
    call move_alloc(new_tau, Sd%tau)
    call move_alloc(new_ityp, Sd%ityp)
    ! call move_alloc(new_fc, fc2d%fc)
    Sd%nat = nR*S%nat
    deallocate(map)
  elseif (n_add > 0) then !> INCLUSION ------------------------------------
    allocate(R_map(Sd%nat))
    R_map = map_sc2uc(S, Sd, sc_grid, "R")
    allocate(map(S%nat, nR))
    map = map_uc2sc(S, Sd, sc_grid)
    allocate(i_map(Sd%nat))
    i_map = map_sc2uc(S, Sd, sc_grid, "nat")
    !
    allocate(D_nx_real(3*S%nat, 3*n_add, nR))
    allocate(inclusion_isc(n_add))
    new_it = 1
    do concurrent(isc=1:Sd%nat, R_map(isc) == -1)
      inclusion_isc(new_it) = isc
      new_it = new_it + 1
    enddo
    !
    do concurrent(R1=1:nR, na1=1:S%nat, j1=1:3, j2=1:3, new_it=1:n_add)
      D_nx_real(j1 + 3*(na1 -1), j2 + 3*(new_it-1), R1) = &
        fc2d%fc(j1 + 3*(map(na1, R1)-1), j2 + 3*(inclusion_isc(new_it)-1), 1)
    enddo
    !
    allocate(inclusion_U(3*n_add, 3*n_add))
    do concurrent(na1=1:n_add, na2=1:n_add, j1=1:3, j2=1:3)
      inclusion_U(j1 + 3*(na1-1), j2 + 3*(na2-1)) = &
        fc2d%fc(j1 + 3*(inclusion_isc(na1)-1), j2 + 3*(inclusion_isc(na2)-1), 1)
    enddo
    !
    allocate(fc2_sc%inclusion_Dnx_out(3*S%nat, 3*n_add, out_grid%nqtot))
    allocate(fc2_sc%inclusion_Dnx_in(3*S%nat, 3*n_add, in_grid%nqtot))
    fc2_sc%inclusion_Dnx_out = 0._dp
    fc2_sc%inclusion_Dnx_in = 0._dp
    do concurrent(R1=1:nR, iq=1:out_grid%nqtot)
      fc2_sc%inclusion_Dnx_out(:,:,iq) = fc2_sc%inclusion_Dnx_out(:,:,iq) + &
        D_nx_real(:,:,R1) * e_iqr(-out_grid%xq(:,iq), index2v_cart(R1, sc_grid, S%at))
    enddo !> It's probably possible to center it, but for now I'll leave it like this
    !> I also don't know if I need to put nR as normalization
    do concurrent(R1=1:nR, iq=1:in_grid%nqtot)
      fc2_sc%inclusion_Dnx_in(:,:,iq) = fc2_sc%inclusion_Dnx_in(:,:,iq) + &
        D_nx_real(:,:,R1) * e_iqr(-in_grid%xq(:,iq), index2v_cart(R1, sc_grid, S%at))
    enddo !> It's probably possible to center it, but for now I'll leave it like this
    !> I also don't know if I need to put nR as normalization
    !
    allocate(fc2_sc%inclusion_eig(3*n_add))
    call mat2_diag(3*n_add, inclusion_U, fc2_sc%inclusion_eig)
    do iq = 1, out_grid%nqtot
      fc2_sc%inclusion_Dnx_out(:,:,iq) = matmul(fc2_sc%inclusion_Dnx_out(:,:,iq), inclusion_U) / nR
    enddo
    do iq = 1, in_grid%nqtot
      fc2_sc%inclusion_Dnx_in(:,:,iq) = matmul(fc2_sc%inclusion_Dnx_in(:,:,iq), inclusion_U)
    enddo
    !
    allocate(new_fc(S%nat3, S%nat3, nR, nR))
    do concurrent(j1=1:3, j2=1:3, isc1=1:Sd%nat, isc2=1:Sd%nat, R_map(isc1) /= -1 .and. R_map(isc2) /= -1)
      new_fc(j1 + 3*(i_map(isc1)-1), j2 + 3*(i_map(isc2)-1), R_map(isc1), R_map(isc2)) = &
        fc2d%fc(j1 + 3*(isc1-1), j2 + 3*(isc2-1), 1)
    enddo
    fc2_sc%fc = new_fc - fc_uc2RR(fc2_periodic)
    deallocate(inclusion_isc, inclusion_U, D_nx_real, map, i_map, R_map)
  else
    fc2_sc%fc = fc_sc2RR(sc_grid, S, Sd, fc2d%fc) - fc_uc2RR(fc2_treated)
  endif ! ---------------------------------------------------------------------------------------
  print"(A,E15.4)", "average fc2 diff", sum(abs(fc2_sc%fc)) / size(fc2_sc%fc)
  call write_fc2_sc_pixels(fc2_sc, "pixels-"//trim(input%file_mat3(11:))//".dat")
  !
  ! call in_grid%destroy()
  ! call impose_asr2('simple', S%nat, fc2_periodic)
  ! call div_mass_fc2(S, fc2_periodic)

  ! CALL out_grid%destroy()
  ! CALL setup_grid(input%grid_type, S%bg, input%nk(1), &
  !   input%nk(2), input%nk(3),&
  !   out_grid, scatter=.false., xq0=input%xk0)
  ! ! do iq = 1, out_grid%nq
  ! !   out_grid%xq(:,iq) = out_grid%xq(:,iq) / 40
  ! ! enddo
  ! if(input%use_symm) call out_grid%symmetrize(S)

  ! out_grid%nqtot = 2
  ! out_grid%nq = 2
  ! out_grid%xq(:,1) = [0._dp, 0.0_dp, 0.1_dp]
  ! out_grid%xq(:,2) = [0._dp, -0.5_dp, -0.5_dp]
  ! call cryst_to_cart(1, out_grid%xq, S%bg, 1)
  ! call cryst_to_cart(out_grid%nq, out_grid%xq, S%at, -1)
  ! do iq = 1, out_grid%nq
  !   ! print"(3F8.2)", out_grid%xq(:,iq)
  !   do i = 1, 3
  !     if(out_grid%xq(i,iq) < 0._dp) &
  !       out_grid%xq(i,iq) = out_grid%xq(i,iq) + 1._dp
  !   enddo
  ! enddo
  ! call cryst_to_cart(out_grid%nq, out_grid%xq, S%bg, 1)

  ! call out_grid%scatter()
  !
  ! call in_grid%symmetrize(S)
  ! call in_grid%scatter()

  ! call impose_asr2("simple", Sd%nat, fc2d)

  ! allocate(dyn(3,3,Sd%nat, Sd%nat), zeu(3,3,Sd%nat))
  ! zeu = 0.0
  ! do i = 1, 3
  !   do j = 1, 3
  !     do na1 = 1, Sd%nat
  !       do na2 = 1, Sd%nat
  !         dyn(i,j,na1,na2) = fc2d%FC(i+3*(na1-1),j+3*(na2-1),1)
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! call set_asr("crystal   ", 3, Sd%nat, Sd%tau, dyn, zeu)

  ! allocate(freqs_sc(Sd%nat3))
  ! call freq_phq_safe([0._dp, 0._dp, 0._dp], Sd, fc2d_centered, freqs_sc)

  ! open(10, file="freqs_sc.dat", status='replace')
  ! do i = 1, size(freqs_sc)
  !   write(10, "(E15.5)") freqs_sc(i)
  ! enddo
  ! close(10)


  ! fc2d%fc(:,:,1) = fc2d%fc(:,:,1) - fc_uc2sc(S, S_sc, sc_grid, fc2_periodic%fc)
  ! call fc2_recenter(S_sc, fc2d, fc2d_centered, 2)
  ! write(filename, "(A, I1, A)") "distance-diff", sc_grid(1), "p.dat"
  ! call distance_uc(fc2d_centered, S_sc, filename)
  ! ! call write_fc2("mat2R_sc_pristine_333", S_sc, fc2d)
  ! call distance_uc(fc2d_centered, S_sc, "distance-pristine.dat")
  ! write(filename, "(A, I1, A)") "distance-uc", sc_grid(1), ".dat"
  ! call distance_uc(fc2_centered, S, filename)
  ! call fftinterp_mat2()
  ! fc2_sc%fc = 0._dp
  ! ! fc2_sc%fc = (fc_sc2RR(sc_grid, S, Sd, fc2d%fc) - fc_uc2RR(fc2_periodic)) * S%sqrtmm1(1)**2
  ! do i = 1, 3
  !   fc2_sc%fc(i,i,1,1) = 0.1028_dp !* fc2_periodic%fc(i,i,1)
  ! enddo


  ! Sd%sqrtmm1(1:3) = Sd%sqrtmm1(4:6)

  ! S_sc%ityp(1) = 2
  ! alpha = 0._dp
  ! fc2_sc%fc = (1 - alpha) * fc_uc2RR(fc2) + alpha * fc_sc2RR(sc_grid, S, Sd, fc2d%fc)
  ! CALL green_inversion(S, fc2_centered, fc2d, fc2_sc, out_grid)

  ! call frobenius_triangle(fc2_sc, S, 'center3.dat')

  ! call fc2_sc%deallocate()
  ! CALL fc2_sc%allocate(S, sc_grid)
  ! fc2_sc%fc = fc_sc2RR(sc_grid, S, Sd, fc2d%fc) - fc_uc2RR(sc_grid, S, fc2%fc)

  ! CALL center2(fc2_sc, sc_grid, S, Sd)
  ! call frobenius_triangle(fc2_sc, S, 'center2.dat')

  input%n_omega = input%n_omega - 1
  CALL fc2_sc%center(sc_grid, S)
  ! fc2_sc%fc = 0._dp
  ! do i = 1,3
  !   fc2_sc%fc(i,i,1,1) = 1._dp
  ! enddo
  !
  ! write(filename, "(A, I1, A)") "distance-cent", sc_grid(1), "p.dat"
  ! max_norm = 0._dp
  ! open(10, file=filename, status='replace')
  ! do R2 = 1, fc2_sc%n_r2
  !   do na2 = 1, S%nat
  !     d2 = S%tau(:,na2) + fc2_sc%xR2(:,R2)
  !     do R1 = 1, fc2_sc%n_r1(R2)
  !       do na1 = 1, S%nat
  !         d1 = S%tau(:,na1) + fc2_sc%xR1(:,R1,R2)
  !         max_norm = max(max_norm, norm2(d1 - d2))
  !         do j1 = 1, 3
  !           do j2 = 1, 3
  !             write(10, "(2E20.8)") norm2(d1 - d2), &
  !               fc2_sc%fc(j1+3*(na1-1), j2+3*(na2-1), R1, R2)
  !           enddo
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(10)
  ! call print_message("Distances written to "//TRIM(filename))
  ! print*, "Maximum distance:", max_norm
  !

  ! call check_derivative_swap(fc2_sc, S%nat3)
  CALL main_defect(S, fc2_centered, fc2_sc, in_grid, out_grid, input)
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
              write(10, "(3E20.8, 3I5)") dist, &
                fc2_centered%fc(j1+3*(na1-1), j2+3*(na2-1), iR), &
                dist2, iR, na1, na2
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
end program
