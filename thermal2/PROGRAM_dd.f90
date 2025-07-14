program defectp
  use defect
  USE fc3_interpolate,  ONLY : forceconst3
  use q_grids, only: setup_grid
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2
  use asr2_module, only: impose_asr2
  use thutils, only: v2index
  use quter_defect
  use test_print
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd, S_
  type(forceconst2_grid) :: fc2d, fc2_centered, fc2_periodic
  type(forceconst2_sc) :: fc2_sc
  type(q_grid) :: in_grid, out_grid, grid_
  type(code_input_type) :: input, input_
  class(forceconst3), pointer :: fc3, fc3_
  real(dp) :: alpha
  integer :: sc_grid(3), i, R, R1, R2
  real(dp), allocatable :: p(:)
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

  if(all(input%sc_grid == -1)) then
    sc_grid = fc2_periodic%nq
  else
    sc_grid = input%sc_grid
    print*, "Using sc_grid = ", input%sc_grid
  end if
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
  CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
    input%nk_in(2), input%nk_in(3),&
    in_grid, scatter=.true., xq0=input%xk0_in)
  ! call in_grid%symmetrize(S)
  ! call in_grid%scatter()

  CALL read_fc2(input%file_mat3, Sd, fc2d)
  CALL aux_system(Sd)
  call impose_asr2(input%asr3, Sd%nat, fc2d)
  call div_mass_fc2(Sd, fc2d)
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
  call print_message("ASR applied to fc2d")
  CALL fc2_sc%allocate(S, sc_grid)
  fc2_sc%fc = fc_sc2RR(sc_grid, S, Sd, fc2d%fc) - fc_uc2RR(fc2_periodic)
  print*, SUM(abs(fc_uc2RR(fc2_periodic)))
  print*, SUM(abs(fc_sc2RR(sc_grid, S, Sd, fc2d%fc)))

  open(10, file="fc2_sc.dat", status='replace')
  allocate(p(product(fc2_centered%nq)))
  do R2 = 1, fc2_sc%n_r2
    do R1 = 1, fc2_sc%n_r1(R2)
      p(r1) = sum(abs(fc2_sc%fc(:,:,R1,R2)))
    enddo
    write(10, "(1000E15.5)") p
  enddo
  close(10)

  ! call fftinterp_mat2()
  ! fc2_sc%fc = 0._dp
  ! ! fc2_sc%fc = (fc_sc2RR(sc_grid, S, Sd, fc2d%fc) - fc_uc2RR(fc2_periodic)) * S%sqrtmm1(1)**2
  ! do i = 1, 3
  !   fc2_sc%fc(i,i,1,1) = 0.1028_dp !* fc2_periodic%fc(i,i,1)
  ! enddo


  ! Sd%sqrtmm1(1:3) = Sd%sqrtmm1(4:6)

  ! call S_uc2sc(S, sc_grid, S_sc)
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
  CALL fc2_sc%center_full(sc_grid, S, Sd)
  ! call check_derivative_swap(fc2_sc, S%nat3)
  CALL main_defect(S, fc2_centered, fc2_sc, in_grid, out_grid, input)
  CALL stop_mpi()

end program
