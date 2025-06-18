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
  type(forceconst2_grid) :: fc2, fc2d, fc2_, fc2_centered
  type(forceconst2_sc) :: fc2_sc
  type(q_grid) :: in_grid, out_grid, grid_
  type(code_input_type) :: input, input_
  class(forceconst3), pointer :: fc3, fc3_
  real(dp) :: alpha
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
  CALL READ_INPUT("TK", input, grid_, S, fc2, fc3)
  ! call in_grid%destroy()

  ! CALL out_grid%destroy()
  CALL setup_grid(input%grid_type, S%bg, input%nk(1), &
    input%nk(2), input%nk(3),&
    out_grid, scatter=.false., xq0=input%xk0)
  ! do iq = 1, out_grid%nq
  !   out_grid%xq(:,iq) = out_grid%xq(:,iq) / 40
  ! enddo
  if(input%use_symm) call out_grid%symmetrize(S)

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


  CALL read_fc2(input%file_mat2_final, Sd, fc2d)
  CALL aux_system(Sd)
  call impose_asr2("diff", Sd%nat, fc2d)
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
  call div_mass_fc2(Sd, fc2d)
  CALL fc2_recenter(S, fc2, fc2_centered, 2)

  ! call S_uc2sc(S, fc2%nq, S_sc)
  ! S_sc%ityp(1) = 2
  CALL fc2_sc%allocate(S, fc2%nq)
  ! alpha = 0._dp
  ! fc2_sc%fc = (1 - alpha) * fc_uc2RR(fc2) + alpha * fc_sc2RR(fc2%nq, S, Sd, fc2d%fc)
  ! CALL fc2_sc%center(fc2%nq, S, Sd)
  ! CALL green_inversion(S, fc2_centered, fc2d, fc2_sc, out_grid)

  ! call check_derivative_swap(fc2_sc, S%nat3)
  ! call frobenius_triangle(fc2_sc, S, 'center3.dat')

  ! call fc2_sc%deallocate()
  ! CALL fc2_sc%allocate(S, fc2%nq)
  ! fc2_sc%fc = fc_sc2RR(fc2%nq, S, Sd, fc2d%fc) - fc_uc2RR(fc2%nq, S, fc2%fc)

  ! CALL center2(fc2_sc, fc2%nq, S, Sd)
  ! call frobenius_triangle(fc2_sc, S, 'center2.dat')

  fc2_sc%fc = fc_sc2RR(fc2%nq, S, Sd, fc2d%fc) ! - fc_uc2RR(fc2)
  CALL main_defect(S, Sd, fc2_centered, fc2_sc, in_grid, out_grid, input)
  CALL stop_mpi()

end program
