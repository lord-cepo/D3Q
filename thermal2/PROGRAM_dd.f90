program defectp
  use defect
  USE fc3_interpolate,  ONLY : forceconst3
  use q_grids, only: setup_grid
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2
  use asr2_module, only: impose_asr2
  use thutils, only: v2index
  use quter_defect, only : map_uc2sc
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd
  type(forceconst2_grid) :: fc2, fc2d
  type(q_grid) :: in_grid, out_grid
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  ! integer :: wait_for_debugger
  integer, allocatable :: atoms(:,:)
  integer :: na1, na2, j1, j2, na1_sc, na2_sc, jn1, jn2, R1, R2, nR
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

  CALL READ_INPUT("TK", input, out_grid, S, fc2, fc3)
  CALL out_grid%destroy()

  CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), &
    input%nk_in(2), input%nk_in(3),&
    in_grid, scatter=.true., xq0=input%xk0_in)

  CALL setup_grid(input%grid_type, S%bg, input%nk(1), &
    input%nk(2), input%nk(3),&
    out_grid, scatter=.true., xq0=input%xk0)

  CALL read_fc2(input%file_mat2_final, Sd, fc2d)
  CALL aux_system(Sd)

  !> probably one should also divide by mass of the UC, I'm not sure #TOFIX
  nR = product(fc2%nq)
  allocate(atoms(S%nat, nR))
  atoms = map_uc2sc(S, Sd, fc2%nq)
  do R1 = 1, nR
    do R2 = 1, nR
      do na1 = 1, S%nat
        na1_sc = atoms(na1,R1)
        do na2 = 1, S%nat
          na2_sc = atoms(na2,R2)
          do j1 = 1, 3
            jn1 = j1 + 3*(na1_sc-1)
            do j2 = 1, 3
              jn2 = j2 + 3*(na2_sc-1)
              fc2d%FC(jn1, jn2, 1) = &
                fc2d%FC(jn1, jn2, 1) * &
                S%sqrtmm1(j1 + 3*(na1-1)) * S%sqrtmm1(j2 + 3*(na2-1)) ! UC here
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  CALL main_defect(S, Sd, fc2, fc2d, in_grid, out_grid, input)
  CALL stop_mpi()

end program
