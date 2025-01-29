program fc2_test
  use defect
  USE fc3_interpolate,  ONLY : forceconst3
  use fc2_interpolate, only : mat2_diag
  use q_grids, only: setup_grid
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2
  use asr2_module, only: impose_asr2
  use thutils, only: v2index
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd
  type(forceconst2_grid) :: fc2, fc2d
  type(q_grid) :: out_grid
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  ! integer :: wait_for_debugger
  ! integer :: na1, na2, j1, j2, na1_sc, na2_sc, jn1, jn2, R1, R2, prod
  integer :: R1, R2, nR, j1, j2, jn1, jn2, na1, na2, na1_sc, na2_sc
  integer :: R(3)
  integer, allocatable :: atoms_sc(:,:)
  real(dp), allocatable :: VKR(:,:), freq_sc(:)
  complex(dp), allocatable :: VKRc(:,:)
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

  CALL read_fc2(input%file_mat2_final, Sd, fc2d)

  !> probably one should also divide by mass of the UC, I'm not sure #TOFIX
  nR = product(fc2%nq)
  allocate(atoms_sc(S%nat, nR))
  allocate(VKR(S%nat3*nR, S%nat3*nR))
  allocate(VKRc(S%nat3*nR, S%nat3*nR))
  atoms_sc = sc2uc(S, Sd, fc2%nq)
  do R1 = 1, nR
    do R2 = 1, nR
      R = index2v(R1, fc2%nq) - index2v(R2, fc2%nq)
      do j1 = 1, 3
        if (R(j1) < 0) R(j1) = R(j1) + fc2%nq(j1)
      enddo
      do na1 = 1, S%nat
        do na2 = 1, S%nat
          do j1 = 1, 3
            jn1 = j1 + 3*(na1-1)
            do j2 = 1, 3
              jn2 = j2 + 3*(na2-1)
              na1_sc = atoms_sc(na1, R1)
              na2_sc = atoms_sc(na2, R2)

              VKR(j1 + 3*(na1_sc-1), j2 + 3*(na2_sc-1)) = &
                fc2%FC(jn1, jn2, v2index(R, fc2%nq))
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo


  allocate(freq_sc(S%nat3*nR))
  VKRc = cmplx(VKR, 0._dp, dp)
  CALL mat2_diag(S%nat3*nR, VKRc, freq_sc)

  freq_sc = SQRT(freq_sc)

  open(10, file='freq_sc.txt', status='unknown')
  do j1 = 1, S%nat3*nR
    write(10,*) freq_sc(j1)
  enddo
  CALL stop_mpi()

end program
