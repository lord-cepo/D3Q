program bands
  use defect
  USE fc3_interpolate,  ONLY : forceconst3
  use fc2_interpolate,  ONLY : add_rgd_blk_d3
  use q_grids, only: setup_grid
  use code_input, only: READ_INPUT
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2
  use asr2_module, only: impose_asr2
  use thutils, only: v2index
  ! use quter_defect, only : map_uc2sc, fc_uc2RR
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd
  type(forceconst2_grid) :: fc2, fc2d, fc2_periodic
  type(forceconst2_sc) :: fc2_sc
  type(q_grid) :: path
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  integer :: iq, ibnd
  real(dp), allocatable :: freq0(:,:), freq(:,:), fc_uc(:,:,:,:)
  complex(dp), allocatable :: D(:,:), U(:,:), D1(:,:)
  ! integer :: wait_for_debugger
  !
  CALL start_mpi()
  !

  CALL READ_INPUT("LW", input, path, S, fc2_periodic, fc3)
  allocate(freq(S%nat3,path%nq), freq0(S%nat3,path%nq))
  allocate(D(S%nat3,S%nat3))
  allocate(U, D1, source=D)
  allocate(fc_uc(S%nat3,S%nat3,product(fc2_periodic%nq),product(fc2_periodic%nq)))

  CALL read_fc2(input%prefix, Sd, fc2d)
  CALL aux_system(Sd)
  call impose_asr2('simple', Sd%nat, fc2d)
  CALL div_mass_fc2(Sd, fc2d)

  CALL fc2_sc%allocate(S, fc2_periodic%nq)

  ! fc2_sc%fc = fc_sc2RR(fc2_periodic%nq, S, Sd, fc2d%fc)! - &
  fc2_sc%fc = fc_uc2RR(fc2_periodic%nq, S, fc2_periodic%fc)

  call fc2_recenter(S, fc2_periodic, fc2, 2)
  call fc2_sc%center(fc2%nq, S, Sd)

  D1 = 0._dp

  ! call setup_grid('simple', S%bg, fc2%nq(1), fc2%nq(2), fc2%nq(3), grid)

  do iq = 1, path%nq
    call freq_phq_safe(path%xq(:,iq), S, fc2, freq0(:,iq), U)
    call fc2_sc%interpolate(path%xq(:,iq), S)
    call fc2_sc%interpolate(path%xq(:,iq), S, D)
    ! D = fc2_sc%mix(:,:,1)
    ! call mat2_diag(S%nat3, D, freq(:,iq))
    ! freq(:,iq) = SQRT(freq(:,iq))
    ! CALL add_rgd_blk_d3(path%xq(:,iq), S, D1)
    do ibnd = 1, S%nat3
      freq(ibnd,iq) = SQRT(REAL(braket(U(:,ibnd), D),dp))
    enddo
    ! call mat2_diag(S%nat3, U, freq(:,iq))
    ! freq(:,iq) = SQRT(freq(:,iq))
  enddo

  open(10, file='bands.dat', status='unknown')
  do iq = 1, path%nq
    write(10, *) path%xq(:,iq), freq0(:,iq), freq(:,iq)
  enddo
  close(10)
  !
  CALL stop_mpi()

end program
