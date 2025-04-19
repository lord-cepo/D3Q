program bands
  use kinds,  ONLY : dp
  USE fc3_interpolate,  ONLY : forceconst3
  use fc2_interpolate,  ONLY : add_rgd_blk_d3, forceconst2_grid, fc2_recenter, freq_phq_safe, fftinterp_mat2
  use q_grids, only: setup_grid, q_grid
  use code_input, only: READ_INPUT, code_input_type
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2, ph_system_info
  use asr2_module, only: impose_asr2
  use thutils, only: v2index, braket, print_message
  use quter_defect!, only : forceconst2_sc, fc_sc2RR, fc_uc2RR, S_uc2sc, center3, freq_phq_degen
  ! use quter_defect, only : map_uc2sc, fc_uc2RR
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd, S_sc
  type(forceconst2_grid) :: fc2, fc2d, fc2_periodic
  type(forceconst2_sc) :: fc2_sc, fc2_uc
  type(q_grid) :: path, grid
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  integer :: iq, ibnd, jq, jbnd, i, j, k
  real(dp), allocatable :: freq0(:,:), freq(:,:,:,:), fc_uc(:,:,:,:), freq1(:,:,:,:)
  complex(dp), allocatable :: D(:,:), Ui(:,:), Uj(:,:), D1(:,:), TU(:,:)
  integer :: ir1, ir2, jn1, jn2
  real(dp), allocatable :: distance(:,:,:,:), distance0(:,:,:,:), distance_uc(:,:,:), dummy(:)
  ! integer :: wait_for_debugger
  real(dp) :: mean, max_diff_fc2d, max_diff_rel
  !
  CALL start_mpi()
  !
  CALL READ_INPUT("LW", input, grid, S, fc2_periodic, fc3)
  ! call setup_grid('simple', S%bg, input%nk(1),input%nk(2),input%nk(3), grid)
  ! grid%xq = grid%xq / input%e0

  allocate(dummy(S%nat3))
  allocate(freq0(S%nat3,grid%nq))
  allocate(freq1(S%nat3,S%nat3,grid%nq,grid%nq))
  allocate(freq(S%nat3,S%nat3,grid%nq,grid%nq))
  ! allocate(freq0, freq1, source=freq)
  allocate(D(S%nat3,S%nat3))
  allocate(Ui, Uj, D1, TU, source=D)
  allocate(fc_uc(S%nat3,S%nat3,product(fc2_periodic%nq),product(fc2_periodic%nq)))

  CALL read_fc2(input%prefix, Sd, fc2d)
  CALL aux_system(Sd)
  call impose_asr2(input%asr2, Sd%nat, fc2d)
  CALL div_mass_fc2(Sd, fc2d)
  max_diff_fc2d = 0._dp
  do ibnd = 1, Sd%nat3
    do jbnd = ibnd+1, Sd%nat3
      mean = (fc2d%fc(ibnd,jbnd,1) + fc2d%fc(jbnd,ibnd,1)) / 2.0_dp
      if (ABS(fc2d%fc(ibnd,jbnd,1) - fc2d%fc(jbnd,ibnd,1)) > max_diff_fc2d) then
        max_diff_fc2d = ABS(fc2d%fc(ibnd,jbnd,1) - fc2d%fc(jbnd,ibnd,1))
        max_diff_rel = max_diff_fc2d / ABS(mean)
      endif
      fc2d%fc(ibnd,jbnd,1) = mean
      fc2d%fc(jbnd,ibnd,1) = mean
    enddo
  enddo
  print*, "max_diff_fc2d", max_diff_fc2d
  print*, "max_diff_rel", max_diff_rel



  CALL fc2_sc%allocate(S, fc2_periodic%nq)
  call fc2_uc%allocate(S, fc2_periodic%nq)

  ! fc2_sc%fc = fc_sc2RR(fc2_periodic%nq, S, Sd, fc2d%fc)
  fc2_sc%fc = fc_sc2RR(fc2_periodic%nq, S, Sd, fc2d%fc) - fc_uc2RR(fc2_periodic%nq, S, fc2_periodic%fc)
  fc2_uc%fc = fc_uc2RR(fc2_periodic%nq, S, fc2_periodic%fc)

  call fc2_recenter(S, fc2_periodic, fc2, 2, distance_uc)
  call S_uc2sc(S, fc2_periodic%nq, S_sc)
  S_sc%ityp(1) = 2
  call center2(fc2_sc, fc2%nq, S, Sd)
  call center3(fc2_uc, fc2%nq, S, Sd)
  ! call center2(fc2_sc, fc2%nq, S, Sd)

  ! D1 = 0._dp

  ! print*, "fc2_uc", size(fc2_uc%fc)
  ! do i = 1, 6
  !   do j = 1, 6
  !     do k = 1, 3
  !       ! print"(2e15.3)", fc2_uc%fc(i,j,k,k+1), fc2_sc%fc(i,j,k,k+1)
  !       if (ABS(fc2_uc%fc(i,j,k,1)-fc2_sc%fc(i,j,k,1)) > 1e-10) then
  !         print"(2e15.3)", fc2_uc%fc(i,j,k,1), fc2_sc%fc(i,j,k,1)
  !         print*, i, j, k
  !         if(any(fc2_uc%yR2(:,k,1) /= fc2_sc%yR2(:,k,1))) print*, "not equal"
  !       endif
  !     enddo
  !   enddo
  ! enddo
  ! do i = 1, product(fc2%nq)
  !   print*, fc2_uc%n_R2(i), fc2_sc%n_R2(i)
  ! enddo

  freq = 0._dp
  freq1 = 0._dp
  do iq = 1, grid%nq
    call freq_phq_safe(grid%xq(:,iq), S, fc2, freq0(:,iq), Ui)
    call fc2_uc%interpolate(grid%xq(:,iq), S)
    call fc2_sc%interpolate(grid%xq(:,iq), S)
    ! do ibnd = 1, S%nat3
    !   freq(ibnd,1,iq,1) = SQRT(REAL(braket(Ui(:,ibnd), D),DP))
    ! enddo

    ! call fc2_sc%interpolate(grid%xq(:,iq), S)
    ! call fc2_sc%interpolate(grid%xq(:,iq), S, D)
    ! do ibnd = 1, S%nat3
    !   freq1(ibnd,1,iq,1) = SQRT(REAL(braket(Ui(:,ibnd), D),DP))
    ! enddo


    ! freq1(:,1,iq,1) = AIMAG(braket(Ui(:,ibnd), D))
    TU = CONJG(TRANSPOSE(Ui(:,:)))

    do jq = 1, grid%nq
      if(iq /= jq .and. input%q_summed) cycle
      call freq_phq_safe(grid%xq(:,jq), S, fc2, dummy, Uj)

      call fc2_sc%interpolate(grid%xq(:,jq), S, D)
      call freq_phq_degen(grid%xq(:,jq), S, fc2, dummy, D, Uj)

      freq(:,:,jq,iq) = ABS(matmul(CONJG(TRANSPOSE(Uj)), matmul(D, Uj)))

      call fc2_uc%interpolate(grid%xq(:,jq), S, D)

      ! CALL fftinterp_mat2(xq, S, fc2, D)

      freq1(:,:,jq,iq) = ABS(matmul(TU, matmul(D, Uj)))
    enddo
  enddo

  print*, SUM(ABS(freq)) / grid%nq / S%nat3

  ! do iq = 1, path%nq
  !   call freq_phq_safe(path%xq(:,iq), S, fc2, freq0(:,iq), U)
  !   call fc2_sc%interpolate(path%xq(:,iq), S)
  !   call fc2_sc%interpolate(path%xq(:,iq), S, D)
  !   ! D = fc2_sc%mix(:,:,1)
  !   ! call mat2_diag(S%nat3, D, freq(:,iq))
  !   ! freq(:,iq) = SQRT(freq(:,iq))
  !   ! CALL add_rgd_blk_d3(path%xq(:,iq), S, D1)
  !   do ibnd = 1, S%nat3
  !     freq(ibnd,iq) = REAL(braket(U(:,ibnd), D),dp)/ 2 / freq0(ibnd,iq)
  !   enddo
  !   ! call mat2_diag(S%nat3, U, freq(:,iq))
  !   ! freq(:,iq) = SQRT(freq(:,iq))
  ! enddo

  ! open(10, file='fc.dat', status='unknown')
  ! ! do iq = 1, path%nq
  ! !   write(10, *) path%xq(:,iq), freq0(:,iq), freq(:,iq)
  ! ! enddo
  ! do ir1 = 1, fc2_sc%n_R1
  !   do ir2 = 1, fc2_sc%n_R2(ir1)
  !     do jn1 = 1, S%nat3
  !       do jn2 = 1, S%nat3
  !         write(10,*) distance(jn1, jn2, ir2, ir1), fc2_sc%fc(jn1,jn2,ir2,ir1)
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(10)
  ! !
  ! open(20, file='fc0.dat', status='unknown')
  ! do ir1 = 1, fc2_uc%n_R1
  !   do ir2 = 1, fc2_uc%n_R2(ir1)
  !     do jn1 = 1, S%nat3
  !       do jn2 = 1, S%nat3
  !         write(20,*) distance0(jn1, jn2, ir2, ir1), fc2_uc%fc(jn1,jn2,ir2,ir1)
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(20)
  ! !
  ! open(30, file='fcuc.dat', status='unknown')
  ! do ir1 = 1, fc2%n_R
  !   do jn1 = 1, S%nat3
  !     do jn2 = 1, S%nat3
  !       write(30,*) distance_uc(jn1, jn2, ir1), fc2%fc(jn1,jn2,ir1)
  !     enddo
  !   enddo
  ! enddo
  ! close(30)
  !
  ! open(40, file='freq.dat', status='unknown')
  ! do iq = 1, 1
  !   do jq = 1, grid%nq
  !     ! if (input%q_summed .and. iq /= jq) cycle
  !     do ibnd = 1, 1
  !       do jbnd = 1, S%nat3
  !         ! if (input%q_resolved .and. ibnd /= jbnd) cycle
  !         write(40, *) freq0(ibnd,iq), freq(jbnd,ibnd,jq,iq), freq1(jbnd,ibnd,jq,iq)
  !       enddo
  !     enddo
  !   enddo
  ! enddo

  open(40, file='f2.dat', status='unknown')
  do iq = 1, grid%nq
    do jq = 1, grid%nq
      if (input%q_summed .and. iq /= jq) cycle
      write(40, '(6E13.4)', advance='no') freq0(:,iq)
      do ibnd = 1, S%nat3
        do jbnd = 1, S%nat3
          if (input%q_resolved .and. ibnd /= jbnd) cycle
          write(40, '(E13.4)', advance='no') freq0(ibnd,iq) + freq(jbnd,ibnd,jq,iq)/freq0(ibnd,iq)/2
        enddo
      enddo
      write(40, *) ''
    enddo
  enddo
  close(40)
  !
  CALL stop_mpi()

end program
