program bands
  use kinds,  ONLY : dp
  USE fc3_interpolate,  ONLY : forceconst3
  use fc2_interpolate,  ONLY : &
    add_rgd_blk_d3, forceconst2_grid, fc2_recenter, &
    freq_phq_safe, fftinterp_mat2, mat2_diag, freq_phq
  use q_grids, only: setup_grid, q_grid
  use code_input, only: READ_INPUT, code_input_type
  use mpi_thermal, only: start_mpi, stop_mpi
  use input_fc, only: read_fc2, aux_system, div_mass_fc2, &
    ph_system_info, deallocate_fc2_grid, allocate_fc2_grid
  use asr2_module, only: impose_asr2
  use thutils, only: v2index, braket, print_message
  use quter_defect!, only : forceconst2_sc, fc_sc2RR, fc_uc2RR, S_uc2sc, center3, freq_phq_degen
  ! use quter_defect, only : map_uc2sc, fc_uc2RR
  IMPLICIT NONE
  !
  type(ph_system_info) :: S, Sd, S_sc
  type(forceconst2_grid) :: fc2, fc2d, fc2_periodic, fc2d_centered, fc20_periodic, fc20_centered
  type(forceconst2_sc) :: fc2_sc, fc2_scd
  type(q_grid) :: path, grid
  type(code_input_type) :: input
  class(forceconst3), pointer :: fc3
  integer :: iq, ibnd, jq, jbnd, j, k, i1, i2, na1, na2
  real(dp), allocatable :: freq0(:,:), freq1(:,:), freq2(:,:), fc_uc(:,:,:,:)
  complex(dp), allocatable :: D(:,:), Ui(:,:), Uj(:,:), D1(:,:), D2(:,:), TU(:,:)
  integer :: ir1, ir2, jn1, jn2, ir1t, ir2t, iconc
  real(dp), allocatable :: distance(:,:,:,:), distance0(:,:,:,:), distance_uc(:,:,:), dummy(:)
  ! integer :: wait_for_debugger
  real(dp) :: mean, max_diff_fc2d, max_diff_rel
  character(len=100) :: filename
  real(dp) :: conc(6)
  real(dp), allocatable :: pixel(:,:)
  integer :: indices(4)
  real(dp) :: ssomma
  integer :: sc_grid(3)
  complex(dp), allocatable :: Dd(:,:)
  CHARACTER(len=6), EXTERNAL :: int_to_char

  !
  CALL start_mpi()
  !
  CALL READ_INPUT("DEF", input, grid, S, fc2_periodic)
  !
  if(all(input%sc_grid == -1)) then
    sc_grid = fc2_periodic%nq
  else
    sc_grid = input%sc_grid
    print*, "Using sc_grid = ", input%sc_grid
  end if
  !
  allocate(freq0(S%nat3,grid%nq))
  allocate(freq1(S%nat3,grid%nq))
  allocate(freq2(S%nat3,grid%nq))
  allocate(D(S%nat3,S%nat3))
  allocate(D1(S%nat3,S%nat3))

  allocate(Ui(S%nat3,S%nat3))

  print*, S%at(:,2)
  print*, S%at(:,3)
  CALL read_fc2(input%file_mat3, Sd, fc2d)
  CALL aux_system(Sd)
  allocate(Dd(Sd%nat3,Sd%nat3))
  call impose_asr2(input%asr3, Sd%nat, fc2d)
  CALL div_mass_fc2(Sd, fc2d)
  ! call S_uc2sc(S, Sd, fc2_periodic%nq, S_sc)
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

  call fc2_recenter(Sd, fc2d, fc2d_centered, 2)
  call fc2_recenter(S, fc2_periodic, fc2, 2)
  !
  conc = [0._dp, 0.05_dp, 0.07_dp, 0.11_dp, 0.16_dp, 0.30_dp] * product(sc_grid)
  do iconc = 1, size(conc)
    CALL fc2_sc%allocate(S, fc2_periodic%nq)
    fc2_sc%fc = (1-conc(iconc)) * fc_uc2RR(fc2_periodic) / product(fc2_periodic%nq)
    CALL fc2_sc%center(fc2%nq, S)
    !
    call fc2_scd%allocate(S, sc_grid)
    fc2_scd%fc = conc(iconc) * fc_sc2RR(sc_grid, S, Sd, fc2d%fc) / product(sc_grid)
    call fc2_scd%center(sc_grid, S)
    !
    ! call cut_distance_rr(S, fc2_sc, 1.35_dp)
    ! call cut_fc_rr(S, fc2_sc, 1e-9_dp)
    ! call equate_distance_rr(S, fc2_sc, 1.0_dp)
    ! call write_distance_uc(S_sc, fc2d_centered, "distance-sc-"//trim(input%file_mat3(11:))//".dat", plus=.true.)
    ! call write_distance_rr(S, fc2_sc, "distance-diff-"//trim(input%file_mat3(11:))//".dat", plus=.false.)
    ! call write_distance_rr(S, fc2_scd, "distance-"//trim(input%file_mat3(:5))//".dat", plus=.false.)
    ! call write_distance_uc(S, fc2, "distance-"//trim(input%file_mat2)//".dat", plus=.false.)
    ! call write_distance_rr(S, fc2_sc, "distance-diff-"//trim(input%file_mat3(11:))//".dat", plus=.true.)
    ! call write_distance_rr(S, fc2_scd, "distance-"//trim(input%file_mat3(11:))//".dat", plus=.true.)
    ! call write_distance_uc(S, fc2, "distance-"//trim(input%file_mat2(11:))//".dat", plus=.true.)
    !
    ! call allocate_fc2_grid(1, Sd%nat, fc20_periodic)
    ! fc20_periodic%nq = [1,1,1]
    ! fc20_periodic%fc(:,:,1) = fc_uc2sc(S, S_sc, fc2%nq, fc2_periodic%fc)
    ! call fc2_recenter(S_sc, fc20_periodic, fc20_centered, 2)
    !
    freq1 = 0._dp
    do iq = 1, grid%nq
      call freq_phq(grid%xq(:,iq), S, fc2, freq0(:,iq), Ui)
      ! call freq_phq_safe(grid%xq(:,iq), S_sc, fc20_centered, freq1(:,iq))
      ! call freq_phq_safe(grid%xq(:,iq), Sd, fc2d_centered, freq2(:,iq))
      ! call freq_phq_safe(grid%xq(:,iq), Sd, fc2d_centered, freq2(:,iq))
      !
      ! call fc2_sc%r2q(grid%xq(:,iq))
      ! call fc2_sc%r2q(grid%xq(:,iq), D)
      ! call mat2_diag(S%nat3, D, freq0(:,iq))

      ! call fftinterp_mat2(grid%xq(:,iq), S, fc2, D)

      call fc2_scd%r2q(grid%xq(:,iq))
      call fc2_scd%r2q(grid%xq(:,iq), D1)

      call fc2_sc%r2q(grid%xq(:,iq))
      call fc2_sc%r2q(grid%xq(:,iq), D)
      ! CALL add_rgd_blk_d3(grid%xq(:,iq), S, D1)

      ! print*, SUM(ABS(D - D1 / product(fc2%nq))**2) / sum(abs(D)**2)

      ! call add_rgd_blk_d3(grid%xq(:,iq), S, D)
      ! do ibnd = 1, S%nat3
      !   freq1(ibnd,iq) = freq0(ibnd,iq) + real(braket(Ui(:,ibnd), D)) / 2 / freq0(ibnd,iq)
      !   ! freq0(ibnd,iq) = REAL(D(ibnd,ibnd), dp)
      ! enddo
      D = D1 + D
      call mat2_diag(S%nat3, D, freq1(:,iq))

      ! freq2(:,iq) = freq0(:,iq) + freq2(:,iq) / 2 / freq0(:,iq)
      !
      ! call fftinterp_mat2(grid%xq(:,iq), S, fc2(i), D)
      ! freq1(3,iq) = real(D(1,1), dp)
      ! freq1(4,iq) = aimag(D(1,1))
    enddo
    ! freq0 = SQRT(freq0/product(fc2%nq))
    freq1 = SQRT(freq1)
    !
    ! open(10, file="dd.dat")
    ! call fftinterp_mat2([0.1_dp,0.2_dp,0.3_dp], S_sc, fc20_centered, Dd)
    ! do ibnd = 1, Sd%nat3
    !   write(10, "(2000E20.8)") Dd(:,ibnd)
    ! enddo
    ! close(10)
    !
    ! call write_freq(freq0, "f"//trim(int_to_char(NINT(conc(iconc)*100)))//".dat")
    call write_freq(freq1, "f"//trim(int_to_char(NINT(conc(iconc)*100/product(sc_grid))))//".dat")
    call fc2_sc%deallocate()
    call fc2_scd%deallocate()
  enddo
  ! call write_freq(freq2, "freq2-"//trim(input%file_mat2(11:))//".dat")
  ! call write_freq(freq1, "freq-sc-"//trim(input%file_mat2(11:))//".dat")
  ! call write_freq(freq2, "freq-"//trim(input%file_mat3(11:))//".dat")
  !
  CALL stop_mpi()
contains
  subroutine write_distance_uc(S, fc2, filename, plus)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    character(len=*), intent(in) :: filename
    logical, intent(in) :: plus
    !
    character(len=100) :: filename_
    integer :: na1, na2, iR, j1, j2, sign
    real(dp) :: ss
    !
    filename_ = filename
    if (plus) then
      sign = 1
      filename_ = "sum-"//trim(filename)
    else
      sign = -1
    endif
    !
    open(10, file=filename_)
    do na1 = 1, S%nat
      do na2 = 1, S%nat
        do iR = 1, fc2%n_R
          do j1 = 1, 3
            do j2 = 1, 3
              ss = fc2%fc(3*(na1-1)+j1, 3*(na2-1)+j2, iR)
              if(abs(ss) > 1e-16_dp) then
                write(10, "(2E20.8,5I4)") &
                  norm2(S%tau(:,na1) + fc2%xr(:,iR) + sign * S%tau(:,na2)), abs(ss), na1, na2, iR, j1, j2
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine
  !
  subroutine write_distance_rr(S, fc2, filename, plus)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(in) :: fc2
    character(len=*), intent(in) :: filename
    logical, intent(in) :: plus
    character(len=100) :: filename_
    !
    integer :: na1, na2, R1, R2, j1, j2, sign
    real(dp) :: ss
    !
    filename_ = filename
    if (plus) then
      sign = 1
      filename_ = "sum-"//trim(filename)
    else
      sign = -1
    endif
    !
    open(10, file=filename_)
    do na2 = 1, S%nat
      do na1 = 1, S%nat
        do R2 = 1, fc2%n_R2
          do R1 = 1, fc2%n_R1(R2)
            do j1 = 1, 3
              do j2 = 1, 3
                ss = fc2%fc(3*(na1-1)+j1, 3*(na2-1)+j2, R1, R2)
                if(abs(ss) > 1e-16_dp) then
                  write(10, "(2E20.8,6I4)") &
                    norm2(S%tau(:,na1) + fc2%xr1(:,R1,R2) + sign*(S%tau(:,na2) + fc2%xr2(:,R2))), abs(ss), na1, na2, R1, R2, j1, j2
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    close(10)
  end subroutine
  !
  subroutine cut_distance_rr(S, fc2, max_dist)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(inout) :: fc2
    real(dp), intent(in) :: max_dist
    !
    integer :: na1, na2, R1, R2
    real(dp) :: ss
    !
    do na2 = 1, S%nat
      do na1 = 1, S%nat
        do R2 = 1, fc2%n_R2
          do R1 = 1, fc2%n_R1(R2)
            if(norm2(S%tau(:,na1) + fc2%xr1(:,R1,R2) - S%tau(:,na2) - fc2%xr2(:,R2)) > max_dist) &
              fc2%fc(3*(na1-1)+1:3*na1, 3*(na2-1)+1:3*na2, R1, R2) = 0._dp
          enddo
        enddo
      enddo
    enddo
  end subroutine
  !
  subroutine equate_distance_rr(S, fc2, min_fc)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(inout) :: fc2
    real(dp), intent(in) :: min_fc
    !
    integer :: na1, na2, R1, R2
    real(dp) :: dist
    !
    do na2 = 1, S%nat
      do na1 = 1, S%nat
        do R2 = 1, fc2%n_R2
          do R1 = 1, fc2%n_R1(R2)
            dist = norm2(S%tau(:,na1) + fc2%xr1(:,R1,R2) - S%tau(:,na2) - fc2%xr2(:,R2))
            if (dist < min_fc) &
              fc2%fc(3*(na1-1)+1:3*na1, 3*(na2-1)+1:3*na2, R1, R2) = 1e-5_dp * exp(-dist)
          enddo
        enddo
      enddo
    enddo
  end subroutine
  !
  subroutine cut_fc_rr(S, fc2, min_fc)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(inout) :: fc2
    real(dp), intent(in) :: min_fc
    !
    integer :: na1, na2, R1, R2, j1, j2
    real(dp) :: dist
    !
    do na2 = 1, S%nat
      do na1 = 1, S%nat
        do R2 = 1, fc2%n_R2
          do R1 = 1, fc2%n_R1(R2)
            do j1 = 1, 3
              do j2 = 1, 3
                if (fc2%fc(3*(na1-1)+j1, 3*(na2-1)+j2, R1, R2) < min_fc) &
                  fc2%fc(3*(na1-1)+j1, 3*(na2-1)+j2, R1, R2) = 0._dp
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine
  !
  subroutine write_freq(freq, filename)
    real(dp), intent(in) :: freq(:,:)
    character(*), intent(in) :: filename
    !
    open(40, file=TRIM(filename))
    do iq = 1, size(freq,2)
      write(40, '(1000E13.4)') freq(:,iq)
    enddo
    close(40)
  end subroutine
  !
end program
