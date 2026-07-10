program vcam
  use kinds, only : dp
  use constants, only : tpi
  use code_input, only : code_input_type, read_input
  use mpi_thermal, only: start_mpi, stop_mpi, num_procs
  use input_fc, only: read_fc2, aux_system, ph_system_info, &
    forceconst2_grid, allocate_fc2_grid
  use q_grids, only : q_grid, setup_grid, q_grid_copy
  use fc2_interpolate, only : fc2_recenter, freq_phq
  use asr2_module, only : impose_asr2
  use quter_defect, only : forceconst2_sc, fc_sc2RR, &
    allocate_fc2_sc, fc_uc2RR, div_mass0_RR, asr3
  use defect, only : full_born_center, write_self, &
  interpolate_self_on_shell, write_lw, write_freq
  use constants, only : RY_TO_CMM1
  use thtetra, only : set_wg, tetra_output
  use mpi_thermal, only: my_id, num_procs, mpi_bsum
  use thutils, only: id_mat
  use functions, only: invzmat
  use ph_velocity, only: velocity
  implicit none
  !
  type(tetra_output) :: wg_out
  type(code_input_type) :: input
  type(ph_system_info) :: S(4), S_, S_vca
  type(forceconst2_grid) :: fc2_, fc2(4), fc2_vca, fc2c_vca
  type(q_grid) :: out_grid, in_grid, sym_grid
  type(forceconst2_sc) :: fc2_sc(2)
  character(len=256) :: input_files(4)
  integer :: i, sc_grid(3), nR, iw, iq, ibnd
  real(dp) :: c, m_vca, def_mass, p(6)
  real(dp), allocatable :: DRR(:,:,:,:)
  complex(dp), allocatable :: Tq(:,:,:,:,:)
  complex(dp), allocatable :: T(:,:,:,:)
  complex(dp), allocatable :: A(:,:)
  complex(dp), allocatable :: self_energy(:,:,:), self_lw(:,:)
  character(len=256) :: out_file
  real(dp), allocatable :: velsq(:,:)
  !
  call start_mpi()
  !
  CALL READ_INPUT("DEF", input, out_grid, S_, fc2_)
  input_files = ["mat2R-si  ", "mat2R-ge4 ", "mat2D-sige", "mat2D-gesi"]
  ! A_(1-x) B_(x)
  do i = 1, 4
    call read_fc2(trim(input_files(i)), S(i), fc2(i))
    call aux_system(S(i))
    if(i < 3) then
      call impose_asr2("simple", S(i)%nat, fc2(i))
    else
      fc2(i)%fc(:,:,1) = (fc2(i)%fc(:,:,1) + transpose(fc2(i)%fc(:,:,1))) / 2._dp
    endif
    ! if(i < 3) call fc2_recenter(S(i), fc2(i), fc2c(i), 2)
  enddo
  !
  c = input%conc

  S_vca = S(1)
  if(allocated(S_vca%sqrtmm1)) deallocate(S_vca%sqrtmm1)
  S_vca%celldm(1) = (1-c) * S(1)%celldm(1) + c * S(2)%celldm(1)
  S_vca%omega = S(1)%omega / S(1)%celldm(1)**3 * S_vca%celldm(1)**3
  S_vca%alat = S_vca%celldm(1)
  S_vca%tpiba = tpi / S_vca%alat
  S_vca%amass(1:S_vca%ntyp) = &
    (1._dp-c) * S(1)%amass(1:S(1)%ntyp) + c * S(2)%amass(1:S(2)%ntyp)
  call aux_system(S_vca)

  call allocate_fc2_grid(fc2(1)%n_R, S(1)%nat, fc2_vca)
  fc2_vca%yR = fc2(1)%yR
  fc2_vca%xR = fc2(1)%xR
  fc2_vca%i_0 = fc2(1)%i_0
  fc2_vca%nq = fc2(1)%nq
  m_vca = S_vca%amass(1)
  fc2_vca%fc = ((1-c) * fc2(1)%fc + c * fc2(2)%fc) / m_vca
  call fc2_recenter(S_vca, fc2_vca, fc2c_vca, 2)
  !
  CALL setup_grid(input%grid_type_in, S_vca%bg, input%nk_in(1), &
    input%nk_in(2), input%nk_in(3),&
    in_grid, scatter=.false., xq0=input%xk0_in)
  call q_grid_copy(in_grid, sym_grid)
  call sym_grid%symmetrize(S_vca)
  !
  sc_grid = input%sc_grid
  nR = product(sc_grid)
  allocate(DRR(S(1)%nat3, S(1)%nat3, nR, nR))
  do i = 1, 2
    call fc2_sc(i)%allocate(S(i), S(i+2), sc_grid)
    DRR = fc_sc2RR(sc_grid, S(i), S(i+2), fc2(i+2)%fc)
    call asr3(DRR)
    fc2_sc(i)%fc = ( DRR - &
      fc_uc2RR(fc2(i))) / m_vca
    fc2_sc(i)%eps = 1._dp - S(3-i)%amass(1) / m_vca
    call fc2_sc(i)%center(sc_grid, S(i))
  enddo
  DRR = (1-c) * fc2_sc(1)%fc - c * fc2_sc(2)%fc
  fc2_sc(1)%fc = (1-c) * DRR
  fc2_sc(2)%fc = -c * DRR
  deallocate(DRR)
  !
  allocate(Tq(S_vca%nat3, S_vca%nat3, out_grid%nqtot, input%n_omega, 2))
  S_vca%atm = S(1)%atm
  call full_born_center(S_vca, input, fc2c_vca, fc2_sc(1), &
    in_grid, sym_grid, out_grid, Tq(:,:,:,:,1))
  S_vca%atm = S(2)%atm
  call full_born_center(S_vca, input, fc2c_vca, fc2_sc(2), &
    in_grid, sym_grid, out_grid, Tq(:,:,:,:,2))
  !
  allocate(T(S_vca%nat3, S_vca%nat3, out_grid%nqtot, input%n_omega))
  T = c * Tq(:,:,:,:,1) + (1-c) * Tq(:,:,:,:,2)
  deallocate(Tq)
  !
  call set_wg(S_vca, fc2c_vca, out_grid, input, wg_out)
  !
  allocate(A(S_vca%nat3, S_vca%nat3))
  allocate(self_energy(S_vca%nat3, out_grid%nqtot, input%n_omega))
  do iw = 1+my_id, input%n_omega, num_procs
    do iq = 1, out_grid%nqtot
      A = 0._dp
      do ibnd = 1, S_vca%nat3
        A(:,ibnd) = wg_out%w(:,wg_out%e(iq),iw) * T(:,ibnd,iq,iw)
      enddo
      A = id_mat(S_vca%nat3) + A
      call invzmat(S_vca%nat3, A)
      T(:,:,iq,iw) = matmul(T(:,:,iq,iw), A)
      do ibnd = 1, S_vca%nat3
        self_energy(ibnd,iq,iw) = T(ibnd,ibnd,iq,iw)
      enddo
    enddo
  enddo
  call mpi_bsum(S_vca%nat3, out_grid%nqtot, input%n_omega, self_energy)
  !
  select case(input%calculation)
   case('lw')
    allocate(self_lw(S_vca%nat3, out_grid%nqtot))
    call interpolate_self_on_shell(wg_out%en, wg_out%f, self_energy, self_lw)
    write(out_file, "(A,F3.1,A)") "lw-vca-", input%conc, ".dat"
    call write_lw(wg_out%f, self_lw, out_file)
    write(out_file, "(A,F3.1,A)") "f-vca-", input%conc, ".dat"
    call write_freq(out_file, out_grid, wg_out%f)
    S_vca%lrigid = .true.
    allocate(velsq(S_vca%nat3, out_grid%nqtot))
    do iq = 1, out_grid%nqtot
      velsq(:,iq) = sum(velocity(S_vca, fc2c_vca, out_grid%xq(:,iq))**2,1)
    enddo
    write(out_file, "(A,F3.1,A)") "vel-vca-", input%conc, ".dat"
    call write_freq(out_file, out_grid, velsq)
   case('self')
    call write_self("self-fb.dat", wg_out%en, self_energy)
  end select
  call stop_mpi()
end program
