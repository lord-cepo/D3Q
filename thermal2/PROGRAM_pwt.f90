program pwt
  use kinds, only : dp
  use mpi_thermal
  use asr2_module, only : impose_asr2
  use q_grids, only : q_grid, setup_grid
  use code_input, only : READ_INPUT
  use fc2_interpolate, only : freq_phq_safe, fftinterp_mat2, mat2_diag
  use input_fc, only: forceconst2_grid, read_fc2, div_mass_fc2
  use ph_system, only : read_system, ph_system_info, aux_system
  use parameters, only : npk
  use constants, only: pi
  use thtetra, only : tetra_init, tetra_weights_delta_sym, tetra_init_sym, &
    tetra_weights_green
  use simtet, only : tetra_init_sym_cmplx, tetra_weights_green_cmplx
  use thutils, only: print_message
  implicit none
  !
  type(ph_system_info) :: S
  type(forceconst2_grid) :: fc2
  type(q_grid) :: grid
  integer, parameter :: N = 25
  integer, parameter :: nef = 520
  ! external kpoint_grid
  integer :: iq, ief, i
  real(dp) :: ef, maxfreq
  complex(dp), allocatable :: wg(:,:,:)
  complex(dp) :: dos(nef)
  real(dp), allocatable :: ek_sym(:,:)
  complex(dp), allocatable :: ek_sym_cmplx(:,:)
  complex(dp), allocatable:: D(:,:)
  real(dp),allocatable  :: e(:)
  real(dp), parameter :: MULT = 1e5_dp
  integer :: iR, a, b, j
  real(dp) :: ssum
  integer :: ibnd, jbnd, irm1
  real(dp) :: mean, max_diff_fc2d, max_diff_rel
  logical :: ok
  !
  call start_mpi()
  !
  call read_fc2('reference/mat2R444', S, fc2)
  call aux_system(S)
  call impose_asr2("simple", S%nat, fc2)
  allocate(D(S%nat3, S%nat3))
  allocate(e(S%nat3))

  ! do i = 1, S%nat
  !   do a = 1, 3
  !     do b = 1, 3
  !       ssum = 0._dp
  !       do j = 1, S%nat
  !         do iR = 1, fc2%n_R
  !           ssum = ssum + fc2%FC(3*(i-1)+a, 3*(j-1)+b, iR)
  !         enddo
  !       enddo
  !       print*, ssum
  !     enddo
  !   enddo
  ! enddo

  do ibnd = 1, S%nat3
    do jbnd = 1, S%nat3
      do iR = 1, fc2%n_R
        do irm1 = 1, fc2%n_R
          if (ALL(fc2%yR(:,irm1) == -fc2%yR(:,iR))) then
            ok = .true.
            exit
          endif
        enddo
        if (.not. ok) call errore("minus R not found")
        mean = (fc2%fc(ibnd,jbnd,ir) + fc2%fc(jbnd,ibnd,irm1)) / 2.0_dp
        if (ABS(fc2%fc(ibnd,jbnd,ir) - fc2%fc(jbnd,ibnd,irm1)) > max_diff_fc2d) then
          max_diff_fc2d = ABS(fc2%fc(ibnd,jbnd,ir) - fc2%fc(jbnd,ibnd,irm1))
          max_diff_rel = max_diff_fc2d / ABS(mean)
          if (max_diff_rel > 1.97285) print*, ibnd, jbnd, ir, irm1
        endif
      enddo
    enddo
  enddo
  print*, "max difference", max_diff_fc2d
  print*, "max rel diff  ", max_diff_rel

  call div_mass_fc2(S, fc2)
  !
  call fftinterp_mat2([0._dp, 0._dp, 0._dp], S, fc2, D)
  call mat2_diag(S%nat3, D, e)
  print"(6E15.4)", e(:6)
  ! call setup_grid('simple', S%bg, N, N, N, grid)
  ! call grid%symmetrize(S)
  ! !
  ! nirr = grid%nqtot
  ! allocate(freqs(S%nat3,nirr), wg(S%nat3,nirr,nef))
  ! do iq = 1, nirr
  !   call freq_phq_safe(xq(:,iq), S, fc2, freqs(:,iq))
  ! enddo
  ! maxfreq = maxval(freqs)
  ! print*, maxfreq
  ! !

  !#################################################
  ! here's the thtetra part
  !#################################################
  !
  ! deallocate(wg)
  do i = 1, 1
    call setup_grid('simple', S%bg, N, N, N, grid)
    if (i == 1) call grid%symmetrize(S)
    allocate(ek_sym(S%nat3, grid%nqtot))
    allocate(wg(S%nat3,N**3,nef))
    ! allocate(wg    (S%nat3, grid%nqtot))
    do iq = 1, grid%nqtot
      call freq_phq_safe(grid%xq(:,iq), S, fc2, ek_sym(:,iq))
    enddo
    maxfreq = maxval(ek_sym) * 1.01_dp
    !
    call tetra_init_sym(grid, S, ek_sym, .false.)
    !
    ! print*, grid%w
    dos = 0._dp
    do ief = 1, nef
      ef = maxfreq / real(nef, dp) * ief
      wg(:,:,ief) = - tetra_weights_green(ef) / pi
      dos(ief) = sum(wg(:,:,ief))
      ! if(ief == 50) print*, wg(3,:)
      ! dos(ief) = sum(wg)
    enddo
    !
    if (i == 1) then
      open(11, file='noopt_sym.dat', status='unknown')
    else
      open(11, file='noopt_nosym.dat', status='unknown')
    endif
    do ief = 1, nef
      ! do iq = 1, grid%nqtot
      !   do ibnd = 1, S%nat3
      write(11,"(3E15.4)") ief*maxfreq/real(nef, dp), dos(ief)
      !   enddo
      ! enddo
    enddo
    close(11)
    !
    call grid%destroy()
    deallocate(ek_sym)
    deallocate(wg)
    call print_message('cycle done')
  enddo
  !
  !#################################################
  ! SimTet part
  !#################################################
  !
  ! call setup_grid('simple', S%bg, N, N, N, grid)
  ! call grid%symmetrize(S)
  ! allocate(ek_sym_cmplx(S%nat3, grid%nqtot))
  ! allocate(ek_sym(S%nat3, grid%nqtot))
  ! allocate(wg(S%nat3,grid%nqtot,nef))
  ! ! allocate(wg    (S%nat3, grid%nqtot))
  ! do iq = 1, grid%nqtot
  !   call freq_phq_safe(grid%xq(:,iq), S, fc2, ek_sym(:,iq))
  ! enddo
  ! maxfreq = maxval(ek_sym)
  ! !
  ! ek_sym_cmplx = cmplx(MULT * ek_sym**2, 0._dp, dp)
  ! call tetra_init_sym_cmplx(grid, S, ek_sym_cmplx)
  ! !
  ! ! print*, grid%w
  ! dos = 0._dp
  ! do ief = 1, nef
  !   ef = maxfreq / real(nef, dp) * ief
  !   wg(:,:,ief) = -2/pi*tetra_weights_green_cmplx(MULT * ef**2)
  !   dos(ief) = sum(AIMAG(matmul(wg(:,:,ief), grid%w))) * MULT
  !   ! if(ief == 50) print*, wg(3,:)
  !   ! dos(ief) = sum(wg)
  ! enddo
  ! !
  ! open(11, file='simtet.dat', status='unknown')
  ! do ief = 1, nef
  !   ! do iq = 1, grid%nqtot
  !   !   do ibnd = 1, S%nat3
  !   write(11,*) ief*maxfreq/real(nef, dp), dos(ief)
  !   !   enddo
  !   ! enddo
  ! enddo
  ! close(11)
  ! !
  ! call grid%destroy()
  ! deallocate(ek_sym)
  ! deallocate(wg)
  call stop_mpi()
end program
