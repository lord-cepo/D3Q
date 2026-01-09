module defect_proj
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  ! use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
  use code_input, only: code_input_type
  use functions, only: invzmat
  use constants, only: tpi
  use quter_defect
  use functions, only: f_gauss
  use fc3_interpolate, only: forceconst3, sparse, d3_mixed, sum_R3
  use merge_degenerate, only: merge_degen
  use test_print, only: allclose
  use simtet, only: tetra_init_sym_cmplx, tetra_weights_green_cmplx
  use ph_velocity, only : velocity
  use constants, only : RY_TO_CMM1
  !
contains
  subroutine project(S, Sd, fc2, fc2d)
    type(ph_system_info), intent(in) :: S, Sd
    type(forceconst2_grid), intent(in) :: fc2, fc2d
    !
    integer :: iq, N, i, na, j
    real(dp), allocatable :: xq(:,:)
    real(dp), allocatable :: R(:,:)
    real(dp), allocatable :: freqs(:,:)
    real(dp), allocatable :: F(:)
    complex(dp), allocatable :: Us(:,:,:)
    complex(dp), allocatable :: v0(:)
    complex(dp), allocatable :: U(:,:), C(:,:), U0(:,:)
    integer, allocatable :: R_map(:), i_map(:)
    !
    N = product(fc2%nq)
    allocate(freqs(S%nat3,N))
    allocate(Us(S%nat3,S%nat3,N))
    allocate(xq(3,N))
    allocate(R(3,N))
    allocate(U(S%nat3*N,S%nat3*N))
    allocate(C(S%nat3*N,S%nat3*N))
    allocate(U0(S%nat3*N,S%nat3*N))
    allocate(F(S%nat3*N))
    allocate(v0(S%nat3*N))
    !
    R_map = map_sc2uc(S, Sd, fc2%nq, "R")
    i_map = map_sc2uc(S, Sd, fc2%nq, "nat")
    xq = grid_vec_cart(fc2%nq, S%bg)
    do iq = 1, size(xq,2)
      xq(:,iq) = xq(:,iq) / fc2%nq
    end do
    !
    R = grid_vec_cart(fc2%nq, S%at)
    !
    do iq = 1, N
      call freq_phq_safe(xq(:,iq), S, fc2, freqs(:,iq), Us(:,:,iq))
    enddo
    do i = 1, 3
      do na = 1, S%nat
        Us(i+(na-1)*3,i,1) = 1 / sqrt(real(S%nat,dp))  ! set acoustic modes at Gamma to uniform translation
      end do
    end do
    ! print*, dot_product(Us(:,1,1), Us(:,1,1))
    !
    call freq_phq_safe([0.0_dp, 0.0_dp, 0.0_dp], Sd, fc2d, F, U)
    do i = 1, 3
      do na = 1, Sd%nat
        U(i+(na-1)*3,i) = 1 / sqrt(real(Sd%nat,dp))  ! set acoustic modes at Gamma to uniform translation
      end do
    end do
    ! print*, dot_product(U(:,1), U(:,1))
    !
    do iq = 1, N
      do i = 1, S%nat3
        do j = 1, Sd%nat
          U0( 3*(j-1)+1:3*j, (iq-1)*S%nat3 + i ) = &
            e_iqr(xq(:,iq), R(:,R_map(j))) * &
            Us(3*(i_map(j)-1)+1:3*i_map(j),i,iq) / &
            sqrt(real(N,dp))
        enddo
      enddo
    enddo

    ! do i = 1, Sd%nat
    !   v0(3*(i-1)+1:3*i) = Us(3*(i_map(i)-1)+1:3*i_map(i),4,1) / sqrt(real(N,dp))
    ! enddo
    !
    C = matmul(conjg(transpose(U)), U0)
    ! print*, abs(C(24,4))**2
    ! print*, abs(C(27,4))**2
    ! print*, freqs(4,1)
    ! print*, abs(dot_product(U(:,N*S%nat3-2), U0(:,4)))**2
    ! print*, F(24)
    ! print*, F(27)
    ! do i = 1, S%nat3*N
    !   print"(2E20.8)", real(U(i,N*S%nat3)), real(v0(i))
    ! enddo
    !
    open(unit=99, file='defect_proj.dat', status='replace')
    do i = 1, size(C,2)
      write(99,*) abs(C(:,i))**2
    end do
    close(99)
    !
    !
    ! open(unit=100, file='dos_uc_sc.dat', status='replace')
    ! do iq = 1, N
    !   do i = 1, S%nat3
    !     write(100,"(2E20.8)") F(i + (iq-1)*S%nat3), freqs(i,iq)
    !   enddo
    ! enddo
    !
    close(100)
    deallocate(freqs, Us, xq, R, U, C, U0, F)
  end subroutine

end module
