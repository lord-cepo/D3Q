module full_born
  use kinds, only: dp
  use thtetra, only: tetra_init_sym, tetra_init, tetra_weights_green, &
    deallocate_tetra, tetra_output, set_wg
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe, &
    fc2_recenter, fftinterp_mat2, mat2_diag
  use thutils
  ! only: outer_product, freq_in_grid, interp1_matrix, &
  ! id_mat, braket, index2v, v2index, e_iqr, interp1_tns4, grid_vec
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  use mpi_thermal, only: mpi_bsum, ionode, num_procs, my_id, ierr
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
  use symm_q_mat, only: apply_sym_q_full
  !
contains
  !
  function g0_p(D0, E)
    real(dp), intent(in) :: D0(:,:)
    real(dp), intent(in) :: E
    !
    complex(dp), allocatable :: g0_p(:,:)
    complex(dp) :: IE(size(D0,1),size(D0,2))
    !
    IE = id_mat(size(D0,1)) * (E**2 + 1e-6_dp*(0.0_dp, 1.0_dp))
    g0_p = IE-D0
    call invzmat(size(D0,1), g0_p)
    !
  end function
  !
  !   function g0_p1(D0, iw, wg, fc2)
  !   real(dp), intent(in) :: D0(:,:)
  !   real(dp), intent(in) :: E
  !   !
  !   complex(dp), allocatable :: g0_p1(:,:)
  !   complex(dp) :: IE(size(D0,1),size(D0,2))
  !   !
  !   IE = id_mat(size(D0,1)) * (E**2 + 1e-6_dp*(0.0_dp, 1.0_dp))
  !   g0_p1 = IE-D0
  !   call invzmat(size(D0,1), g0_p1)
  !   !
  ! end function
  !
  subroutine full_born_p(input, S, Sd, fc2, D0, V, grid)
    type(code_input_type), intent(in) :: input
    type(ph_system_info), intent(in) :: S, Sd
    type(forceconst2_grid), intent(in) :: fc2
    real(dp), intent(in) :: D0(:,:)
    complex(dp), intent(in) :: V(:,:)
    type(q_grid), intent(in) :: grid
    !
    real(dp) :: x, dx
    integer :: x0
    complex(dp) :: lws(S%nat3)
    real(dp) :: freqs(S%nat3, grid%nqtot)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    !
    complex(dp), allocatable :: g0(:,:), VgV(:,:), Vg(:,:)
    complex(dp), allocatable :: TR(:,:,:)
    ! type(forceconst2_sc) :: T_sc
    integer :: iw, N, nR, iR, na1, na2, iq, j, ibnd
    !
    character(20) :: filename
    complex(dp) :: lwsq(S%nat3, grid%nqtot)
    !
    real(dp), allocatable :: diffs(:,:,:,:,:), weights(:,:,:,:)
    integer, allocatable:: n_weights(:,:,:)
    complex(dp), allocatable :: Tq(:,:,:,:)
    real(dp) :: omega, omega_max
    !
    call freq_in_grid(S, fc2, grid, freqs, Us)
    N = size(D0,1)
    nR = size(V,1) / S%nat3
    allocate(g0(N,N), Vg(N,N), VgV(N,N))
    allocate(TR(S%nat3, S%nat3, nR))
    allocate(Tq(S%nat3, S%nat3, grid%nqtot, input%n_omega))
    !
    ! call print_matrix(nR, V, "V")
    ! call print_matrix(nR, cmplx(D0, 0._dp, dp), "D0")
    !
    omega_max = maxval(freqs) * 1.2_dp
    !
    call minimal_image_2(S, input%sc_grid, diffs, n_weights, weights)
    Tq = 0.0_dp
    do iw = 1, input%n_omega
      print*, iw
      omega = (iw-1) * omega_max / input%n_omega
      VgV = 0.0_dp
      Vg = 0.0_dp
      !
      g0 = g0_p(D0, omega)
      write(filename, "(A2,I2.2)") "g0", iw
      call print_matrix(nR, g0, filename)
      ! VgV = matmul(V, matmul(g0, V))
      call zgemm_N(N, V, g0, Vg)
      call zgemm_N(N, Vg, V, VgV)
      !
      ! write(filename, "(A,I4.4)") "VgV_", iw
      ! call print_matrix(nR, VgV, filename)
      !
      call symmetrize_t(input%sc_grid, fc_sc2RR_cmplx(input%sc_grid, S, Sd, VgV), TR)
      !
      do iq = 1, grid%nqtot
        do iR = 1, nR
          do na1 = 1, S%nat
            do na2 = 1, S%nat
              do j = 1, n_weights(iR,na1,na2)
                Tq((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq,iw) = &
                  Tq((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iq,iw) + &
                  TR( (na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,iR) + &
                  e_iqr(grid%xq(:,iq), diffs(:,j,iR,na1,na2)) * weights(j,iR,na1,na2)
              enddo
            enddo
          enddo
        enddo
        Tq(:,:,iq,iw) = matmul(conjg(transpose(Us(:,:,iq))), matmul(Tq(:,:,iq,iw), Us(:,:,iq)))
      enddo
      !
    enddo
    !
    open(123, file="fbp.dat")
    ! do iw = 1, input%n_omega
    !   write(123, "(2E20.8)") wg%en(iw), A(iw)
    ! enddo
    ! do iw = 1, input%n_omega
    !   dos = 0._dp
    !   do iq = 1, sym_grid%nqtot
    !     do ibnd = 1, S%nat3
    !       dos = dos + G_q(ibnd,ibnd,iq,iw) * wg%qw(iq)
    !     enddo
    !   enddo
    !   write(123, "(3E20.8)") wg%en(iw), dos
    ! enddo
    do iq = 1, grid%nqtot
      do ibnd = 1, S%nat3
        x = freqs(ibnd,iq)*input%n_omega/omega_max+ 1
        x0 = INT(x)
        dx = x - x0
        lws(ibnd) = (1.0_dp - dx) * Tq(ibnd,ibnd,iq,x0) + dx * Tq(ibnd,ibnd,iq,x0+1)
      enddo
      call merge_degen(S%nat3, lws, freqs(:,iq))
      do ibnd = 1, S%nat3
        write(123, "(3E20.8)") freqs(ibnd,iq), lws(ibnd)!/sym_freqs(ibnd,iq)
      enddo
    enddo
    close(123)
    !
    open(103, file="self.dat")
    do iw = 1, input%n_omega
      do iq = 1, grid%nqtot
        do ibnd = 1, 1
          lwsq(ibnd,iq) = Tq(ibnd,ibnd,iq,iw)
        enddo
      enddo
      write(103, "(1200E20.8)") aimag(lwsq(1,:))
    enddo
    close(103)
    !
    deallocate(g0)
  contains
    subroutine print_matrix(nR, A, fn)
      integer, intent(in) :: nR
      complex(dp), intent(in) :: A(:,:)
      character(len=*), intent(in) :: fn
      !
      integer :: iR, jR, nat3
      complex(dp) :: Tij
      !
      nat3 = size(A,1) / nR
      !
      open(10, file=trim(fn)//".re")
      open(11, file=trim(fn)//".im")
      do iR = 1, nR
        do jR = 1, nR
          Tij = sum(A((iR-1)*nat3+1:iR*nat3,(jR-1)*nat3+1:jR*nat3))
          write(10, "(E20.8)") real(Tij, dp)
          write(11, "(E20.8)") aimag(Tij)
        enddo
      enddo
      close(10)
      close(11)
    end subroutine
    !
  end subroutine
  !
  subroutine symmetrize_t(sc_grid, TRR, TR)
    integer, intent(in) :: sc_grid(3)
    complex(dp), intent(in) :: TRR(:,:,:,:)
    complex(dp), intent(out) :: TR(size(TRR,1), size(TRR,2), size(TRR,4))
    !
    integer :: iR, jR, nR
    integer :: Ri(3), Rj(3)
    !
    nR = product(sc_grid)
    !
    TR = 0.0_dp
    do iR = 1, nR
      Ri = index2v(iR, sc_grid)
      do jR = 1, nR
        Rj = index2v(jR, sc_grid)
        TR(:,:,iR) = TR(:,:,iR) + TRR(:,:,v2index(mod(Ri + Rj, sc_grid), sc_grid), jR)
      enddo
    enddo
    !
    TR = TR / real(nR, dp)
    !
  end subroutine
  !
  subroutine full_born_analytical(input, S, fc2, grid, sym_grid, out_grid, epsilon)
    use constants, only: RY_TO_CMM1
    !
    type(code_input_type), intent(in) :: input
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(q_grid), intent(in) :: grid, sym_grid, out_grid
    real(dp), intent(in) :: epsilon
    !
    real(dp) :: freqs(S%nat3, grid%nqtot)
    complex(dp) :: Us(S%nat3, S%nat3, grid%nqtot)
    type(tetra_output) :: wg
    complex(dp) :: g0(S%nat3,S%nat3)
    integer :: iq, jq, ibnd, jbnd, iw
    real(dp) :: Vb(S%nat3,S%nat3), V(S%nat3,S%nat3)
    !
    complex(dp) :: I_gV(S%nat3,S%nat3), T(S%nat3,S%nat3)
    complex(dp) :: green_weight(S%nat3,sym_grid%nqtot), self_on_shell
    real(dp) :: fwhm(S%nat3,out_grid%nqtot)
    complex(dp) :: self(input%n_omega)
    complex(dp) :: outer_products(S%nat3,S%nat3,S%nat3,grid%nqtot)
    !
    real(dp)    :: out_freqs(S%nat3,out_grid%nqtot)
    complex(dp) :: out_Us(S%nat3,S%nat3,out_grid%nqtot)
    !
    call freq_in_grid(S, fc2, out_grid, out_freqs, out_Us)
    call freq_in_grid(S, fc2, grid, freqs, Us)
    call set_wg(S, fc2, sym_grid, input, wg, mult = 1.05_dp)
    if(input%n_omega <= 1) &
      call errore("full_born_analytical", "n_omega must be > 1 to interpolate T", 1)
    !
    print*, "max freq", wg%en(input%n_omega) * RY_TO_CMM1
    Vb = 0.0_dp
    do ibnd = 1, 3
      Vb(ibnd,ibnd) = 1._dp
    enddo
    !
    do iq = 1, grid%nqtot
      do ibnd = 1, S%nat3
        outer_products(:,:,ibnd,iq) = outer_product(Us(:,ibnd,iq))
      enddo
    enddo
    self = 0._dp
    do iw = 1+my_id, input%n_omega, num_procs
      g0 = 0.0_dp
      do jq = 1, grid%nqtot
        do jbnd = 1, S%nat3
          g0 = g0 + outer_products(:,:,jbnd,jq) * wg%w(jbnd,wg%e(jq),iw)
        enddo
      enddo
      !
      V = Vb * epsilon * wg%en(iw)**2
      I_gV = id_mat(S%nat3) - (1-input%conc)*matmul(g0, V) !+ Vb * (0._dp, 1e-6_dp)
      call invzmat(S%nat3, I_gV )
      T = 2 * input%conc * matmul(V, I_gV)
      !
      ! In the analytical Si/Ge mass case T is scalar in phonon space, so
      ! store only the scalar self-energy and interpolate it on shell below.
      self(iw) = T(1,1)
    enddo
    call mpi_bsum(input%n_omega, self)
    !
    fwhm = 0._dp
    do iq = 1, out_grid%nqtot
      do ibnd = 1, S%nat3
        if(out_freqs(ibnd,iq) <= 0._dp) cycle
        self_on_shell = interp1_scl(self, out_freqs(ibnd,iq)*input%n_omega/wg%max_f)
        fwhm(ibnd,iq) = -aimag(self_on_shell) / out_freqs(ibnd,iq) * RY_TO_CMM1
      enddo
    enddo
    !
    if(ionode) then
      open(124, file="fba.dat")
      do iq = 1, out_grid%nqtot
        write(124, "(I6,100E20.8)") iq, out_freqs(:,iq)*RY_TO_CMM1, fwhm(:,iq)
      enddo
      close(124)
    endif
    !
  end subroutine
end module
