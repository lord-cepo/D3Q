! This module is rewritten from the tetra.f90 in PW/src
!
MODULE thtetra
  !
  ! Tetrahedron method, linear and optimized. opt is better for all purpose
  ! weights for delta integration are obtained by multiplying gi with Iik in
  ! https://iopscience.iop.org/article/10.1088/0022-3719/12/15/008
  !
  ! another useful link, for delta integration where the integrand is 1 (DOS)
  ! http://staff.ustc.edu.cn/~zqj/posts/LinearTetrahedronMethod/#fn:tet_weight
  !
  ! optimized tetrahedron method, used in QE for theta integration in opt_tetra_weights
  ! https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.094515
  ! they multiply ni with Jik in the first article, then they transform (fit) through wlsm matrices
  USE kinds, ONLY: DP
  USE mpi_thermal, ONLY: my_id, num_procs, mpi_bsum, ionode
  use input_fc, only: ph_system_info, allocate_fc2_grid
  use q_grids, only: q_grid, q_grid_copy, q_grid_symmetrize
  use fc2_interpolate, only: forceconst2_grid, freq_phq_safe
  use thutils, only : freq_in_grid
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  INTEGER :: ntetra = 0
  !! number of tetrahedra
  INTEGER :: nntetra
  !! k-points per tetrahedron used to compute weights.
  !! 4 for linear / 20 for optimized tetrahedron method
  INTEGER, ALLOCATABLE :: tetra(:,:)
  !! index of k-points in a given tetrahedron shape (nntetra,ntetra)
  REAL(DP), ALLOCATABLE :: wlsm(:,:)
  !! Weights for the optimized tetrahedron method
  INTEGER :: nqs
  !! number of q-points in IRREDUCIBLE BZ
  INTEGER :: nqtot
  !! number of q-points in the whole BZ
  INTEGER :: nbnd
  !! number of bands
  INTEGER, ALLOCATABLE :: itetra(:,:,:)
  !! order index of vertices of each tetrahedron (4,nbnd,ntetra)
  REAL(DP), ALLOCATABLE :: ek_sort(:,:,:)
  !! sorted energies for each tetrahedron (4,nbnd,ntetra)
  REAL(DP), ALLOCATABLE :: ek_in(:,:)
  !! band structure of the phonons (nbnd,nqtot)
  integer, allocatable :: ii_tetra(:,:) ! (nntetra, ntetra)
  !! gives the VALID k point index for each VALID tetrahedron
  integer, allocatable :: iisize_tetra(:) ! (ntetra)
  !! number of VALID neighbouring k points in each VALID tetrahedron
  integer, allocatable :: nt_tetra(:) ! (ntetra)
  !! index of VALID tetrahedra
  integer, allocatable :: equiv(:)
  !! equivalence points in symmetrized grid
  integer :: nvalid
  !! number of VALID tetrahedra
  real(dp) :: MIN_DISTANCE
  !! if symmetry is used
  logical :: symmetry
  !! it helps the convergence of the real part bringing everything close to the unity
  real(dp) :: MULTIPLIER

  ! INTEGER, allocatable :: which_tetra(:,:,:)
  !! inverse of tetra: given a q point, it gives all the tetrahedra that contain it

  REAL(DP), PARAMETER :: tet_cutoff = 1.0E-3_DP
  REAL(DP), PARAMETER :: min_relative_distance = 1.0E-9_DP
  LOGICAL :: opt_flag
  !
  PUBLIC :: nntetra, tetra_weights_green, tetra_init_sym
  PUBLIC :: tetra_init, deallocate_tetra, tetra_weights_delta
  PUBLIC :: tetra_weights_delta_sym, rm_degen_vertices
  PUBLIC :: equiv_grid, ek_sort, nqtot, tetra_output
  PUBLIC :: set_wg, deallocate_tetra_output

  EXTERNAL :: errore, hpsort

  !
  type tetra_output
    complex(dp), allocatable :: w(:,:,:)
    !! tetrahedron weights for each q-point (ibnd, iq)
    integer :: ntot
    !! number of q-points in the whole BZ
    integer :: nsym
    !! number of q-points in the IRREDUCIBLE BZ
    integer, allocatable :: e(:)
    !! equivalence points in symmetrized grid
    real(dp), allocatable :: qw(:)
    !! q point weight, copied from the symmetrized grid
    real(dp), allocatable :: f(:,:)
    !! symmetrized frequencies
    real(dp), allocatable :: en(:)
    !! energies in the real axis
    real(dp) :: max_f
  end type tetra_output

CONTAINS
  !
  subroutine deallocate_tetra_output(wg)
    type(tetra_output), intent(inout) :: wg
    if (allocated(wg%w)) deallocate(wg%w)
    if (allocated(wg%e)) deallocate(wg%e)
    if (allocated(wg%qw)) deallocate(wg%qw)
    if (allocated(wg%f)) deallocate(wg%f)
    if (allocated(wg%en)) deallocate(wg%en)
  end subroutine
!
  subroutine tetra_init_grid_sym(grid, S, fc2, wg, n_omega, mult)
    type(q_grid), intent(in) :: grid
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(tetra_output), intent(out) :: wg
    integer, intent(in) :: n_omega
    real(dp), intent(in) :: mult
    ! type(q_grid), intent(out), optional :: grid_sym_
    !
    real(dp), allocatable :: freqs_sym(:,:)
    integer :: iw
    !
    allocate(freqs_sym(S%nat3, grid%nqtot))
    call freq_in_grid(S, fc2, grid, freqs_sym)
    call tetra_init_sym(grid, S, freqs_sym**2, .false., wg)
    allocate(wg%w(S%nat3, wg%nsym, n_omega))
    call move_alloc(freqs_sym, wg%f)
    wg%max_f = maxval(wg%f) * mult
    allocate(wg%en(n_omega))
    do iw = 1, n_omega
      wg%en(iw) = (iw-1) * wg%max_f / REAL(n_omega, dp)
    enddo
    !
  end subroutine
  !
  subroutine set_wg(S, fc2, grid, n_omega, wg, mult)
    use thutils, only : freq_in_grid
    use merge_degenerate, only: merge_degen
    !
    type(ph_system_info), intent(in) :: S
    type(forceconst2_grid), intent(in) :: fc2
    type(q_grid), intent(in) :: grid
    integer, intent(in) :: n_omega
    type(tetra_output), intent(out) :: wg
    real(dp), intent(in), optional :: mult
    !
    real(dp) :: mult_
    integer :: iq, ibnd, iw
    real(dp) :: freqs(S%nat3,grid%nqtot)
    !
    if(present(mult)) then
      mult_ = mult
    else
      mult_ = 1.2_dp
    end if
    !
    call tetra_init_grid_sym(grid, S, fc2, wg, n_omega, mult_)
    !
    do iw = 1, n_omega
      wg%w(:,:,iw) = tetra_weights_green(wg%en(iw)**2)
      do iq = 1, grid%nqtot
        do ibnd = 1, S%nat3
          if(isnan(ABS(wg%w(ibnd,wg%e(iq),iw)))) wg%w(ibnd,wg%e(iq),iw) = 0._dp
        enddo
      enddo
      do iq = 1, grid%nqtot
        call merge_degen(S%nat3, wg%w(:,iq,iw), wg%f(:,iq))
      enddo
    enddo
  end subroutine
  !
  subroutine equiv_grid(grid, S, equiv_, first_point)
    use symm_base, only : symms => s, nsym, time_reversal, t_rev
    USE noncollin_module,   ONLY : colin_mag
    use q_grids, only : q_grid, setup_grid
    use ph_system, only : ph_system_info
    !
    type(q_grid), intent(in) :: grid
    type(ph_system_info), intent(in) :: S
    integer, allocatable, intent(out) :: equiv_(:)
    integer, allocatable, optional, intent(out) :: first_point(:)
    !
    type(q_grid) :: full_grid
    integer :: ik, jk, isym, jkp
    real(dp), dimension(3) :: xkr, deltap, deltam
    REAL(DP), PARAMETER :: eps = 1e-5_dp
    !
    if (present(first_point)) then
      allocate(first_point(grid%nqtot))
      first_point = 0
    endif
    !
    if(grid%symmetrized) then
      call setup_grid(grid%type, S%bg, grid%n(1), grid%n(2), grid%n(3), &
        full_grid, xq0=grid%xq0, scatter = .false., quiet = .true.)
      nqtot = full_grid%nqtot
      allocate(equiv_(nqtot))
      equiv_ = 0
      CALL cryst_to_cart( grid%nq, grid%xq, S%at, -1 )
      call cryst_to_cart( nqtot, full_grid%xq, S%at, -1 )

      DO ik = 1, nqtot
        ! if (ik == 11) print"(3F8.2)", full_grid%xq(:,ik)
        DO jk = 1, grid%nq
          jkp = jk + grid%iq0
          DO isym = 1, nsym
            !
            xkr(1:3) = MATMUL(REAL(symms(:,:,isym), dp), grid%xq(:,jk))
            IF (t_rev(isym) == 1 .AND. colin_mag < 2) xkr(1:3) = - xkr(1:3)
            !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
            deltap = xkr - full_grid%xq(:,ik)
            deltap = deltap - NINT(deltap)
            deltam = xkr + full_grid%xq(:,ik)
            deltam = deltam - NINT(deltam)
            !  deltap is the difference vector, brought back in the first BZ
            !  deltam is the same but with k => -k (for time reversal)
            IF ( norm2(deltap) < eps .OR. ( time_reversal .AND. &
              norm2(deltam) < eps ) ) THEN
              !  equivalent irreducible k-point found
              equiv_(ik) = jkp
              if (present(first_point)) then
                if (first_point(jkp) == 0) first_point(jkp) = ik
              endif
              GOTO 15
            ENDIF
            !
          ENDDO
        ENDDO
        !  equivalent irreducible k-point found - something wrong
        if(.not. grid%scattered) CALL errore( 'opt_tetra_init', 'cannot locate  k point', ik )
        !
15      CONTINUE
        !
      ENDDO
      !
      DO jk = 1, grid%nq
        jkp = jk + grid%iq0
        DO ik = 1, product(grid%n)
          IF (equiv_(ik) == jkp) GOTO 20
        ENDDO
        !  this failure of the algorithm may indicate that the displaced grid
        !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
        print"(3F8.2)", grid%xq(:,jk)
        print*, grid%nqtot
        CALL errore( 'opt_tetra_init', 'cannot remap grid on k-point list', jkp )
        !
20      CONTINUE
      ENDDO
      !
      !  bring irreducible k-points back to cartesian axis
      !
      CALL cryst_to_cart( grid%nq, grid%xq, S%bg, 1 )
    else
      nqtot = grid%nqtot
      allocate(equiv_(grid%nqtot))
      equiv_ = 0
      do jk = 1, grid%nq
        jkp = jk + grid%iq0
        equiv_(jkp) = jkp
        if (present(first_point)) first_point(jkp) = jkp
      enddo
    endif
    if(grid%scattered) call mpi_bsum(nqtot, equiv_)
    if(grid%scattered) call mpi_bsum(grid%nqtot, first_point)
  end subroutine
  !--------------------------------------------------------------------------
  SUBROUTINE tetra_init_sym(grid, S, ek, opt, wg)
    !-----------------------------------------------------------------------------
    !! This rouotine sets the corners and additional points for each tetrahedron.
    !
    use ph_system, only : ph_system_info
    use q_grids,   only : q_grid, setup_grid
    IMPLICIT NONE
    !
    type(ph_system_info), intent(in) :: S
    !! usual system info
    REAL(DP), INTENT(IN) :: ek(:,:)
    !! energy in the form ek(ibnd, iq)
    type(q_grid), intent(in) :: grid
    !! accepts symmetrized grids
    logical, intent(in), optional :: opt
    !! if .true., uses opt_tetra methods
    type(tetra_output), intent(out) :: wg
    !! tetrahedron weights infos
    REAL(DP), PARAMETER :: eps = 1e-5_dp
    !
    INTEGER :: i1, i2, i3, itet, itettot, ii, ik,  rest, &
      ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3), ibnd, &
      itvalid, count, iks(20)
    ! integer :: tetra_ik(nq(1) * nq(2) * nq(3))
    !
    REAL(DP) :: l(4), bvec2(3,3), bvec3(3,4) !xkg(3, product(nq))
    integer, allocatable :: first_point(:)
    !
    IF(ntetra /= 0) CALL deallocate_tetra()
    !
    nbnd = SIZE(ek,1)
    nqs =  size(ek,2)
    !
    ! IF(nqs /= SIZE(ek,2)) then
    !   print*, "size of q grid in freq_", size(ek,2)
    !   print*, "size of original grid", nq
    !   CALL errore("tetra_init", "n(1) * n(2) * n(3) /= SIZE(ek,2)", SIZE(ek,2))
    ! ENDIF
    ntetra  = 6*product(grid%n)
    !
    ALLOCATE(ek_sort(4,nbnd,ntetra))
    ek_sort = 0.0_dp
    ALLOCATE(ek_in(nbnd,nqs))
    ek_in = ek
    ALLOCATE(itetra (4,nbnd,ntetra))
    allocate(iisize_tetra(ntetra))
    allocate(nt_tetra(ntetra))
    ! ALLOCATE(which_tetra (2,nqs,24))
    !
    opt_flag = .true.
    if (present(opt)) opt_flag = opt
    ! if(PRESENT(opt)) opt_flag = opt
    !
    ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
    !
    bvec2(1:3,1) = S%bg(:,1) / REAL(grid%n(1), dp)
    bvec2(1:3,2) = S%bg(:,2) / REAL(grid%n(2), dp)
    bvec2(1:3,3) = S%bg(:,3) / REAL(grid%n(3), dp)
    !
    bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
    bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
    bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
    bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
    !
    DO ii = 1, 4
      l(ii) = DOT_PRODUCT(bvec3(1:3, ii), bvec3(1:3, ii))
    ENDDO
    !
    ii = MINLOC(l(1:4),1)
    !
    ivvec0(1:4) = (/ 0, 0, 0, 0 /)
    !
    divvec(1:4,1) = (/ 1, 0, 0, 0 /)
    divvec(1:4,2) = (/ 0, 1, 0, 0 /)
    divvec(1:4,3) = (/ 0, 0, 1, 0 /)
    divvec(1:4,4) = (/ 0, 0, 0, 1 /)
    !
    ivvec0(ii) = 1
    divvec(ii, ii) = - 1
    !
    ! Divide a subcell into 6 tetrahedra
    !
    itet = 0
    DO i1 = 1, 3
      DO i2 = 1, 3
        IF(i2 == i1) CYCLE
        DO i3 = 1, 3
          IF(i3 == i1 .OR. i3 == i2) CYCLE
          !
          itet = itet + 1
          !
          ivvec(1:3,1,itet) = ivvec0(1:3)
          ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
          ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
          ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
          !
        ENDDO
      ENDDO
    ENDDO
    !
    ! Additional points surrounding the tetrahedron
    !
    ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
    ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
    ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
    ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
    !
    ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
    ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
    ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
    ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
    !
    ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
    ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
    ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
    ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
    !
    ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
    ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
    ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
    ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
    !
    ! Set the weight for the each tetrahedron method
    !
    ! WRITE(stdout,*) "    [opt_tetra]  Optimized tetrahedron method is used."
    !
    IF (opt_flag) THEN
      !
      nntetra = 20
      allocate(ii_tetra(nntetra, ntetra))
      IF (.NOT. ALLOCATED(tetra)) ALLOCATE( tetra(nntetra,ntetra) )
      IF (.NOT. ALLOCATED(wlsm))  ALLOCATE( wlsm(4,nntetra) )
      !
      wlsm(1, 1: 4) = REAL((/1440,    0,   30,    0/), dp)
      wlsm(2, 1: 4) = REAL((/   0, 1440,    0,   30/), dp)
      wlsm(3, 1: 4) = REAL((/  30,    0, 1440,    0/), dp)
      wlsm(4, 1: 4) = REAL((/   0,   30,    0, 1440/), dp)
      !
      wlsm(1, 5: 8) = REAL((/ -38,    7,   17,  -28/), dp)
      wlsm(2, 5: 8) = REAL((/ -28,  -38,    7,   17/), dp)
      wlsm(3, 5: 8) = REAL((/  17,  -28,  -38,    7/), dp)
      wlsm(4, 5: 8) = REAL((/   7,   17,  -28,  -38/), dp)
      !
      wlsm(1, 9:12) = REAL((/ -56,    9,  -46,    9/), dp)
      wlsm(2, 9:12) = REAL((/   9,  -56,    9,  -46/), dp)
      wlsm(3, 9:12) = REAL((/ -46,    9,  -56,    9/), dp)
      wlsm(4, 9:12) = REAL((/   9,  -46,    9,  -56/), dp)
      !
      wlsm(1,13:16) = REAL((/ -38,  -28,   17,    7/), dp)
      wlsm(2,13:16) = REAL((/   7,  -38,  -28,   17/), dp)
      wlsm(3,13:16) = REAL((/  17,    7,  -38,  -28/), dp)
      wlsm(4,13:16) = REAL((/ -28,   17,    7,  -38/), dp)
      !
      wlsm(1,17:20) = REAL((/ -18,  -18,   12,  -18/), dp)
      wlsm(2,17:20) = REAL((/ -18,  -18,  -18,   12/), dp)
      wlsm(3,17:20) = REAL((/  12,  -18,  -18,  -18/), dp)
      wlsm(4,17:20) = REAL((/ -18,   12,  -18,  -18/), dp)
      !
      wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260.0_dp
      !
    ELSE
      !
      nntetra = 4
      allocate(ii_tetra(nntetra, ntetra))
      IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra(nntetra,ntetra) )
      IF(.NOT. ALLOCATED(wlsm))  ALLOCATE ( wlsm(4,nntetra) )
      wlsm(:,:) = 0.0_dp
      !
      wlsm(1,1) = 1.0_dp
      wlsm(2,2) = 1.0_dp
      wlsm(3,3) = 1.0_dp
      wlsm(4,4) = 1.0_dp
      !
    ENDIF
    !  locate k-points of the uniform grid in the list of irreducible k-points
    !  that was previously calculated
    !
    !  bring irreducible k-points to crystal axis
    !
    ! nqtot = nqs
    symmetry = grid%symmetrized
    call equiv_grid(grid, S, equiv, first_point)
    itvalid = 0
    itetra = 0
    tetra = 0
    ! tetra_ik = 0
    DO itettot = 1+my_id, ntetra, num_procs
      itet = mod(itettot,6) + 1
      rest = itettot / 6
      i3 = mod(rest,grid%n(3)) + 1
      rest = rest / grid%n(3)
      i2 = mod(rest,grid%n(2)) + 1
      rest = rest / grid%n(2)
      i1 = mod(rest,grid%n(1)) + 1

      count = 0
      ! DO ibnd = 1, nbnd
      !
      DO ii = 1, nntetra
        !
        ikv(1:3) = (/i1, i2, i3/) - 1
        ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
        ikv(1:3) = MODULO(ikv(1:3), (/grid%n(1), grid%n(2), grid%n(3)/))
        !
        ik = ikv(3) + grid%n(3) * (ikv(2) + grid%n(2) * ikv(1)) + 1
        !
        iks(ii) = ik
        if (any(ik == first_point)) then
          count = count + 1
          ii_tetra(count,itvalid+1) = ii
        endif
        !
        ! tetra(ii, itettot) = equiv(ik)
        ! !
        ! ek_sort(:,ibnd,itettot) = ek_sort(:,ibnd,itettot) + wlsm(:,ii) * ek(ibnd,equiv(ik))

      END DO ! ii
      !

      ! ENDDO ! ibnd
      !
      if (count > 0) then
        itvalid = itvalid + 1
        nt_tetra(itvalid) = itettot
        iisize_tetra(itvalid) = count
        do ibnd = 1, nbnd
          do ii = 1, nntetra
            ik = iks(ii)
            tetra(ii, itvalid) = equiv(ik)
            ek_sort(:,ibnd,itvalid) = ek_sort(:,ibnd,itvalid) + &
              wlsm(:,ii) * ek(ibnd,equiv(ik))
          enddo
          itetra(1,ibnd,itvalid) = 0 ! needed to initialize index inside hpsort
          CALL hpsort( 4, ek_sort(:,ibnd,itvalid), itetra(:,ibnd,itvalid))
        enddo
      endif
    ENDDO ! itettot
    !
    MULTIPLIER = SUM(ek_sort)/REAL(SIZE(ek_sort), dp)
    ! print*, "multiplier", MULTIPLIER
    ek_sort = ek_sort / MULTIPLIER
    MIN_DISTANCE = MULTIPLIER * min_relative_distance
    nvalid = itvalid
    ! print*, "number of tetra to be used", nvalid, "out of", ntetra
    ! print*, "tetra", nqtot, nqs, symmetry
    wg%ntot = nqtot
    wg%nsym = nqs
    allocate(wg%e(nqs))
    wg%e = equiv
    allocate(wg%qw(nqs))
    wg%qw = grid%w
    ! print"(A,I4.1,A)", "grid has ", nqs, " q-points in the IRREDUCIBLE BZ"
  END SUBROUTINE
  !
  SUBROUTINE tetra_init(nq, bg, ek, opt)
    !-----------------------------------------------------------------------------
    !! This rouotine sets the corners and additional points for each tetrahedron.
    !
    use ph_system, only : ph_system_info
    use q_grids,   only : q_grid, setup_grid
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nq(3)
    !! number of q-points in each direction
    ! type(ph_system_info), intent(in) :: S
    !! usual system info
    REAL(DP), INTENT(IN) :: ek(:,:)
    !! energy in the form ek(ibnd, iq)
    ! type(q_grid), intent(in) :: grid
    real(dp) :: bg(3,3)
    !! accepts symmetrized grids
    logical, intent(in), optional :: opt

    ! LOGICAL, INTENT(IN) :: is_mpi
    !! if .true., the grid is scattered
    ! LOGICAL, INTENT(IN), OPTIONAL :: opt
    ! !! if .true., uses opt_tetra methods

    REAL(DP), PARAMETER :: eps = 1e-5_dp
    !
    INTEGER :: i1, i2, i3, itet, itettot, ii, ik,  &
      ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3), ibnd, &
      rest
    ! integer :: tetra_ik(nq(1) * nq(2) * nq(3))
    !
    REAL(DP) :: l(4), bvec2(3,3), bvec3(3,4) !xkg(3, product(nq))
    !
    IF(ntetra /= 0) CALL deallocate_tetra()
    !
    nbnd = SIZE(ek,1)
    nqs = product(nq)
    nqtot = nqs
    !
    ! IF(nqs /= SIZE(ek,2)) then
    !   print*, "size of q grid in freq_", size(ek,2)
    !   print*, "size of original grid", nq
    !   CALL errore("tetra_init", "n(1) * n(2) * n(3) /= SIZE(ek,2)", SIZE(ek,2))
    ! ENDIF
    ntetra  = 6*nqs
    !
    ALLOCATE(ek_sort(4,nbnd,ntetra))
    ALLOCATE(ek_in(nbnd,nqs))
    ek_in = ek
    ALLOCATE(itetra (4,nbnd,ntetra))
    allocate(iisize_tetra(ntetra))
    allocate(nt_tetra(ntetra))
    ! ALLOCATE(which_tetra (2,nqs,24))
    !
    opt_flag = .true.
    if (present(opt)) opt_flag = opt
    ! if(PRESENT(opt)) opt_flag = opt
    !
    ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
    !
    bvec2(1:3,1) = bg(:,1) / REAL(nq(1), dp)
    bvec2(1:3,2) = bg(:,2) / REAL(nq(2), dp)
    bvec2(1:3,3) = bg(:,3) / REAL(nq(3), dp)
    !
    bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
    bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
    bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
    bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
    !
    DO ii = 1, 4
      l(ii) = DOT_PRODUCT(bvec3(1:3, ii), bvec3(1:3, ii))
    ENDDO
    !
    ii = MINLOC(l(1:4),1)
    !
    ivvec0(1:4) = (/ 0, 0, 0, 0 /)
    !
    divvec(1:4,1) = (/ 1, 0, 0, 0 /)
    divvec(1:4,2) = (/ 0, 1, 0, 0 /)
    divvec(1:4,3) = (/ 0, 0, 1, 0 /)
    divvec(1:4,4) = (/ 0, 0, 0, 1 /)
    !
    ivvec0(ii) = 1
    divvec(ii, ii) = - 1
    !
    ! Divide a subcell into 6 tetrahedra
    !
    itet = 0
    DO i1 = 1, 3
      DO i2 = 1, 3
        IF(i2 == i1) CYCLE
        DO i3 = 1, 3
          IF(i3 == i1 .OR. i3 == i2) CYCLE
          !
          itet = itet + 1
          !
          ivvec(1:3,1,itet) = ivvec0(1:3)
          ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
          ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
          ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
          !
        ENDDO
      ENDDO
    ENDDO
    !
    ! Additional points surrounding the tetrahedron
    !
    ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
    ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
    ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
    ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
    !
    ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
    ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
    ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
    ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
    !
    ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
    ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
    ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
    ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
    !
    ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
    ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
    ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
    ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
    !
    ! Set the weight for the each tetrahedron method
    !
    ! WRITE(stdout,*) "    [opt_tetra]  Optimized tetrahedron method is used."
    !
    IF (opt_flag) THEN
      !
      nntetra = 20
      IF (.NOT. ALLOCATED(tetra)) ALLOCATE( tetra(nntetra,ntetra) )
      IF (.NOT. ALLOCATED(wlsm))  ALLOCATE( wlsm(4,nntetra) )
      !
      wlsm(1, 1: 4) = REAL((/1440,    0,   30,    0/), dp)
      wlsm(2, 1: 4) = REAL((/   0, 1440,    0,   30/), dp)
      wlsm(3, 1: 4) = REAL((/  30,    0, 1440,    0/), dp)
      wlsm(4, 1: 4) = REAL((/   0,   30,    0, 1440/), dp)
      !
      wlsm(1, 5: 8) = REAL((/ -38,    7,   17,  -28/), dp)
      wlsm(2, 5: 8) = REAL((/ -28,  -38,    7,   17/), dp)
      wlsm(3, 5: 8) = REAL((/  17,  -28,  -38,    7/), dp)
      wlsm(4, 5: 8) = REAL((/   7,   17,  -28,  -38/), dp)
      !
      wlsm(1, 9:12) = REAL((/ -56,    9,  -46,    9/), dp)
      wlsm(2, 9:12) = REAL((/   9,  -56,    9,  -46/), dp)
      wlsm(3, 9:12) = REAL((/ -46,    9,  -56,    9/), dp)
      wlsm(4, 9:12) = REAL((/   9,  -46,    9,  -56/), dp)
      !
      wlsm(1,13:16) = REAL((/ -38,  -28,   17,    7/), dp)
      wlsm(2,13:16) = REAL((/   7,  -38,  -28,   17/), dp)
      wlsm(3,13:16) = REAL((/  17,    7,  -38,  -28/), dp)
      wlsm(4,13:16) = REAL((/ -28,   17,    7,  -38/), dp)
      !
      wlsm(1,17:20) = REAL((/ -18,  -18,   12,  -18/), dp)
      wlsm(2,17:20) = REAL((/ -18,  -18,  -18,   12/), dp)
      wlsm(3,17:20) = REAL((/  12,  -18,  -18,  -18/), dp)
      wlsm(4,17:20) = REAL((/ -18,   12,  -18,  -18/), dp)
      !
      wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260.0_dp
      !
    ELSE
      !
      nntetra = 4
      IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra(nntetra,ntetra) )
      IF(.NOT. ALLOCATED(wlsm))  ALLOCATE ( wlsm(4,nntetra) )
      wlsm(:,:) = 0.0_dp
      !
      wlsm(1,1) = 1.0_dp
      wlsm(2,2) = 1.0_dp
      wlsm(3,3) = 1.0_dp
      wlsm(4,4) = 1.0_dp
      !
    ENDIF
    ek_sort = 0._dp
    itetra = 0
    tetra = 0
    DO itettot = 1+my_id, ntetra, num_procs
      itet = mod(itettot,6) + 1
      rest = itettot / 6
      i3 = mod(rest,nq(3)) + 1
      rest = rest / nq(3)
      i2 = mod(rest,nq(2)) + 1
      rest = rest / nq(2)
      i1 = mod(rest,nq(1)) + 1
      !
      DO ibnd = 1, nbnd
        DO ii = 1, nntetra
          !
          ikv(1:3) = (/i1, i2, i3/) - 1
          ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
          ikv(1:3) = MODULO(ikv(1:3), nq)
          !
          ik = ikv(3) + nq(3) * (ikv(2) + nq(2) * ikv(1)) + 1
          !
          tetra(ii, itettot) = ik
          !
          ek_sort(:,ibnd,itettot) = ek_sort(:,ibnd,itettot) + wlsm(:,ii) * ek(ibnd,ik)
        enddo
        itetra(1,ibnd,itettot) = 0 ! needed to initialize index inside hpsort
        CALL hpsort( 4, ek_sort(:,ibnd,itettot), itetra(:,ibnd,itettot))
      END DO
      !
    ENDDO !
    call mpi_bsum(nntetra, ntetra, tetra)
    call mpi_bsum(4, nbnd, ntetra, ek_sort)
    call mpi_bsum(4, nbnd, ntetra, itetra)
    !
    MULTIPLIER = SUM(ek_sort)/REAL(SIZE(ek_sort), dp)
    ! print*, "multiplier", MULTIPLIER
    ek_sort = ek_sort / MULTIPLIER
    MIN_DISTANCE = MULTIPLIER * min_relative_distance
  END SUBROUTINE
  !
  subroutine tetra_weights_delta(ef, wI)
    USE constants, ONLY : pi
    !-----------------------------------------------------------------------------------
    !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
    !! The resulting wg can be used as sum(Ak * wk)
    !-----------------------------------------------------------------------------------
    REAL(DP), intent(out) :: wI(nbnd, nqs)
    !! COMPLEX Integration weight of each k
    REAL(DP), INTENT(IN) :: ef
    !! The Fermi energy
    INTEGER :: ik, ibnd, ii, it, jbnd, kbnd
    REAL(DP) :: e(4), wI0(4), wg1
    !
    wI = 0._dp
    !
    DO it = 1+my_id, ntetra, num_procs
      !
      ! nt = nt_tetra(it)
      !
      DO ibnd = 1, nbnd
        !
        e = ek_sort(:,ibnd,it)
        ! print"(4E12.4)", e, ef
        ! D = -e
        ! call rm_degen_vertices(ef, e)
        ! e = -D
        wI0 = delta_vertices(ef, e)
        ! if(ef < maxval(e) .and. ef > minval(e)) print*, "non va"
        ! if(any(wi0 > 0)) print"(4E12.4)", e, ef
        !
        DO ii = 1, nntetra
          !
          ik = tetra(ii, it)
          ! IF(opt_flag) THEN
          ! if(ik > nqs) cycle
          wI(ibnd,ik) = wI(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,it),ii), wI0(:))
          ! ELSE
          !   ik_s = tetra(itetra(ii,ibnd,nt), nt)
          !   wI(ibnd,ik_s) = wI(ibnd,ik_s) + wI0(ii)
          !   wR(ibnd,ik_s) = wR(ibnd,ik_s) + wR0(ii)
          ! ENDIF
        ENDDO
        !
      ENDDO ! ibnd
      !
    ENDDO ! nt
    ! wg = wg / REAL(ntetra, dp)
    wI = wI / ntetra
    !
    ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
    CALL mpi_bsum(nbnd, nqs, wI)
    !
    DO ik = 1, nqs
      DO ibnd = 1, nbnd
        !
        wg1 = wI(ibnd,ik)
        !
        DO jbnd = ibnd + 1, nbnd
          !
          IF (ABS(ek_in(ibnd,ik) - ek_in(jbnd,ik)) < MIN_DISTANCE) THEN
            wg1 = wg1 + wI(jbnd,ik)
          ELSE
            !
            DO kbnd = ibnd, jbnd - 1
              wI(kbnd,ik) = wg1 / REAL(jbnd - ibnd, dp)
            ENDDO
            !
            EXIT
          ENDIF
          !
        ENDDO
        !
      ENDDO
    ENDDO
    !
  END subroutine
  !
  subroutine tetra_weights_delta_sym(ef, wI)
    USE constants, ONLY : pi
    !-----------------------------------------------------------------------------------
    !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
    !! The resulting wg can be used as sum(Ak * wk)
    !-----------------------------------------------------------------------------------
    REAL(DP), intent(out) :: wI(nbnd, nqs)
    !! COMPLEX Integration weight of each k
    REAL(DP), INTENT(IN) :: ef
    !! The Fermi energy
    INTEGER :: ik, ibnd, ii_, ii, it
    REAL(DP) :: e(4), wI0(4)
    !
    wI = 0._dp
    !
    DO it = 1+my_id, nvalid, num_procs
      !
      ! nt = nt_tetra(it)
      !
      DO ibnd = 1, nbnd
        !
        e = ek_sort(:,ibnd,it)
        ! print"(4E12.4)", e, ef
        wI0 = delta_vertices(ef, e)
        ! if(any(wi0 > 0)) print"(4E12.4)", e, ef
        !
        DO ii_ = 1, iisize_tetra(it)
          !
          ii = ii_tetra(ii_, it)
          ! ii = ii_
          ik = tetra(ii, it)
          ! IF(opt_flag) THEN
          ! if(ik > nqs) cycle
          wI(ibnd,ik) = wI(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,it),ii), wI0(:))
          ! ELSE
          !   ik_s = tetra(itetra(ii,ibnd,nt), nt)
          !   wI(ibnd,ik_s) = wI(ibnd,ik_s) + wI0(ii)
          !   wR(ibnd,ik_s) = wR(ibnd,ik_s) + wR0(ii)
          ! ENDIF
        ENDDO
        !
      ENDDO ! ibnd
      !
    ENDDO ! nt
    ! wg = wg / REAL(ntetra, dp)
    wI = wI / (6.0_dp * nqs)
    !
    ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
    CALL mpi_bsum(nbnd, nqs, wI)
    !
  END subroutine
  !
  subroutine rm_degen_vertices(hwe, D)
    real(dp), INTENT(IN) :: hwe
    real(dp), INTENT(INOUT) :: D(4)
    !
    real(dp) :: DAV, D_small_prev, D_large_prev, hw
    !
    hw = hwe
    if(any(D + hwe == 0._dp)) hw = hw + 1e-10_dp
    DAV = (D(2) + D(3))/2.0_dp
    if ((D(3) - D(2))/(DAV + hw) < tet_cutoff) then
      D_small_prev = D(2); D_large_prev = D(3)
      D(3) = DAV + 0.5_dp*abs(DAV + hw)*tet_cutoff
      D(2) = DAV - 0.5_dp*abs(DAV + hw)*tet_cutoff
      if (D(1) > D(2)) D(1) = D(1) + (D(2) - D_small_prev)
      if (D(3) > D(4)) D(4) = D(4) + (D(3) - D_large_prev)
    endif
    DAV = (D(1) + D(2))/2.0_dp
    if ((D(2) - D(1))/(DAV + hw) < tet_cutoff) then
      D(1) = DAV - 0.5_dp*abs(DAV + hw)*tet_cutoff
    endif
    DAV = (D(3) + D(4))/2.0_dp
    if ((D(4) - D(3))/(DAV + hw) < tet_cutoff) then
      D(4) = DAV + 0.5_dp*abs(DAV + hw)*tet_cutoff
    endif
  end subroutine
  !
  FUNCTION real_vertices(hwe, D) result(wR0)
    !
    real(dp), INTENT(IN) :: hwe
    real(dp), INTENT(IN) :: D(4)
    !
    real(dp) :: wR0(4), hw
    real(dp) :: dd(3), ll(3), ff, bb(4), cc(4, 3)
    integer :: a, b, c, i
    !
    ! intermediate variables for case 1 and 3
    !
    hw = hwe
    if(any(D + hw == 0._dp)) hw = 1e-5_dp
    do i = 1, 3
      dd(i) = (D(4) - D(i))/(D(i) + hw)
      ll(i) = tetrahedron_log1p(dd(i))
    enddo
    !
    ff = 1.0_dp
    do i = 1, 3
      a = i
      b = mod(i, 3) + 1
      c = mod(i + 1, 3) + 1
      cc(a, a) = -(1.0_dp + dd(a))*(3.0_dp*dd(a)**2 - 2.0_dp*(dd(b) + dd(c))*dd(a) + dd(b)*dd(c)) &
        *((dd(b) - dd(c))*dd(b)*dd(c))**2
      cc(b, a) = -dd(a)*(1.0_dp + dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
      cc(c, a) = dd(a)*(1.0_dp + dd(c))*(dd(a) - dd(b))*((dd(b) - dd(c))*dd(b)*dd(c))**2
      cc(4, a) = -(dd(a) - dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
      bb(a) = cc(4, a)*dd(a)
      ff = ff*(1.0_dp + dd(a))/(dd(a)*(dd(a) - dd(b)))**2
    enddo
    bb(4) = -dd(1)*dd(2)*dd(3)*((dd(1) - dd(2))*(dd(2) - dd(3))*(dd(3) - dd(1)))**2
    ff = -ff/(D(4) + hw)
    !
    do i = 1, 4
      wR0(i) = (cc(i, 1)*ll(1) + cc(i, 2)*ll(2) + cc(i, 3)*ll(3) + bb(i))
    enddo
    !
    wR0 = wR0 * ff
    if(any(abs(wR0)>1e10_dp)) then
      print*, hw, D
      print*, wR0
      print*, dd
      print*, ll
      print*, bb
      print*, cc
      print*, ff
      call errore("real_vertices", "weight too large", 1)
    endif
  END FUNCTION

  FUNCTION delta_vertices(ef, e) result(wI0)
    !
    real(dp), INTENT(IN) :: ef
    real(dp) :: e(4)
    !
    real(dp) :: wI0(4)
    !
    real(dp) :: C, a(4,4)
    !
    integer :: i, ii
    logical :: near(4,4)
    !
    !
    !
    IF(ef > e(4) .or. ef < e(1)) THEN
      wI0 = 0.0_dp
      RETURN
    ENDIF
    !
    DO ii = 1, 4
      DO i = 1, 4
        near(ii,i) = ABS(e(i)-e(ii)) < MIN_DISTANCE
        IF ( e(ii) == e(i) ) then
          a(ii,i) = 0.0_dp
        else
          a(ii,i) = ( ef - e(i) ) / (e(ii) - e(i) )
        ENDIF
      ENDDO
    ENDDO
    !
    !
    ! if((near(1,2) .and. near(3,4)) .or. &
    !   (near(1,2) .and. near(2,3)) .or. &
    !   (near(2,3) .and. near(3,4))) then
    !   call rm_degen_vertices(ef, e)
    ! endif
    !
    IF( e(1) <= ef .AND. ef <= e(2) ) THEN
      !
      C = a(2,1) * a(3,1)
      wI0(1) = a(1,2) + a(1,3) + a(1,4)
      wI0(2:4) = a(2:4,1)

      wI0 = wI0 * C
      IF (near(1,2)) wI0 = 0.0_dp
      !
    ELSEIF( e(2) < ef .AND. ef <= e(3)) THEN
      !
      C = a(2,3) * a(3,1) + a(3,2) * a(2,4)
      !
      wI0(1) = a(1,4) * C + a(1,3) * a(3,1) * a(2,3)
      wI0(2) = a(2,3) * C + a(2,4)**2 * a(3,2)
      wI0(3) = a(3,2) * C + a(3,1)**2 * a(2,3)
      wI0(4) = a(4,1) * C + a(4,2) * a(2,4) * a(3,2)
      !
      if (near(2,3)) then
        wI0 = 0.0_dp
        wI0(2:3) = 1._dp
      endif
    ELSEIF ( e(3) < ef .AND. ef <= e(4)) THEN
      !
      C = a(2,4) * a(3,4)
      !
      wI0(1:3) = a(1:3,4)
      wI0(4) = a(4,1) + a(4,2) + a(4,3)
      !
      wI0 = wI0 * C
      IF (near(3,4)) wI0 = 0.0_dp
      !
    ENDIF

    wI0 = wI0 / (e(4) - e(1))
    !
  END FUNCTION
  !
  !
  !
  FUNCTION tetra_weights_green(ef) RESULT(wg_sym)
    USE constants, ONLY : pi
    !-----------------------------------------------------------------------------------
    !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
    !! The resulting wg can be used as sum(Ak * wk)
    !-----------------------------------------------------------------------------------
    !! COMPLEX Integration weight of each k
    REAL(DP), INTENT(IN) :: ef
    !! The Fermi energy
    !
    ! ... local variables
    !
    real(dp) :: D(4)
    !! end wannier90 tetra

    REAL(DP) :: wI(nbnd,nqs), wR(nbnd, nqs)
    complex(dp) :: wg_sym(nbnd, nqs)
    real(dp) :: ef_mult

    INTEGER :: ik, nt, ibnd, ii, ntmax, iimax, ii_
    REAL(DP) :: e(4), wI0(4), wR0(4)

    ! for real part calc
    ! REAL(DP) :: wR0(4), ef_e(4), log_ef_e(4), prod_a(4), sum_a(4), second_term(4)
    ! INTEGER :: i3, j3
    !
    ef_mult = ef / MULTIPLIER
    wg_sym = 0._dp
    wI = 0._dp
    wR = 0._dp
    !
    if (symmetry) then
      ntmax = nvalid
    else
      ntmax = ntetra
    endif
    !
    DO nt = 1, ntmax
      !
      DO ibnd = 1, nbnd
        !
        e = ek_sort(:,ibnd,nt)
        CALL rm_degen_vertices(ef_mult, e)
        !
        wI0 = delta_vertices(ef_mult, e)
        D = -e
        wR0 = real_vertices(ef_mult, D)
        !
        !
        if(symmetry) then
          iimax = iisize_tetra(nt)
        else
          iimax = nntetra
        endif
        !
        DO ii_ = 1, iimax
          !
          if(symmetry) then
            ii = ii_tetra(ii_, nt)
          else
            ii = ii_
          endif
          ! ii = ii_
          ik = tetra(ii, nt)
          ! IF(opt_flag) THEN
          wI(ibnd,ik) = wI(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,nt),ii), wI0(1:4))
          wR(ibnd,ik) = wR(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,nt),ii), WR0(1:4))
          ! ELSE
          !   ik_s = tetra(itetra(ii,ibnd,nt), nt)
          !   wI(ibnd,ik_s) = wI(ibnd,ik_s) + wI0(ii)
          !   wR(ibnd,ik_s) = wR(ibnd,ik_s) + wR0(ii)
          ! ENDIF
        ENDDO
        !
      ENDDO ! ibnd
      !
    ENDDO ! nt
    ! wg = wg / REAL(ntetra, dp)
    wg_sym = CMPLX(wR, -pi*wI, kind=DP) / (ntetra * MULTIPLIER)
    !
    ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
    CALL mpi_bsum(nbnd, nqs, wg_sym)
    !
  END FUNCTION

  ! subroutine sort(list)
  !   !! Swap sort list of reals

  !   real(DP), intent(inout) :: list(:)
  !   real(DP) :: aux, tmp
  !   integer :: i, j, n

  !   n = size(list)

  !   do i = 1, n
  !     aux = list(i)
  !     do j = i + 1, n
  !       if (aux > list(j)) then
  !         tmp = list(j)
  !         list(j) = aux
  !         list(i) = tmp
  !         aux = tmp
  !       end if
  !     end do
  !   end do
  ! end subroutine sort
  !
  !--------------------------------------------------------------------
  SUBROUTINE deallocate_tetra( )
    !--------------------------------------------------------
    !! Deallocate tetra and wlsm
    !
    ntetra = 0
    nntetra = 0
    nqs = 0
    nbnd = 0
    nvalid = 0
    nqtot = 0
    IF (ALLOCATED(tetra  ))      DEALLOCATE (tetra  )
    IF (ALLOCATED(wlsm   ))      DEALLOCATE (wlsm   )
    IF (ALLOCATED(ek_sort))      DEALLOCATE (ek_sort)
    IF (ALLOCATED(itetra ))      DEALLOCATE (itetra )
    if (allocated(ii_tetra))     deallocate (ii_tetra)
    if (allocated(ek_in))        deallocate (ek_in   )
    if (allocated(nt_tetra))     deallocate (nt_tetra)
    if (allocated(iisize_tetra)) deallocate (iisize_tetra)
    if (allocated(equiv))        deallocate (equiv)
    !
  END SUBROUTINE deallocate_tetra
  !
  !
  ! FUNCTION tetra_weights_theta( nks, nbnd, et, ef) RESULT(wg)
  !   !--------------------------------------------------------------------
  !   !! Calculates weights with the tetrahedron method (P.E.Bloechl).
  !   !! Fermi energy has to be calculated in previous step.
  !   !! Generalization to noncollinear case courtesy of Iurii Timrov.
  !   !! @Note (P. Delugas 8/10/2019) Needs to be called only after initializations,
  !   !!       stops the program with an error call otherwise.
  !   !!
  !   !
  !   USE kinds
  !   !
  !   IMPLICIT NONE
  !   !
  !   INTEGER, INTENT(IN) :: nks
  !   !! Total # of k in irreducible BZ
  !   INTEGER, INTENT(IN) :: nbnd
  !   !! number of bands
  !   REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !   !! eigenvalues of the hamiltonian
  !   REAL(DP) :: wg(nbnd,nks)
  !   !! the weight of each k point and band
  !   ! wg must be (inout) and not (out) because if is/=0 only terms for
  !   ! spin=is are initialized; the remaining terms should be kept, not lost.
  !   REAL(DP), INTENT(IN) :: ef
  !   !! Fermi energy

  !   ! ... local variables
  !   !
  !   REAL(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra(4), dosef
  !   INTEGER :: ibnd, nt, nk, i, kp1, kp2, kp3, kp4, itetra(4)
  !   !
  !   nk = 0
  !   wg = 0._dp
  !   !
  !   DO nt = 1, ntetra
  !     DO ibnd = 1, nbnd
  !       !
  !       ! etetra are the energies at the vertexes of the nt-th tetrahedron
  !       !
  !       DO i = 1, 4
  !         etetra(i) = et (ibnd, tetra(i,nt) + nk)
  !       ENDDO
  !       itetra (1) = 0
  !       CALL hpsort( 4, etetra, itetra )
  !       !
  !       ! ...sort in ascending order: e1 < e2 < e3 < e4
  !       !
  !       e1 = etetra(1)
  !       e2 = etetra(2)
  !       e3 = etetra(3)
  !       e4 = etetra(4)
  !       !
  !       ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
  !       !
  !       kp1 = tetra(itetra(1), nt) + nk
  !       kp2 = tetra(itetra(2), nt) + nk
  !       kp3 = tetra(itetra(3), nt) + nk
  !       kp4 = tetra(itetra(4), nt) + nk
  !       !
  !       ! calculate weights wg
  !       !
  !       IF (ef>=e4) THEN
  !         !
  !         wg(ibnd, kp1) = wg(ibnd, kp1) + 0.25d0 / ntetra
  !         wg(ibnd, kp2) = wg(ibnd, kp2) + 0.25d0 / ntetra
  !         wg(ibnd, kp3) = wg(ibnd, kp3) + 0.25d0 / ntetra
  !         wg(ibnd, kp4) = wg(ibnd, kp4) + 0.25d0 / ntetra
  !         !
  !       ELSEIF (ef<e4 .AND. ef>=e3) THEN
  !         !
  !         c4 = 0.25d0 / ntetra * (e4 - ef)**3 / (e4 - e1) / (e4 - e2) &
  !           / (e4 - e3)
  !         dosef = 3.d0 / ntetra * (e4 - ef)**2 / (e4 - e1) / (e4 - e2) &
  !           / (e4 - e3)
  !         wg(ibnd,kp1) = wg(ibnd,kp1) + 0.25d0 / ntetra - c4 * &
  !           (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
  !           (ibnd, kp1) ) / 40.d0
  !         wg(ibnd,kp2) = wg(ibnd,kp2) + 0.25d0 / ntetra - c4 * &
  !           (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
  !           (ibnd, kp2) ) / 40.d0
  !         wg(ibnd,kp3) = wg(ibnd,kp3) + 0.25d0 / ntetra - c4 * &
  !           (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
  !           (ibnd, kp3) ) / 40.d0
  !         wg(ibnd,kp4) = wg(ibnd,kp4) + 0.25d0 / ntetra - c4 * &
  !           (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
  !           + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
  !           et(ibnd,kp4) ) / 40.d0
  !         !
  !       ELSEIF (ef<e3 .AND. ef>=e2) THEN
  !         !
  !         c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
  !         c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
  !           / (e4 - e1) / (e3 - e2) / (e3 - e1)
  !         c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
  !           / (e3 - e2) / (e4 - e1)
  !         dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
  !           (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
  !           * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
  !         wg(ibnd, kp1) = wg(ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
  !           / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
  !           (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
  !         wg(ibnd, kp2) = wg(ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
  !           * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
  !           (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
  !         wg(ibnd, kp3) = wg(ibnd, kp3) + (c1 + c2) * (ef - e1) &
  !           / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
  !           (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
  !         wg(ibnd, kp4) = wg(ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
  !           / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
  !           e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
  !         !
  !       ELSEIF (ef<e2 .AND. ef>=e1) THEN
  !         !
  !         c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
  !           / (e4 - e1)
  !         dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
  !           / (e4 - e1)
  !         wg(ibnd, kp1) = wg(ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
  !           * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
  !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
  !         wg(ibnd, kp2) = wg(ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
  !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
  !         wg(ibnd, kp3) = wg(ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
  !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
  !         wg(ibnd, kp4) = wg(ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
  !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0

  !         ! c4 = (ef-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)/4/ntetra
  !         ! wg(ibnd,kp1) = wg(ibnd,kp1) + (1 + (ef-e2)/(e1-e2) + (ef-e3)/(e1-e3) + (ef-e4)/(e1-e4))*c4
  !         ! wg(ibnd,kp2) = wg(ibnd,kp2) + (ef-e1)/(e2-e1)*c4
  !         ! wg(ibnd,kp3) = wg(ibnd,kp3) + (ef-e1)/(e3-e1)*c4
  !         ! wg(ibnd,kp4) = wg(ibnd,kp4) + (ef-e1)/(e4-e1)*c4
  !       ENDIF
  !       !
  !     ENDDO
  !   ENDDO
  !   !
  ! END FUNCTION tetra_weights_theta

  PURE function tetrahedron_log1p(x)
    implicit none

    real(dp) :: tetrahedron_log1p
    real(dp), intent(in) :: x
    real(dp) :: y, z

    if (ABS(x) > 0.5_dp) then
      tetrahedron_log1p = LOG(ABS(1.0_dp + x))
    else
      y = 1.0_dp + x
      z = y - 1.0_dp
      if (z == 0) then
        tetrahedron_log1p = x
      else
        tetrahedron_log1p = x*LOG(y)/z
      endif
    endif

  end function tetrahedron_log1p

END MODULE thtetra
