module simtet
  use q_grids, only : q_grid
  use ph_system, only : ph_system_info
  use thtetra, only : equiv_grid
  use kinds, only : DP
  use mpi_thermal, only : mpi_bsum, num_procs, my_id
  !
  INTEGER :: ntetra = 0
  !! number of tetrahedra
  INTEGER :: nntetra
  !! k-points per tetrahedron used to compute weights.
  !! 4 for linear / 20 for optimized tetrahedron method
  INTEGER, ALLOCATABLE :: tetra(:,:)
  !! index of k-points in a given tetrahedron shape (nntetra,ntetra)
  INTEGER :: nqs
  !! number of q-points in IRREDUCIBLE BZ
  INTEGER :: nqtot
  !! number of q-points in the whole BZ
  INTEGER :: nbnd
  !! number of bands
  INTEGER, ALLOCATABLE :: itetra(:,:,:)
  !! order index of vertices of each tetrahedron (4,nbnd,ntetra)
  complex(DP), ALLOCATABLE :: ek_sort(:,:,:)
  !! sorted energies for each tetrahedron (4,nbnd,ntetra)
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
  real(dp) :: TETRA_MULTIPLIER
  !! brings ek to unity to avoid float overflow
  integer :: NUM_PROCS_TETRA, MY_ID_TETRA
  !
  PUBLIC
contains
  SUBROUTINE tetra_init_sym_nosym(grid, S, ek, ek_tetra, mpi)
    !-----------------------------------------------------------------------------
    !! This rouotine sets the corners and additional points for each tetrahedron.
    !
    use ph_system, only : ph_system_info
    use q_grids,   only : q_grid, setup_grid
    IMPLICIT NONE
    !
    type(ph_system_info), intent(in) :: S
    !! usual system info
    complex(DP), INTENT(IN) :: ek(:,:)
    !! energy in the form ek(ibnd, iq)
    type(q_grid), intent(in) :: grid
    !! accepts symmetrized grids
    complex(dp), intent(out) :: ek_tetra(4,size(ek, 1),product(grid%n)*6)
    !! size (4,nbnd,ntetra)
    logical, intent(in), optional :: mpi

    REAL(DP), PARAMETER :: eps = 1e-5_dp
    !
    INTEGER :: i1, i2, i3, itet, itettot, ii, ik,  rest, &
      ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3), ibnd
    ! integer :: tetra_ik(nq(1) * nq(2) * nq(3))
    !
    REAL(DP) :: l(4), bvec2(3,3), bvec3(3,4) !xkg(3, product(nq))
    external :: SIM0ONEI
    !
    NUM_PROCS_TETRA = num_procs
    MY_ID_TETRA = my_id
    if(present(mpi)) then
      if(.not.mpi) then
        NUM_PROCS_TETRA = 1
        MY_ID_TETRA = 0
      endif
    endif
    !
    nbnd = SIZE(ek,1)
    nqs =  size(ek,2)
    !
    ntetra  = 6*product(grid%n)
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
    !
    call equiv_grid(grid, S, equiv)

    ! tetra_ik = 0
    DO itettot = 1+MY_ID_TETRA, ntetra, NUM_PROCS_TETRA
      itet = mod(itettot,6) + 1
      rest = itettot / 6
      i3 = mod(rest,grid%n(3)) + 1
      rest = rest / grid%n(3)
      i2 = mod(rest,grid%n(2)) + 1
      rest = rest / grid%n(2)
      i1 = mod(rest,grid%n(1)) + 1

      do ibnd = 1, nbnd
        DO ii = 1, 4
          !
          ikv(1:3) = (/i1, i2, i3/) - 1
          ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
          ikv(1:3) = MODULO(ikv(1:3), (/grid%n(1), grid%n(2), grid%n(3)/))
          !
          ik = ikv(3) + grid%n(3) * (ikv(2) + grid%n(2) * ikv(1)) + 1
          !
          ek_tetra(ii,ibnd,itettot) = ek(ibnd,equiv(ik))
        END DO ! ii
        !
      enddo
    ENDDO ! itettot
    !
    ! TETRA_MULTIPLIER = real(size(ek), dp) / sum(ek)
    ! ek_sort = ek_sort * TETRA_MULTIPLIER
  END SUBROUTINE
  !
  function bz_integral(e, S, grid, ek)
    real(dp), intent(in) :: e
    type(ph_system_info), intent(in) :: S
    type(q_grid) :: grid
    complex(dp), intent(in) :: ek(:,:)
    complex(dp) :: bz_integral
    !
    complex(dp) :: ek_tetra(4,size(ek, 1),product(grid%n)*6)
    integer :: it, ibnd
    complex(dp) :: dummy, tetra_integral
    !
    call tetra_init_sym_nosym(grid, S, ek, ek_tetra)
    !
    bz_integral = 0.0_dp
    do it = 1, 6*product(grid%n)
      do ibnd = 1, size(ek, 1)
        call SIM0ONEI(tetra_integral, dummy, e - ek_tetra(:, ibnd, it))
        bz_integral = bz_integral + tetra_integral
      enddo
    enddo
  end function
  !
  SUBROUTINE tetra_init_sym_cmplx(grid, S, ek, mpi)
    !-----------------------------------------------------------------------------
    !! This rouotine sets the corners and additional points for each tetrahedron.
    !
    use ph_system, only : ph_system_info
    use q_grids,   only : q_grid, setup_grid
    IMPLICIT NONE
    !
    type(ph_system_info), intent(in) :: S
    !! usual system info
    complex(DP), INTENT(IN) :: ek(:,:)
    !! energy in the form ek(ibnd, iq)
    type(q_grid), intent(in) :: grid
    !! accepts symmetrized grids
    logical, intent(in), optional :: mpi

    ! LOGICAL, INTENT(IN) :: is_mpi
    !! if .true., the grid is scattered
    ! LOGICAL, INTENT(IN), OPTIONAL :: opt
    ! !! if .true., uses opt_tetra methods

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
    NUM_PROCS_TETRA = num_procs
    MY_ID_TETRA = my_id
    if(present(mpi)) then
      if(.not.mpi) then
        NUM_PROCS_TETRA = 1
        MY_ID_TETRA = 0
      endif
    endif
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
    ALLOCATE(itetra (4,nbnd,ntetra))
    allocate(iisize_tetra(ntetra))
    allocate(nt_tetra(ntetra))
    ! ALLOCATE(which_tetra (2,nqs,24))
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
    !
    !
    nntetra = 4
    allocate(ii_tetra(nntetra, ntetra))
    IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra(nntetra,ntetra) )
    !  locate k-points of the uniform grid in the list of irreducible k-points
    !  that was previously calculated
    !
    !  bring irreducible k-points to crystal axis
    !
    ! nqtot = nqs
    call equiv_grid(grid, S, equiv, first_point)
    itvalid = 0
    itetra = 0
    tetra = 0
    ! tetra_ik = 0
    DO itettot = 1+MY_ID_TETRA, ntetra, NUM_PROCS_TETRA
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
            ek_sort(ii,ibnd,itvalid) = ek(ibnd,equiv(ik))
          enddo
        enddo
      endif
    ENDDO ! itettot
    !
    TETRA_MULTIPLIER = real(size(ek), dp) / ABS(sum(ek))
    ek_sort = ek_sort * TETRA_MULTIPLIER
    nvalid = itvalid
  END SUBROUTINE
  !
  FUNCTION tetra_weights_green_cmplx(ef) RESULT(wg_sym)
    USE constants, ONLY : pi
    use thtetra, only : rm_degen_vertices
    !-----------------------------------------------------------------------------------
    !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
    !! The resulting wg can be used as sum(Ak * wk)
    !-----------------------------------------------------------------------------------
    real(DP), INTENT(IN) :: ef
    !! The Fermi energy
    !
    ! ... local variables
    !
    complex(dp) :: wg_sym(nbnd, nqs)
    real(dp) :: r(4), i(4)
    INTEGER :: ik, nt, ibnd, ii,  ii_
    complex(DP) :: e(4)
    complex(dp) :: w(4), dummy(4)
    real(dp) :: ef_

    ! for real part calc
    ! REAL(DP) :: wR0(4), ef_e(4), log_ef_e(4), prod_a(4), sum_a(4), second_term(4)
    ! INTEGER :: i3, j3
    !
    ef_ = ef * TETRA_MULTIPLIER
    wg_sym = 0._dp
    !
    DO nt = 1, nvalid
      !
      DO ibnd = 1, nbnd
        !
        e = ek_sort(:,ibnd,nt)
        ! r = real(e, dp)
        ! i = -aimag(e)
        ! call rm_degen_vertices(ef_, r)
        ! call rm_degen_vertices(ef_, i)
        ! e = cmplx(r, -i, dp)
        call SIM0TWOI(w, dummy, ef_ - e)
        !
        DO ii_ = 1, iisize_tetra(nt)
          !
          ii = ii_tetra(ii_, nt)
          ik = tetra(ii, nt)
          wg_sym(ibnd,ik) = wg_sym(ibnd,ik) + w(ii)
        ENDDO
        !
      ENDDO ! ibnd
      !
    ENDDO ! nt
    ! wg = wg / REAL(ntetra, dp)
    wg_sym = wg_sym / ntetra * TETRA_MULTIPLIER
    !
    ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
    if (NUM_PROCS_TETRA > 1) &
      CALL mpi_bsum(nbnd, nqs, wg_sym)
    !
  END FUNCTION
  !
  SUBROUTINE deallocate_tetra( )
    !--------------------------------------------------------
    !! Deallocate tetra and wlsm
    !
    ntetra = 0
    nntetra = 0
    nqs = 0
    nbnd = 0
    nvalid = 0
    IF (ALLOCATED(tetra  ))      DEALLOCATE (tetra  )
    IF (ALLOCATED(ek_sort))      DEALLOCATE (ek_sort)
    IF (ALLOCATED(itetra ))      DEALLOCATE (itetra )
    if (allocated(ii_tetra))     deallocate (ii_tetra)
    if (allocated(nt_tetra))     deallocate (nt_tetra)
    if (allocated(iisize_tetra)) deallocate (iisize_tetra)
    if (allocated(equiv))        deallocate (equiv)
    !
  END SUBROUTINE deallocate_tetra
end module
