module quter_defect
  use kinds,              only : dp
  use fc2_interpolate,    only : forceconst2_grid
  use ph_system,          only : ph_system_info
  use thutils,            only : v2index, index2v
  implicit none
  !
  integer, parameter :: nfar = 2
  !
  ! type forceconst2_sc
  !   INTEGER :: n_R = 0, i_0 = -1
  !   COMPLEX(dp), allocatable :: FC(:,:,:) ! (jn1,jn2,iR,jR)
  !   integer, allocatable  :: yR(:,:)    ! (3,iR,3,jR)
  !   real(dp), allocatable :: xR(:,:)    ! (3,iR,3,jR)
  !   integer               :: nq(3)          ! sc size
  !   logical :: centered = .false.           ! centered or not
  ! contains
  !   procedure :: allocate => allocate_fc2_sc
  ! end type
  !
contains
  !
  function sc2uc(S, S_sc, sc_grid, which)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    character(*) :: which
    !
    integer :: sc2uc(S_sc%nat)
    integer :: i, isc, iR, el
    real(dp) :: r_cryst(3), tau(3), tau_sc(3)
    sc2uc = -1
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        tau_sc = S_sc%tau(:,isc)
        call cryst_to_cart(1, tau_sc, S_sc%bg, -1)
        tau = S%tau(:,i)
        call cryst_to_cart(1, tau, S%bg, -1)
        r_cryst = tau_sc*sc_grid - tau
        if (NORM2(r_cryst - NINT(r_cryst))<1e-2) then
          iR = v2index(NINT(r_cryst), sc_grid)
          if (which == "R") then
            el = iR
          elseif (which == "nat") then
            el = i
          else
            CALL errore("sc2uc", "which is not R or nat", -1)
          endif
          if (sc2uc(el) > -1) CALL errore("sc2uc", "found same atom in same R", -1)
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("sc2uc", "R is out of bound", -ABS(iR))
          sc2uc(isc) = el
        endif
        !
      enddo
    enddo
    if (ANY(sc2uc == -1)) CALL errore("map_sc_atm", "some atoms are not mapped", -1)
  end function
  !
  ! SUBROUTINE allocate_fc2_sc(fc, n_R, nat)
  !   IMPLICIT NONE
  !   INTEGER,INTENT(in) :: n_R, nat
  !   CLASS(forceconst2_sc),INTENT(inout) :: fc
  !   CHARACTER(len=16),PARAMETER :: sub = "allocate_fc2_sc"
  !   !
  !   IF(allocated(fc%yR) .or. allocated(fc%xR) .or. allocated(fc%FC)) &
  !     CALL errore(sub, 'some element is already allocated', 1)
  !   !
  !   ALLOCATE(fc%yR(6,n_R))
  !   ALLOCATE(fc%xR(6,n_R))
  !   ALLOCATE(fc%FC(3*nat,3*nat,n_R))
  !   fc%n_R = n_R
  !   !
  ! END SUBROUTINE
  !
  subroutine center_sc(fc, grid, S, S_sc)
    class(forceconst2_grid), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    !
    ! type(forceconst2_grid) :: fsc
    integer :: nRbig, l1, l2, l3, big_grid(3), new_yR(6), small_Ri(3), small_Rj(3)
    real(dp) :: Rbig(3, (2*nfar+1)**3), wg_tot, big_yR(3)
    real(dp), allocatable :: new_fc(:,:,:)
    real(dp) :: dist(3)
    integer :: na1, na2, j1, j2, ixR, nR, jR_big
    integer, dimension(S_sc%nat) ::  map_R, map_nat
    integer :: R_list(PRODUCT(grid), PRODUCT(grid)*(2*nfar+1)**3)
    integer :: index_i, index_j, nxR
    integer, allocatable :: new_yR_list(:,:)
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: wg, rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,S_sc%at)
    !
    ! Construct a big enough lattice of Supercell **supercell** vectors
    fc%centered = .true.
    nR = PRODUCT(grid)
    big_grid = grid * (2*nfar+1)
    allocate(new_yR_list(6, nR**2*100))
    allocate(new_fc(S%nat3, S%nat3, nR**2*100))
    nRbig=0
    DO l1=-nfar, nfar
      DO l2=-nfar, nfar
        DO l3=-nfar, nfar
          nRbig=nRbig+1
          Rbig(:, nRbig) = S_sc%at(:,1)*l1 +S_sc%at(:,2)*l2 +S_sc%at(:,3)*l3
        END DO
      END DO
    END DO
    IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)
    !
    map_R = sc2uc(S, S_sc, grid, "R")
    map_nat = sc2uc(S, S_sc, grid, "nat")
    !
    nxR = 0

    R_list = -1
    do na1 = 1, S_sc%nat
      do na2 = 1, S_sc%nat
        small_Ri = index2v(map_R(na1), grid)
        ! call cryst_to_cart(1, small_Ri, S%at, 1)
        small_Rj = index2v(map_R(na2), grid)
        ! call cryst_to_cart(1, small_Rj, S%at, 1)
        wg_tot = 0._dp
        do jR_big = 1, nRbig
          dist = S_sc%tau(:,na1) - Rbig(:,jR_big) - S_sc%tau(:,na2)
          wg = wsweight(dist,rws,nrws)
          wg_tot = wg_tot + wg
          ! wg = wg/DFLOAT(nqt)
          if (nfar == 0) wg = 1._dp
          if (wg > 1e-6) then
            new_yR(1:3) = small_Ri
            big_yR = Rbig(:,jR_big)
            call cryst_to_cart(1, big_yR, S_sc%bg, -1)
            new_yR(4:6) = NINT(big_yR*grid + small_Rj)
            index_i = v2index(new_yR(1:3), grid)
            index_j = v2index(new_yR(4:6)+grid*nfar, big_grid)
            if (R_list(index_i, index_j) == -1) then
              nxR = nxR + 1
              ixR = nxR
              R_list(index_i, index_j) = ixR
              new_yR_list(:, ixR) = new_yR
            else
              ixR = R_list(index_i, index_j)
            endif
            if(norm2(REAL(new_yR, DP))<1e-6*S%alat) fc%i_0 = ixR
            do j1 = 1, 3
              do j2 = 1, 3
                new_fc(j1+3*(map_nat(na1)-1),j2+3*(map_nat(na2)-1),ixR) = &
                  fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, 1) * wg
                ! if(map_R(na1) == 1 .and. map_R(na2) == 3) print*, fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, 1) * wg
              enddo
            enddo
          endif
        enddo
        if (ABS(wg_tot -1)<1e-6) CALL errore("center_sc", "sum of weights is not 1", -1)
      enddo
    enddo
    !
    deallocate(fc%yR, fc%xR, fc%FC)
    ALLOCATE(fc%yR(6,nxR))
    ALLOCATE(fc%xR(6,nxR))
    ALLOCATE(fc%FC(S%nat3,S%nat3,nxR))
    fc%n_R = nxR
    fc%nq = grid
    fc%yR = new_yR_list(:,1:nxR)
    fc%xR = REAL(new_yR_list(:,1:nxR), DP)
    CALL cryst_to_cart(nxR, fc%xR(1:3,:), S%at, 1)
    CALL cryst_to_cart(nxR, fc%xR(4:6,:), S%at, 1)
    fc%FC = new_fc(:,:,1:nxR)
  end subroutine
  !
  function fc_gamma2RR(grid, S, S_sc, fc)
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    real(dp), intent(in) :: FC(S_sc%nat3,S_sc%nat3)
    !
    real(dp) :: fc_gamma2RR(S%nat3, S%nat3, product(grid), product(grid))
    integer :: sc_na1, sc_na2, j1, j2
    integer, dimension(S_sc%nat) :: map_R, map_nat
    !
    map_R = sc2uc(S, S_sc, grid, "R")
    map_nat = sc2uc(S, S_sc, grid, "nat")
    fc_gamma2RR = 0._dp
    do sc_na1 = 1, S_sc%nat
      do sc_na2 = 1, S_sc%nat
        do j1 = 1, 3
          do j2 = 1, 3
            fc_gamma2RR(j1+(map_nat(sc_na1)-1)*3, j2+(map_nat(sc_na2)-1)*3, map_R(sc_na1), map_R(sc_na2)) = &
              FC(j1+(sc_na1-1)*3, j2+(sc_na2-1)*3)
          enddo
        enddo
      enddo
    enddo
  end function
end module
