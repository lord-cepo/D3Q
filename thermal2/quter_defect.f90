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
  function map_sc_atm(S, S_sc, sc_grid)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    !
    integer :: map_sc_atm(S_sc%nat)
    integer :: i, isc, iR
    real(dp) :: r_cryst(3), tau(3), tau_sc(3)
    map_sc_atm = -1
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        tau_sc = S_sc%tau(:,isc)
        call cryst_to_cart(1, tau_sc, S_sc%bg, -1)
        tau = S%tau(:,i)
        call cryst_to_cart(1, tau, S%bg, -1)
        r_cryst = tau_sc*sc_grid - tau
        if (NORM2(r_cryst - NINT(r_cryst))<1e-2) then
          iR = v2index(NINT(r_cryst), sc_grid)
          if (map_sc_atm(i+(iR-1)*S%nat) > -1) CALL errore("map_sc_atm", "found same atom in same R", -1)
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("map_sc_atm", "R is out of bound", -ABS(iR))
          map_sc_atm(isc) = i+(iR-1)*S%nat
        endif
        !
      enddo
    enddo
    if (ANY(map_sc_atm == -1)) CALL errore("map_sc_atm", "some atoms are not mapped", -1)
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
    integer :: nRbig, l1, l2, l3
    real(dp) :: Rbig(3, (2*nfar+1)**3)
    real(dp) :: new_xR(6, PRODUCT(grid)*(2*nfar+1)**3)
    real(dp) :: new_fc(S%nat3, S%nat3, PRODUCT(grid)*(2*nfar+1)**3)
    real(dp) :: dist(3)
    integer :: na1, na2, j1, j2, iR, jR, small_Ri(3), small_Rj(3), ixR
    integer, dimension(S_sc%nat) :: map, map_R, map_nat
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: wg, rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    !
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,S_sc%at)
    !
    ! Construct a big enough lattice of Supercell **supercell** vectors
    fc%centered = .true.
    nRbig = (2*nfar+1)**3
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
    map = map_sc_atm(S, S_sc, grid)
    map_R = (map-1) / S%nat + 1
    map_nat = MOD(map-1, S%nat) + 1
    !
    ixR = 0
    do na1 = 1, S%nat
      do na2 = 1, S%nat
        do j1 = 1, 3
          do j2 = 1, 3
            do iR = 1, nRbig
              do jR = 1, nRbig
                small_Ri = index2v(map_R(na1), grid)
                small_Rj = index2v(map_R(na2), grid)
                dist = Rbig(:,iR) + S_sc%tau(:,na1) - Rbig(:,jR) - S_sc%tau(:,na2)
                wg = wsweight(dist,rws,nrws)
                ! wg = wg/DFLOAT(nqt)
                if (wg /= 0) then
                  ixR = ixR + 1
                  new_xR(1:3,ixR) = Rbig(:,iR) + small_Ri
                  new_xR(4:6,ixR) = Rbig(:,jR) + small_Rj
                  if(norm2(new_xR(:,ixR))<1e-6*S%alat) fc%i_0 = ixR
                  new_fc(map_nat(na1),map_nat(na2),ixR) = fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, 1) * wg
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !
    deallocate(fc%yR, fc%xR, fc%FC)
    ALLOCATE(fc%yR(6,ixR))
    ALLOCATE(fc%xR(6,ixR))
    ALLOCATE(fc%FC(S%nat3,S%nat3,ixR))
    fc%n_R = ixR
    fc%nq = grid
    fc%xR = new_xR(:,1:ixR)
    CALL cryst_to_cart(ixR, new_xR(1:3,1:ixR), S%bg, -1)
    CALL cryst_to_cart(ixR, new_xR(4:6,1:ixR), S%bg, -1)
    fc%yR = NINT(new_xR(:, 1:ixR))
    fc%FC = new_fc(:,:,1:ixR)
  end subroutine
end module
