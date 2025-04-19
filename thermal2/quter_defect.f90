module quter_defect
  use kinds,              only : dp
  use fc2_interpolate,    only : forceconst2_grid
  use ph_system,          only : ph_system_info
  use thutils,            only : v2index, index2v
  implicit none
  !
  integer, parameter :: nfar = 2
  !
  type forceconst2_sc
    INTEGER :: n_R1 = 0
    INTEGER, allocatable :: n_R2(:)
    ! INTEGER :: i_0(2) = 0
    real(dp), allocatable :: FC(:,:,:,:) ! (jn1,jn2,iR,jR)
    complex(dp), allocatable :: mix(:,:,:) ! (jn1,jn2,iR)
    integer,  allocatable :: yR2(:,:,:)    ! (3,iR,which)
    real(dp), allocatable :: xR2(:,:,:)    ! (3,iR,which)
    integer,  allocatable :: yR1(:,:)    ! (3,iR)
    real(dp), allocatable :: xR1(:,:)    ! (3,iR)
    integer               :: nq(3)       ! sc size
    integer :: stage = -1                ! centered or not

  contains
    procedure :: allocate => allocate_fc2_sc
    procedure :: center => center_sc
    procedure :: cryst => construct_cryst
    procedure :: cart => construct_cart
    generic   :: interpolate => interp_1st_step, interp_2nd_step
    PROCEDURE, PASS :: interp_1st_step
    PROCEDURE, PASS :: interp_2nd_step
  end type
  !
contains
  !
  !> map_atm_sc(inat, iR) = inat_sc
  function map_uc2sc(S, S_sc, sc_grid, idef)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    integer, intent(out), optional :: idef
    !
    integer :: map_uc2sc(S%nat, PRODUCT(sc_grid))
    integer :: i, isc, iR
    real(dp) :: r_cryst(3)
    map_uc2sc = -1
    if (present(idef)) idef = -1
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        r_cryst = S_sc%tau(:,isc)*sc_grid-S%tau(:,i)
        call cryst_to_cart(1, r_cryst, S_sc%bg, -1)
        if (NORM2(r_cryst - NINT(r_cryst))<1e-2) then
          iR = v2index(NINT(r_cryst), sc_grid)
          if (map_uc2sc(i, iR) > -1) CALL errore("map_uc2sc", "found same atom in same R", 1)
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("map_uc2sc", "R is out of bound", ABS(iR))
          map_uc2sc(i, iR) = isc
          if(present(idef) .and. S%ityp(i) /= S_sc%ityp(isc)) then
            if(idef == -1) then
              idef = isc
            else
              call errore("map_uc2sc", "Found two defects", 1)
            endif
          endif
          !
        endif
      enddo
    enddo
    if(present(idef)) then
      if(idef == -1) call errore("map_uc2sc", "defect not found", 1)
    endif
    if (ANY(map_uc2sc == -1)) CALL errore("map_uc2sc", "some atoms are not mapped", 1)
  end function
!
  function map_sc2uc(S, S_sc, sc_grid, which)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    character(*) :: which
    !
    integer :: map_sc2uc(S_sc%nat)
    integer :: i, isc, iR, el
    real(dp) :: r_cryst(3), tau(3), tau_sc(3)
    map_sc2uc = -1
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
            CALL errore("map_sc2uc", "which is not R or nat", 1)
          endif
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("map_sc2uc", "R is out of bound", ABS(iR))
          map_sc2uc(isc) = el
        endif
        !
      enddo
    enddo
    if (ANY(map_sc2uc == -1)) CALL errore("map_sc_atm", "some atoms are not mapped", 1)
  end function
  !
  SUBROUTINE allocate_fc2_sc(fc, S, grid)
    use thutils, only : grid_vec
    IMPLICIT NONE
    INTEGER,INTENT(in) :: grid(3)
    type(ph_system_info),INTENT(in) :: S
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    CHARACTER(len=16),PARAMETER :: sub = "allocate_fc2_sc"
    integer :: n_R, iR
    !
    IF(allocated(fc%yR2) .or. allocated(fc%FC) &
      .or. allocated(fc%xR2) .or. allocated(fc%mix)) &
      CALL errore(sub, 'some element is already allocated', 1)
    !
    n_R = PRODUCT(grid)
    fc%n_R1 = n_R
    allocate(fc%n_R2(n_R))
    fc%n_R2 = n_R
    !
    ALLOCATE(fc%yR2(3,n_R,n_R), fc%yR1(3,n_R))
    fc%yR1 = grid_vec(grid)
    do iR = 1, n_R
      fc%yR2(:,:,iR) = fc%yR1
    enddo
    !
    ALLOCATE(fc%xR2(3,n_R,n_R), fc%xR1(3,n_R))
    call fc%cart(S, 1)
    call fc%cart(S, 2)
    ALLOCATE(fc%FC(S%nat3,S%nat3,n_R,n_R))
    fc%nq = grid
    fc%stage = -1
    !
  END SUBROUTINE
  !
  subroutine construct_cart(fc, S, which)
    class(forceconst2_sc), intent(inout) :: fc
    type(ph_system_info), intent(in) :: S
    integer, intent(in) :: which
    !
    real(dp), allocatable :: R(:,:)
    integer :: i
    !
    if(which == 2) then
      do i = 1, fc%n_R1
        allocate(R(3,fc%n_R2(i)))
        R = REAL(fc%yR2(:,:,i), DP)
        call cryst_to_cart(fc%n_R2(i), R, S%at, 1)
        fc%xR2(:,:,i) = R
        deallocate(R)
      enddo
      !
    elseif(which == 1) then
      allocate(R(3,fc%n_R1))
      R = REAL(fc%yR1, DP)
      call cryst_to_cart(fc%n_R1, R, S%at, 1)
      fc%xR1 = R
      deallocate(R)
    else
      call errore("construct_cart", "which is not 1 or 2", 1)
      !
    endif
  end subroutine
  !
  subroutine construct_cryst(fc, S, which)
    class(forceconst2_sc), intent(inout) :: fc
    type(ph_system_info), intent(in) :: S
    integer, intent(in) :: which
    !
    real(dp), allocatable :: R(:,:)
    integer :: i
    !
    if (which == 2) then
      do i = 1, fc%n_R1
        allocate(R(3,fc%n_R2(i)))
        R = fc%xR2(:,:,i)
        call cryst_to_cart(fc%n_R2(i), R, S%bg, -1)
        if (ALL(ABS(R - NINT(R)) < 1e-6)) then
          fc%yR2(:,:,i) = NINT(R)
        else
          call errore("construct_cryst", "R is not integer", 1)
        endif
        deallocate(R)
      enddo
      !
    elseif(which == 1) then
      allocate(R(3,fc%n_R1))
      R = fc%xR1
      call cryst_to_cart(fc%n_R1, R, S%bg, -1)
      if (ALL(ABS(R - NINT(R)) < 1e-6)) then
        fc%yR1 = NINT(R)
      else
        call errore("construct_cryst", "R is not integer", 1)
      endif
      deallocate(R)
    else
      call errore("construct_cryst", "which is not 1 or 2", 1)
    endif
  end subroutine
  !
  subroutine S_uc2sc(S, grid, S_sc)
    use thutils, only : grid_vec
    use ph_system, only : aux_system
    type(ph_system_info), intent(in) :: S
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(out) :: S_sc
    !
    integer :: nR, iR, i, t
    integer, allocatable :: sc_grid(:,:)
    real(dp) :: tau(3,S%nat)
    !
    nR = PRODUCT(grid)
    S_sc%ntyp = S%ntyp
    S_sc%amass = S%amass
    S_sc%amass_variance = S%amass_variance
    S_sc%atm = S%atm
    S_sc%nat = S%nat * nR
    allocate(S_sc%ityp(S_sc%ntyp))
    S_sc%ityp = reshape(spread(S%ityp, dim=1, ncopies=nR), [S_sc%nat])
    S_sc%ibrav = S%ibrav
    S_sc%symm_type = S%symm_type
    S_sc%celldm = S%celldm
    S_sc%celldm(1) = S%celldm(1) * grid(1) ! bad coding here
    S_sc%at = S%at
    S_sc%bg = S%bg
    S_sc%omega = S%omega * nR
    S_sc%alat = S%alat * grid(1) ! also here
    S_sc%tpiba = S%tpiba / grid(1)
    S_sc%epsil = S%epsil
    S_sc%lrigid = S%lrigid
    allocate(S_sc%tau(3, S_sc%nat))
    tau = S%tau
    call cryst_to_cart(S%nat, tau, S%bg, -1)
    sc_grid = grid_vec(grid)
    do ir = 1, nR
      do t = 1, S%nat
        do i = 1, 3
          S_sc%tau(i,t + (ir-1)*S%nat) = (tau(i,t) + sc_grid(i,ir) ) / grid(i)
        enddo
      enddo
    enddo
    call cryst_to_cart(S_sc%nat, S_sc%tau, S_sc%at, 1)
    ! no zeu and dzeu !
    call aux_system(S_sc)
  end subroutine
  !
  subroutine center_sc(fc, grid, S, S_sc, distance)
    use functions, only : refold_bz
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    real(dp), allocatable, optional, intent(out) :: distance(:,:,:,:)
    !
    ! type(forceconst2_grid) :: fsc
    real(dp) :: dist(3), Rbig(3), wg_tot
    !
    integer :: na1, na2, j1, j2, nR, jR_big, R1, R2
    integer :: R_list(PRODUCT(grid)*(2*nfar+1)**3)
    integer :: map_sc(S%nat, PRODUCT(grid))
    integer :: nRbig, R_vec(3)
    integer :: ind, ixR
    integer :: nxR(PRODUCT(grid))
    integer, dimension(3) :: far_grid, Rbig_from_0, Rbig_shift
    !
    integer, allocatable :: new_yR_list(:,:,:)
    real(dp), allocatable :: new_fc(:,:,:,:)
    integer :: counter
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: wg, rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,S_sc%at)
    !
    fc%stage = 0
    nR = PRODUCT(grid)
    if (nfar == 0) return
    far_grid = 2*nfar+1
    nRbig = PRODUCT(far_grid)
    allocate(new_yR_list(3, nR*nRbig,nR))
    allocate(new_fc(S%nat3, S%nat3, nR*nRbig, nR))
    !
    if(present(distance)) allocate(distance, source=new_fc)
    map_sc = map_uc2sc(S, S_sc, grid)
    !
    nxR = 0
    new_fc = 0._dp
    !
    counter = 0
    do R1 = 1, nR
      R_list = -1
      do R2 = 1, nR
        do na1 = 1, S%nat
          do na2 = 1, S%nat
            wg_tot = 0._dp
            do jR_big = 1, nRbig
              !> I create a normal [0:N-1]^3 grid
              Rbig_from_0 = index2v(jR_big, far_grid)
              !> Then I shift it so that the center in [-N/2:N/2]
              Rbig_shift = Rbig_from_0 - nfar
              !> I work in S_sc%at alat units
              Rbig = REAL(Rbig_shift, DP)
              call cryst_to_cart(1, Rbig, S_sc%at, 1)
              !> Rbig is the R vector in the supercell, so we don't need to multiply by grid,
              !> cause tau is in [0,1] in crystal units
              dist = S_sc%tau(:,map_sc(na1,R1)) - Rbig - S_sc%tau(:,map_sc(na2,R2))
              !> Compute the Wigner-Seitz weight
              wg = wsweight(dist,rws,nrws)
              wg_tot = wg_tot + wg
              ! if (nfar == 0) wg = 1._dp
              if (wg /= 0) then
                !> R2 is in the unit cell, so we multiply Rbig by grid to transform it in a
                !> supercell vector of the unit cell
                R_vec = - index2v(R2, grid) - grid * Rbig_shift
                !> I use the 0-indexing to populate R_list at the correct non-negative integer
                ind = v2index(index2v(R2, grid) + grid * Rbig_from_0, grid * far_grid)
                !> R_list contains the ixR or -1. The nxR is refreshed at each step
                if (R_list(ind) == -1) then
                  nxR(R1) = nxR(R1) + 1
                  ixR = nxR(R1)
                  R_list(ind) = ixR
                  new_yR_list(:,ixR,R1) = R_vec
                else
                  ixR = R_list(ind)
                endif
                !> I find Gamma Gamma
                ! if(ALL(R_vec == 0)) fc%i_0 = ixR
                !> and populate force constants with usual index wrapping
                do j1 = 1, 3
                  do j2 = 1, 3
                    new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR, R1) = &
                      fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, R2, R1) * wg
                    if (present(distance)) &
                      distance(j1+3*(na1-1), j2+3*(na2-1), ixR, R1) = norm2(dist)
                  enddo
                enddo
              endif
            enddo
            if (ABS(wg_tot -1) > 1e-6) then
              print"(A,F14.6)", "wg_tot is", wg_tot
              CALL errore("center_sc", "sum of weights is not 1", 1)
            endif
          enddo
        enddo
      enddo
    enddo
    !
    deallocate(fc%yR2, fc%xR2, fc%FC)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR2(3,maxval(nxR),fc%n_R1))
    ALLOCATE(fc%xR2(3,maxval(nxR),fc%n_R1))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR),nR))
    fc%n_R2 = nxR
    fc%nq = grid
    do R1 = 1, nR
      ! fc%xr1(:,R1) = refold_bz(fc%xR1(:,R1), S%at)
      fc%yR2(:,:nxR(R1),R1) = new_yR_list(:,:nxR(R1),R1)
    enddo
    call fc%cart(S_sc, 2)
    ! call fc%cryst(S_sc, 1)
    fc%FC = new_fc(:,:,:maxval(nxR),:)
  end subroutine
  !
  subroutine center2(fc, grid, S, S_sc)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    !
    ! type(forceconst2_grid) :: fsc
    real(dp), dimension(3) :: taudef, d1, d2
    real(dp) :: perix, peri_min
    real(dp), parameter :: eps_peri = 1e-3
    integer, parameter :: nperix = (2*nfar+1)**3
    integer :: nperi
    integer :: SAFE_ALLOCATION
    !
    integer :: na1, na2, j1, j2, nR, R1, R2, R1_big, R2_big
    integer :: map_sc(S%nat, PRODUCT(grid))
    integer :: far_grid_cryst(3,(2*nfar+1)**3)
    real(dp) :: far_grid_cart(3,(2*nfar+1)**3)
    integer :: nRbig, counter
    integer :: idef, iperi, nxR1, ixR1, ixR2
    integer, allocatable :: R1_list(:), R2_list(:,:), yR1_list(:,:), yR2_list(:,:,:), nxR2(:)
    integer, dimension(3) :: far_mesh
    integer :: farx_list(3,2,nperix), ind(2,nperix)
    !
    real(dp), allocatable :: new_fc(:,:,:,:)
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,S_sc%at)
    !
    fc%stage = 0
    nR = PRODUCT(grid)
    if (nfar == 0) return
    far_mesh = 2*nfar+1
    nRbig = PRODUCT(far_mesh)
    !
    far_grid_cryst = grid_vec_cryst(far_mesh, -nfar)
    far_grid_cart = grid_vec_cart(far_mesh, S_sc%at, -nfar)
    !
    SAFE_ALLOCATION = 10 * nR
    allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(R1_list(nR*nRbig))
    allocate(R2_list(nR*nRbig,nR*nRbig))
    allocate(yR1_list(3, nR*nRbig))
    allocate(yR2_list(3, nR*nRbig, nR*nRbig))
    allocate(nxR2(nR*nRbig))
    !
    map_sc = map_uc2sc(S, S_sc, grid, idef)
    taudef = S_sc%tau(:,idef)
    !
    new_fc = 0._dp
    R1_list = -1
    R2_list = -1
    nxR2 = 0
    nxR1 = 0
    counter = 0
    !
    do R1 = 1, nR
      do na1 = 1, S%nat
        do R2 = 1, nR
          do na2 = 1, S%nat
            nperi = 0
            peri_min = 0._dp
            do R1_big = 1, nRbig
              d1 = S_sc%tau(:,map_sc(na1,R1)) + far_grid_cart(:,R1_big)
              do R2_big = 1, nRbig
                d2 = S_sc%tau(:,map_sc(na2,R2)) + far_grid_cart(:,R2_big)
                perix = (norm2(d1 - taudef) + norm2(d2 - taudef)) + norm2(d1 - d2)
                IF (perix < peri_min-eps_peri .or. nperi==0 ) THEN
                  nperi = 1
                  farx_list = 0
                  ind = 0
                  peri_min = perix
                  farx_list(:,1,nperi) = index2v(R1, grid) + grid * far_grid_cryst(:,R1_big)
                  farx_list(:,2,nperi) = - index2v(R2, grid) - grid * far_grid_cryst(:,R2_big)
                  ind(1,nperi) = R1 + nR*(R1_big-1)
                  ind(2,nperi) = R2 + nR*(R2_big-1)
                ELSE IF ( ABS(perix-peri_min) <= eps_peri ) THEN
                  nperi = nperi + 1
                  IF(nperi > nperix) CALL errore("center2", "nperix is too small", nperi)
                  peri_min = (peri_min*(nperi-1)+perix)/DBLE(nperi)
                  farx_list(:,1,nperi) = index2v(R1, grid) + grid * far_grid_cryst(:,R1_big)
                  farx_list(:,2,nperi) = - index2v(R2, grid) - grid * far_grid_cryst(:,R2_big)
                  ind(1,nperi) = R1 + nR*(R1_big-1)
                  ind(2,nperi) = R2 + nR*(R2_big-1)
                END IF
              enddo
            enddo
            if (nperi > 1) counter = counter + 1
            !
            do iperi = 1, nperi
              call add_ind(R1_list, ind(1,iperi), nxR1, ixR1)
              yR1_list(:,ixR1) = farx_list(:,1,iperi)
              !
              call add_ind(R2_list(:,ixR1), ind(2,iperi), nxR2(ixR1), ixR2)
              yR2_list(:,ixR2,ixR1) = farx_list(:,2,iperi)
              !> and populate force constants with usual index wrapping
              do j1 = 1, 3
                do j2 = 1, 3
                  new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR2, ixR1) = &
                    fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, R2, R1) / nperi
                  ! if (present(distance)) &
                  !   distance(j1+3*(na1-1), j2+3*(na2-1), ixR2, ixR1) = norm2(dist)
                enddo
              enddo
            enddo
            !
          enddo
        enddo
      enddo
    enddo
    !
    deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R2)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR1(3,nxR1))
    ALLOCATE(fc%xR1(3,nxR1))
    ALLOCATE(fc%yR2(3,maxval(nxR2),nxR1))
    ALLOCATE(fc%xR2(3,maxval(nxR2),nxR1))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR2),nxR1))
    ALLOCATE(fc%n_R2(nxR1))
    fc%n_R2 = nxR2(:nxR1)
    fc%n_R1 = nxR1
    fc%nq = grid
    do ixR1 = 1, nxR1
      fc%yR1(:,ixR1) = yR1_list(:,ixR1)
      fc%yR2(:,:nxR2(ixR1),ixR1) = yR2_list(:,:nxR2(ixR1),ixR1)
      fc%FC(:,:,:nxR2(ixR1),ixR1) = new_fc(:,:,:nxR2(ixR1),ixR1)
    enddo
    call fc%cart(S_sc, 1)
    call fc%cart(S_sc, 2)

    call print_message("fine del centering")

    print*, "------------------------"
    print*, "number of R1", nxR1
    print*, "number of R1,R2", sum(nxR2)
    ! print*, "number of R2", nxR2(:nxR1)
    ! do R1 = 1, nxR1
    !   print"(3I4,3F10.3)", fc%yR1(:,R1), fc%xR1(:,R1)
    ! enddo
    print*, "counter is  ", counter, "over", nR**2*S%nat**2
    print*, "------------------------"

    deallocate(new_fc, R1_list, R2_list, yR1_list, yR2_list, nxR2)
  end subroutine
  !
  subroutine center3(fc, grid, S, S_sc)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    !
    ! type(forceconst2_grid) :: fsc
    real(dp), dimension(3) :: taudef, d1, d2
    real(dp), parameter :: eps_peri = 1e-3
    integer, parameter :: nperix = (2*nfar+1)**3
    integer :: SAFE_ALLOCATION
    !
    integer :: na1, na2, j1, j2, nR, R1, R2, R1_big, R2_big
    integer :: map_sc(S%nat, PRODUCT(grid))
    integer :: far_grid_cryst(3,(2*nfar+1)**3)
    real(dp) :: far_grid_cart(3,(2*nfar+1)**3)
    integer :: nRbig
    integer :: idef, nxR1, ixR1, ixR2
    integer, allocatable :: R1_list(:), R2_list(:,:), yR1_list(:,:), yR2_list(:,:,:), nxR2(:)
    real(dp), allocatable :: weights1(:), weights2(:)
    integer, allocatable :: inds1(:), inds2(:)
    integer, dimension(3) :: far_mesh
    !
    real(dp), allocatable :: new_fc(:,:,:,:)
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,S_sc%at)
    !
    fc%stage = 0
    nR = PRODUCT(grid)
    if (nfar == 0) return
    far_mesh = 2*nfar+1
    nRbig = PRODUCT(far_mesh)
    !
    far_grid_cryst = grid_vec_cryst(far_mesh, -nfar)
    far_grid_cart = grid_vec_cart(far_mesh, S_sc%at, -nfar)
    !
    SAFE_ALLOCATION = 10 * nR
    allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(R1_list(nR*nRbig))
    allocate(R2_list(nR*nRbig, nR*nRbig))
    allocate(yR1_list(3, SAFE_ALLOCATION))
    allocate(yR2_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(nxR2(SAFE_ALLOCATION))
    !
    map_sc = map_uc2sc(S, S_sc, grid, idef)
    taudef = S_sc%tau(:,idef)
    !
    new_fc = 0._dp
    R1_list = -1
    R2_list = -1
    nxR2 = 0
    nxR1 = 0
    !
    do R1 = 1, nR
      do na1 = 1, S%nat
        d1 = S_sc%tau(:,map_sc(na1,R1))-taudef
        call inside_ws(far_grid_cart, d1, nrws, rws, weights1, inds1)
        do R1_big = 1, size(weights1)
          call add_ind(R1_list, R1 + nR*(inds1(R1_big)-1), nxR1, ixR1)
          yR1_list(:,ixR1) = index2v(R1, grid) + grid * far_grid_cryst(:,inds1(R1_big))
          do R2 = 1, nR
            do na2 = 1, S%nat
              d2 = S_sc%tau(:,map_sc(na2,R2)) - S_sc%tau(:,map_sc(na1,R1)) - far_grid_cart(:,inds1(R1_big))
              call inside_ws(far_grid_cart, d2, nrws, rws, weights2, inds2)
              do R2_big = 1, size(weights2)
                call add_ind(R2_list(:,ixR1), R2 + nR*(inds2(R2_big)-1), nxR2(ixR1), ixR2)
                yR2_list(:,ixR2,ixR1) = - index2v(R2, grid) - grid * far_grid_cryst(:,inds2(R2_big))
                do j1 = 1, 3
                  do j2 = 1, 3
                    new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR2, ixR1) = &
                      fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, R2, R1) * weights1(R1_big) * weights2(R2_big)
                  enddo
                enddo
              enddo
            enddo
            !
          enddo
        enddo
      enddo
    enddo
    !
    deallocate(weights1, weights2, inds1, inds2)
    deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R2)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR1(3,nxR1))
    ALLOCATE(fc%xR1(3,nxR1))
    ALLOCATE(fc%yR2(3,maxval(nxR2),nxR1))
    ALLOCATE(fc%xR2(3,maxval(nxR2),nxR1))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR2),nxR1))
    ALLOCATE(fc%n_R2(nxR1))
    fc%n_R2 = nxR2(:nxR1)
    fc%n_R1 = nxR1
    fc%nq = grid
    fc%yR2 = -1000
    do ixR1 = 1, nxR1
      fc%yR1(:,ixR1) = yR1_list(:,ixR1)
      fc%yR2(:,:nxR2(ixR1),ixR1) = yR2_list(:,:nxR2(ixR1),ixR1)
      fc%FC(:,:,:nxR2(ixR1),ixR1) = new_fc(:,:,:nxR2(ixR1),ixR1)
    enddo
    call fc%cart(S_sc, 1)
    call fc%cart(S_sc, 2)

    call print_message("fine del centering")

    deallocate(new_fc, R1_list, R2_list, yR1_list, yR2_list, nxR2)
  end subroutine
  !
  subroutine inside_ws(R_list, r, nrws, rws, weights, ind_out)
    real(dp), intent(in) :: R_list(:,:)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: rws(:,:)
    !
    real(dp), allocatable, intent(out) :: weights(:)
    integer, allocatable, intent(out) :: ind_out(:)
    !
    real(dp), allocatable :: weights_(:)
    integer, allocatable :: ind_out_(:)
    integer :: i, nweights, nrws
    real(dp) :: wg, wg_tot
    real(dp), external :: wsweight
    !
    allocate(weights_(size(R_list, 2)))
    allocate(ind_out_(size(R_list, 2)))
    wg_tot = 0._dp
    nweights = 0
    do i = 1, size(R_list, 2)
      wg = wsweight(r+R_list(:,i),rws,nrws)
      if (wg > 1e-8) then
        wg_tot = wg_tot + wg
        nweights = nweights + 1
        weights_(nweights) = wg
        ind_out_(nweights) = i
      endif
    enddo
    if (ABS(wg_tot-1)>1e-8) then
      print*, "wg_tot is", wg_tot
      call errore("inside_ws", "sum of weights is not 1", 1)
    endif
    allocate(weights(nweights), ind_out(nweights))
    weights = weights_(:nweights)
    ind_out = ind_out_(:nweights)
    deallocate(weights_, ind_out_)
  end subroutine
  !
  subroutine add_ind(ind_list, ind, nR, iR)
    integer, intent(inout) :: ind_list(:)
    integer, intent(in) :: ind
    integer, intent(inout) :: nR
    integer, intent(out) :: iR
    !
    if(ind_list(ind) == -1) then
      nR = nR + 1
      iR = nR
      ind_list(ind) = iR
    else
      iR = ind_list(ind)
    endif
  end subroutine
  !
  subroutine center_grid_sc(fc, grid, S, S_sc)
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
    map_R = map_sc2uc(S, S_sc, grid, "R")
    map_nat = map_sc2uc(S, S_sc, grid, "nat")
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
  function fc_sc2RR(grid, S, S_sc, fc)
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    real(dp), intent(in) :: FC(S_sc%nat3,S_sc%nat3)
    !
    real(dp) :: fc_sc2RR(S%nat3, S%nat3, product(grid), product(grid))
    integer :: sc_na1, sc_na2, j1, j2
    integer, dimension(S_sc%nat) :: map_R, map_nat
    !
    map_R = map_sc2uc(S, S_sc, grid, "R")
    map_nat = map_sc2uc(S, S_sc, grid, "nat")
    fc_sc2RR = 0._dp
    do sc_na1 = 1, S_sc%nat
      do sc_na2 = 1, S_sc%nat
        do j1 = 1, 3
          do j2 = 1, 3
            fc_sc2RR(j1+(map_nat(sc_na1)-1)*3, j2+(map_nat(sc_na2)-1)*3, map_R(sc_na2), map_R(sc_na1)) = &
              FC(j1+(sc_na1-1)*3, j2+(sc_na2-1)*3)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function fc_uc2sc(S, S_sc, sc_grid, fc)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    real(dp), intent(in) :: fc(S%nat3,S%nat3,PRODUCT(sc_grid))
    integer :: R1, R2, nR, j1, j2, jn1, jn2, na1, na2, na1_sc, na2_sc, R(3)
    real(dp) :: fc_uc2sc(S_sc%nat3, S_sc%nat3)
    integer :: atoms_sc(S%nat, product(sc_grid))
    !
    nR = product(sc_grid)
    atoms_sc = map_uc2sc(S, S_sc, sc_grid)
    do R1 = 1, nR
      do R2 = 1, nR
        R = index2v(R2, sc_grid) - index2v(R1, sc_grid)
        do j1 = 1, 3
          if (R(j1) < 0) R(j1) = R(j1) + sc_grid(j1)
        enddo
        do na1 = 1, S%nat
          do na2 = 1, S%nat
            do j1 = 1, 3
              jn1 = j1 + 3*(na1-1)
              do j2 = 1, 3
                jn2 = j2 + 3*(na2-1)
                na1_sc = atoms_sc(na1, R1)
                na2_sc = atoms_sc(na2, R2)

                fc_uc2sc(j1 + 3*(na1_sc-1), j2 + 3*(na2_sc-1)) = &
                  fc(jn1, jn2, v2index(R, sc_grid))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function fc_uc2RR(sc_grid, S, fc)
    TYPE(ph_system_info), intent(in) :: S
    integer, intent(in) :: sc_grid(3)
    real(dp), intent(in) :: fc(S%nat3,S%nat3,PRODUCT(sc_grid))
    integer :: R1, R2, nR, j1, R(3)
    real(dp) :: fc_uc2RR(S%nat3, S%nat3, product(sc_grid), product(sc_grid))
    !
    nR = product(sc_grid)
    do R1 = 1, nR
      do R2 = 1, nR
        R = index2v(R1, sc_grid) - index2v(R2, sc_grid)
        do j1 = 1, 3
          if (R(j1) < 0) R(j1) = R(j1) + sc_grid(j1)
        enddo
        fc_uc2RR(:,:,R2,R1) = fc(:,:,v2index(R, sc_grid))
      enddo
    enddo
  end function
  !
  SUBROUTINE interp_1st_step(fc, xq, S)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    TYPE(ph_system_info),INTENT(in) :: S
    !
    INTEGER :: R1, R2
    REAL(DP), allocatable :: varg(:), vcos(:), vsin(:)
    COMPLEX(DP), allocatable :: vphase(:)
    !
    if(allocated(fc%mix)) deallocate(fc%mix)
    !
    allocate(fc%mix(S%nat3, S%nat3, fc%n_R1))
    fc%mix = 0._dp
    !
    do R1 = 1, fc%n_R1
      allocate(varg(fc%n_R2(R1)), vcos(fc%n_R2(R1)), &
        vsin(fc%n_R2(R1)), vphase(fc%n_R2(R1)))
      FORALL(R2=1:fc%n_R2(R1)) varg(R2) =  tpi * &
        dot_product(xq, fc%xR2(:,R2,R1))
      !
      ! Pre-compute phase to use the vectorized MKL subroutines
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R2(R1), varg, vcos)
      CALL vdSin(fc%n_R2(R1), varg, vsin)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      !
      do R2 = 1, fc%n_R2(R1)
        fc%mix(:,:,R1) = fc%mix(:,:,R1) + vphase(R2) * fc%FC(:,:,R2,R1)
      enddo
      ! fc%mix(:,:,R1) = SUM(vphase * fc%FC(:,:,:,R1), dim=3)
      !
      deallocate(varg, vcos, vsin, vphase)
    enddo
    !
  END SUBROUTINE
  !
  SUBROUTINE interp_2nd_step(fc, xq, S, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    TYPE(ph_system_info),INTENT(in) :: S
    complex(dp), intent(out) :: D(S%nat3, S%nat3)
    !
    INTEGER :: i
    REAL(DP), dimension(fc%n_R1) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(fc%n_R1)
    !
    !
    FORALL(i=1:fc%n_R1) varg(i) =  tpi * &
      dot_product(xq, fc%xR1(:,i))
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
    CALL vdCos(fc%n_R1, varg, vcos)
    CALL vdSin(fc%n_R1, varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    D = 0._dp
    do i = 1, fc%n_R1
      D = D + vphase(i) * fc%mix(:,:,i)
    enddo
    !
    D = D / PRODUCT(fc%nq)
    !
  END SUBROUTINE
  !
  SUBROUTINE interp_at_once(fc, xq1, xq2, nat3, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc), INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq1(3), xq2(3)
    integer, intent(in) :: nat3
    complex(dp), intent(out) :: D(nat3, nat3)
    !
    INTEGER :: i, j
    REAL(DP), dimension(fc%n_R1,fc%n_R2(1)) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(fc%n_R1,fc%n_R2(1))
    !
    FORALL(i=1:fc%n_R1, j=1:fc%n_R2(1)) varg(i,j) = &
      tpi * (dot_product(xq1, fc%xR1(:,i)) + &
      dot_product(xq2, fc%xR2(:,j,1)))
    !
    vcos = DCOS(varg)
    vsin = DSIN(varg)
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    D = 0._dp
    do j = 1, fc%n_R2(1)
      do i = 1, fc%n_R1
        D = D + vphase(i,j) * fc%fc(:,:,i,j)
      enddo
    enddo
    !
  END SUBROUTINE
  !
  subroutine build_mass_ratios(S, Sd, grid, mass_def, iR_def, na_def)
    type(ph_system_info), intent(in):: S, Sd
    integer, intent(in) :: grid(3)
    real(dp), intent(out) :: mass_def
    integer, intent(out) :: iR_def, na_def
    !
    integer :: which_iR(S%nat)
    real(dp):: mass_ratio(S%nat)
    integer :: minority, majority, tot_defects
    integer :: map(S%nat, product(grid)), typ(product(grid))
    integer :: na, iR, nR
    !
    nR = product(grid)
    map = map_uc2sc(S, Sd, grid)
    mass_ratio = 0._dp
    do na = 1, S%nat
      do iR = 1, nR
        typ(iR) = Sd%ityp(map(na, iR))
      enddo
      if (all(typ == typ(1))) then
        which_iR(na) = -1
        cycle
      endif
      minority = typ(1)
      which_iR(na) = 1
      do iR = 2, nR-1
        if (Sd%ityp(map(na, iR)) == minority) then
          minority = typ(iR+1)
          which_iR(na) = iR+1
        endif
      enddo
      majority = (SUM(typ)-minority)/(nR-1)
      mass_ratio(na) = (Sd%amass(minority) - Sd%amass(majority)) / &
        (Sd%amass(majority))
    enddo
    !
    tot_defects = 0
    do na = 1, S%nat
      if (which_iR(na) < 0) cycle
      tot_defects = tot_defects + 1
      iR_def = which_iR(na)
      na_def = na
      mass_def = mass_ratio(na)
    enddo
    if (tot_defects /= 1) &
      CALL errore("build_mass_ratios", "there should be only one defect", ABS(tot_defects))
  end subroutine
  !
  SUBROUTINE div_mass_fcsc(S, Sd,fc)
    USE kinds, only : DP
    use thutils, only: bz2simple
    IMPLICIT NONE
    TYPE(forceconst2_sc) :: fc
    TYPE(ph_system_info)   :: S, Sd
    !
    INTEGER :: i, j, iR1, iR2, map(S%nat,product(fc%nq)), si, sj
    !
    IF(.not.ALLOCATED(Sd%sqrtmm1)) &
      call errore('div_mass_fc2', 'missing sqrtmm1, call aux_system first', 1)

    map = map_uc2sc(S, Sd, fc%nq)
    DO iR1 = 1, fc%n_R1
      si = v2index(bz2simple(fc%yR1(:,iR1), fc%nq), fc%nq)
      do iR2 = 1, fc%n_R2(iR1)
        sj = v2index(bz2simple(fc%yR2(:,iR2,iR1), fc%nq), fc%nq)
        DO j = 1, S%nat3
          DO i = 1, S%nat3
            fc%FC(i, j, iR2, iR1) = fc%FC(i, j, iR2, iR1) * Sd%sqrtmm1(map((i-1)/3+1,si))*Sd%sqrtmm1(map((j-1)/3+1,sj))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE
  !
  subroutine matd2RR(fc3, S, fcsc)
    use fc3_interpolate, only : sparse
    type(sparse), intent(in) :: fc3
    type(ph_system_info), intent(in) :: S
    integer :: mesh_cube, len2_new, i2_new, i3_new, i, j
    integer :: minimum, maximum, mesh(3)
    integer, dimension(size(fc3%yR2, 2)) :: iR2, iR3
    integer, allocatable :: ind2(:), ind3(:,:), len3_new(:)
    type(forceconst2_sc), intent(out) :: fcsc
    !
    minimum = min(minval(fc3%yR2), minval(fc3%yR3))
    maximum = max(maxval(fc3%yR2), maxval(fc3%yR3))
    mesh = maximum - minimum + 1
    mesh_cube = product(mesh)
    !
    do i = 1, fc3%n_R
      iR2(i) = v2index(fc3%yR2(:,i) - minimum, mesh)
      iR3(i) = v2index(fc3%yR3(:,i) - minimum, mesh)
    end do
    !
    allocate(ind2(mesh_cube), ind3(mesh_cube, mesh_cube))
    allocate(len3_new(mesh_cube))
    !
    ind2 = -1
    ind3 = -1
    len2_new = 0
    len3_new = 0
    do i = 1, fc3%n_R
      call add_ind(ind2, iR2(i), len2_new, i2_new)
      call add_ind(ind3(:,i2_new), iR3(i), len3_new(i2_new), i3_new)
    enddo
    !
    allocate(fcsc%n_R2(len2_new))
    fcsc%n_R2 = len3_new(:len2_new)
    fcsc%n_R1 = len2_new
    !
    allocate(fcsc%yR1(3,len2_new))
    allocate(fcsc%xR1(3,len2_new))
    allocate(fcsc%yR2(3,maxval(len3_new),len2_new))
    allocate(fcsc%xR2(3,maxval(len3_new),len2_new))
    allocate(fcsc%fc(S%nat3,S%nat3,maxval(len3_new),len2_new))
    !
    ind2 = -1
    ind3 = -1
    len2_new = 0
    len3_new = 0
    do i = 1, fc3%n_R
      call add_ind(ind2, iR2(i), len2_new, i2_new)
      fcsc%yR1(:,i2_new) = fc3%yR2(:,i)
      call add_ind(ind3(:,i2_new), iR3(i), len3_new(i2_new), i3_new)
      fcsc%yR2(:,i3_new,i2_new) = fc3%yR3(:,i)
      do j = 1, fc3%n_terms(i)
        fcsc%fc(fc3%dat(i)%idx(2,j),fc3%dat(i)%idx(3,j),i3_new,i2_new) = &
          fc3%dat(i)%fc(j)
      enddo
    enddo
    !
    call fcsc%cart(S, 1)
    call fcsc%cart(S, 2)
    fcsc%nq = fc3%nq
    !
    deallocate(ind2, ind3, len3_new)
  end subroutine
  !
  SUBROUTINE freq_phq_degen(xq, S, fc2, freq, V1, U)
    use fc2_interpolate, only : mat2_diag, fftinterp_mat2
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants,          ONLY : RY_TO_CMM1
    use thutils, only: near
    IMPLICIT NONE
    REAL(DP),INTENT(in)               :: xq(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    complex(dp), intent(in)           :: V1(S%nat3,S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    REAL(DP),PARAMETER :: epsq = 1.e-8_dp
    complex(dp), allocatable :: Udeg(:,:), Vdeg(:,:)
    real(dp), allocatable :: dummy(:)
    REAL(DP) :: cq(3), chk(3)
    LOGICAL :: gamma
    integer :: i,j, dim_deg
    !
    ! RAF
    !U = CONJG(U)
    CALL fftinterp_mat2(xq, S, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    cq = xq
    CALL cryst_to_cart(1,cq,S%at,-1)
    gamma = ALL( ABS(cq-NINT(cq))<epsq)
    IF( gamma )THEN
      freq(1:3) = 0._dp
      U(:,1:3) = (0._dp, 0._dp)
    ENDIF

    !> I lift the degeneracy by diagonalizing the perturbation
    !> we obtain new kets and frequencies, needed for 2nd order (FGR)
    do i = 1, S%nat3-1
      do j = i+1, S%nat3
        if (abs(freq(i)-freq(j)) > 1.d-8) exit
      enddo
      dim_deg = j-i
      if(dim_deg > 1) then
        allocate(Udeg(S%nat3,  dim_deg))
        allocate(Vdeg(dim_deg, dim_deg))
        allocate(dummy(dim_deg))
        do j = 1, dim_deg
          Udeg(:,j) = U(:,i+j-1)
        enddo

        Vdeg = matmul(TRANSPOSE(CONJG(Udeg)), matmul(V1, Udeg))
        call mat2_diag(dim_deg, Vdeg, dummy) !freq(i:i+dim_deg-1))
        U(:,i:i+dim_deg-1) = matmul(Udeg, Vdeg)
        if(all(abs(dummy -dummy(1)) < 1e-6 * dummy(1))) then

        endif
        ! freq(i:i+dim_deg-1) = freq(i:i+dim_deg-1) + dummy
        deallocate(Udeg, Vdeg, dummy)
      endif
    enddo


    chk(:) = cq(:)*fc2%nq(:)
    chk=chk-NINT(chk(:))
    IF(fc2%periodic .AND. SUM(ABS(chk)) > 1.d-8) CALL errore("interp","cannot interpolate with periodic matrices",1 )

    IF(ANY(freq<0._dp)) THEN
      WRITE(*,*) gamma
      WRITE(*,"(e12.5,e12.5,e12.5)") cq
      WRITE(*,"(e12.5,e12.5,e12.5)") xq
      WHERE    (freq >  0.0)
        freq = DSQRT(freq)
      ELSEWHERE(freq < 0.0)
        freq = -DSQRT(-freq)
      ENDWHERE
      WRITE(*,"('negative freq = ',12e12.4)")freq*RY_TO_CMM1
      CALL errore("freq_phq_safe", "cannot continue with negative frequencies",1)
    ELSE
      !
      freq = DSQRT(freq)
    END IF
    !
  END SUBROUTINE
  !
  subroutine freq_in_grid_degen(S, fc2, fc2sc, grid, freqs, Us)
    use q_grids, only : q_grid
    use mpi_thermal, only : mpi_bsum
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc) :: fc2sc
    !
    complex(dp) :: V1(S%nat3,S%nat3)
    real(dp), intent(out):: freqs(S%nat3, grid%nqtot)
    complex(dp), intent(out), optional :: Us(S%nat3, S%nat3, grid%nqtot)
    !
    integer :: iq, iqp
    !
    freqs = 0.0_dp
    DO iq = 1, grid%nq
      iqp = iq + grid%iq0
      call fc2sc%interpolate(grid%xq(:,iq), S)
      call fc2sc%interpolate(grid%xq(:,iq), S, V1)
      CALL freq_phq_degen(grid%xq(:,iq), S, fc2, freqs(:,iqp), V1, Us(:,:,iqp))
    END DO
    if (grid%scattered) CALL mpi_bsum(S%nat3, grid%nqtot, freqs)
  end subroutine
end module
