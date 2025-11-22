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
    INTEGER :: n_R2 = 0
    INTEGER, allocatable :: n_R1(:)
    !! indices : (R2)
    real(dp), allocatable :: FC(:,:,:,:)
    !! indices : (jn1,jn2,R1,R2)
    complex(dp), allocatable :: mix(:,:,:)
    !! indices : (jn1,jn2,R2)
    integer,  allocatable :: yR1(:,:,:)
    !! indices : (3,R1,R2)
    real(dp), allocatable :: xR1(:,:,:)
    !! indices : (3,R1,R2)
    integer,  allocatable :: yR2(:,:)
    !! indices : (3,R1)
    real(dp), allocatable :: xR2(:,:)
    !! indices : (3,R1)
    integer               :: nq(3)
    !! sc size
    integer :: stage = -1
    !! centered or not
    real(dp) :: taudef(3) = 0._dp
    !! tau of the defect atom
    integer :: iR1 = 0
    !! used for q2r
    integer, allocatable :: defects(:,:)
    !! indices : (3, n_defects) (nat, R, isc)
    !!
    real(dp), allocatable :: mass_ratios(:)
    !
    complex(dp), allocatable :: inclusion_Dnx_out(:,:,:)
    !! only for inclusion, they are Dnx in (Dnn Dnx Dxn Dxx), where x is inclusion and n is normal atom.
    !! indices : (3*nat, 3*n_add, q_out)
    complex(dp), allocatable :: inclusion_Dnx_in(:,:,:)
    !! only for inclusion, they are Dnx in (Dnn Dnx Dxn Dxx), where x is inclusion and n is normal atom.
    !! indices : (3*nat, 3*n_add, q_in)
    real(dp), allocatable :: inclusion_eig(:)
    !! only for inclusion, they are eigenvalues of Dxx
  contains
    procedure :: allocate => allocate_fc2_sc
    procedure :: deallocate => deallocate_fc2_sc
    ! procedure :: center => center_sc
    procedure :: cryst => construct_cryst
    procedure :: center_full => center2
    procedure :: center => center5
    procedure :: cart => construct_cart
    generic   :: r2q => r2q_1st_step, r2q_2nd_step, r2q_at_once
    ! generic   :: q2r => q2r_1st_step, q2r_2nd_step
    ! PROCEDURE, PASS :: q2r_1st_step
    ! PROCEDURE, PASS :: q2r_2nd_step
    PROCEDURE, PASS :: r2q_1st_step
    PROCEDURE, PASS :: r2q_2nd_step
    PROCEDURE, PASS :: r2q_at_once
  end type
  !
contains
  !
  !> map_atm_sc(inat, iR) = inat_sc
  function map_uc2sc(S, S_sc, sc_grid)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    ! real(dp), intent(out), optional :: taudef(3)
    integer :: ndef
    !
    logical :: idef(S_sc%nat)
    integer :: n_vac
    integer :: map_uc2sc(S%nat, PRODUCT(sc_grid))
    integer :: i, isc, iR
    real(dp) :: r_cryst(3)
    map_uc2sc = -1
    idef = .false.
    ndef = 0
    ! if(present(taudef)) taudef = 0._dp
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        r_cryst = S_sc%tau(:,isc)*sc_grid-S%tau(:,i)
        call cryst_to_cart(1, r_cryst, S_sc%bg, -1)
        if (NORM2(r_cryst - NINT(r_cryst))<1e-1_dp) then
          iR = v2index(NINT(r_cryst), sc_grid)
          if (map_uc2sc(i, iR) > -1) CALL errore("map_uc2sc", "found same atom in same R", 1)
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("map_uc2sc", "R is out of bound", ABS(iR))
          map_uc2sc(i, iR) = isc
          ! if(S%atm(S%ityp(i)) /= S_sc%atm(S_sc%ityp(isc))) idef(isc) = .true.
          ! if(present(taudef)) taudef = taudef + S_sc%tau(:,isc)
          ! ndef = ndef + 1
        endif
      enddo
    enddo
    ! if(present(taudef)) taudef = taudef / REAL(ndef, DP)
    ! print"(A,I3.3,A)", "found ", ndef, " defect atoms in the supercell"
    n_vac = COUNT(map_uc2sc == -1)
    if (n_vac + S_sc%nat /= S%nat * product(sc_grid)) then
      print*, map_uc2sc
      CALL errore("map_uc2sc", "some atoms are not mapped", 1)
    endif

  end function
!
  subroutine find_defects(S, S_sc, fc, sites)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    CLASS(forceconst2_sc), intent(inout) :: fc
    integer, allocatable, intent(in) :: sites(:)
    !
    integer :: map(S%nat, product(fc%nq))
    integer :: i, iR, ndef
    integer, allocatable :: defects_tmp(:,:)
    real(dp), allocatable :: mass_ratios_tmp(:)
    !
    map = map_uc2sc(S, S_sc, fc%nq)
    !
    allocate(fc%defects(3, S%nat * product(fc%nq)))
    allocate(fc%mass_ratios(S%nat * product(fc%nq)))
    !
    ndef = 0
    do concurrent(iR=1:product(fc%nq), i=1:S%nat, &
      S_sc%amass(S_sc%ityp(map(i,iR))) /= S%amass(S%ityp(i)))
      ndef = ndef + 1
      fc%defects(:, ndef) = [i, iR, map(i,iR)]
      fc%mass_ratios(ndef) = - (S_sc%amass(S_sc%ityp(map(i,iR))) - S%amass(S%ityp(i))) / S%amass(S%ityp(i))
    enddo
    !
    if(allocated(sites) .and. ndef == 1) then
      do i = 1, size(sites)
        ndef = ndef + 1
        fc%defects(:, ndef) = [sites(i), fc%defects(2,1), map(sites(i), fc%defects(2,1))]
        fc%mass_ratios(ndef) = fc%mass_ratios(1)
      enddo
    endif
    !
    allocate(defects_tmp(3, ndef))
    allocate(mass_ratios_tmp(ndef))
    defects_tmp = fc%defects(:, :ndef)
    mass_ratios_tmp = fc%mass_ratios(:ndef)
    call move_alloc(defects_tmp, fc%defects)
    call move_alloc(mass_ratios_tmp, fc%mass_ratios)
    !
    fc%taudef = 0._dp
    do i = 1, ndef
      fc%taudef = fc%taudef + S%tau(:, fc%defects(1,i))
    enddo
    fc%taudef = fc%taudef / REAL(ndef, DP)
    !
  end subroutine
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
        if (NORM2(r_cryst - NINT(r_cryst))<1e-1) then
          if (map_sc2uc(isc) /= -1) &
            call errore("map_sc2uc", "We already have found this defect", isc)
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
    if (COUNT(map_sc2uc == -1) + S%nat * product(sc_grid) /= S_sc%nat) then
      print*, "grid is", sc_grid
      print"(100I5.2)", map_sc2uc
      CALL errore("map_sc2uc", "some atoms are not mapped", 1)
    endif
  end function
  !
  SUBROUTINE allocate_fc2_nodef(fc, S, grid)
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
    fc%n_R2 = n_R
    allocate(fc%n_R1(n_R))
    fc%n_R1 = n_R
    !
    ALLOCATE(fc%yR1(3,n_R,n_R), fc%yR2(3,n_R))
    fc%yR2 = grid_vec(grid)
    do iR = 1, n_R
      fc%yR1(:,:,iR) = fc%yR2
    enddo
    !
    ALLOCATE(fc%xR1(3,n_R,n_R), fc%xR2(3,n_R))
    call fc%cart(S, 1)
    call fc%cart(S, 2)
    ALLOCATE(fc%FC(S%nat3,S%nat3,n_R,n_R))
    fc%nq = grid
    fc%stage = -1
    fc%fc = 0._dp
    !
  END SUBROUTINE
  !
  SUBROUTINE allocate_fc2_sc(fc, S, S_sc, grid, sites)
    use thutils, only : grid_vec
    IMPLICIT NONE
    INTEGER,INTENT(in) :: grid(3)
    type(ph_system_info),INTENT(in) :: S
    type(ph_system_info),INTENT(in) :: S_sc
    integer, allocatable, intent(in) :: sites(:)
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    CHARACTER(len=16),PARAMETER :: sub = "allocate_fc2_sc"
    integer :: n_R, iR
    !
    IF(allocated(fc%yR2) .or. allocated(fc%FC) &
      .or. allocated(fc%xR2) .or. allocated(fc%mix)) &
      CALL errore(sub, 'some element is already allocated', 1)
    !
    n_R = PRODUCT(grid)
    fc%n_R2 = n_R
    allocate(fc%n_R1(n_R))
    fc%n_R1 = n_R
    !
    ALLOCATE(fc%yR1(3,n_R,n_R), fc%yR2(3,n_R))
    fc%yR2 = grid_vec(grid)
    do iR = 1, n_R
      fc%yR1(:,:,iR) = fc%yR2
    enddo
    !
    ALLOCATE(fc%xR1(3,n_R,n_R), fc%xR2(3,n_R))
    call fc%cart(S, 1)
    call fc%cart(S, 2)
    ALLOCATE(fc%FC(S%nat3,S%nat3,n_R,n_R))
    fc%nq = grid
    fc%stage = -1
    fc%fc = 0._dp
    !
    call find_defects(S, S_sc, fc, sites)
    !
  END SUBROUTINE
  !
  subroutine deallocate_fc2_sc(fc)
    class(forceconst2_sc), intent(inout) :: fc
    !
    if (allocated(fc%yR2)) deallocate(fc%yR2)
    if (allocated(fc%xR2)) deallocate(fc%xR2)
    if (allocated(fc%FC)) deallocate(fc%FC)
    if (allocated(fc%mix)) deallocate(fc%mix)
    if (allocated(fc%yR1)) deallocate(fc%yR1)
    if (allocated(fc%xR1)) deallocate(fc%xR1)
    if (allocated(fc%n_R1)) deallocate(fc%n_R1)
    fc%n_R2 = 0
    fc%taudef = 0._dp
    fc%stage = -1
    fc%nq = 0
  end subroutine
  !
  subroutine construct_cart(fc, S, which)
    class(forceconst2_sc), intent(inout) :: fc
    type(ph_system_info), intent(in) :: S
    integer, intent(in) :: which
    !
    real(dp), allocatable :: R(:,:)
    integer :: i
    !
    if(which == 1) then
      do i = 1, fc%n_R2
        allocate(R(3,fc%n_R1(i)))
        R = REAL(fc%yR1(:,:fc%n_R1(i),i), DP)
        call cryst_to_cart(fc%n_R1(i), R, S%at, 1)
        fc%xR1(:,:fc%n_R1(i),i) = R
        deallocate(R)
      enddo
      !
    elseif(which == 2) then
      allocate(R(3,fc%n_R2))
      R = REAL(fc%yR2, DP)
      call cryst_to_cart(fc%n_R2, R, S%at, 1)
      fc%xR2 = R
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
    if (which == 1) then
      do i = 1, fc%n_R2
        allocate(R(3,fc%n_R1(i)))
        R = fc%xR1(:,:fc%n_R1(i),i)
        call cryst_to_cart(fc%n_R1(i), R, S%bg, -1)
        if (ALL(ABS(R - NINT(R)) < 1e-6)) then
          fc%yR1(:,:fc%n_R1(i),i) = NINT(R)
        else
          call errore("construct_cryst", "R is not integer", 1)
        endif
        deallocate(R)
      enddo
      !
    elseif(which == 2) then
      allocate(R(3,fc%n_R2))
      R = fc%xR2
      call cryst_to_cart(fc%n_R2, R, S%bg, -1)
      if (ALL(ABS(R - NINT(R)) < 1e-6)) then
        fc%yR2 = NINT(R)
      else
        call errore("construct_cryst", "R is not integer", 1)
      endif
      deallocate(R)
    else
      call errore("construct_cryst", "which is not 1 or 2", 1)
    endif
  end subroutine
  !
  subroutine S_uc2sc(S, Sd, grid, S_sc)
    use thutils, only : grid_vec
    use ph_system, only : aux_system
    type(ph_system_info), intent(in) :: S, Sd
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(out) :: S_sc
    !
    integer :: nR, iR, i, t
    integer, allocatable :: sc_grid(:,:)
    real(dp) :: tau(3,S%nat)
    integer :: map(S%nat, product(grid)), na
    !
    nR = PRODUCT(grid)
    S_sc%ntyp = S%ntyp
    S_sc%amass = S%amass
    S_sc%amass_variance = S%amass_variance
    S_sc%atm = S%atm
    S_sc%nat = S%nat * nR
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
    S_sc%lrigid = .false. ! I'm not caring about effective charges up to now S%lrigid
    allocate(S_sc%tau(3, S_sc%nat))
    tau = S%tau
    call cryst_to_cart(S%nat, tau, S%bg, -1)
    !
    map = map_uc2sc(S, Sd, grid)
    do iR = 1, nR
      do na = 1, S%nat
        S_sc%tau(:,map(na, iR)) = (tau(:,na) + index2v(iR, grid))/grid
      enddo
    enddo
    allocate(S_sc%ityp(S_sc%nat))
    S_sc%ityp = Sd%ityp
    call cryst_to_cart(S_sc%nat, S_sc%tau, S%at, 1)
    !
    call aux_system(S_sc)
  end subroutine
  !
  ! subroutine center_sc(fc, grid, S, S_sc, distance)
  !   use functions, only : refold_bz
  !   class(forceconst2_sc), intent(inout) :: fc
  !   integer, intent(in) :: grid(3)
  !   type(ph_system_info), intent(in) :: S, S_sc
  !   real(dp), allocatable, optional, intent(out) :: distance(:,:,:,:)
  !   !
  !   ! type(forceconst2_grid) :: fsc
  !   real(dp) :: dist(3), Rbig(3), wg_tot
  !   !
  !   integer :: na1, na2, j1, j2, nR, jR_big, R1, R2
  !   integer :: R_list(PRODUCT(grid)*(2*nfar+1)**3)
  !   integer :: map_sc(S%nat, PRODUCT(grid))
  !   integer :: nRbig, R_vec(3)
  !   integer :: ind, ixR
  !   integer :: nxR(PRODUCT(grid))
  !   integer, dimension(3) :: far_grid, Rbig_from_0, Rbig_shift
  !   !
  !   integer, allocatable :: new_yR_list(:,:,:)
  !   real(dp), allocatable :: new_fc(:,:,:,:)
  !   integer :: counter
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws
  !   REAL(DP) :: wg, rws(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   ! initialize WS r-vectors
  !   CALL wsinit(rws,nrwsx,nrws,S_sc%at)
  !   !
  !   fc%stage = 0
  !   nR = PRODUCT(grid)
  !   if (nfar == 0) return
  !   far_grid = 2*nfar+1
  !   nRbig = PRODUCT(far_grid)
  !   allocate(new_yR_list(3, nR*nRbig,nR))
  !   allocate(new_fc(S%nat3, S%nat3, nR*nRbig, nR))
  !   !
  !   ! if(present(distance)) allocate(distance, source=new_fc)
  !   map_sc = map_uc2sc(S, S_sc, grid)
  !   !
  !   nxR = 0
  !   new_fc = 0._dp
  !   !
  !   counter = 0
  !   do R1 = 1, nR
  !     R_list = -1
  !     do R2 = 1, nR
  !       do na1 = 1, S%nat
  !         do na2 = 1, S%nat
  !           wg_tot = 0._dp
  !           do jR_big = 1, nRbig
  !             !> I create a normal [0:N-1]^3 grid
  !             Rbig_from_0 = index2v(jR_big, far_grid)
  !             !> Then I shift it so that the center in [-N/2:N/2]
  !             Rbig_shift = Rbig_from_0 - nfar
  !             !> I work in S_sc%at alat units
  !             Rbig = REAL(Rbig_shift, DP)
  !             call cryst_to_cart(1, Rbig, S_sc%at, 1)
  !             !> Rbig is the R vector in the supercell, so we don't need to multiply by grid,
  !             !> cause tau is in [0,1] in crystal units
  !             dist = S_sc%tau(:,map_sc(na1,R1)) - Rbig - S_sc%tau(:,map_sc(na2,R2))
  !             !> Compute the Wigner-Seitz weight
  !             wg = wsweight(dist,rws,nrws)
  !             wg_tot = wg_tot + wg
  !             ! if (nfar == 0) wg = 1._dp
  !             if (wg /= 0) then
  !               !> R2 is in the unit cell, so we multiply Rbig by grid to transform it in a
  !               !> supercell vector of the unit cell
  !               R_vec = - index2v(R2, grid) - grid * Rbig_shift
  !               !> I use the 0-indexing to populate R_list at the correct non-negative integer
  !               ind = v2index(index2v(R2, grid) + grid * Rbig_from_0, grid * far_grid)
  !               !> R_list contains the ixR or -1. The nxR is refreshed at each step
  !               if (R_list(ind) == -1) then
  !                 nxR(R1) = nxR(R1) + 1
  !                 ixR = nxR(R1)
  !                 R_list(ind) = ixR
  !                 new_yR_list(:,ixR,R1) = R_vec
  !               else
  !                 ixR = R_list(ind)
  !               endif
  !               !> I find Gamma Gamma
  !               ! if(ALL(R_vec == 0)) fc%i_0 = ixR
  !               !> and populate force constants with usual index wrapping
  !               do j1 = 1, 3
  !                 do j2 = 1, 3
  !                   new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR, R1) = &
  !                     fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, R2, R1) * wg
  !                   ! if (present(distance)) &
  !                   !   distance(j1+3*(na1-1), j2+3*(na2-1), ixR, R1) = norm2(dist)
  !                 enddo
  !               enddo
  !             endif
  !           enddo
  !           if (ABS(wg_tot -1) > 1e-6) then
  !             print"(A,F14.6)", "wg_tot is", wg_tot
  !             CALL errore("center_sc", "sum of weights is not 1", 1)
  !           endif
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   !
  !   deallocate(fc%yR2, fc%xR2, fc%FC)
  !   !> the yR list is populated until nxR(which), but the size can be larger.
  !   !> the values after nxR(which) are not even initialized (they are garbage).
  !   ALLOCATE(fc%yR2(3,maxval(nxR),fc%n_R1))
  !   ALLOCATE(fc%xR2(3,maxval(nxR),fc%n_R1))
  !   ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR),nR))
  !   fc%n_R2 = nxR
  !   fc%nq = grid
  !   do R1 = 1, nR
  !     ! fc%xr1(:,R1) = refold_bz(fc%xR1(:,R1), S%at)
  !     fc%yR2(:,:nxR(R1),R1) = new_yR_list(:,:nxR(R1),R1)
  !   enddo
  !   call fc%cart(S_sc, 2)
  !   ! call fc%cryst(S_sc, 1)
  !   fc%FC = new_fc(:,:,:maxval(nxR),:)
  ! end subroutine
  !
  subroutine center2(fc, grid, S, S_sc)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    !
    ! type(forceconst2_grid) :: fsc
    real(dp), dimension(3) :: d2, d1
    real(dp) :: perix, peri_min, max_norm
    real(dp), parameter :: eps_peri = 1e-4_dp
    integer, parameter :: nperix = (2*nfar+1)**3
    integer :: nperi
    integer :: SAFE_ALLOCATION
    !
    integer :: na2, na1, j1, j2, nR, R2, R1, R2_big, R1_big
    integer :: map_sc(S%nat, PRODUCT(grid))
    integer, allocatable :: far_grid_cryst(:,:)
    real(dp), allocatable :: far_grid_cart(:,:)
    integer :: nRbig, counter
    integer :: iperi, nxR2, ixR2, ixR1
    integer, allocatable :: R2_list(:), R1_list(:,:), yR2_list(:,:), yR1_list(:,:,:), nxR1(:)
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
    far_grid_cryst = grid_vec_cryst(far_mesh, .true.)
    far_grid_cart = grid_vec_cart(far_mesh, S_sc%at, .true.)
    !
    SAFE_ALLOCATION = 10 * nR
    allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(R2_list(nR*nRbig))
    allocate(R1_list(nR*nRbig,nR*nRbig))
    allocate(yR2_list(3, nR*nRbig))
    allocate(yR1_list(3, nR*nRbig, nR*nRbig))
    allocate(nxR1(nR*nRbig))
    !
    map_sc = map_uc2sc(S, S_sc, grid)
    !
    new_fc = 0._dp
    R2_list = -1
    R1_list = -1
    nxR1 = 0
    nxR2 = 0
    counter = 0
    !
    ! open(10, file='distance-defect.dat', status='replace')
    do R2 = 1, nR
      do na2 = 1, S%nat
        do R1 = 1, nR
          do na1 = 1, S%nat
            nperi = 0
            peri_min = 0._dp
            do R2_big = 1, nRbig
              d2 = S_sc%tau(:,map_sc(na2,R2)) + far_grid_cart(:,R2_big)
              do R1_big = 1, nRbig
                d1 = S_sc%tau(:,map_sc(na1,R1)) + far_grid_cart(:,R1_big)
                perix = norm2(d2 - d1) + ((norm2(d2 + d1 - fc%taudef))) * 1e-3_dp
                IF (perix < peri_min-eps_peri .or. nperi==0) THEN
                  nperi = 1
                  farx_list = 0
                  ind = 0
                  peri_min = perix
                  farx_list(:,1,nperi) = index2v(R1, grid) + grid * far_grid_cryst(:,R1_big)
                  farx_list(:,2,nperi) = index2v(R2, grid) + grid * far_grid_cryst(:,R2_big)
                  ind(1,nperi) = R1 + nR*(R1_big-1)
                  ind(2,nperi) = R2 + nR*(R2_big-1)
                ELSE IF ( ABS(perix-peri_min) <= eps_peri ) THEN
                  nperi = nperi + 1
                  IF(nperi > nperix) CALL errore("center2", "nperix is too small", nperi)
                  peri_min = (peri_min*(nperi-1)+perix)/DBLE(nperi)
                  farx_list(:,1,nperi) = index2v(R1, grid) + grid * far_grid_cryst(:,R1_big)
                  farx_list(:,2,nperi) = index2v(R2, grid) + grid * far_grid_cryst(:,R2_big)
                  ind(2,nperi) = R2 + nR*(R2_big-1)
                  ind(1,nperi) = R1 + nR*(R1_big-1)
                END IF
              enddo
            enddo
            if (nperi > 1) counter = counter + 1
            !
            do iperi = 1, nperi
              call add_ind(R2_list, ind(2,iperi), nxR2, ixR2)
              yR2_list(:,ixR2) = farx_list(:,2,iperi)
              !
              call add_ind(R1_list(:,ixR2), ind(1,iperi), nxR1(ixR2), ixR1)
              yR1_list(:,ixR1,ixR2) = farx_list(:,1,iperi)
              !> and populate force constants with usual index wrapping
              R1_big = ind(1,iperi)/ nR + 1
              R2_big = ind(2,iperi)/ nR + 1
              d2 = S_sc%tau(:,map_sc(na2,R2)) + far_grid_cart(:,R2_big)
              d1 = S_sc%tau(:,map_sc(na1,R1)) + far_grid_cart(:,R1_big)
              ! max_norm = MAX(max_norm, norm2(d1-d2))
              if (all(fc%fc(3*(na1-1)+1:3*na1, 3*(na2-1)+1:3*na2, R1, R2) == 0._dp)) &
                call errore("center2", "force constants are zero in a point", 1)
              do j1 = 1, 3
                do j2 = 1, 3
                  new_fc(j1+(na1-1)*3, j2+(na2-1)*3, ixR1, ixR2) = &
                    fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, R1, R2) / nperi
                  ! write(10, "(4E15.5)") norm2(d1-fc%taudef), norm2(d2-fc%taudef), &
                  !   norm2(d1-d2), new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR1, ixR2)
                enddo
              enddo
            enddo
            !
          enddo
        enddo
      enddo
    enddo
    ! close(10)
    !
    deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR2(3,nxR2))
    ALLOCATE(fc%xR2(3,nxR2))
    ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
    ALLOCATE(fc%n_R1(nxR2))
    fc%n_R1 = nxR1(:nxR2)
    fc%n_R2 = nxR2
    fc%nq = grid
    do ixR2 = 1, nxR2
      fc%yR2(:,ixR2) = yR2_list(:,ixR2)
      fc%yR1(:,:nxR1(ixR2),ixR2) = yR1_list(:,:nxR1(ixR2),ixR2)
      fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
    enddo
    call fc%cart(S_sc, 1)
    call fc%cart(S_sc, 2)

    call print_message("end of centering")

    deallocate(new_fc, R2_list, R1_list, yR2_list, yR1_list, nxR1)
  end subroutine
  !
  subroutine center4(fc, grid, S)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S
    !
    ! type(forceconst2_grid) :: fsc
    real(dp), dimension(3) :: d1, d2
    real(dp), parameter :: eps_peri = 1e-3
    integer, parameter :: nperix = (2*nfar+1)**3
    integer :: SAFE_ALLOCATION
    !
    integer :: na1, na2, j1, j2, nR, iR1, iR2, iR1_big, iR2_big, jn1, jn2
    real(dp), allocatable :: far_grid_cart(:,:)
    integer, allocatable :: far_grid_cryst(:,:)
    integer :: nRbig
    integer :: nxR2, ixR2, ixR1
    integer, allocatable :: R1_list(:,:), R2_list(:), nxR1(:)
    real(dp), allocatable :: xR2_list(:,:), xR1_list(:,:,:)
    real(dp), allocatable :: weights1(:), weights2(:)
    integer, allocatable :: inds1(:), inds2(:)
    integer, dimension(3) :: far_mesh
    real(dp) :: big_at(3,3), max_norm
    !
    real(dp), allocatable :: new_fc(:,:,:,:)
    character(100) :: filename
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    ! initialize WS r-vectors
    forall(j1 = 1:3) big_at(:,j1) = S%at(:,j1) * grid(j1)
    CALL wsinit(rws,nrwsx,nrws,big_at)
    !
    fc%stage = 0
    nR = PRODUCT(grid)
    if (nfar == 0) return
    far_mesh = 2*nfar+1
    nRbig = PRODUCT(far_mesh)
    !
    allocate(far_grid_cart(3, nRbig*nR))
    allocate(far_grid_cryst(3, nrbig*nr))
    far_grid_cart = grid_vec_cart(far_mesh*grid, S%at, center=.true.)
    far_grid_cryst = grid_vec_cryst(far_mesh*grid, center=.true.)
    !
    SAFE_ALLOCATION = 20 * nR
    allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(R2_list(nR*nRbig))
    allocate(R1_list(nR*nRbig, nR*nRbig))
    allocate(xR1_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(xR2_list(3, SAFE_ALLOCATION))
    allocate(nxR1(SAFE_ALLOCATION))
    !
    new_fc = 0._dp
    R1_list = -1
    R2_list = -1
    nxR1 = 0
    nxR2 = 0
    max_norm = 0._dp
    !
    write(filename, "(A,I1,A)") 'distance-center4-', grid(1), '.dat'
    open(10, file=trim(filename), status='replace')
    do na2 = 1, S%nat
      d2 = S%tau(:,na2) - fc%taudef
      call inside_ws(far_grid_cart, d2, nrws, rws, weights2, inds2)
      do iR2_big = 1, size(weights2)
        iR2 = inside_periodic(inds2(iR2_big), grid)
        call add_ind(R2_list, inds2(iR2_big), nxR2, ixR2)
        xR2_list(:,ixR2) = far_grid_cart(:,inds2(iR2_big))
        do na1 = 1, S%nat
          d1 = S%tau(:,na1) - fc%taudef
          call inside_ws(far_grid_cart, d1, nrws, rws, weights1, inds1)
          do iR1_big = 1, size(weights1)
            iR1 = inside_periodic(inds1(iR1_big), grid)
            call add_ind(R1_list(:,ixR2), inds1(iR1_big), nxR1(ixR2), ixR1)
            xR1_list(:,ixR1,ixR2) = far_grid_cart(:,inds1(iR1_big))
            ! if (iR1 == iR2 .and. any(xR2_list(:,ixR2) /= xR1_list(:,ixR1,ixR2))) then
            !   call cryst_to_cart(1, xR1_list(:,ixR1,ixR2), S%bg, -1)
            !   call cryst_to_cart(1, xR2_list(:,ixR2), S%bg, -1)
            !   print*, na1, na2, iR1, iR2, xR1_list(:,ixR1,ixR2), xR2_list(:,ixR2)
            !   print*, inds1(iR1_big)
            !   print*, inside_periodic(inds1(iR1_big), grid)
            !   print*, far_grid_cart(:,inds1(iR1_big))
            !   call errore("center4", "iR1 == iR2 but xR1_list /= xR2_list", 1)
            ! endif
            ! max_norm = max(max_norm, norm2(d1 + xR1_list(:,ixR1,ixR2) - d2 - xR2_list(:,ixr2)))
            ! if(norm2(d1 + xR1_list(:,ixR1,ixR2) - d2 - xR2_list(:,ixr2)) < 2._dp) then
            do j1 = 1, 3
              jn1 = 3*(na1-1) + j1
              do j2 = 1, 3
                jn2 = 3*(na2-1) + j2
                new_fc(jn1, jn2, ixR1, ixR2) = &
                  fc%FC(jn1, jn2, iR1, iR2) * &
                  weights1(iR1_big) * weights2(iR2_big)
                write(10, "(3E20.8)") &
                  norm2(d1 + xR1_list(:,ixR1,ixR2)), &
                  norm2(d2 + xR2_list(:,ixr2)), &
                  new_fc(jn1, jn2, ixR1, ixR2)
              enddo
            enddo
            ! endif
            !
          enddo
        enddo
      enddo
    enddo
    print*, "max_norm:", max_norm
    close(10)
    !
    deallocate(weights1, weights2)
    deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR2(3,nxR2))
    ALLOCATE(fc%xR2(3,nxR2))
    ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
    ALLOCATE(fc%n_R1(nxR2))
    fc%n_R1 = nxR1(:nxR2)
    fc%n_R2 = nxR2
    fc%nq = grid
    do ixR2 = 1, nxR2
      fc%xR2(:,ixR2) = xR2_list(:,ixR2)
      fc%xR1(:,:nxR1(ixR2),ixR2) = xR1_list(:,:nxR1(ixR2),ixR2)
      fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
    enddo
    call fc%cryst(S, 1)
    call fc%cryst(S, 2)
    !
    call print_message("end of centering R1 && R2")
    deallocate(new_fc, xR2_list, xR1_list, nxR1)
  contains
    function inside_periodic(iR, sc_grid) result(iRp)
      integer, intent(in) :: iR
      integer, intent(in) :: sc_grid(3)
      !
      integer :: iRp
      integer, dimension(3) :: R, Rp, grid_
      !> This function returns the periodic image of R in the [-nfar:nfar]^3 box
      !> It is used to compute the periodic image of the force constants
      !> in the supercell
      grid_ = (2*nfar + 1) * sc_grid
      R = index2v(iR, grid_)
      Rp = mod(R-grid_/2 + sc_grid*10, sc_grid)
      iRp = v2index(Rp, sc_grid)
      !
    end function inside_periodic
  end subroutine
  !
  subroutine center5(fc, grid, S)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S
    !
    ! type(forceconst2_grid) :: fsc
    real(dp), dimension(3) :: d1, d2
    real(dp), parameter :: eps_peri = 1e-3
    integer, parameter :: nperix = (2*nfar+1)**3
    integer :: SAFE_ALLOCATION
    !
    integer :: na1, na2, j1, j2, nR, iR1, iR2, iR1_big, iR2_big, jn1, jn2
    real(dp), allocatable :: sum_grid_cart(:,:), diff_grid_cart(:,:)
    integer, allocatable :: far_grid_cryst(:,:)
    integer :: nRbig
    integer :: nxR2, ixR2, ixR1
    integer, allocatable :: R1_list(:,:), R2_list(:), nxR1(:)
    integer, allocatable :: yR2_list(:,:), yR1_list(:,:,:)
    real(dp), allocatable :: weights1(:), weights2(:)
    integer, allocatable :: inds1(:), inds2(:)
    integer, dimension(3) :: far_mesh
    real(dp) :: big_at(3,3)
    real(dp), dimension(3) :: R1, R2
    !
    real(dp), allocatable :: new_fc(:,:,:,:)
    character(100) :: filename
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    ! initialize WS r-vectors
    forall(j1 = 1:3) big_at(:,j1) = S%at(:,j1) * grid(j1)
    CALL wsinit(rws,nrwsx,nrws,big_at)
    !
    fc%stage = 0
    nR = PRODUCT(grid)
    if (nfar == 0) return
    far_mesh = 2*nfar+1
    nRbig = PRODUCT(far_mesh)
    !
    ! allocate(sum_grid_cart(3, 8*nRbig*nR))
    ! allocate(diff_grid_cart(3, nRbig*nR))
    ! allocate(far_grid_cryst(3, nRbig))
    sum_grid_cart = grid_vec_cart(far_mesh*grid*2, S%at, center=.true.)/2
    diff_grid_cart = grid_vec_cart(far_mesh*grid, S%at, center=.true.)
    ! far_grid_cryst = grid_vec_cryst(far_mesh, center=.true.)
    !
    SAFE_ALLOCATION = 20 * nR
    allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(R2_list(nR*nRbig))
    allocate(R1_list(nR*nRbig, nR*nRbig))
    allocate(yR1_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(yR2_list(3, SAFE_ALLOCATION))
    allocate(nxR1(SAFE_ALLOCATION))
    !
    new_fc = 0._dp
    R1_list = -1
    R2_list = -1
    nxR1 = 0
    nxR2 = 0
    !
    do na1 = 1, S%nat
      do na2 = 1, S%nat
        d2 = S%tau(:,na1) - S%tau(:,na2)
        d1 = (S%tau(:,na1) + S%tau(:,na2)) / 2._dp - fc%taudef
        call inside_ws(diff_grid_cart, d2, nrws, rws, weights2, inds2)
        call inside_ws(sum_grid_cart, d1, nrws, rws, weights1, inds1)
        do iR2_big = 1, size(weights2)
          do iR1_big = 1, size(weights1)
            R1 = sum_grid_cart(:,inds1(iR1_big)) + diff_grid_cart(:,inds2(iR2_big))/2._dp
            R2 = sum_grid_cart(:,inds1(iR1_big)) - diff_grid_cart(:,inds2(iR2_big))/2._dp
            call cryst_to_cart(1, R1, S%bg, -1)
            call cryst_to_cart(1, R2, S%bg, -1)
            if(any(abs((NINT(R1*2) - R1*2)) > 1e-10_dp)) &
              call errore("center4", "R1_cart is not integer", 1)
            if(any(abs((NINT(R2*2) - R2*2)) > 1e-10_dp)) &
              call errore("center4", "R2_cart is not integer", 1)
            if(any(abs(NINT(R1) - R1) > 1e-10_dp)) cycle
            if(any(abs(NINT(R2) - R2) > 1e-10_dp)) cycle
            iR1 = iR_of(NINT(R1), grid)
            iR2 = iR_of(NINT(R2), grid)
            if (any(abs(nint(R1)) >= grid-1) .or. any(abs(nint(R2)) >= grid-1)) cycle
            call add_ind(R2_list, iR_of(NINT(R2), grid*far_mesh), nxR2, ixR2)
            yR2_list(:,ixR2) = NINT(R2)
            call add_ind(R1_list(:,ixR2), iR_of(NINT(R1), grid*far_mesh), nxR1(ixR2), ixR1)
            yR1_list(:,ixR1,ixR2) = NINT(R1)
            !
            do j1 = 1, 3
              jn1 = 3*(na1-1) + j1
              do j2 = 1, 3
                jn2 = 3*(na2-1) + j2
                new_fc(jn1, jn2, ixR1, ixR2) = &
                  fc%FC(jn1, jn2, iR1, iR2) * &
                  weights1(iR1_big) * weights2(iR2_big)
              enddo
            enddo
            !
          enddo
        enddo
      enddo
    enddo
    !
    deallocate(weights1, weights2)
    deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR2(3,nxR2))
    ALLOCATE(fc%xR2(3,nxR2))
    ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
    ALLOCATE(fc%n_R1(nxR2))
    fc%n_R1 = nxR1(:nxR2)
    fc%n_R2 = nxR2
    fc%nq = grid
    do ixR2 = 1, nxR2
      fc%yR2(:,ixR2) = yR2_list(:,ixR2)
      fc%yR1(:,:nxR1(ixR2),ixR2) = yR1_list(:,:nxR1(ixR2),ixR2)
      fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
    enddo
    call fc%cart(S, 1)
    call fc%cart(S, 2)
    !
    call print_message("end of centering  R1-R2 && R1+R2")
    deallocate(new_fc, yR2_list, yR1_list, nxR1)
  end subroutine
  !
  function iR_of(R, sc_grid) result(iR)
    integer, intent(in) :: R(3)
    integer, intent(in) :: sc_grid(3)
    !
    integer :: i, iR
    integer :: Rplus(3)
    !
    Rplus = R
    do i = 1, 3
      do while (Rplus(i) < 0)
        Rplus(i) = Rplus(i) + sc_grid(i)
      enddo
      do while (Rplus(i) >= sc_grid(i))
        Rplus(i) = Rplus(i) - sc_grid(i)
      enddo
    enddo
    !
    iR = v2index(Rplus, sc_grid)
  end function iR_of
  !
  SUBROUTINE minimal_image(S, sc_grid, diffs, n_weights, weights)
    use thutils, only: e_iqr, grid_vec_cart, cryst2cart
    type(ph_system_info), intent(in) :: S
    integer, intent(in) :: sc_grid(3)
    !
    real(dp) :: big_at(3,3)
    integer :: iR, j, na1, na2
    real(dp), allocatable :: big_R(:,:), R_grid(:,:)
    integer :: far_grid(3)
    real(dp) :: diff(3)
    real(dp), allocatable :: ws_weights(:)
    integer, allocatable :: inds(:)
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    !
    real(dp), allocatable, intent(out) :: diffs(:,:,:,:,:), weights(:,:,:,:)
    integer, allocatable, intent(out):: n_weights(:,:,:)
    !
    allocate(diffs(3,10,product(sc_grid),S%nat,S%nat))
    allocate(n_weights(product(sc_grid),S%nat,S%nat))
    allocate(weights(10,product(sc_grid),S%nat,S%nat))
    !
    far_grid = 2*nfar+1
    ! initialize WS r-vectors
    forall(j = 1:3) big_at(:,j) = S%at(:,j) * sc_grid(j)
    CALL wsinit(rws,nrwsx,nrws,big_at)
    !
    big_R = grid_vec_cart(far_grid, big_at, center=.true.)
    R_grid = grid_vec_cart(sc_grid, S%at)
    do iR = 1, product(sc_grid)
      do na2 = 1, S%nat
        do na1 = 1, S%nat
          diff = R_grid(:,iR) + S%tau(:,na1) - S%tau(:,na2)
          call inside_ws(big_R, diff, nrws, rws, ws_weights, inds)
          n_weights(iR,na1,na2) = size(ws_weights)
          do j = 1, size(ws_weights)
            diffs(:,j,ir,na1,na2) = - R_grid(:,iR) - big_R(:,inds(j))
            ! print*, cryst2cart(diffs(:,j,na1,na2,ir), S%bg, -1)
            weights(j,ir,na1,na2) = ws_weights(j)
          enddo
        enddo
      enddo
    enddo
  end subroutine
  !
  SUBROUTINE minimal_image_2(S, sc_grid, diffs, n_weights, weights)
    use thutils, only: e_iqr, grid_vec_cart, cryst2cart
    type(ph_system_info), intent(in) :: S
    integer, intent(in) :: sc_grid(3)
    !
    real(dp) :: big_at(3,3)
    integer :: iR, j, na1, na2
    real(dp), allocatable :: big_R(:,:), R_grid(:,:)
    integer :: far_grid(3)
    real(dp) :: diff(3)
    real(dp), allocatable :: ws_weights(:)
    integer, allocatable :: inds(:)
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=2000
    INTEGER :: nrws
    REAL(DP) :: rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    integer :: n_ir(product(sc_grid))
    integer :: R_int(3)
    !
    real(dp), allocatable, intent(out) :: diffs(:,:,:,:,:), weights(:,:,:,:)
    integer, allocatable, intent(out):: n_weights(:,:,:)
    !
    allocate(diffs(3,10,product(sc_grid),S%nat,S%nat))
    allocate(n_weights(product(sc_grid), S%nat,S%nat))
    allocate(weights(10,product(sc_grid),S%nat,S%nat))
    !
    far_grid = 2*nfar+1
    ! initialize WS r-vectors
    forall(j = 1:3) big_at(:,j) = S%at(:,j) * sc_grid(j)
    CALL wsinit(rws,nrwsx,nrws,big_at)
    !
    R_grid = grid_vec_cart(far_grid*sc_grid, S%at, center=.true.)
    do na2 = 1, S%nat
      do na1 = 1, S%nat
        diff = S%tau(:,na1) - S%tau(:,na2)
        call inside_ws(R_grid, diff, nrws, rws, ws_weights, inds)
        n_ir = 0
        do j = 1, size(ws_weights)
          R_int = NINT(cryst2cart(R_grid(:,inds(j)), S%bg, -1))
          iR = iR_of(R_int, sc_grid)
          n_iR(iR) = n_iR(iR) + 1
          diffs(:,n_ir(ir),ir,na1,na2) = - R_grid(:,inds(j))
          ! print*, cryst2cart(diffs(:,j,na1,na2,ir), S%bg, -1)
          weights(n_ir(ir),ir,na1,na2) = ws_weights(j)
        enddo
        n_weights(:,na1,na2) = n_ir
      enddo
    enddo
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
    real(dp), dimension(3) :: d_diff, d_sum
    real(dp), parameter :: eps_peri = 1e-3
    integer, parameter :: nperix = (2*nfar+1)**3
    integer :: SAFE_ALLOCATION
    !
    integer :: na1, na2, j1, j2, nR, iR1, iR2, iR_diff, iR_sum
    integer :: map_sc(S%nat, PRODUCT(grid))
    integer, allocatable :: far_grid_cryst(:,:)
    real(dp), allocatable :: far_grid_cart(:,:), grid_cart(:,:)
    integer :: nRbig
    integer :: nxR2, ixR2, ixR1
    integer, allocatable :: yR1_list(:,:,:), yR2_list(:,:), yR_diff(:,:), yR_sum(:,:,:), nxR1(:)
    real(dp), allocatable :: weights_diff(:), weights_sum(:)
    integer, allocatable :: inds_diff(:), inds_sum(:)
    integer, dimension(3) :: far_mesh, R_diff, R_sum, R1, R2
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
    far_grid_cryst = grid_vec_cryst(far_mesh, .true.)
    far_grid_cart = grid_vec_cart(far_mesh, S_sc%at, .true.)
    grid_cart = grid_vec_cart(grid, S_sc%at)
    !
    SAFE_ALLOCATION = 20 * nR
    allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(yR2_list(3, SAFE_ALLOCATION))
    allocate(yR1_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
    allocate(nxR1(SAFE_ALLOCATION))
    !
    map_sc = map_uc2sc(S, S_sc, grid)
    !
    new_fc = 0._dp
    yR1_list = -1
    yR2_list = -1
    nxR1 = 0
    nxR2 = 0
    !

    do na2 = 1, S%nat
      do na1 = 1, S%nat
        do iR2 = 1, nR
          do iR1 = 1, nR
            d_diff = S_sc%tau(:,map_sc(na1,iR1)) - S_sc%tau(:,map_sc(na2,iR2))
            call inside_ws(far_grid_cart, d_diff, nrws, rws, weights_diff, inds_diff)
            do iR_diff = 1, size(weights_diff)
              R_diff = index2v(iR1, grid) - index2v(iR2, grid) + far_grid_cryst(:,inds_diff(iR_diff))
              ! call add_ind(yR1_list, inds_diff(iR_diff), nxR_diff, ixR_diff)
              ! yR_diff(:,ixR_diff) = index2v(R1, grid) + grid * far_grid_cryst(:,inds_diff(iR_diff))
              d_sum = S_sc%tau(:,map_sc(na1,iR1)) + S_sc%tau(:,map_sc(na2,iR2)) - fc%taudef
              call inside_ws(far_grid_cart, d_sum, nrws, rws, weights_sum, inds_sum)
              do iR_sum = 1, size(weights_sum)
                R_sum = index2v(iR1, grid) + index2v(iR2, grid) + grid * far_grid_cryst(:,inds_sum(iR_sum))
                R2 = R_sum - R_diff
                call add_R(yR2_list, R2, nxR2, ixR2)
                R1 = R_sum + R_diff
                call add_R(yR1_list(:,:,ixR2), R1, nxR1(ixR2), ixR1)
                do j1 = 1, 3
                  do j2 = 1, 3
                    new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR1, ixR2) = &
                      fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, iR1, iR2) * &
                      weights_diff(iR_diff) * weights_sum(iR_sum)
                    ! write(10, '(2E15.5)') norm2(fc2_sc%xR1(:,R1,R2) - fc2_sc%xR2(:,R2)), &
                    !   SUM(ABS(fc2_sc%fc((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,R1,R2)))
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
    deallocate(weights_diff, weights_sum)
    deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
    !> the yR list is populated until nxR(which), but the size can be larger.
    !> the values after nxR(which) are not even initialized (they are garbage).
    ALLOCATE(fc%yR2(3,nxR2))
    ALLOCATE(fc%xR2(3,nxR2))
    ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
    ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
    ALLOCATE(fc%n_R1(nxR2))
    fc%n_R1 = nxR1(:nxR2)
    fc%n_R2 = nxR2
    fc%nq = grid
    do ixR2 = 1, nxR2
      fc%yR2(:,ixR2) = yR2_list(:,ixR2)
      fc%yR1(:,:nxR1(ixR2),ixR2) = yR1_list(:,:nxR1(ixR2),ixR2)
      fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
    enddo
    call fc%cart(S_sc, 1)
    call fc%cart(S_sc, 2)

    ! do ixR2 = 1, fc%n_R2
    !   fc%xR2(:,ixR2)= fc%xR2(:,ixR2) / grid
    !   do ixR1 = 1, fc%n_R1(ixR2)
    !     fc%xR1(:,ixR1, ixR2) = fc%xR1(:,ixR1, ixR2) / grid
    !   enddo
    ! enddo
    !
    call print_message("end of centering")
    deallocate(new_fc, yR2_list, yR1_list, nxR1)
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
    integer :: i, nweights, nrws, nR
    real(dp) :: wg, wg_tot
    real(dp), external :: wsweight
    !
    nR = size(R_list, 2) / (nfar*2 + 1)**3
    allocate(weights_(size(R_list, 2)))
    allocate(ind_out_(size(R_list, 2)))
    if (allocated(ind_out)) deallocate(ind_out)
    if (allocated(weights)) deallocate(weights)
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
    !
    wg_tot = wg_tot / REAL(nR, DP)
    !
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
  subroutine add_R(R_list, R, nR, iR)
    integer, intent(inout) :: R_list(:,:)
    integer, intent(in) :: R(3)
    integer, intent(inout) :: nR
    integer, intent(out) :: iR
    !
    integer :: i
    !
    do i = 1, nR
      if(all(R_list(:,i) == R)) then
        iR = i
        return
      endif
    enddo
    !
    nR = nR + 1
    iR = nR
    if(nR > size(R_list, 2)) call errore("add_R", "SAFE ALLOCATION is too short", 1)
    R_list(:,iR) = R
  end subroutine
  ! subroutine center_grid_sc(fc, grid, S, S_sc)
  !   class(forceconst2_grid), intent(inout) :: fc
  !   integer, intent(in) :: grid(3)
  !   type(ph_system_info), intent(in) :: S, S_sc
  !   !
  !   ! type(forceconst2_grid) :: fsc
  !   integer :: nRbig, l1, l2, l3, big_grid(3), new_yR(6), small_Ri(3), small_Rj(3)
  !   real(dp) :: Rbig(3, (2*nfar+1)**3), wg_tot, big_yR(3)
  !   real(dp), allocatable :: new_fc(:,:,:)
  !   real(dp) :: dist(3)
  !   integer :: na1, na2, j1, j2, ixR, nR, jR_big
  !   integer, dimension(S_sc%nat) ::  map_R, map_nat
  !   integer :: R_list(PRODUCT(grid), PRODUCT(grid)*(2*nfar+1)**3)
  !   integer :: index_i, index_j, nxR
  !   integer, allocatable :: new_yR_list(:,:)
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws
  !   REAL(DP) :: wg, rws(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   ! initialize WS r-vectors
  !   CALL wsinit(rws,nrwsx,nrws,S_sc%at)
  !   !
  !   ! Construct a big enough lattice of Supercell **supercell** vectors
  !   fc%centered = .true.
  !   nR = PRODUCT(grid)
  !   big_grid = grid * (2*nfar+1)
  !   allocate(new_yR_list(6, nR**2*100))
  !   allocate(new_fc(S%nat3, S%nat3, nR**2*100))
  !   nRbig=0
  !   DO l1=-nfar, nfar
  !     DO l2=-nfar, nfar
  !       DO l3=-nfar, nfar
  !         nRbig=nRbig+1
  !         Rbig(:, nRbig) = S_sc%at(:,1)*l1 +S_sc%at(:,2)*l2 +S_sc%at(:,3)*l3
  !       END DO
  !     END DO
  !   END DO
  !   IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)
  !   !
  !   map_R = map_sc2uc(S, S_sc, grid, "R")
  !   map_nat = map_sc2uc(S, S_sc, grid, "nat")
  !   !
  !   nxR = 0

  !   R_list = -1
  !   do na1 = 1, S_sc%nat
  !     do na2 = 1, S_sc%nat
  !       small_Ri = index2v(map_R(na1), grid)
  !       ! call cryst_to_cart(1, small_Ri, S%at, 1)
  !       small_Rj = index2v(map_R(na2), grid)
  !       ! call cryst_to_cart(1, small_Rj, S%at, 1)
  !       wg_tot = 0._dp
  !       do jR_big = 1, nRbig
  !         dist = S_sc%tau(:,na1) - Rbig(:,jR_big) - S_sc%tau(:,na2)
  !         wg = wsweight(dist,rws,nrws)
  !         wg_tot = wg_tot + wg
  !         ! wg = wg/DFLOAT(nqt)
  !         if (nfar == 0) wg = 1._dp
  !         if (wg > 1e-6) then
  !           new_yR(1:3) = small_Ri
  !           big_yR = Rbig(:,jR_big)
  !           call cryst_to_cart(1, big_yR, S_sc%bg, -1)
  !           new_yR(4:6) = NINT(big_yR*grid + small_Rj)
  !           index_i = v2index(new_yR(1:3), grid)
  !           index_j = v2index(new_yR(4:6)+grid*nfar, big_grid)
  !           if (R_list(index_i, index_j) == -1) then
  !             nxR = nxR + 1
  !             ixR = nxR
  !             R_list(index_i, index_j) = ixR
  !             new_yR_list(:, ixR) = new_yR
  !           else
  !             ixR = R_list(index_i, index_j)
  !           endif
  !           if(norm2(REAL(new_yR, DP))<1e-6*S%alat) fc%i_0 = ixR
  !           do j1 = 1, 3
  !             do j2 = 1, 3
  !               new_fc(j1+3*(map_nat(na1)-1),j2+3*(map_nat(na2)-1),ixR) = &
  !                 fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, 1) * wg
  !             enddo
  !           enddo
  !         endif
  !       enddo
  !       if (ABS(wg_tot -1)<1e-6) CALL errore("center_sc", "sum of weights is not 1", -1)
  !     enddo
  !   enddo
  !   !
  !   deallocate(fc%yR, fc%xR, fc%FC)
  !   ALLOCATE(fc%yR(6,nxR))
  !   ALLOCATE(fc%xR(6,nxR))
  !   ALLOCATE(fc%FC(S%nat3,S%nat3,nxR))
  !   fc%n_R = nxR
  !   fc%nq = grid
  !   fc%yR = new_yR_list(:,1:nxR)
  !   fc%xR = REAL(new_yR_list(:,1:nxR), DP)
  !   CALL cryst_to_cart(nxR, fc%xR(1:3,:), S%at, 1)
  !   CALL cryst_to_cart(nxR, fc%xR(4:6,:), S%at, 1)
  !   fc%FC = new_fc(:,:,1:nxR)
  ! end subroutine
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
            fc_sc2RR(j1+(map_nat(sc_na1)-1)*3, j2+(map_nat(sc_na2)-1)*3, map_R(sc_na1), map_R(sc_na2)) = &
              FC(j1+(sc_na1-1)*3, j2+(sc_na2-1)*3)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function fc_sc2RR_cmplx(grid, S, S_sc, fc)
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    complex(dp), intent(in) :: FC(S_sc%nat3,S_sc%nat3)
    !
    complex(dp) :: fc_sc2RR_cmplx(S%nat3, S%nat3, product(grid), product(grid))
    integer :: sc_na1, sc_na2, j1, j2
    integer, dimension(S_sc%nat) :: map_R, map_nat
    !
    map_R = map_sc2uc(S, S_sc, grid, "R")
    map_nat = map_sc2uc(S, S_sc, grid, "nat")
    fc_sc2RR_cmplx = 0._dp
    do sc_na1 = 1, S_sc%nat
      do sc_na2 = 1, S_sc%nat
        do j1 = 1, 3
          do j2 = 1, 3
            fc_sc2RR_cmplx(j1+(map_nat(sc_na1)-1)*3, j2+(map_nat(sc_na2)-1)*3, map_R(sc_na1), map_R(sc_na2)) = &
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
    fc_uc2sc = 0._dp !> needed for inclusion defect
    nR = product(sc_grid)
    atoms_sc = map_uc2sc(S, S_sc, sc_grid)
    do R1 = 1, nR
      do R2 = 1, nR
        R = index2v(R1, sc_grid) - index2v(R2, sc_grid)
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
  function fc_uc2RR(fc)
    type(forceconst2_grid), intent(in) :: fc
    integer :: R1, R2, nR, j1, R(3)
    real(dp), allocatable :: fc_uc2RR(:,:,:,:)
    integer :: nat3
    !
    nR = product(fc%nq)
    nat3 = size(fc%fc, 1)
    allocate(fc_uc2RR(nat3, nat3, nR, nR))
    do R2 = 1, nR
      do R1 = 1, nR
        R = index2v(R1, fc%nq) - index2v(R2, fc%nq)
        do j1 = 1, 3
          if (R(j1) < 0) R(j1) = R(j1) + fc%nq(j1)
        enddo
        fc_uc2RR(:,:,R1,R2) = fc%fc(:,:,v2index(R, fc%nq))
      enddo
    enddo
  end function
  !
!   subroutine q2r_1st_step(fc, iR1, grid, V)
!     use q_grids, only : q_grid
!     use constants, only : tpi
!     CLASS(forceconst2_sc),INTENT(inout) :: fc
!     INTEGER, INTENT(in) :: iR1
!     complex(DP),INTENT(in) :: V(:,:,:,:)
!     type(q_grid), intent(in) :: grid
!     !
!     INTEGER :: qj, Nq, qi, nat3
!     REAL(DP), allocatable :: varg(:), vcos(:), vsin(:)
!     COMPLEX(DP), allocatable :: vphase(:)
!     !
!     Nq = product(grid%n)
!     nat3 = size(fc%FC, 1)
!     if(allocated(fc%mix)) deallocate(fc%mix)
!     fc%iR1 = iR1
!     !
!     allocate(fc%mix(nat3, nat3, Nq))
!     fc%mix = 0._dp
!     !
!     allocate(varg(Nq), vcos(Nq), &
!       vsin(Nq), vphase(Nq))
!     do qi = 1, Nq
!       FORALL(qj=1:Nq) varg(qj) =  tpi * &
!         dot_product(fc%xR1(:,iR1), grid%xq(:,qj))
!       !
!       ! Pre-compute phase to use the vectorized MKL subroutines
! #if defined(__INTEL) && defined(__HASVTRIG)
! !dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
!       CALL vdCos(Nq, varg, vcos)
!       CALL vdSin(Nq, varg, vsin)
! #else
!       vcos = DCOS(varg)
!       vsin = DSIN(varg)
! #endif
!       vphase =  CMPLX( vcos, -vsin, kind=DP  )
!       !
!       do qj = 1, Nq
!         fc%mix(:,:,qi) = fc%mix(:,:,qi) + vphase(qj) * V(:,:,qi,qj)
!       enddo
!       ! fc%mix(:,:,R1) = SUM(vphase * fc%FC(:,:,:,R1), dim=3)
!       !
!     enddo
!     deallocate(varg, vcos, vsin, vphase)
!   end subroutine
!   !
!   subroutine q2r_2nd_step(fc, iR2, grid)
!     use q_grids, only : q_grid
!     use constants, only : tpi
!     CLASS(forceconst2_sc),INTENT(inout) :: fc
!     integer, INTENT(in) :: iR2
!     type(q_grid), intent(in) :: grid
!     !
!     INTEGER :: qi, Nq
!     REAL(DP), allocatable, dimension(:) :: varg, vcos, vsin
!     COMPLEX(DP), allocatable :: vphase(:)
!     !
!     Nq = product(grid%n)
!     allocate(varg(Nq), vcos(Nq), vsin(Nq), vphase(Nq))
!     !
!     FORALL(qi=1:Nq) varg(qi) =  tpi * &
!       dot_product(fc%xR2(:,iR2, fc%iR1), grid%xq(:,qi))
!     !
!     ! Pre-compute phase to use the vectorized MKL subroutines
! #if defined(__INTEL) && defined(__HASVTRIG)
! !dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
!     CALL vdCos(Nq, varg, vcos)
!     CALL vdSin(Nq, varg, vsin)
! #else
!     vcos = DCOS(varg)
!     vsin = DSIN(varg)
! #endif
!     vphase =  CMPLX( vcos, vsin, kind=DP  )
!     !
!     fc%fc(:,:,:,fc%iR1) = 0._dp
!     do qi = 1, Nq
!       fc%fc(:,:,iR2,fc%iR1) = fc%fc(:,:,iR2,fc%iR1) + vphase(qi) * fc%mix(:,:,qi)
!     enddo
!     !
!   END SUBROUTINE

  !
  SUBROUTINE r2q_1st_step(fc, xq)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    !
    INTEGER :: R2, R1, nat3
    REAL(DP), allocatable :: varg(:), vcos(:), vsin(:)
    COMPLEX(DP), allocatable :: vphase(:)
    !
    nat3 = size(fc%FC, 1)
    if(allocated(fc%mix)) deallocate(fc%mix)
    !
    allocate(fc%mix(nat3, nat3, fc%n_R2))
    fc%mix = 0._dp
    !
    do R2 = 1, fc%n_R2
      allocate(varg(fc%n_R1(R2)), vcos(fc%n_R1(R2)), &
        vsin(fc%n_R1(R2)), vphase(fc%n_R1(R2)))
      FORALL(R1=1:fc%n_R1(R2)) varg(R1) =  tpi * &
        dot_product(xq, fc%xR1(:,R1,R2))
      !
      ! Pre-compute phase to use the vectorized MKL subroutines
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R1(R2), varg, vcos)
      CALL vdSin(fc%n_R1(R2), varg, vsin)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      !
      do R1 = 1, fc%n_R1(R2)
        fc%mix(:,:,R2) = fc%mix(:,:,R2) + vphase(R1) * fc%FC(:,:,R1,R2)
      enddo
      ! fc%mix(:,:,R1) = SUM(vphase * fc%FC(:,:,:,R1), dim=3)
      !
      deallocate(varg, vcos, vsin, vphase)
    enddo
    !
  END SUBROUTINE
  !
  SUBROUTINE r2q_2nd_step(fc, xq, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    complex(dp), intent(out) :: D(size(fc%fc, 1), size(fc%fc, 1))
    !
    INTEGER :: i
    REAL(DP), dimension(fc%n_R2) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(fc%n_R2)
    !
    !
    FORALL(i=1:fc%n_R2) varg(i) =  tpi * &
      dot_product(xq, fc%xR2(:,i))
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
    CALL vdCos(fc%n_R2, varg, vcos)
    CALL vdSin(fc%n_R2, varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, vsin, kind=DP  )
    !
    D = 0._dp
    do i = 1, fc%n_R2
      D = D + vphase(i) * fc%mix(:,:,i)
    enddo
    !
    ! D = D / PRODUCT(fc%nq)
    !
  END SUBROUTINE
  !
  SUBROUTINE r2q_at_once(fc, xq1, xq2, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc), INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq1(3), xq2(3)
    complex(dp), intent(out) :: D(size(fc%fc, 1), size(fc%fc, 1))
    !
    INTEGER :: ir1, ir2
    REAL(DP), dimension(size(fc%fc,3),fc%n_R2) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(size(fc%fc,3),fc%n_R2)
    !
    do ir2 = 1, fc%n_R2
      do ir1 = 1, fc%n_R1(ir2)
        varg(ir1,ir2) = tpi * (-dot_product(xq1, fc%xR1(:,ir1,ir2)) + &
          dot_product(xq2, fc%xR2(:,ir2)))
      enddo
    enddo
    !
    vcos = DCOS(varg)
    vsin = DSIN(varg)
    vphase =  CMPLX( vcos, vsin, kind=DP  )
    !
    D = 0._dp
    do ir2 = 1, fc%n_R2
      do ir1 = 1, fc%n_R1(ir2)
        D = D + vphase(ir1,ir2) * fc%fc(:,:,ir1,ir2)
      enddo
    enddo
    !
    ! D = D / PRODUCT(fc%nq)
    !
  END SUBROUTINE
  !
  function flatten_RR_cmplx(mat) result(mat_flat)
    complex(dp), intent(in) :: mat(:,:,:,:)
    integer :: nat3i, nat3j, nqi, nqj
    !
    complex(dp), allocatable :: mat_flat(:,:)
    integer :: i, j, na1, na2, n1, n2
    !
    nat3i = size(mat,1)
    nat3j = size(mat,2)
    nqi = size(mat,3)
    nqj = size(mat,4)

    allocate(mat_flat(nat3i*nqi, nat3j*nqj))
    do n1 = 1, nqi
      do n2 = 1, nqj
        do na2 = 1, nat3j
          j = na2 + (n2-1)*nat3j
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat_flat(i,j) = mat(na1,na2,n1,n2)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function flatten_RR_real(mat) result(mat_flat)
    real(dp), intent(in) :: mat(:,:,:,:)
    integer :: nat3i, nat3j, nqi, nqj
    !
    real(dp), allocatable :: mat_flat(:,:)
    integer :: i, j, na1, na2, n1, n2
    !
    nat3i = size(mat,1)
    nat3j = size(mat,2)
    nqi = size(mat,3)
    nqj = size(mat,4)

    allocate(mat_flat(nat3i*nqi, nat3j*nqj))
    do n1 = 1, nqi
      do n2 = 1, nqj
        do na2 = 1, nat3j
          j = na2 + (n2-1)*nat3j
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat_flat(i,j) = mat(na1,na2,n1,n2)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function unflatten_RR_cmplx(mat_flat, nqi, nqj) result(mat)
    complex(dp), intent(in) :: mat_flat(:,:)
    integer,  intent(in) :: nqi, nqj
    integer :: nat3i, nat3j
    complex(dp), allocatable :: mat(:,:,:,:)
    integer :: n1, n2, na1, na2, i, j

    nat3i = size(mat_flat,1) / nqi
    nat3j = size(mat_flat,2) / nqj
    allocate(mat(nat3i, nat3j, nqi, nqj))

    do n2 = 1, nqj
      do na2 = 1, nat3j
        j = na2 + (n2-1)*nat3j
        do n1 = 1, nqi
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat(na1,na2,n1,n2) = mat_flat(i,j)
          end do
        end do
      end do
    end do
  end function

  function unflatten_RR_real(mat_flat, nqi, nqj) result(mat)
    real(dp), intent(in) :: mat_flat(:,:)
    integer,  intent(in) :: nqi, nqj
    integer :: nat3i, nat3j
    real(dp), allocatable :: mat(:,:,:,:)
    integer :: n1, n2, na1, na2, i, j

    nat3i = size(mat_flat,1) / nqi
    nat3j = size(mat_flat,2) / nqj
    allocate(mat(nat3i, nat3j, nqi, nqj))

    do n2 = 1, nqj
      do na2 = 1, nat3j
        j = na2 + (n2-1)*nat3j
        do n1 = 1, nqi
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat(na1,na2,n1,n2) = mat_flat(i,j)
          end do
        end do
      end do
    end do
  end function
  !
  SUBROUTINE div_mass0_fcsc(S, fc)
    USE kinds, only : DP
    use thutils, only: bz2simple
    IMPLICIT NONE
    TYPE(forceconst2_sc) :: fc
    TYPE(ph_system_info)   :: S
    !
    INTEGER :: i, j, iR1, iR2
    !
    DO iR2 = 1, fc%n_R2
      do iR1 = 1, fc%n_R1(iR2)
        do concurrent(i=1:S%nat3, j=1:S%nat3)
          fc%FC(i, j, iR1, iR2) = fc%FC(i, j, iR1, iR2) * S%sqrtmm1(i) * S%sqrtmm1(j)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE
  !
  ! SUBROUTINE div_mass_fcsc(S, Sd,fc)
  !   USE kinds, only : DP
  !   use thutils, only: bz2simple
  !   IMPLICIT NONE
  !   TYPE(forceconst2_sc) :: fc
  !   TYPE(ph_system_info)   :: S, Sd
  !   !
  !   INTEGER :: i, j, iR1, iR2, map(S%nat,product(fc%nq)), si, sj
  !   !
  !   IF(.not.ALLOCATED(Sd%sqrtmm1)) &
  !     call errore('div_mass_fc2', 'missing sqrtmm1, call aux_system first', 1)

  !   map = map_uc2sc(S, Sd, fc%nq)
  !   DO iR1 = 1, fc%n_R1
  !     si = v2index(bz2simple(fc%yR1(:,iR1), fc%nq), fc%nq)
  !     do iR2 = 1, fc%n_R2(iR1)
  !       sj = v2index(bz2simple(fc%yR2(:,iR2,iR1), fc%nq), fc%nq)
  !       DO j = 1, S%nat3
  !         DO i = 1, S%nat3
  !           fc%FC(i, j, iR2, iR1) = fc%FC(i, j, iR2, iR1) * Sd%sqrtmm1(map((i-1)/3+1,si))*Sd%sqrtmm1(map((j-1)/3+1,sj))
  !         ENDDO
  !       ENDDO
  !     ENDDO
  !   ENDDO
  !   !
  ! END SUBROUTINE
  !
  ! subroutine matd2RR(fc3, S, fcsc)
  !   use fc3_interpolate, only : sparse
  !   type(sparse), intent(in) :: fc3
  !   type(ph_system_info), intent(in) :: S
  !   integer :: mesh_cube, len2_new, i2_new, i3_new, i, j
  !   integer :: minimum, maximum, mesh(3)
  !   integer, dimension(size(fc3%yR2, 2)) :: iR2, iR3
  !   integer, allocatable :: ind2(:), ind3(:,:), len3_new(:)
  !   type(forceconst2_sc), intent(out) :: fcsc
  !   !
  !   minimum = min(minval(fc3%yR2), minval(fc3%yR3))
  !   maximum = max(maxval(fc3%yR2), maxval(fc3%yR3))
  !   mesh = maximum - minimum + 1
  !   mesh_cube = product(mesh)
  !   !
  !   do i = 1, fc3%n_R
  !     iR2(i) = v2index(fc3%yR2(:,i) - minimum, mesh)
  !     iR3(i) = v2index(fc3%yR3(:,i) - minimum, mesh)
  !   end do
  !   !
  !   allocate(ind2(mesh_cube), ind3(mesh_cube, mesh_cube))
  !   allocate(len3_new(mesh_cube))
  !   !
  !   ind2 = -1
  !   ind3 = -1
  !   len2_new = 0
  !   len3_new = 0
  !   do i = 1, fc3%n_R
  !     call add_ind(ind2, iR2(i), len2_new, i2_new)
  !     call add_ind(ind3(:,i2_new), iR3(i), len3_new(i2_new), i3_new)
  !   enddo
  !   !
  !   allocate(fcsc%n_R2(len2_new))
  !   fcsc%n_R2 = len3_new(:len2_new)
  !   fcsc%n_R1 = len2_new
  !   !
  !   allocate(fcsc%yR1(3,len2_new))
  !   allocate(fcsc%xR1(3,len2_new))
  !   allocate(fcsc%yR2(3,maxval(len3_new),len2_new))
  !   allocate(fcsc%xR2(3,maxval(len3_new),len2_new))
  !   allocate(fcsc%fc(S%nat3,S%nat3,maxval(len3_new),len2_new))
  !   !
  !   ind2 = -1
  !   ind3 = -1
  !   len2_new = 0
  !   len3_new = 0
  !   do i = 1, fc3%n_R
  !     call add_ind(ind2, iR2(i), len2_new, i2_new)
  !     fcsc%yR1(:,i2_new) = fc3%yR2(:,i)
  !     call add_ind(ind3(:,i2_new), iR3(i), len3_new(i2_new), i3_new)
  !     fcsc%yR2(:,i3_new,i2_new) = fc3%yR3(:,i)
  !     do j = 1, fc3%n_terms(i)
  !       fcsc%fc(fc3%dat(i)%idx(2,j),fc3%dat(i)%idx(3,j),i3_new,i2_new) = &
  !         fc3%dat(i)%fc(j)
  !     enddo
  !   enddo
  !   !
  !   call fcsc%cart(S, 1)
  !   call fcsc%cart(S, 2)
  !   fcsc%nq = fc3%nq
  !   !
  !   deallocate(ind2, ind3, len3_new)
  ! end subroutine
  !
  subroutine diag_degen_cmplx(nat3, dim_deg, V1, Udeg, e1)
    use fc2_interpolate, only : mat2_diag
    !
    integer, intent(in) :: nat3, dim_deg
    complex(dp), intent(in) :: V1(nat3,nat3)
    complex(dp), intent(inout) :: Udeg(nat3,dim_deg)
    complex(dp), intent(out), optional :: e1(dim_deg)
    !
    complex(dp) :: Vdeg(dim_deg,dim_deg)
    !
    Vdeg = matmul(TRANSPOSE(CONJG(Udeg)), matmul(V1, Udeg))
    call mat2_diag(dim_deg, Vdeg, e1) !freq(i:i+dim_deg-1))
    Udeg = matmul(Udeg, Vdeg)
  end subroutine
  !
  subroutine diag_degen_real(nat3, dim_deg, V1, Udeg, e1)
    use fc2_interpolate, only : mat2_diag
    !
    integer, intent(in) :: nat3, dim_deg
    complex(dp), intent(in) :: V1(nat3,nat3)
    complex(dp), intent(inout) :: Udeg(nat3,dim_deg)
    real(dp), intent(out), optional :: e1(dim_deg)
    !
    complex(dp) :: Vdeg(dim_deg,dim_deg)
    !
    Vdeg = matmul(TRANSPOSE(CONJG(Udeg)), matmul(V1, Udeg))
    call mat2_diag(dim_deg, Vdeg, e1) !freq(i:i+dim_deg-1))
    Udeg = matmul(Udeg, Vdeg)
  end subroutine
  !
  subroutine lift_degen(nat3, freq, V1, U, e1)
    use thutils, only: near, id_mat, outer_product, near, braket
    use functions, only: quicksort_idx
    use fc2_interpolate, only : mat2_diag
    !
    integer, intent(in) :: nat3
    real(dp), intent(in) :: freq(nat3)
    complex(dp), intent(in) :: V1(nat3,nat3)
    complex(dp), intent(inout) :: U(nat3,nat3)
    real(dp), intent(out), optional :: e1(nat3)
    !
    REAL(DP),PARAMETER :: epsq = 1.e-8_dp
    integer :: i,j, dim_deg
    real(dp) :: e1_(nat3)
    !
    i = 1
    do ! i cycle
      do j = i+1, nat3
        if (.not. near(freq(i), freq(j))) exit
      enddo
      dim_deg = j-i
      if(dim_deg > 1) then
        call diag_degen_real(nat3, dim_deg, V1, U(:,i:i+dim_deg-1), e1_(i:i+dim_deg-1))
      else
        e1_(i) = REAL(braket(U(:,i), V1), DP)
      endif
      i = i + dim_deg
      if(i > nat3) exit
    enddo
    if (present(e1)) e1 = e1_
  end subroutine
  !
  SUBROUTINE freq_phq_degen(xq, S, fc2, freq, V1, U, e1)
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
    real(dp), intent(out), optional   :: e1(S%nat3)
    REAL(DP),PARAMETER :: epsq = 1.e-8_dp
    REAL(DP) :: cq(3), chk(3)
    LOGICAL :: gamma
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
    if (present(e1)) then
      call lift_degen(S%nat3, freq, V1, U, e1)
    else
      call lift_degen(S%nat3, freq, V1, U)
    endif

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
  subroutine freq_in_grid_degen(S, fc2, fc2sc, grid, freqs, Us, freqs1)
    use q_grids, only : q_grid
    use mpi_thermal, only : mpi_bsum
    use merge_degenerate, only: merge_degen
    !
    type(ph_system_info), intent(in) :: S
    type(q_grid), intent(in) :: grid
    type(forceconst2_grid), intent(in) :: fc2
    type(forceconst2_sc) :: fc2sc
    !
    complex(dp) :: V1(S%nat3,S%nat3)
    real(dp), intent(out):: freqs(S%nat3, grid%nqtot)
    complex(dp), intent(out) :: Us(S%nat3, S%nat3, grid%nqtot)
    real(dp), intent(out)    :: freqs1(S%nat3, grid%nqtot)
    !
    integer :: iq, iqp
    !
    freqs = 0.0_dp
    DO iq = 1, grid%nq
      iqp = iq + grid%iq0
      call fc2sc%r2q(grid%xq(:,iq))
      call fc2sc%r2q(grid%xq(:,iq), V1)
      CALL freq_phq_degen(grid%xq(:,iq), S, fc2, freqs(:,iqp), V1, Us(:,:,iqp), freqs1(:,iqp))
      ! call merge_degen(S%nat3, freqs(:,iqp), freqs(:,iqp))
    END DO
    if (grid%scattered) CALL mpi_bsum(S%nat3, grid%nqtot, freqs)
  end subroutine
end module
