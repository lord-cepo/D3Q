module quter_defect
  use kinds,              only : dp
  use fc2_interpolate,    only : forceconst2_grid
  use ph_system,          only : ph_system_info
  use thutils,            only : v2index, index2v, index2v_cart, print_message
  implicit none
  !
  integer, parameter :: nfar = 2
  !
  type :: list1
    integer, pointer :: list(:,:)
  end type
  !
  type forceconst2_sc
    INTEGER :: n_R2 = 0
    character(len=16) :: def_type
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
    integer, allocatable :: R_map(:)
    !!
    integer, allocatable :: at_map(:)
    !!
    integer, allocatable :: sc_map(:,:)
    !!
    real(dp) :: eps = 0._dp
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
    use thutils, only : cryst2cart
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    ! real(dp), intent(out), optional :: taudef(3)
    integer :: ndef
    !
    logical :: idef(S_sc%nat)
    integer, allocatable :: map_uc2sc(:,:)
    integer :: i, isc, iR, n_add
    real(dp) :: r_cryst(3)
    !
    allocate(map_uc2sc(S%nat, PRODUCT(sc_grid)))
    map_uc2sc = -1
    idef = .false.
    ndef = 0
    ! if(present(taudef)) taudef = 0._dp
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        r_cryst = cryst2cart(S_sc%tau(:,isc), S_sc%bg, -1)*sc_grid - &
          cryst2cart(S%tau(:,i), S%bg, -1)
        ! call cryst_to_cart(1, r_cryst, S_sc%bg, -1)
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
    n_add = S_sc%nat - product(sc_grid)*S%nat
    if (n_add < 0 .and. 0 /= n_add + count(map_uc2sc == -1)) then
      ! print*, map_uc2sc
      print*, n_add, count(map_uc2sc == -1)
      CALL errore("map_uc2sc", "some atoms are not mapped", 1)
    endif

  end function
!
  SUBROUTINE get_equiv_sites( nat, na0, equiv_list )
    USE kinds, ONLY : DP
    use symm_base, only : nsym, irt
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: nat, na0
    INTEGER, allocatable, INTENT(OUT) :: equiv_list(:)

    INTEGER :: isym, na, neq
    LOGICAL :: used(nat)
    integer :: temp_list(nat)

    used(:)       = .FALSE.
    neq           = 0
    if (na0 < 1 .or. na0 > nat) &
      call errore("get_equiv_sites", "defect atom index is out of range", max(1, abs(na0)))
    if (.not. allocated(irt)) &
      call errore("get_equiv_sites", "symmetry atom map is not initialized", 1)
    if (size(irt, 1) < nsym .or. size(irt, 2) < nat) &
      call errore("get_equiv_sites", "symmetry atom map has incompatible dimensions", 1)
    !
    DO isym = 1, nsym
      na = irt(isym, na0)
      if (na < 1 .or. na > nat) &
        call errore("get_equiv_sites", "invalid atom index in symmetry map", max(1, abs(na)))
      IF (.NOT. used(na) .and. na /= na0) THEN
        neq = neq + 1
        temp_list(neq) = na
        used(na) = .TRUE.
      END IF
    END DO
    !
    allocate(equiv_list(neq))
    if (neq > 0) equiv_list = temp_list(:neq)
    !
  END SUBROUTINE
!
  subroutine find_defects(S, S_sc, fc)
    use q_grids, only : symmetrize_system
    TYPE(ph_system_info), intent(in) :: S, S_sc
    CLASS(forceconst2_sc), intent(inout) :: fc
    !
    integer, allocatable :: sites(:)
    integer :: map(S%nat, product(fc%nq))
    integer, allocatable :: map_sc(:)
    integer :: i, iR, isc, ndef, ndef_pos
    integer, allocatable :: defects_tmp(:,:)
    real(dp), allocatable :: mass_ratios_tmp(:)
    !
    map = map_uc2sc(S, S_sc, fc%nq)
    !
    allocate(fc%defects(3, S%nat * product(fc%nq)))
    allocate(fc%mass_ratios(S%nat * product(fc%nq)))
    !
    ndef = 0
    do iR = 1, product(fc%nq)
      do i = 1, S%nat
        if(S_sc%atm(S_sc%ityp(map(i,iR))) /= S%atm(S%ityp(i))) then
          ndef = ndef + 1
          fc%defects(:, ndef) = [i, iR, map(i,iR)]
          ! fc%mass_ratios(ndef) = - (S_sc%amass(S_sc%ityp(map(i,iR))) - S%amass(S%ityp(i))) / &
          !   S%amass(S%ityp(i))
        endif
      enddo
    enddo
    ! do isc = 1, S_sc%nat
    !   if(map_sc(isc) /= -1) cycle
    !   ndef = ndef + 1
    !   fc%defects(:, ndef) = [0, 0, isc]
    !   ! fc%mass_ratios(ndef) = 0._dp
    ! enddo
    !
    !
    if(ndef == 1 .and. fc%defects(1,1) > 0) then
      fc%taudef = S%tau(:, fc%defects(1,1)) + index2v_cart(fc%defects(2,1), fc%nq, S%at)
      call symmetrize_system(S)
      call get_equiv_sites(S%nat, fc%defects(1,1), sites)
      do i = 1, size(sites)
        ndef = ndef + 1
        fc%defects(:, ndef) = [sites(i), 0, 0]
        fc%mass_ratios(ndef) = fc%mass_ratios(1)
      enddo
    else
      call errore("find_defects", "more than one defect atom is not supported yet", 1)
      ! fc%taudef = 0._dp
      ! ndef_pos = 0
      ! do i = 1, ndef
      !   if(fc%defects(1,i) > 0 .and. fc%defects(2,i) > 0) then
      !     fc%taudef = fc%taudef + S%tau(:, fc%defects(1,i)) + index2v_cart(fc%defects(2,i), fc%nq, S%at)
      !     ndef_pos = ndef_pos + 1
      !   elseif(fc%defects(3,i) > 0) then
      !     fc%taudef = fc%taudef + S_sc%tau(:, fc%defects(3,i))
      !     ndef_pos = ndef_pos + 1
      !   endif
      ! enddo
      ! if(ndef_pos > 0) fc%taudef = fc%taudef / REAL(ndef_pos, DP)
    endif
    !
    allocate(defects_tmp(3, ndef))
    allocate(mass_ratios_tmp(ndef))
    defects_tmp = fc%defects(:, :ndef)
    mass_ratios_tmp = fc%mass_ratios(:ndef)
    call move_alloc(defects_tmp, fc%defects)
    call move_alloc(mass_ratios_tmp, fc%mass_ratios)
    !
  end subroutine
  !
  function map_sc2uc(S, S_sc, sc_grid, which)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    character(*) :: which
    !
    integer, allocatable :: map_sc2uc(:)
    integer :: i, isc, iR, el
    real(dp) :: r_cryst(3), tau(3), tau_sc(3)
    !
    allocate(map_sc2uc(S_sc%nat))
    map_sc2uc = -1
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        tau_sc = S_sc%tau(:,isc)
        call cryst_to_cart(1, tau_sc, S_sc%bg, -1)
        tau = S%tau(:,i)
        call cryst_to_cart(1, tau, S%bg, -1)
        r_cryst = tau_sc*sc_grid - tau
        if (NORM2(r_cryst - NINT(r_cryst))<1e-1) then
          if (map_sc2uc(isc) /= -1) then
            print*, "i, isc, r_cryst", i, isc, map_sc2uc(isc)
            call errore("map_sc2uc", "We already have found this defect", isc)
          endif
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
      print*, map_sc2uc
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
  SUBROUTINE allocate_fc2_sc(fc, S, S_sc, grid) !, sites)
    use thutils, only : grid_vec
    IMPLICIT NONE
    INTEGER,INTENT(in) :: grid(3)
    type(ph_system_info),INTENT(in) :: S
    type(ph_system_info),INTENT(in) :: S_sc
    ! integer, allocatable, intent(in) :: sites(:)
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
    fc%R_map = map_sc2uc(S, S_sc, grid, "R")
    fc%at_map = map_sc2uc(S, S_sc, grid, "nat")
    fc%sc_map = map_uc2sc(S, S_sc, grid)
    call find_defects(S, S_sc, fc)
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
    if (allocated(fc%defects)) deallocate(fc%defects)
    if (allocated(fc%mass_ratios)) deallocate(fc%mass_ratios)
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
    integer :: nR, iR
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
  subroutine center2(fc, grid, S, S_sc)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    !
    ! type(forceconst2_grid) :: fsc
    real(dp), dimension(3) :: d2, d1
    real(dp) :: perix, peri_min
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
  ! subroutine center4(fc, grid, S)
  !   use functions, only : refold_bz
  !   use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
  !   class(forceconst2_sc), intent(inout) :: fc
  !   integer, intent(in) :: grid(3)
  !   type(ph_system_info), intent(in) :: S
  !   !
  !   ! type(forceconst2_grid) :: fsc
  !   real(dp), dimension(3) :: d1, d2
  !   real(dp), parameter :: eps_peri = 1e-3
  !   integer, parameter :: nperix = (2*nfar+1)**3
  !   integer :: SAFE_ALLOCATION
  !   !
  !   integer :: na1, na2, j1, j2, nR, iR1, iR2, iR1_big, iR2_big, jn1, jn2
  !   real(dp), allocatable :: far_grid_cart(:,:)
  !   integer, allocatable :: far_grid_cryst(:,:)
  !   integer :: nRbig
  !   integer :: nxR2, ixR2, ixR1
  !   integer, allocatable :: R1_list(:,:), R2_list(:), nxR1(:)
  !   real(dp), allocatable :: xR2_list(:,:), xR1_list(:,:,:)
  !   real(dp), allocatable :: weights1(:), weights2(:)
  !   integer, allocatable :: inds1(:), inds2(:)
  !   integer, dimension(3) :: far_mesh
  !   real(dp) :: big_at(3,3), max_norm
  !   !
  !   real(dp), allocatable :: new_fc(:,:,:,:)
  !   character(100) :: filename
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws
  !   REAL(DP) :: rws(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   ! initialize WS r-vectors
  !   forall(j1 = 1:3) big_at(:,j1) = S%at(:,j1) * grid(j1)
  !   CALL wsinit(rws,nrwsx,nrws,big_at)
  !   !
  !   fc%stage = 0
  !   nR = PRODUCT(grid)
  !   if (nfar == 0) return
  !   far_mesh = 2*nfar+1
  !   nRbig = PRODUCT(far_mesh)
  !   !
  !   allocate(far_grid_cart(3, nRbig*nR))
  !   allocate(far_grid_cryst(3, nrbig*nr))
  !   far_grid_cart = grid_vec_cart(far_mesh*grid, S%at, center=.true.)
  !   far_grid_cryst = grid_vec_cryst(far_mesh*grid, center=.true.)
  !   !
  !   SAFE_ALLOCATION = 20 * nR
  !   allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
  !   allocate(R2_list(nR*nRbig))
  !   allocate(R1_list(nR*nRbig, nR*nRbig))
  !   allocate(xR1_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
  !   allocate(xR2_list(3, SAFE_ALLOCATION))
  !   allocate(nxR1(SAFE_ALLOCATION))
  !   !
  !   new_fc = 0._dp
  !   R1_list = -1
  !   R2_list = -1
  !   nxR1 = 0
  !   nxR2 = 0
  !   max_norm = 0._dp
  !   !
  !   write(filename, "(A,I1,A)") 'distance-center4-', grid(1), '.dat'
  !   open(10, file=trim(filename), status='replace')
  !   do na2 = 1, S%nat
  !     d2 = S%tau(:,na2) - fc%taudef
  !     call inside_ws(far_grid_cart, d2, nrws, rws, weights2, inds2)
  !     do iR2_big = 1, size(weights2)
  !       iR2 = inside_periodic(inds2(iR2_big), grid)
  !       call add_ind(R2_list, inds2(iR2_big), nxR2, ixR2)
  !       xR2_list(:,ixR2) = far_grid_cart(:,inds2(iR2_big))
  !       do na1 = 1, S%nat
  !         d1 = S%tau(:,na1) - fc%taudef
  !         call inside_ws(far_grid_cart, d1, nrws, rws, weights1, inds1)
  !         do iR1_big = 1, size(weights1)
  !           iR1 = inside_periodic(inds1(iR1_big), grid)
  !           call add_ind(R1_list(:,ixR2), inds1(iR1_big), nxR1(ixR2), ixR1)
  !           xR1_list(:,ixR1,ixR2) = far_grid_cart(:,inds1(iR1_big))
  !           ! if (iR1 == iR2 .and. any(xR2_list(:,ixR2) /= xR1_list(:,ixR1,ixR2))) then
  !           !   call cryst_to_cart(1, xR1_list(:,ixR1,ixR2), S%bg, -1)
  !           !   call cryst_to_cart(1, xR2_list(:,ixR2), S%bg, -1)
  !           !   print*, na1, na2, iR1, iR2, xR1_list(:,ixR1,ixR2), xR2_list(:,ixR2)
  !           !   print*, inds1(iR1_big)
  !           !   print*, inside_periodic(inds1(iR1_big), grid)
  !           !   print*, far_grid_cart(:,inds1(iR1_big))
  !           !   call errore("center4", "iR1 == iR2 but xR1_list /= xR2_list", 1)
  !           ! endif
  !           ! max_norm = max(max_norm, norm2(d1 + xR1_list(:,ixR1,ixR2) - d2 - xR2_list(:,ixr2)))
  !           ! if(norm2(d1 + xR1_list(:,ixR1,ixR2) - d2 - xR2_list(:,ixr2)) < 2._dp) then
  !           do j1 = 1, 3
  !             jn1 = 3*(na1-1) + j1
  !             do j2 = 1, 3
  !               jn2 = 3*(na2-1) + j2
  !               new_fc(jn1, jn2, ixR1, ixR2) = &
  !                 fc%FC(jn1, jn2, iR1, iR2) * &
  !                 weights1(iR1_big) * weights2(iR2_big)
  !               write(10, "(3E20.8)") &
  !                 norm2(d1 + xR1_list(:,ixR1,ixR2)), &
  !                 norm2(d2 + xR2_list(:,ixr2)), &
  !                 new_fc(jn1, jn2, ixR1, ixR2)
  !             enddo
  !           enddo
  !           ! endif
  !           !
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   print*, "max_norm:", max_norm
  !   close(10)
  !   !
  !   deallocate(weights1, weights2)
  !   deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
  !   !> the yR list is populated until nxR(which), but the size can be larger.
  !   !> the values after nxR(which) are not even initialized (they are garbage).
  !   ALLOCATE(fc%yR2(3,nxR2))
  !   ALLOCATE(fc%xR2(3,nxR2))
  !   ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%n_R1(nxR2))
  !   fc%n_R1 = nxR1(:nxR2)
  !   fc%n_R2 = nxR2
  !   fc%nq = grid
  !   do ixR2 = 1, nxR2
  !     fc%xR2(:,ixR2) = xR2_list(:,ixR2)
  !     fc%xR1(:,:nxR1(ixR2),ixR2) = xR1_list(:,:nxR1(ixR2),ixR2)
  !     fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
  !   enddo
  !   call fc%cryst(S, 1)
  !   call fc%cryst(S, 2)
  !   !
  !   call print_message("end of centering R1 && R2")
  !   deallocate(new_fc, xR2_list, xR1_list, nxR1)
  ! contains
  !   function inside_periodic(iR, sc_grid) result(iRp)
  !     integer, intent(in) :: iR
  !     integer, intent(in) :: sc_grid(3)
  !     !
  !     integer :: iRp
  !     integer, dimension(3) :: R, Rp, grid_
  !     !> This function returns the periodic image of R in the [-nfar:nfar]^3 box
  !     !> It is used to compute the periodic image of the force constants
  !     !> in the supercell
  !     grid_ = (2*nfar + 1) * sc_grid
  !     R = index2v(iR, grid_)
  !     Rp = mod(R-grid_/2 + sc_grid*10, sc_grid)
  !     iRp = v2index(Rp, sc_grid)
  !     !
  !   end function inside_periodic
  ! end subroutine
  !
  subroutine center5(fc, grid, S)
    use functions, only : refold_bz
    use thutils, only: grid_vec_cart, grid_vec_cryst, print_message, cryst2cart
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
    integer :: nRbig
    integer :: nxR2, ixR2, ixR1
    integer, allocatable :: R1_list(:,:), R2_list(:), nxR1(:)
    integer, allocatable :: yR2_list(:,:), yR1_list(:,:,:)
    real(dp), allocatable :: weights1(:), weights2(:)
    integer, allocatable :: inds1(:), inds2(:)
    integer, dimension(3) :: far_mesh
    real(dp) :: big_at(3,3)
    real(dp), dimension(3) :: R1_cart, R2_cart, R1, R2
    real(dp) :: block_norm
    !
    real(dp), allocatable :: new_fc(:,:,:,:)
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
        d2 = (S%tau(:,na1) - S%tau(:,na2))
        d1 = (S%tau(:,na1) + S%tau(:,na2))/2 - fc%taudef
        call inside_ws(diff_grid_cart, d2, nrws, rws, weights2, inds2)
        call inside_ws(sum_grid_cart, d1, nrws, rws, weights1, inds1)
        do iR2_big = 1, size(weights2)
          do iR1_big = 1, size(weights1)
            R1_cart = sum_grid_cart(:,inds1(iR1_big)) + diff_grid_cart(:,inds2(iR2_big))/2
            R2_cart = sum_grid_cart(:,inds1(iR1_big)) - diff_grid_cart(:,inds2(iR2_big))/2
            R1 = cryst2cart(R1_cart, S%bg, -1)
            R2 = cryst2cart(R2_cart, S%bg, -1)
            if(any(abs((NINT(R1*2) - R1*2)) > 1e-10_dp)) &
              call errore("center4", "R1_cart is not integer", 1)
            if(any(abs((NINT(R2*2) - R2*2)) > 1e-10_dp)) &
              call errore("center4", "R2_cart is not integer", 1)
            if(any(abs(NINT(R1) - R1) > 1e-10_dp)) cycle
            if(any(abs(NINT(R2) - R2) > 1e-10_dp)) cycle
            iR1 = iR_of(NINT(R1), grid)
            iR2 = iR_of(NINT(R2), grid)
            ! if (any(abs(nint(R1)) >= 2) .or. any(abs(nint(R2)) >= 2)) cycle
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
  ! subroutine center6(fc, grid, S)
  !   use functions, only : refold_bz
  !   use quter_module, only: R_list_idx
  !   use thutils, only: grid_vec_cart, grid_vec_cryst, print_message, cryst2cart
  !   class(forceconst2_sc), intent(inout) :: fc
  !   integer, intent(in) :: grid(3)
  !   type(ph_system_info), intent(in) :: S
  !   !
  !   ! type(forceconst2_grid) :: fsc
  !   real(dp), dimension(3) :: d_plus, d_minus
  !   real(dp), parameter :: eps_peri = 1e-3
  !   integer, parameter :: nperix = (2*nfar+1)**3
  !   integer :: SAFE_ALLOCATION
  !   !
  !   integer :: na1, na2, j1, j2, nR, iR1, iR2, iR_plus, iR_minus, jn1, jn2
  !   real(dp), allocatable :: far_grid_cart(:,:), R_cart(:,:)
  !   integer, allocatable :: far_grid_cryst(:,:)
  !   integer :: nRbig
  !   integer :: nxR2, ixR2, ixR1
  !   integer, pointer :: R2_list(:,:)
  !   type(list1), allocatable :: R1_list(:)
  !   integer, allocatable :: nxR1(:)
  !   integer, allocatable :: yR2_list(:,:), yR1_list(:,:,:)
  !   real(dp), allocatable :: weights_minus(:), weights_plus(:)
  !   integer, allocatable :: inds_plus(:), inds_minus(:)
  !   integer, dimension(3) :: far_mesh
  !   real(dp) :: big_at(3,3)
  !   real(dp), dimension(3) :: T1, T2, T1_try
  !   integer, dimension(3) :: R1, R2
  !   !
  !   real(dp), allocatable :: new_fc(:,:,:,:)
  !   character(100) :: filename
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws, nrws2
  !   REAL(DP) :: rws(0:3,nrwsx), rws2(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   ! initialize WS r-vectors
  !   forall(j1 = 1:3) big_at(:,j1) = S%at(:,j1) * grid(j1)
  !   CALL wsinit(rws,nrwsx,nrws,big_at)
  !   CALL wsinit(rws2,nrwsx,nrws2,big_at*2)
  !   !
  !   fc%stage = 0
  !   nR = PRODUCT(grid)
  !   if (nfar == 0) return
  !   far_mesh = 2*nfar+1
  !   nRbig = PRODUCT(far_mesh)
  !   !
  !   ! allocate(sum_grid_cart(3, 8*nRbig*nR))
  !   ! allocate(diff_grid_cart(3, nRbig*nR))
  !   ! allocate(far_grid_cryst(3, nRbig))
  !   far_grid_cart = grid_vec_cart(far_mesh, big_at, center=.true.)
  !   R_cart = grid_vec_cart(grid, S%at)
  !   ! far_grid_cryst = grid_vec_cryst(far_mesh, center=.true.)
  !   !
  !   SAFE_ALLOCATION = 20 * nR
  !   allocate(R1_list(SAFE_ALLOCATION))
  !   allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
  !   allocate(yR1_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
  !   allocate(yR2_list(3, SAFE_ALLOCATION))
  !   allocate(nxR1(SAFE_ALLOCATION))
  !   !
  !   new_fc = 0._dp
  !   nxR1 = 0
  !   nxR2 = 0
  !   !
  !   open(110, file='fc2_distance-6.dat')
  !   do na1 = 1, S%nat
  !     do na2 = 1, S%nat
  !       do iR1 = 1, nR
  !         do iR2 = 1, nR
  !           d_minus = S%tau(:,na1) + R_cart(:,iR1) - S%tau(:,na2) - R_cart(:,iR2)
  !           d_plus = S%tau(:,na1) + R_cart(:,iR1) + S%tau(:,na2) + R_cart(:,iR2)
  !           call inside_ws(far_grid_cart, d_minus, nrws, rws, weights_minus, inds_minus)
  !           call inside_ws(far_grid_cart, d_plus, nrws2, rws2, weights_plus, inds_plus)
  !           do iR_minus = 1, size(weights_minus)
  !             do iR_plus = 1, size(weights_plus)
  !               T1 = (far_grid_cart(:,inds_plus(iR_plus)) + far_grid_cart(:,inds_minus(iR_minus)))/2._dp
  !               T2 = (far_grid_cart(:,inds_plus(iR_plus)) - far_grid_cart(:,inds_minus(iR_minus)))/2._dp
  !               T1_try = cryst2cart(T1, S%bg, -1)
  !               if (any(abs(NINT(T1_try) - T1_try) > 1e-5_dp)) cycle
  !               R1 = NINT(cryst2cart(R_cart(:,iR1) + T1, S%bg, -1))
  !               R2 = NINT(cryst2cart(R_cart(:,iR2) + T2, S%bg, -1))
  !               ixR2 = R_list_idx_int(nxR2, R2_list, R2)
  !               yR2_list(:,ixR2) = R2
  !               ixR1 = R_list_idx_int(nxR1(ixR2), R1_list(ixR2)%list, R1)
  !               yR1_list(:,ixR1,ixR2) = R1
  !               !
  !               do j1 = 1, 3
  !                 jn1 = 3*(na1-1) + j1
  !                 do j2 = 1, 3
  !                   jn2 = 3*(na2-1) + j2
  !                   new_fc(jn1, jn2, ixR1, ixR2) = &
  !                     fc%FC(jn1, jn2, iR1, iR2) * &
  !                     weights_plus(iR_plus) * weights_minus(iR_minus)
  !                   ! if(abs(norm2(d_plus + T1 + T2) - 0.5_dp) < 1e-3_dp) print*, d_plus, T1, T2
  !                   write(110, "(3E20.8)") norm2(d_minus + T1 - T2), norm2(d_plus + T1 + T2), new_fc(jn1, jn2, ixR1, ixR2)
  !                 enddo
  !               enddo
  !               !
  !             enddo
  !           enddo
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   close(110)
  !   !
  !   deallocate(weights_plus, weights_minus)
  !   deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
  !   !> the yR list is populated until nxR(which), but the size can be larger.
  !   !> the values after nxR(which) are not even initialized (they are garbage).
  !   ALLOCATE(fc%yR2(3,nxR2))
  !   ALLOCATE(fc%xR2(3,nxR2))
  !   ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%n_R1(nxR2))
  !   fc%n_R1 = nxR1(:nxR2)
  !   fc%n_R2 = nxR2
  !   fc%nq = grid
  !   do ixR2 = 1, nxR2
  !     fc%yR2(:,ixR2) = yR2_list(:,ixR2)
  !     fc%yR1(:,:nxR1(ixR2),ixR2) = yR1_list(:,:nxR1(ixR2),ixR2)
  !     fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
  !   enddo
  !   call fc%cart(S, 1)
  !   call fc%cart(S, 2)
  !   !
  !   call print_message("end of centering  R1-R2 && R1+R2")
  !   deallocate(new_fc, yR2_list, yR1_list, nxR1)
  ! end subroutine
  !
  INTEGER FUNCTION R_list_idx_int(NR, listR, R) RESULT(idx)
    IMPLICIT NONE
    INTEGER, INTENT(inout)  :: NR
    integer, INTENT(inout),POINTER :: listR(:,:)
    integer, INTENT(in)    :: R(3)
    !
    INTEGER :: i
    integer, POINTER :: tmpR(:,:)
    !
    IF(NR == 0)THEN
!       print*, "new list of R initialized"
      ALLOCATE(listR(3,1))
      idx=1
      listR(:,1) = R
      NR = 1
      RETURN
    ENDIF
    !
    DO i = 1,NR
      IF(all(R-listR(:,i) == 0))THEN
        idx = i
        RETURN
      ENDIF
    ENDDO
    !
    tmpR => listR
    ALLOCATE(listR(3,NR+1))
    listr(1:3,1:NR) = tmpR(1:3,1:NR)
    DEALLOCATE(tmpR)
    NR = NR + 1
    idx = NR
    listR(1:3,idx)   = R
  END FUNCTION
  !
  ! SUBROUTINE minimal_image(S, sc_grid, diffs, n_weights, weights)
  !   use thutils, only: e_iqr, grid_vec_cart, cryst2cart
  !   type(ph_system_info), intent(in) :: S
  !   integer, intent(in) :: sc_grid(3)
  !   !
  !   real(dp) :: big_at(3,3)
  !   integer :: iR, j, na1, na2
  !   real(dp), allocatable :: big_R(:,:), R_grid(:,:)
  !   integer :: far_grid(3)
  !   real(dp) :: diff(3)
  !   real(dp), allocatable :: ws_weights(:)
  !   integer, allocatable :: inds(:)
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws
  !   REAL(DP) :: rws(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   !
  !   real(dp), allocatable, intent(out) :: diffs(:,:,:,:,:), weights(:,:,:,:)
  !   integer, allocatable, intent(out):: n_weights(:,:,:)
  !   !
  !   allocate(diffs(3,10,product(sc_grid),S%nat,S%nat))
  !   allocate(n_weights(product(sc_grid),S%nat,S%nat))
  !   allocate(weights(10,product(sc_grid),S%nat,S%nat))
  !   !
  !   far_grid = 2*nfar+1
  !   ! initialize WS r-vectors
  !   forall(j = 1:3) big_at(:,j) = S%at(:,j) * sc_grid(j)
  !   CALL wsinit(rws,nrwsx,nrws,big_at)
  !   !
  !   big_R = grid_vec_cart(far_grid, big_at, center=.true.)
  !   R_grid = grid_vec_cart(sc_grid, S%at)
  !   do iR = 1, product(sc_grid)
  !     do na2 = 1, S%nat
  !       do na1 = 1, S%nat
  !         diff = R_grid(:,iR) + S%tau(:,na1) - S%tau(:,na2)
  !         call inside_ws(big_R, diff, nrws, rws, ws_weights, inds)
  !         n_weights(iR,na1,na2) = size(ws_weights)
  !         do j = 1, size(ws_weights)
  !           diffs(:,j,ir,na1,na2) = - R_grid(:,iR) - big_R(:,inds(j))
  !           ! print*, cryst2cart(diffs(:,j,na1,na2,ir), S%bg, -1)
  !           weights(j,ir,na1,na2) = ws_weights(j)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! end subroutine
  ! !
  ! SUBROUTINE minimal_image_2(S, sc_grid, diffs, n_weights, weights)
  !   use thutils, only: e_iqr, grid_vec_cart, cryst2cart
  !   type(ph_system_info), intent(in) :: S
  !   integer, intent(in) :: sc_grid(3)
  !   !
  !   real(dp) :: big_at(3,3)
  !   integer :: iR, j, na1, na2
  !   real(dp), allocatable :: big_R(:,:), R_grid(:,:)
  !   integer :: far_grid(3)
  !   real(dp) :: diff(3)
  !   real(dp), allocatable :: ws_weights(:)
  !   integer, allocatable :: inds(:)
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws
  !   REAL(DP) :: rws(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   integer :: n_ir(product(sc_grid))
  !   integer :: R_int(3)
  !   !
  !   real(dp), allocatable, intent(out) :: diffs(:,:,:,:,:), weights(:,:,:,:)
  !   integer, allocatable, intent(out):: n_weights(:,:,:)
  !   !
  !   allocate(diffs(3,10,product(sc_grid),S%nat,S%nat))
  !   allocate(n_weights(product(sc_grid), S%nat,S%nat))
  !   allocate(weights(10,product(sc_grid),S%nat,S%nat))
  !   !
  !   far_grid = 2*nfar+1
  !   ! initialize WS r-vectors
  !   forall(j = 1:3) big_at(:,j) = S%at(:,j) * sc_grid(j)
  !   CALL wsinit(rws,nrwsx,nrws,big_at)
  !   !
  !   R_grid = grid_vec_cart(far_grid*sc_grid, S%at, center=.true.)
  !   do na2 = 1, S%nat
  !     do na1 = 1, S%nat
  !       diff = S%tau(:,na1) - S%tau(:,na2)
  !       call inside_ws(R_grid, diff, nrws, rws, ws_weights, inds)
  !       n_ir = 0
  !       do j = 1, size(ws_weights)
  !         R_int = NINT(cryst2cart(R_grid(:,inds(j)), S%bg, -1))
  !         iR = iR_of(R_int, sc_grid)
  !         n_iR(iR) = n_iR(iR) + 1
  !         diffs(:,n_ir(ir),ir,na1,na2) = - R_grid(:,inds(j))
  !         ! print*, cryst2cart(diffs(:,j,na1,na2,ir), S%bg, -1)
  !         weights(n_ir(ir),ir,na1,na2) = ws_weights(j)
  !       enddo
  !       n_weights(:,na1,na2) = n_ir
  !     enddo
  !   enddo
  ! end subroutine
  !
  ! subroutine center3(fc, grid, S, S_sc)
  !   use functions, only : refold_bz
  !   use thutils, only: grid_vec_cart, grid_vec_cryst, print_message
  !   class(forceconst2_sc), intent(inout) :: fc
  !   integer, intent(in) :: grid(3)
  !   type(ph_system_info), intent(in) :: S, S_sc
  !   !
  !   ! type(forceconst2_grid) :: fsc
  !   real(dp), dimension(3) :: d_diff, d_sum
  !   real(dp), parameter :: eps_peri = 1e-3
  !   integer, parameter :: nperix = (2*nfar+1)**3
  !   integer :: SAFE_ALLOCATION
  !   !
  !   integer :: na1, na2, j1, j2, nR, iR1, iR2, iR_diff, iR_sum
  !   integer :: map_sc(S%nat, PRODUCT(grid))
  !   integer, allocatable :: far_grid_cryst(:,:)
  !   real(dp), allocatable :: far_grid_cart(:,:), grid_cart(:,:)
  !   integer :: nRbig
  !   integer :: nxR2, ixR2, ixR1
  !   integer, allocatable :: yR1_list(:,:,:), yR2_list(:,:), yR_diff(:,:), yR_sum(:,:,:), nxR1(:)
  !   real(dp), allocatable :: weights_diff(:), weights_sum(:)
  !   integer, allocatable :: inds_diff(:), inds_sum(:)
  !   integer, dimension(3) :: far_mesh, R_diff, R_sum, R1, R2
  !   !
  !   real(dp), allocatable :: new_fc(:,:,:,:)
  !   !
  !   ! Stuff used to compute Wigner-Seitz weights:
  !   INTEGER, PARAMETER:: nrwsx=2000
  !   INTEGER :: nrws
  !   REAL(DP) :: rws(0:3,nrwsx)
  !   REAL(DP),EXTERNAL :: wsweight
  !   ! initialize WS r-vectors
  !   CALL wsinit(rws,nrwsx,nrws,S_sc%at)
  !   !
  !   fc%stage = 0
  !   nR = PRODUCT(grid)
  !   if (nfar == 0) return
  !   far_mesh = 2*nfar+1
  !   nRbig = PRODUCT(far_mesh)
  !   !
  !   far_grid_cryst = grid_vec_cryst(far_mesh, center=.true.)
  !   far_grid_cart = grid_vec_cart(far_mesh, S_sc%at, center=.true.)
  !   grid_cart = grid_vec_cart(grid, S_sc%at)
  !   !
  !   SAFE_ALLOCATION = 20 * nR
  !   allocate(new_fc(S%nat3, S%nat3, SAFE_ALLOCATION, SAFE_ALLOCATION))
  !   allocate(yR2_list(3, SAFE_ALLOCATION))
  !   allocate(yR1_list(3, SAFE_ALLOCATION, SAFE_ALLOCATION))
  !   allocate(nxR1(SAFE_ALLOCATION))
  !   !
  !   new_fc = 0._dp
  !   yR1_list = -1
  !   yR2_list = -1
  !   nxR1 = 0
  !   nxR2 = 0
  !   !
  !   open(110, file='fc2_distance-c.dat')
  !   do na2 = 1, S%nat
  !     do na1 = 1, S%nat
  !       do iR2 = 1, nR
  !         do iR1 = 1, nR
  !           d_diff = S_sc%tau(:,fc%sc_map(na1,iR1)) - S_sc%tau(:,fc%sc_map(na2,iR2))
  !           call inside_ws(far_grid_cart, d_diff, nrws, rws, weights_diff, inds_diff)
  !           do iR_diff = 1, size(weights_diff)
  !             R_diff = index2v(iR1, grid) - index2v(iR2, grid) + far_grid_cryst(:,inds_diff(iR_diff))
  !             ! call add_ind(yR1_list, inds_diff(iR_diff), nxR_diff, ixR_diff)
  !             ! yR_diff(:,ixR_diff) = index2v(R1, grid) + grid * far_grid_cryst(:,inds_diff(iR_diff))
  !             d_sum = S_sc%tau(:,fc%sc_map(na1,iR1)) + S_sc%tau(:,fc%sc_map(na2,iR2)) - fc%taudef
  !             call inside_ws(far_grid_cart, d_sum, nrws, rws, weights_sum, inds_sum)
  !             do iR_sum = 1, size(weights_sum)
  !               R_sum = index2v(iR1, grid) + index2v(iR2, grid) + grid * far_grid_cryst(:,inds_sum(iR_sum))
  !               R2 = R_sum - R_diff
  !               call add_R(yR2_list, R2, nxR2, ixR2)
  !               R1 = R_sum + R_diff
  !               call add_R(yR1_list(:,:,ixR2), R1, nxR1(ixR2), ixR1)
  !               do j1 = 1, 3
  !                 do j2 = 1, 3
  !                   new_fc(j1+3*(na1-1), j2+3*(na2-1), ixR1, ixR2) = &
  !                     fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, iR1, iR2) * &
  !                     weights_diff(iR_diff) * weights_sum(iR_sum)
  !                   ! write(10, '(2E15.5)') norm2(fc2_sc%xR1(:,R1,R2) - fc2_sc%xR2(:,R2)), &
  !                   !   SUM(ABS(fc2_sc%fc((na1-1)*3+1:na1*3,(na2-1)*3+1:na2*3,R1,R2)))
  !                 enddo
  !               enddo
  !             enddo
  !           enddo
  !           !
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  !   !
  !   deallocate(weights_diff, weights_sum)
  !   deallocate(fc%yR2, fc%xR2, fc%FC, fc%xR1, fc%yR1, fc%n_R1)
  !   !> the yR list is populated until nxR(which), but the size can be larger.
  !   !> the values after nxR(which) are not even initialized (they are garbage).
  !   ALLOCATE(fc%yR2(3,nxR2))
  !   ALLOCATE(fc%xR2(3,nxR2))
  !   ALLOCATE(fc%yR1(3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%xR1(3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%FC(S%nat3,S%nat3,maxval(nxR1),nxR2))
  !   ALLOCATE(fc%n_R1(nxR2))
  !   fc%n_R1 = nxR1(:nxR2)
  !   fc%n_R2 = nxR2
  !   fc%nq = grid
  !   do ixR2 = 1, nxR2
  !     fc%yR2(:,ixR2) = yR2_list(:,ixR2)
  !     fc%yR1(:,:nxR1(ixR2),ixR2) = yR1_list(:,:nxR1(ixR2),ixR2)
  !     fc%FC(:,:,:nxR1(ixR2),ixR2) = new_fc(:,:,:nxR1(ixR2),ixR2)
  !   enddo
  !   call fc%cart(S_sc, 1)
  !   call fc%cart(S_sc, 2)

  !   ! do ixR2 = 1, fc%n_R2
  !   !   fc%xR2(:,ixR2)= fc%xR2(:,ixR2) / grid
  !   !   do ixR1 = 1, fc%n_R1(ixR2)
  !   !     fc%xR1(:,ixR1, ixR2) = fc%xR1(:,ixR1, ixR2) / grid
  !   !   enddo
  !   ! enddo
  !   !
  !   call print_message("end of centering")
  !   deallocate(new_fc, yR2_list, yR1_list, nxR1)
  ! end subroutine
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
  function fc_RR2sc(grid, S, S_sc, fc)
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    real(dp), intent(in) :: fc(S%nat3,S%nat3, product(grid), product(grid))
    !
    real(dp) :: fc_RR2sc(S_sc%nat3, S_sc%nat3)
    integer :: sc_na1, sc_na2, j1, j2
    integer, dimension(S_sc%nat) :: map_R, map_nat
    !
    map_R = map_sc2uc(S, S_sc, grid, "R")
    map_nat = map_sc2uc(S, S_sc, grid, "nat")
    fc_RR2sc = 0._dp
    do sc_na1 = 1, S_sc%nat
      do sc_na2 = 1, S_sc%nat
        do j1 = 1, 3
          do j2 = 1, 3
            fc_RR2sc(j1+(sc_na1-1)*3, j2+(sc_na2-1)*3) = &
              fc(j1+(map_nat(sc_na1)-1)*3, j2+(map_nat(sc_na2)-1)*3, map_R(sc_na1), map_R(sc_na2))
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
  subroutine asr3_div_mass_RR(asr, S, Sd, grid, fc2RR)
    character(*), intent(in) :: asr
    type(ph_system_info), intent(in) :: S, Sd
    integer, intent(in) :: grid(3)
    real(dp), intent(inout) :: fc2RR(:,:,:,:)
    !
    if(trim(asr) /= 'no') call asr3(fc2RR, asr)
    call div_mass_RR(S, Sd, grid, fc2RR)
    if(trim(asr) /= 'no') call print_message("ASR applied to DRR")
  end subroutine
  !
  subroutine div_mass_RR(S, Sd, grid, fc2RR)
    type(ph_system_info), intent(in) :: S, Sd
    integer, intent(in) :: grid(3)
    real(dp), intent(inout) :: fc2RR(:,:,:,:)
    !
    integer, allocatable :: sc_map_rr(:,:)
    integer :: na1, na2, j1, j2, jn1, jn2
    integer :: iR1, iR2, isc1, isc2, nR
    !
    if(.not. allocated(Sd%sqrtmm1)) &
      call errore('div_mass_RR', 'missing sqrtmm1 in Sd, call aux_system first', 1)
    if(size(fc2RR, 1) /= S%nat3 .or. size(fc2RR, 2) /= S%nat3) &
      call errore('div_mass_RR', 'RR atom dimensions do not match S', 1)
    !
    nR = product(grid)
    if(size(fc2RR, 3) /= nR .or. size(fc2RR, 4) /= nR) &
      call errore('div_mass_RR', 'RR grid dimensions do not match grid', 1)
    !
    sc_map_rr = map_uc2sc(S, Sd, grid)
    if(any(sc_map_rr < 1)) &
      call errore('div_mass_RR', 'some RR atoms are not mapped to Sd', count(sc_map_rr < 1))
    !
    do iR2 = 1, nR
      do na2 = 1, S%nat
        isc2 = sc_map_rr(na2, iR2)
        do j2 = 1, 3
          jn2 = j2 + 3*(na2-1)
          do iR1 = 1, nR
            do na1 = 1, S%nat
              isc1 = sc_map_rr(na1, iR1)
              do j1 = 1, 3
                jn1 = j1 + 3*(na1-1)
                fc2RR(jn1, jn2, iR1, iR2) = fc2RR(jn1, jn2, iR1, iR2) * &
                  Sd%sqrtmm1(j1 + 3*(isc1-1)) * Sd%sqrtmm1(j2 + 3*(isc2-1))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !
  end subroutine
  !
  subroutine div_mass0_RR(S, Sd, grid, fc2RR)
    type(ph_system_info), intent(in) :: S, Sd
    integer, intent(in) :: grid(3)
    real(dp), intent(inout) :: fc2RR(:,:,:,:)
    !
    integer :: na1, na2, j1, j2, jn1, jn2
    integer :: iR1, iR2, nR
    !
    if(size(fc2RR, 1) /= S%nat3 .or. size(fc2RR, 2) /= S%nat3) &
      call errore('div_mass0_RR', 'RR atom dimensions do not match S', 1)
    !
    nR = product(grid)
    if(size(fc2RR, 3) /= nR .or. size(fc2RR, 4) /= nR) &
      call errore('div_mass0_RR', 'RR grid dimensions do not match grid', 1)
    !
    do iR2 = 1, nR
      do na2 = 1, S%nat
        do j2 = 1, 3
          jn2 = j2 + 3*(na2-1)
          do iR1 = 1, nR
            do na1 = 1, S%nat
              do j1 = 1, 3
                jn1 = j1 + 3*(na1-1)
                fc2RR(jn1, jn2, iR1, iR2) = fc2RR(jn1, jn2, iR1, iR2) * &
                  S%sqrtmm1(j1 + 3*(na1-1)) * S%sqrtmm1(j2 + 3*(na2-1))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !
  end subroutine
  !
  subroutine freq_in_grid_degen(S, fc2, fc2sc, grid, freqs, Us, freqs1)
    use q_grids, only : q_grid
    use mpi_thermal, only : mpi_bsum
    use merge_degenerate, only: merge_degen
    use defutils, only : freq_phq_degen
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
  !
  ! subroutine asr2(fc2)
  !   real(dp), intent(inout) :: fc2(:,:,:)
  !   integer :: na1, na2, j1, j2, jn1, jn2
  !   integer :: iR, nR, nat
  !   real(dp) :: sum
  !   !
  !   nR = size(fc2, 3)
  !   nat = size(fc2, 1)/3
  !   !
  !   do na1 = 1, nat
  !     do j1 = 1, 3
  !       jn1 = j1 + 3*(na1-1)
  !       do j2 = 1, 3
  !         sum = 0._dp
  !         do iR = 1, nR
  !           do na2 = 1, nat
  !             jn2 = j2 + 3*(na2-1)
  !             sum = sum + fc2(jn1, jn2, iR) + fc2(j2 + 3*(na1-1), j1 + 3*(na2-1), iR)
  !           enddo
  !         enddo
  !       enddo
  !       fc2(jn1, jn1, fc2%i_0) = fc2(jn1, jn1, fc2%i_0) - sum/2
  !     enddo
  !   enddo
  !   !
  ! end subroutine
  !
  subroutine asr3(fc2RR, method)
    ! Enforce translational ASR on a defect FC matrix while preserving the
    ! exchange symmetry Phi(I,J)=Phi(J,I).  I=(R,atom,Cartesian).
    !
    ! method='local' applies the least-Frobenius-norm symmetric correction on
    ! the existing FC support.  method='diff' retains the previous diagonal
    ! correction for comparison; it does not generally satisfy ASR exactly.
    ! method='project' applies the reference dense P*Phi*P projection, where
    ! P removes the three uniform translations.  Both give Phi*T=0.
    use mpi_thermal, only : ionode
    real(dp), intent(inout) :: fc2RR(:,:,:,:)
    character(*), intent(in), optional :: method
    character(16) :: selected_method
    integer :: nR, nat, n, iR1, iR2, i, j, a, b, row, col
    real(dp) :: max_asr_before, max_asr_after, max_exchange
    real(dp), allocatable :: matrix(:,:), asr_before(:,:,:,:), asr_after(:,:,:,:)
    !
    nR = size(fc2RR, 4)
    nat = size(fc2RR, 1) / 3
    n = 3 * nat * nR
    selected_method = 'local'
    if (present(method)) selected_method = adjustl(trim(method))
    if (selected_method == 'legacy') selected_method = 'diff'
    !
    allocate(asr_before(3,3,nat,nR), asr_after(3,3,nat,nR), matrix(n,n))
    call asr3_residual(fc2RR, asr_before)
    max_asr_before = maxval(abs(asr_before))
    !
    do iR1 = 1, nR
      do i = 1, nat
        do a = 1, 3
          row = a + 3*(i-1) + 3*nat*(iR1-1)
          do iR2 = 1, nR
            do j = 1, nat
              do b = 1, 3
                col = b + 3*(j-1) + 3*nat*(iR2-1)
                matrix(row,col) = fc2RR(a + 3*(i-1), b + 3*(j-1), iR1, iR2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    select case (trim(selected_method))
    case ('local')
      ! Symmetrize once before applying the constrained correction, so both
      ! exchange symmetry and ASR are imposed simultaneously.
      matrix = 0.5_dp * (matrix + transpose(matrix))
      call asr3_local_projection(matrix)
    case ('project')
      matrix = 0.5_dp * (matrix + transpose(matrix))
      call asr3_global_projection(matrix)
    case ('diff')
      call asr3_diff_correction(matrix, nat, nR)
    case default
      call errore('asr3', 'unknown ASR method: '//trim(selected_method), 1)
    end select
    !
    do iR1 = 1, nR
      do i = 1, nat
        do a = 1, 3
          row = a + 3*(i-1) + 3*nat*(iR1-1)
          do iR2 = 1, nR
            do j = 1, nat
              do b = 1, 3
                col = b + 3*(j-1) + 3*nat*(iR2-1)
                fc2RR(a + 3*(i-1), b + 3*(j-1), iR1, iR2) = matrix(row,col)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !
    call asr3_residual(fc2RR, asr_after)
    max_asr_after = maxval(abs(asr_after))
    max_exchange = maxval(abs(matrix - transpose(matrix)))
    if (ionode) then
      write(*,"(A,A)") 'ASR3 method: ', trim(selected_method)
      write(*,"(A,E20.8)") 'max exchange residual after ASR3: ', max_exchange
      write(*,"(A,E20.8)") 'max ASR residual before ASR3: ', max_asr_before
      write(*,"(A,E20.8)") 'max ASR residual after ASR3:  ', max_asr_after
      open(91, file='asr3_residual.dat', status='replace', action='write')
      write(91,'(A,A)') '# method ', trim(selected_method)
      write(91,'(A)') '# R1 atom cart |ASR residual| before |ASR residual| after'
      do iR1 = 1, nR
        do i = 1, nat
          do a = 1, 3
            write(91,'(3I7,2ES24.14)') iR1, i, a, &
              sqrt(sum(asr_before(a,:,i,iR1)**2)), sqrt(sum(asr_after(a,:,i,iR1)**2))
          enddo
        enddo
      enddo
      close(91)
    endif
  end subroutine asr3
  !
  subroutine asr3_diff_correction(matrix, nat, nR)
    ! Original ASR3 implementation: alter only the same-atom, same-cell
    ! blocks.  It is retained to quantify its residual, not as a strict ASR.
    real(dp), intent(inout) :: matrix(:,:)
    integer, intent(in) :: nat, nR
    integer :: iR1, iR2, i, j, a, b, row, col, row_t, col_t, n3
    real(dp) :: delta
    !
    n3 = 3 * nat
    do iR1 = 1, nR
      do i = 1, nat
        do a = 1, 3
          do b = 1, 3
            delta = 0._dp
            do iR2 = 1, nR
              do j = 1, nat
                row = a + 3*(i-1) + n3*(iR1-1)
                col = b + 3*(j-1) + n3*(iR2-1)
                row_t = b + 3*(i-1) + n3*(iR1-1)
                col_t = a + 3*(j-1) + n3*(iR2-1)
                delta = delta + matrix(row,col) + matrix(row_t,col_t)
              enddo
            enddo
            row = a + 3*(i-1) + n3*(iR1-1)
            col = b + 3*(i-1) + n3*(iR1-1)
            matrix(row,col) = matrix(row,col) - delta / 2._dp
          enddo
        enddo
      enddo
    enddo
  end subroutine asr3_diff_correction
  !
  subroutine asr3_residual(fc2RR, residual)
    real(dp), intent(in) :: fc2RR(:,:,:,:)
    real(dp), intent(out) :: residual(:,:,:,:)
    integer :: nR, nat, iR1, iR2, i, j, a, b
    ! residual(a,b,i,R1) = sum_{R2,j} Phi_ia,jb(R1,R2)
    nR = size(fc2RR, 4)
    nat = size(fc2RR, 1) / 3
    if (any(shape(residual) /= [3, 3, nat, nR])) &
      call errore('asr3_residual', 'incompatible residual shape', 1)
    residual = 0._dp
    do iR1 = 1, nR
      do i = 1, nat
        do a = 1, 3
          do b = 1, 3
            do iR2 = 1, nR
              do j = 1, nat
                residual(a,b,i,iR1) = residual(a,b,i,iR1) + &
                  fc2RR(a + 3*(i-1), b + 3*(j-1), iR1, iR2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine asr3_residual
  !
  subroutine asr3_global_projection(matrix)
    ! K <- P K P, P = I - T (T^T T)^-1 T^T for uniform translations.
    real(dp), intent(inout) :: matrix(:,:)
    integer :: n, n_per_cart, i, j, alpha, beta
    real(dp), allocatable :: row_sum(:,:), col_sum(:,:), total_sum(:,:)
    !
    n = size(matrix, 1)
    if (size(matrix, 2) /= n .or. mod(n, 3) /= 0) &
      call errore('asr3_global_projection', 'matrix must be square with dimension divisible by 3', 1)
    n_per_cart = n / 3
    allocate(row_sum(n,3), col_sum(3,n), total_sum(3,3))
    row_sum = 0._dp
    col_sum = 0._dp
    total_sum = 0._dp
    do i = 1, n
      do j = 1, n
        beta = modulo(j-1, 3) + 1
        alpha = modulo(i-1, 3) + 1
        row_sum(i,beta) = row_sum(i,beta) + matrix(i,j)
        col_sum(alpha,j) = col_sum(alpha,j) + matrix(i,j)
        total_sum(alpha,beta) = total_sum(alpha,beta) + matrix(i,j)
      enddo
    enddo
    do i = 1, n
      alpha = modulo(i-1, 3) + 1
      do j = 1, n
        beta = modulo(j-1, 3) + 1
        matrix(i,j) = matrix(i,j) - row_sum(i,beta) / real(n_per_cart,dp) &
          - col_sum(alpha,j) / real(n_per_cart,dp) &
          + total_sum(alpha,beta) / real(n_per_cart*n_per_cart,dp)
      enddo
    enddo
  end subroutine asr3_global_projection
  !
  subroutine asr3_local_projection(matrix)
    ! Minimum-Frobenius-norm symmetric correction with support restricted to
    ! FC elements already nonzero (plus all on-site 3x3 blocks).  It solves
    ! min ||dK||_F subject to dK=dK^T and (K+dK)T=0.
    real(dp), intent(inout) :: matrix(:,:)
    integer :: n, m, nvar, i, j, c1, c2, info
    integer, allocatable :: ii(:), jj(:), ipiv(:)
    logical, allocatable :: allowed(:,:)
    real(dp) :: threshold, weight_inverse, correction
    real(dp), allocatable :: normal(:,:), rhs(:)
    !
    n = size(matrix, 1)
    if (size(matrix, 2) /= n .or. mod(n, 3) /= 0) &
      call errore('asr3_local_projection', 'matrix must be square with dimension divisible by 3', 1)
    m = 3 * n
    threshold = max(1.e-14_dp, 1.e-12_dp * maxval(abs(matrix)))
    allocate(allowed(n,n))
    allowed = .false.
    do i = 1, n
      allowed(i,i) = .true.
      do j = i + 1, n
        if (abs(matrix(i,j)) > threshold .or. abs(matrix(j,i)) > threshold) then
          allowed(i,j) = .true.
          allowed(j,i) = .true.
        endif
      enddo
    enddo
    nvar = 0
    do i = 1, n
      do j = i, n
        if (allowed(i,j)) nvar = nvar + 1
      enddo
    enddo
    allocate(ii(nvar), jj(nvar), normal(m,m), rhs(m), ipiv(m))
    nvar = 0
    do i = 1, n
      do j = i, n
        if (.not. allowed(i,j)) cycle
        nvar = nvar + 1
        ii(nvar) = i
        jj(nvar) = j
      enddo
    enddo
    ! rhs is -K*T.  Constraint (i,beta) is stored at i+n*(beta-1).
    rhs = 0._dp
    do i = 1, n
      do j = 1, n
        rhs(i + n*(modulo(j-1,3))) = rhs(i + n*(modulo(j-1,3))) - matrix(i,j)
      enddo
    enddo
    ! Assemble C W^-1 C^T.  A diagonal correction has Frobenius weight 1;
    ! an off-diagonal symmetric pair has weight 2.
    normal = 0._dp
    do i = 1, nvar
      c1 = ii(i) + n * modulo(jj(i)-1, 3)
      c2 = jj(i) + n * modulo(ii(i)-1, 3)
      weight_inverse = 1._dp
      if (ii(i) /= jj(i)) weight_inverse = 0.5_dp
      normal(c1,c1) = normal(c1,c1) + weight_inverse
      if (c2 /= c1) then
        normal(c2,c2) = normal(c2,c2) + weight_inverse
        normal(c1,c2) = normal(c1,c2) + weight_inverse
        normal(c2,c1) = normal(c2,c1) + weight_inverse
      endif
    enddo
    call dgesv(m, 1, normal, m, ipiv, rhs, m, info)
    if (info /= 0) call errore('asr3_local_projection', 'singular local ASR constraint system', info)
    do i = 1, nvar
      c1 = ii(i) + n * modulo(jj(i)-1, 3)
      c2 = jj(i) + n * modulo(ii(i)-1, 3)
      weight_inverse = 1._dp
      if (ii(i) /= jj(i)) weight_inverse = 0.5_dp
      correction = weight_inverse * rhs(c1)
      if (c2 /= c1) correction = correction + weight_inverse * rhs(c2)
      matrix(ii(i),jj(i)) = matrix(ii(i),jj(i)) + correction
      if (ii(i) /= jj(i)) matrix(jj(i),ii(i)) = matrix(jj(i),ii(i)) + correction
    enddo
  end subroutine asr3_local_projection
  !
  subroutine write_fc2_sc_RR(S, fc2_sc, filename)
    type(ph_system_info), intent(in) :: S
    type(forceconst2_sc), intent(in) :: fc2_sc
    character(*), intent(in) :: filename
    integer :: na1, na2, iR1, iR2
    real(dp) :: block_norm
    !
    open(113, file=filename, status='replace', action='write')
    do na1 = 1, S%nat
      do na2 = 1, S%nat
        do iR2 = 1, fc2_sc%n_R2
          do iR1 = 1, fc2_sc%n_R1(iR2)
            block_norm = sqrt(sum(fc2_sc%fc(3*(na1-1)+1:3*na1, 3*(na2-1)+1:3*na2, iR1, iR2)**2))
            write(113, "(3E20.8,4I8)") norm2(fc2_sc%xR1(:,iR1,iR2) - fc2_sc%xR2(:,iR2) + S%tau(:,na1) - S%tau(:,na2)), &
              norm2(fc2_sc%xR1(:,iR1,iR2) + fc2_sc%xR2(:,iR2) + S%tau(:,na1) + S%tau(:,na2) - 2*fc2_sc%taudef), &
              block_norm, na1, na2, iR1, iR2
          enddo
        enddo
      enddo
    enddo
    close(113)
  end subroutine
  !
end module
