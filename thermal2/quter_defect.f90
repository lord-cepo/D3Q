module quter_defect
  use kinds,              only : dp
  use fc2_interpolate,    only : forceconst2_grid
  use ph_system,          only : ph_system_info
  use thutils,            only : v2index, index2v
  implicit none
  !
  integer, parameter :: nfar = 0
  !
  type forceconst2_sc
    INTEGER :: n_R(2) = 0, i_0(2) = -1
    real(dp), allocatable :: FC(:,:,:,:) ! (jn1,jn2,iR,jR)
    complex(dp), allocatable :: mix(:,:,:) ! (jn1,jn2,iR)
    integer,  allocatable :: yR(:,:,:)    ! (3,iR,which)
    real(dp), allocatable :: xR(:,:,:)    ! (3,iR,which)
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
  function map_uc2sc(S, S_sc, sc_grid)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    !
    integer :: map_uc2sc(S%nat, PRODUCT(sc_grid))
    integer :: i, isc, iR
    real(dp) :: r_cryst(3)
    map_uc2sc = -1
    do i = 1, S%nat
      do isc = 1, S_sc%nat
        r_cryst = S_sc%tau(:,isc)*sc_grid-S%tau(:,i)
        call cryst_to_cart(1, r_cryst, S_sc%bg, -1)
        if (NORM2(r_cryst - NINT(r_cryst))<1e-2) then
          iR = v2index(NINT(r_cryst), sc_grid)
          if (map_uc2sc(i, iR) > -1) CALL errore("map_uc2sc", "found same atom in same R", 1)
          if (iR < 1 .or. iR > PRODUCT(sc_grid)) CALL errore("map_uc2sc", "R is out of bound", ABS(iR))
          map_uc2sc(i, iR) = isc
        endif
        !
      enddo
    enddo
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
    integer :: n_R
    !
    IF(allocated(fc%yR) .or. allocated(fc%FC) &
      .or. allocated(fc%xR) .or. allocated(fc%mix)) &
      CALL errore(sub, 'some element is already allocated', 1)
    !
    n_R = PRODUCT(grid)
    fc%n_R(:) = n_R
    !
    ALLOCATE(fc%yR(3,n_R,2))
    fc%yR(:,:,1) = grid_vec(grid)
    fc%yR(:,:,2) = grid_vec(grid)
    !
    ALLOCATE(fc%xR(3,n_R,2))
    call fc%cart(S)
    ALLOCATE(fc%FC(S%nat3,S%nat3,n_R,n_R))
    fc%nq = grid
    fc%stage = -1
    fc%i_0(:) = 1
    !
  END SUBROUTINE
  !
  subroutine construct_cart(fc, S)
    class(forceconst2_sc), intent(inout) :: fc
    type(ph_system_info), intent(in) :: S
    !
    real(dp), allocatable :: R(:,:)
    integer :: i
    !
    do i = 1, 2
      allocate(R(3,fc%n_R(i)))
      R = REAL(fc%yR(:,:,i), DP)
      call cryst_to_cart(fc%n_R(i), R, S%at, 1)
      fc%xR(:,:,i) = R
      deallocate(R)
    enddo
    !
  end subroutine
  !
  subroutine construct_cryst(fc, S)
    class(forceconst2_sc), intent(inout) :: fc
    type(ph_system_info), intent(in) :: S
    !
    real(dp), allocatable :: R(:,:)
    integer :: i
    !
    do i = 1, 2
      allocate(R(3,fc%n_R(i)))
      R = fc%xR(:,:,i)
      call cryst_to_cart(fc%n_R(i), R, S%bg, -1)
      if (ALL(ABS(R - NINT(R)) < 1e-6)) then
        fc%yR(:,:,i) = NINT(R)
      else
        call errore("construct_cryst", "R is not integer", 1)
      endif
      deallocate(R)
    enddo
  end subroutine
  !
  subroutine center_sc(fc, grid, S, S_sc)
    class(forceconst2_sc), intent(inout) :: fc
    integer, intent(in) :: grid(3)
    type(ph_system_info), intent(in) :: S, S_sc
    !
    ! type(forceconst2_grid) :: fsc
    real(dp) :: dist(3), Rbig(3), wg_tot
    !
    integer :: na1, na2, j1, j2, nR, jR_big, R1, R2
    integer :: R_list(PRODUCT(grid)*(2*nfar+1)**3, 2)
    integer :: map_sc(S%nat, PRODUCT(grid))
    integer :: nRbig, R_vec(3,2), which
    integer, dimension(2) :: ind, nxR, ixR
    integer, dimension(3) :: far_grid, Rbig_from_0
    !
    integer, parameter :: SAFE_LIMIT = 10
    !
    integer, allocatable :: new_yR_list(:,:,:)
    real(dp), allocatable :: new_fc(:,:,:,:)
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
    if (nfar == 0) then
      call fc%cart(S_sc)
      return
    endif
    far_grid = 2*nfar+1
    nRbig = PRODUCT(far_grid)
    allocate(new_yR_list(3, nR*SAFE_LIMIT,2))
    allocate(new_fc(S%nat3, S%nat3, nR, nR*SAFE_LIMIT))
    !
    map_sc = map_uc2sc(S, S_sc, grid)
    !
    nxR = 0
    new_fc = 0._dp
    !
    R_list = -1
    do na1 = 1, S%nat
      do na2 = 1, S%nat
        do R1 = 1, nR
          R_vec(:,1) = index2v(R1, grid)
          ind(1) = R1
          do R2 = 1, nR
            wg_tot = 0._dp
            do jR_big = 1, nRbig
              Rbig_from_0 = index2v(jR_big, far_grid)
              Rbig = REAL(Rbig_from_0 - nfar, DP)
              call cryst_to_cart(1, Rbig, S_sc%at, 1)
              dist = S_sc%tau(:,map_sc(na1,R1)) - Rbig - S_sc%tau(:,map_sc(na2,R2))
              wg = wsweight(dist,rws,nrws)
              wg_tot = wg_tot + wg
              ! if (nfar == 0) wg = 1._dp
              if (wg > 1e-6) then
                R_vec(:,2) = index2v(R2, grid) + grid * (Rbig_from_0 - nfar)
                ind(2) = v2index(R_vec(:,2) + grid * Rbig_from_0, grid * far_grid)
                do which = 1, 2
                  if (R_list(ind(which),which) == -1) then
                    nxR(which) = nxR(which) + 1
                    ixR(which) = nxR(which)
                    R_list(ind(which),which) = ixR(which)
                    new_yR_list(:,ixR(which),which) = R_vec(:,which)
                  else
                    ixR(which) = R_list(ind(which), which)
                  endif
                enddo
                if(ALL(R_vec == 0)) fc%i_0 = ixR
                do j1 = 1, 3
                  do j2 = 1, 3
                    new_fc(j1+3*(na1-1),j2+3*(na2-1),ixR(1),ixR(2)) = &
                      fc%FC(j1+(na1-1)*3, j2+(na2-1)*3, R1, R2) * wg
                  enddo
                enddo
              endif
            enddo
            if (ABS(wg_tot -1)<1e-6) CALL errore("center_sc", "sum of weights is not 1", 1)
          enddo
        enddo
      enddo
    enddo
    !
    deallocate(fc%yR, fc%xR, fc%FC)
    ALLOCATE(fc%yR(3,maxval(nxR),2))
    ALLOCATE(fc%xR(3,maxval(nxR),2))
    ALLOCATE(fc%FC(S%nat3,S%nat3,nxR(1),nxR(2)))
    fc%n_R(:) = nxR
    fc%nq = grid
    do which = 1, 2
      fc%yR(:,1:nxR(which),which) = new_yR_list(:,1:nxR(which),which)
    enddo
    call fc%cart(S_sc)
    fc%FC = new_fc(:,:,1:nxR(1),1:nxR(2))
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
    map_R = map_sc2uc(S, S_sc, grid, "R")
    map_nat = map_sc2uc(S, S_sc, grid, "nat")
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
  !
  function fc_uc2sc(S, S_sc, sc_grid, fc2)
    TYPE(ph_system_info), intent(in) :: S, S_sc
    integer, intent(in) :: sc_grid(3)
    type(forceconst2_grid), intent(in) :: fc2
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
                  fc2%FC(jn1, jn2, v2index(R, sc_grid))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end function
  !
  SUBROUTINE interp_1st_step(fc, xq, S, which)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc),INTENT(inout) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    integer, intent(in) :: which
    TYPE(ph_system_info),INTENT(in) :: S
    !
    INTEGER :: i, j
    REAL(DP), dimension(fc%n_R(which)) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(fc%n_R(which))
    !
    if(allocated(fc%mix)) deallocate(fc%mix)
    !
    if (which == 1 .or. which == 2) then
      fc%stage = 3-which
      allocate(fc%mix(S%nat3, S%nat3, fc%n_R(fc%stage)))
      FORALL(i=1:fc%n_R(which)) varg(i) =  tpi * &
        dot_product(xq, fc%xR(:,i,which))
    else
      CALL errore("interp_1st_step", "the variable 'which' is not 1 or 2", 1)
    endif
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
    CALL vdCos(fc%n_R(which), varg, vcos)
    CALL vdSin(fc%n_R(which), varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    fc%mix = 0._dp
    do j = 1, fc%n_R(2)
      do i = 1, fc%n_R(1)
        if (which == 1) then
          fc%mix(:,:,j) = fc%mix(:,:,j) + vphase(i) * fc%fc(:,:,i,j)
        else
          fc%mix(:,:,i) = fc%mix(:,:,i) + vphase(j) * fc%fc(:,:,i,j)
        endif
      enddo
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
    REAL(DP), dimension(fc%n_R(fc%stage)) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(fc%n_R(fc%stage))
    !
    if (fc%stage == 1 .or. fc%stage == 2) then
      FORALL(i=1:fc%n_R(fc%stage)) varg(i) =  tpi * &
        dot_product(xq, fc%xR(:,i,fc%stage))
    else
      CALL errore("interp_2nd_step", "the variable 'stage' is not 1 or 2", 1)
    endif
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
    CALL vdCos(fc%n_R(which), varg, vcos)
    CALL vdSin(fc%n_R(which), varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    D = 0._dp
    do i = 1, fc%n_R(fc%stage)
      D = D + vphase(i) * fc%mix(:,:,i)
    enddo
    !
  END SUBROUTINE
  !
  SUBROUTINE interp_at_once(fc, xq1, xq2, S, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    CLASS(forceconst2_sc), INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq1(3), xq2(3)
    TYPE(ph_system_info), INTENT(in) :: S
    complex(dp), intent(out) :: D(S%nat3, S%nat3)
    !
    INTEGER :: i, j
    REAL(DP), dimension(fc%n_R(1),fc%n_R(2)) :: varg, vcos, vsin
    COMPLEX(DP) :: vphase(fc%n_R(1),fc%n_R(2))
    !
    FORALL(i=1:fc%n_R(1), j=1:fc%n_R(2)) varg(i,j) = &
      tpi * (dot_product(xq1, fc%xR(:,i,1)) + &
      dot_product(xq2, fc%xR(:,j,2)))
    !
    vcos = DCOS(varg)
    vsin = DSIN(varg)
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    D = 0._dp
    do j = 1, fc%n_R(2)
      do i = 1, fc%n_R(1)
        D = D + vphase(i,j) * fc%fc(:,:,i,j)
      enddo
    enddo
    !
  END SUBROUTINE
  !
end module
