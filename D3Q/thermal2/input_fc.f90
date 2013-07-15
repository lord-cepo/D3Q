!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE input_fc
  !
  USE kinds,            ONLY : DP
  USE parameters,       ONLY : ntypx
  USE io_global,        ONLY : stdout
  !
  ! \/o\________\\\_________________________________________/^>
  TYPE ph_system_info
    ! atoms
    INTEGER              :: ntyp
    REAL(DP)             :: amass(ntypx)
    CHARACTER(len=3  )   :: atm(ntypx)
    ! atoms basis
    INTEGER              :: nat
    REAL(DP),ALLOCATABLE :: tau(:,:), zeu(:,:,:)
    INTEGER, ALLOCATABLE :: ityp(:)
    ! unit cell, and reciprocal
    INTEGER              :: ibrav
    CHARACTER(len=9)     :: symm_type
    REAL(DP)             :: celldm(6), at(3,3), bg(3,3)
    REAL(DP)             :: omega
    ! phonon switches (mostly unused here)
    REAL(DP)             :: epsil(3,3)
    LOGICAL              :: lrigid
    ! ''''''''''''''''''''''''''''''''''''''''
    ! auxiliary quantities:
    REAL(DP),ALLOCATABLE :: sqrtmm1(:) ! 1/sqrt(amass)
    INTEGER :: nat3, nat32, nat33
  END TYPE
  ! \/o\________\\\_________________________________________/^>
  TYPE forceconst2_grid
    ! q points
    INTEGER :: n_R = 0, i_0 = -1
    INTEGER,ALLOCATABLE  :: yR(:,:) ! crystalline coords  3*n_R
    REAL(DP),ALLOCATABLE :: xR(:,:) ! cartesian coords    3*n_R
    REAL(DP),ALLOCATABLE :: FC(:,:,:) ! 3*nat,3*nat, n_R
  END TYPE forceconst2_grid
  ! \/o\________\\\_________________________________________/^>
  TYPE forceconst3_grid
    ! q points
    INTEGER :: n_R = 0
    INTEGER :: i_00 = -1 ! set to the index where both R2 and R3 are zero
    INTEGER,ALLOCATABLE :: i_R20(:), i_R30(:) ! lists of indexes where R2 (R3) is zero
    !INTEGER,ALLOCATABLE :: iR2eR3(:) ! list of indexes R2(i) =  R3(i_R2eR3(i))
    !INTEGER,ALLOCATABLE :: iR3eR2(:) ! list of indexes R3(i) =  R2(i_R3eR2(i))
    !INTEGER,ALLOCATABLE :: iR2mR3(:) ! list of indexes R2(i) = -R3(i_R2mR3(i))
    !INTEGER,ALLOCATABLE :: iR3mR2(:) ! list of indexes R3(i) = -R2(i_R3mR2(i))
    
    INTEGER,ALLOCATABLE  :: yR2(:,:), yR3(:,:) ! crystalline coords  3*n_R
    REAL(DP),ALLOCATABLE :: xR2(:,:), xR3(:,:) ! cartesian coords    3*n_R
    REAL(DP),ALLOCATABLE :: FC(:,:,:,:) ! 3*nat,3*nat, n_R
  END TYPE forceconst3_grid
  !
  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  
  FUNCTION same_system(S,Z) RESULT (same)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in) :: S,Z
    LOGICAL :: same
    REAL(DP),PARAMETER :: eps = 1.d-6
    !
    ! NOT checking : atm, amass, symm_type
    !
    same = .true.
    same = same .and. (S%ntyp == Z%ntyp)
    same = same .and. (S%nat == Z%nat)
    same = same .and. (S%ibrav == Z%ibrav)
    same = same .and. ALL( S%ityp(1:S%ntyp) == Z%ityp(1:Z%ntyp))
    
    IF(allocated(S%tau).and.allocated(Z%tau)) &
      same = same .and. ALL( ABS(S%tau -Z%tau) < eps)
    IF(allocated(S%zeu).and.allocated(Z%zeu)) &
      same = same .and. ALL( ABS(S%zeu -Z%zeu) < eps)

    same = same .and. ALL( ABS(S%celldm -Z%celldm) < eps)
    same = same .and. ALL( ABS(S%at -Z%at) < eps)
    same = same .and. ALL( ABS(S%bg -Z%bg) < eps)
    same = same .and. ( ABS(S%omega -Z%omega) < eps)

!     same = same .and. (S%lrigid .or. Z%lrigid)
    same = same .and. ALL( ABS(S%epsil -Z%epsil) < eps)
    
  END FUNCTION same_system
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE read_system(unit, S)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    !
    CHARACTER(len=11),PARAMETER :: sub = "read_system"
    !
    INTEGER :: unit, ios, dummy
    !
    INTEGER :: nt, na
    !
    READ(unit,*,iostat=ios) S%ntyp, S%nat, S%ibrav, S%celldm(1:6)
    IF(ios/=0) CALL errore(sub,"reading S%ntyp, S%nat, S%ibrav, S%celldm(1:6)", 1)
    !
    IF(S%ibrav==0)THEN
      READ(unit,*,iostat=ios) S%at(1:3,1:3)
      IF(ios/=0) CALL errore(sub,"reading S%at(1:3,1:3)", 1)
    ENDIF
    !
    ! generate at, bg, volume
    IF (S%ibrav /= 0) THEN
      CALL latgen(S%ibrav, S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega)
      S%at = S%at / S%celldm(1)  !  bring at from bohr to units of alat 
    ENDIF
    CALL volume(S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega)
    CALL recips(S%at(:,1), S%at(:,2), S%at(:,3), S%bg(:,1), S%bg(:,2), S%bg(:,3))
    !
    DO nt = 1, S%ntyp
      READ(unit,*,iostat=ios) dummy, S%atm(nt), S%amass(nt)
      IF(ios/=0) CALL errore(sub,"reading S%ityp(na), S%atm(nt), S%amass(nt)", nt)
    ENDDO
    !
    ALLOCATE(S%ityp(S%nat), S%tau(3,S%nat))
    DO na = 1, S%nat
      READ(unit,*,iostat=ios) dummy, S%ityp(na), S%tau(:,na)
      IF(ios/=0) CALL errore(sub,"reading S%ityp(na), S%atm(nt), S%amass(nt)", nt)
    ENDDO
    !
  END SUBROUTINE read_system
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! compute redundant but useful quantities
  SUBROUTINE aux_system(S)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    INTEGER :: i, na
    S%nat3 = 3*S%nat
    S%nat32 = S%nat3**2
    S%nat33 = S%nat3**3
    
    IF(allocated(S%sqrtmm1)) CALL errore("aux_system","should be called twice",1)
    ALLOCATE(S%sqrtmm1(S%nat3))
    DO i = 1,S%nat3
      na = (i-1)/3 +1
      S%sqrtmm1(i) = 1._dp/SQRT(S%amass(S%ityp(na)))
    ENDDO
    
  END SUBROUTINE aux_system
  ! <<^V^\\=========================================//-//-//========//O\\//
  SUBROUTINE read_fc2(filename, S, fc)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    !
    CHARACTER(len=8),PARAMETER :: sub = "read_fc2"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER :: na1,  na2,  j1,  j2, jn1, jn2, jn3
    INTEGER :: na1_, na2_, j1_, j2_
    INTEGER :: n_R, i
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    CALL read_system(unit, S)
    !
    READ(unit, *) jn1, jn2, jn3
    WRITE(stdout,*) "Original FC2 grid:", jn1, jn2, jn3
    !
    DO na1=1,S%nat
    DO na2=1,S%nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3            
          !
          READ(unit,*) j1_, j2_, na1_, na2_
          IF ( ANY((/na1,na2,j1,j2/) /= (/na1_,na2_,j1_,j2_/)) ) &
            CALL errore(sub,'not matching na1,na2,j1,j2',1)
          !
          READ(unit,*) n_R
          IF ( fc%n_R == 0) CALL allocate_fc2_grid(n_R, S%nat, fc)
          !
          DO i = 1, fc%n_R
            READ(unit,*) fc%yR(:,i), fc%FC(jn1,jn2,i)
            ! also, find the index of R=0
            IF( ALL(fc%yR(:,i)==0) ) fc%i_0 = i
          ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    !
    CLOSE(unit)
    ! Compute list of R in cartesian coords
    fc%xR = REAL(fc%yR, kind=DP)
    CALL cryst_to_cart(n_R,fc%xR,S%at,1)
    !
  END SUBROUTINE read_fc2
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE allocate_fc2_grid(n_R, nat, fc)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n_R, nat
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    CHARACTER(len=16),PARAMETER :: sub = "allocate_fc2_grid"
    !
    IF(allocated(fc%yR) .or. allocated(fc%xR) .or. allocated(fc%FC)) &
      CALL errore(sub, 'some element is already allocated', 1)
    !
    ALLOCATE(fc%yR(3,n_R))
    ALLOCATE(fc%xR(3,n_R))
    ALLOCATE(fc%FC(3*nat,3*nat,n_R))
    fc%n_R = n_R
    !
  END SUBROUTINE allocate_fc2_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE deallocate_fc2_grid(fc)
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    CHARACTER(len=18),PARAMETER :: sub = "deallocate_fc2_grid"
    !
    IF(.NOT.(allocated(fc%yR) .and. allocated(fc%xR) .and. allocated(fc%FC))) &
      CALL errore(sub, 'some element is already deallocated', 1)
    !
    DEALLOCATE(fc%yR)
    DEALLOCATE(fc%xR)
    DEALLOCATE(fc%FC)
    fc%n_R = 0
    fc%i_0 = -1
    !
  END SUBROUTINE deallocate_fc2_grid
  ! <<^V^\\=========================================//-//-//========//O\\//
  SUBROUTINE read_fc3(filename, S, fc)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    TYPE(forceconst3_grid),INTENT(inout) :: fc
    !
    CHARACTER(len=8),PARAMETER :: sub = "read_fc3"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
    INTEGER :: na1_, na2_, na3_, j1_, j2_, j3_
    INTEGER :: n_R, i
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    CALL read_system(unit, S)
    !
    READ(unit, *) jn1, jn2, jn3
    WRITE(stdout,*) "Original FC3 grid:", jn1, jn2, jn3
    !
    DO na1=1,S%nat
    DO na2=1,S%nat 
    DO na3=1,S%nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3
      DO j3=1,3     
      jn3 = j3 + (na3-1)*3
          !
          READ(unit,*) j1_, j2_, j3_, na1_, na2_, na3_
          IF ( ANY((/na1,na2,na3,j1,j2,j3/) /= (/na1_,na2_,na3_,j1_,j2_,j3_/)) ) &
            CALL errore(sub,'not matching na1,na2,na3,j1,j2,j3',1)
          !
          READ(unit,*) n_R
          IF ( fc%n_R == 0) CALL allocate_fc3_grid(n_R, S%nat, fc)
          !
          DO i = 1, fc%n_R
            READ(unit,*) fc%yR2(:,i), fc%yR3(:,i), fc%FC(jn1,jn2,jn3,i)
            ! also, find the index of R=0
            IF( ALL(fc%yR2(:,i)==0) .and. ALL(fc%yR3(:,i)==0) ) fc%i_00 = i
          ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    !
    CLOSE(unit)
    ! Prepare the lists of R in cartesian coords
    fc%xR2 = REAL(fc%yR2, kind=DP)
    CALL cryst_to_cart(n_R,fc%xR2,S%at,1)
    fc%xR3 = REAL(fc%yR3, kind=DP)
    CALL cryst_to_cart(n_R,fc%xR3,S%at,1)
    !
  END SUBROUTINE read_fc3
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE allocate_fc3_grid(n_R, nat, fc)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n_R, nat
    TYPE(forceconst3_grid),INTENT(inout) :: fc
    CHARACTER(len=16),PARAMETER :: sub = "allocate_fc3_grid"
    !
    IF(allocated(fc%yR2) .or. allocated(fc%xR2) .or. allocated(fc%FC) .or. &
       allocated(fc%yR3) .or. allocated(fc%xR3) ) &
      CALL errore(sub, 'some element is already allocated', 1)
    !
    ALLOCATE(fc%yR2(3,n_R), fc%yR3(3,n_R))
    ALLOCATE(fc%xR2(3,n_R), fc%xR3(3,n_R))
    ALLOCATE(fc%FC(3*nat,3*nat,3*nat,n_R))
    fc%n_R = n_R
    !
  END SUBROUTINE allocate_fc3_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE deallocate_fc3_grid(fc)
    IMPLICIT NONE
    TYPE(forceconst3_grid),INTENT(inout) :: fc
    CHARACTER(len=18),PARAMETER :: sub = "deallocate_fc3_grid"
    !
    IF(.NOT.(allocated(fc%yR2) .and. allocated(fc%xR2) .and. allocated(fc%FC) .and. &
             allocated(fc%yR3) .and. allocated(fc%xR3) )) &
      CALL errore(sub, 'some element is already deallocated', 1)
    !
    DEALLOCATE(fc%yR2, fc%yR3)
    DEALLOCATE(fc%xR2, fc%yR3)
    DEALLOCATE(fc%FC)
    fc%n_R = 0
    fc%i_00 = -1
    !
  END SUBROUTINE deallocate_fc3_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE div_mass_fc2 (S,fc)
    USE kinds, only : DP
    IMPLICIT NONE
    TYPE(forceconst2_grid) :: fc
    TYPE(ph_system_info)   :: S
    !
    INTEGER :: i, j, i_R
    !
    IF(.not.ALLOCATED(S%sqrtmm1)) &
      call errore('div_mass_fc2', 'missing sqrtmm1, call aux_system first', 1)
    
    DO i_R = 1, fc%n_R
      DO j = 1, S%nat3
      DO i = 1, S%nat3
        fc%FC(i, j, i_R) = fc%FC(i, j, i_R) * S%sqrtmm1(i)*S%sqrtmm1(j)
      ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE div_mass_fc2
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE div_mass_fc3 (S,fc)
    USE kinds, only : DP
    IMPLICIT NONE
    TYPE(forceconst3_grid) :: fc
    TYPE(ph_system_info)   :: S
    !
    INTEGER :: i, j, k, i_R
    !
    IF(.not.ALLOCATED(S%sqrtmm1)) &
      call errore('div_mass_fc2', 'missing sqrtmm1, call aux_system first', 1)
    
    DO i_R = 1, fc%n_R
      DO k = 1, S%nat3
      DO j = 1, S%nat3
      DO i = 1, S%nat3
        fc%FC(i, j, k, i_R) = fc%FC(i, j, k, i_R) &
                    * S%sqrtmm1(i)*S%sqrtmm1(j)*S%sqrtmm1(k)
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE div_mass_fc3
  ! \/o\________\\\_________________________________________/^>
  
END MODULE input_fc







