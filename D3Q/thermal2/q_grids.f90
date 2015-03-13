!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE q_grids

  USE kinds,     ONLY : DP
  
  TYPE q_grid
    CHARACTER(len=9) :: basis = ''
    INTEGER :: n(3) = -1
    INTEGER :: nq = 0
    REAL(DP),ALLOCATABLE :: xq(:,:)
  END TYPE q_grid
  !
  TYPE  q_basis
    INTEGER :: nbnd  ! number of phonon bands
    INTEGER :: nconf ! number of temperatures
    REAL(DP),ALLOCATABLE :: w(:,:)    ! phonon frequencies
    REAL(DP),ALLOCATABLE :: c(:,:,:)  ! phonon group velocity, \nabla_q w
    REAL(DP),ALLOCATABLE :: be(:,:,:) ! bose-einstein occupation (one per configuration)
    !REAL(DP),ALLOCATABLE :: b(:,:,:)  ! b vector, as in PRB 88, 045430 (2013)
    TYPE(q_grid),POINTER :: grid
    CONTAINS
      procedure :: B  => B_right_hand_side
  END TYPE q_basis
  !
  CONTAINS
!   ! \/o\________\\\_________________________________________/^>
!   !
!   ! Nasty subroutine that sets some global variables of QE.
!
!       NOTE: all the at, bg, nat, tau and so on must be set in 
!             the QE global modules for this to work
!
!   SUBROUTINE setup_symmetry(S)
!     USE input_fc,  ONLY : ph_system_info 
!     USE symm_base, ONLY : nrot, nsym, set_sym_bl, find_sym
!     IMPLICIT NONE
!     TYPE(ph_system_info),INTENT(in)   :: S ! = System
!     REAL(DP) :: m_loc(3,S%nat) ! fake starting magnetisation
!     m_loc = 0._dp
!     ! ######################### symmetry setup #########################
!     ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~ 
!     CALL set_sym_bl ( )
!     WRITE(*, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
!     !
!     ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
!     CALL find_sym ( S%nat, S%tau, S%ityp, 6,6,6, .false., m_loc )
!     WRITE(*, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
!   END SUBROUTINE setup_symmetry
!   ! \/o\________\\\_________________________________________/^>
  SUBROUTINE setup_simple_grid(S, n1,n2,n3, grid, xq0)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    !
    INTEGER :: i,j,k, idx
    grid%n(1) = n1
    grid%n(2) = n2
    grid%n(3) = n3
    grid%nq = n1*n2*n3
    !
    IF(allocated(grid%xq)) CALL errore("setup_simple_grid", "grid is already allocated", 1)
    ALLOCATE(grid%xq(3,grid%nq))
    !
    idx = 0
    DO i = 0, n1-1
    DO j = 0, n2-1
    DO k = 0, n3-1
      !
      idx = idx+1
      grid%xq(1,idx) = REAL(i,kind=DP)/REAL(n1,kind=DP)
      grid%xq(2,idx) = REAL(j,kind=DP)/REAL(n2,kind=DP)
      grid%xq(3,idx) = REAL(k,kind=DP)/REAL(n3,kind=DP)
      !
    ENDDO
    ENDDO
    ENDDO
    !
    CALL cryst_to_cart(grid%nq,grid%xq,S%bg, +1)
    grid%basis = 'cartesian'
    !
    IF(present(xq0)) THEN
      DO idx = 1,grid%nq
        grid%xq(:,idx) = grid%xq(:,idx) + xq0
      ENDDO
    ENDIF
    !
  END SUBROUTINE setup_simple_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE setup_path(xqi, xqf, nq, path)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nq
    REAL(DP),INTENT(in) :: xqi(3), xqf(3)
    TYPE(q_grid),INTENT(inout) :: path
    !
    INTEGER :: i, n0
    REAL(DP),PARAMETER:: eps = 1.d-8
    REAL(DP),ALLOCATABLE :: auxq(:,:)
    REAL(DP) :: dq(3)
    !
    IF(nq==0) RETURN
    !
    IF(path%nq >0)THEN
      ALLOCATE(auxq(3,path%nq))
      auxq = path%xq
      DEALLOCATE(path%xq)
      ! Correctly chain up series of points
      IF(SUM(ABS(auxq(:,path%nq)-xqi))<eps) THEN
        n0 = path%nq 
        path%nq = path%nq + nq -1
      ELSE
        n0 = path%nq +1
        path%nq = path%nq + nq
      ENDIF
      !
      ALLOCATE(path%xq(3,path%nq))
      path%xq(:,1:size(auxq,2)) = auxq(:,1:size(auxq,2))
      DEALLOCATE(auxq)
    ELSE
      n0 = 1
      path%nq = nq
      ALLOCATE(path%xq(3,path%nq))
    ENDIF
    !
    IF(nq>1)THEN
      dq = (xqf-xqi)/(nq-1)
      DO i = n0, path%nq
        path%xq(:,i) = xqi + dq * (i-n0)
      ENDDO
    ELSE
        ! With nq=1 xqf is appended to the list, xqi is ignored
        ! note that repeated points are dropped
        path%xq(:,path%nq) = xqf
    ENDIF
    !
  END SUBROUTINE setup_path
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE prepare_q_basis(qgrid, qbasis, nconf, T, S, fc2)
    USE fc2_interpolate,    ONLY : freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity_proj
    USE input_fc,           ONLY : ph_system_info, forceconst2_grid
    IMPLICIT NONE
    TYPE(q_grid),INTENT(in),TARGET  :: qgrid
    TYPE(q_basis),INTENT(out) :: qbasis
    TYPE(ph_system_info)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    INTEGER             :: nconf
    REAL(DP),INTENT(in) :: T(nconf)
    !
    INTEGER :: iq, it
    COMPLEX(DP),ALLOCATABLE :: U(:,:)
    !
    IF(qgrid%nq<=0) CALL errore("prepare_q_basis", "grid must be set up first", 1)
    !
    qbasis%nbnd = S%nat3
    qbasis%nconf = nconf
    qbasis%grid => qgrid
    !
    ALLOCATE(qbasis%w(qbasis%nbnd, qgrid%nq))
    ALLOCATE(qbasis%c(3, qbasis%nbnd, qgrid%nq))
    ALLOCATE(qbasis%be(qbasis%nbnd, qbasis%nconf, qgrid%nq))
    !
    ALLOCATE(U(qbasis%nbnd,qbasis%nbnd))
    
    DO iq = 1,qgrid%nq
      qbasis%c(:,:,iq) = velocity_proj(S, fc2, qgrid%xq(:,iq))
      CALL  freq_phq_safe(qgrid%xq(:,iq), S, fc2, qbasis%w(:,iq), U)
      DO it = 1, qbasis%nconf
        CALL  bose_phq(T(it), qbasis%nbnd, qbasis%w(:,iq), qbasis%be(:,it,iq))
      ENDDO
    ENDDO
    !
    DEALLOCATE(U)
    
  END SUBROUTINE prepare_q_basis
  ! \/o\________\\\_________________________________________/^>
  REAL(DP) FUNCTION B_right_hand_side(self, ix, it, nu, iq) &
           RESULT (B)
    IMPLICIT NONE
    CLASS(q_basis),INTENT(in) :: self
    INTEGER :: ix, nu, iq, it
    !
    B = -self%c(ix,nu,iq) * self%w(nu,iq) * self%be(it,nu,iq) * (self%be(it,nu,iq)+1)
    !
  END FUNCTION B_right_hand_side

  !
END MODULE q_grids




