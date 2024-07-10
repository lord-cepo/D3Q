!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! References:
! [1] Calandra, Lazzeri, Mauri : Physica C 456 (2007) 38-44
! [2] Li. et.al. PhysRevLett.112.175501
#define timer_CALL CALL

MODULE linewidth

#include "mpi_thermal.h"
  USE kinds,       ONLY : DP
  USE mpi_thermal, ONLY : my_id, num_procs, mpi_bsum, allgather_mat, scatteri_tns, mpi_sum_mat
  USE constants,   ONLY : RY_TO_CMM1
  USE timers

CONTAINS

  FUNCTION linewidth_q_dense(xq0, nconf, T, S, grid, fc2, fc3, quality, freq1, U1, lw_UN)
    USE q_grids,          ONLY : q_grid, setup_grid
    USE constants,        ONLY : pi
    USE input_fc,         ONLY : ph_system_info
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE functions,        ONLY : refold_bz
    USE thtetra,          ONLY : tetra_init, tetra_weights_delta
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    !
    INTEGER, INTENT(IN) :: quality
    INTEGER :: quality3
    INTEGER :: idx_, i_, j_, k_
    REAL(DP)  :: ni(3), xq0_cryst(3), big_ni(3), center_xq(3,quality**3)
    REAL(DP),OPTIONAL,INTENT(out)   :: lw_UN(S%nat3,2,nconf) ! Normal/Umklapp contribution
    REAL(DP),OPTIONAL,INTENT(in)    :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    COMPLEX(DP) :: dummy(S%nat3,S%nat3)
    REAL(DP) :: wgout(S%nat3**2,S%nat3,2)
    ! FUNCTION RESULT:
    REAL(DP) :: linewidth_q_dense(S%nat3,nconf)
    REAL(DP) :: lw(S%nat3,nconf) !aux

    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, it, nu0(3), j_un, ibnd, jbnd, index_double, i
    INTEGER,PARAMETER :: normal=1, umklapp=2
    INTEGER, PARAMETER :: X=1, C=2
    ! REAL(DP) :: threshold, xq_cryst(3)
    !> Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=5000
    ! INTEGER :: nrws
    REAL(DP) :: atws(3,3) ! supercell of size nq1 x nq2 x nq3
    ! REAL(DP) :: wgin(grid%nqtot*quality**3), rws(0:3,nrwsx)
    ! REAL(DP), EXTERNAL :: wsweight
    ! EXTERNAL wsinit
    REAL(DP) :: freqs_doubled_X(S%nat3**2, quality**3), freqs_doubled_C(S%nat3**2, quality**3)
    ! REAL(DP), ALLOCATABLE :: gather_doubled_X(:,:), gather_doubled_C(:,:)
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3), xq_aux(3), aux(S%nat3)
    LOGICAL :: X_flag, C_flag
    ! REAL(DP) :: bose_doubled(S%nat3**2, grid%nq,2,nconf)
    !
    quality3 = quality**3
    !> WS cell initialization: atws describe the grid%xq seen as lattice points
    ! CALL setup_grid(grid%type, S%bg, grid%n(1)*quality, grid%n(2)*quality, grid%n(3)*quality, &
    ! dense_grid, xq0=grid%xq0, scatter=grid%scattered, quiet=.true.)
    DO i = 1,3
      atws(:,i) = S%bg(:,i)/grid%n(i)
    END DO
    CALL tetra_init( (/quality,quality,quality/), atws, .false., .false.)
    !> initialize WS r-vectors
    ! CALL wsinit(rws,nrwsx,nrws,atws)

    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    lw = 0._dp
    IF(present(lw_UN)) lw_UN = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    timer_CALL t_freq%start()
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    IF(present(freq1) .and. present(U1)) THEN
      freq(:,1) = freq1
      U(:,:,1)    = U1
    ELSE
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    ENDIF
    timer_CALL t_freq%stop()

    DO iq = 1, grid%nq
      !
      timer_CALL t_freqd%start()
      xq0_cryst = grid%xq(:,iq)
      CALL cryst_to_cart(1, xq0_cryst, S%at, -1)
      big_ni = REAL(grid%n*quality, DP)
      idx_ = 0

      DO i_ = 0, quality-1
        DO j_ = 0, quality-1
          DO k_ = 0, quality-1
            !
            idx_ = idx_ + 1
            ni = REAL((/i_, j_, k_/), DP)
            center_xq(:,idx_) = (ni + (1-quality)/2)/big_ni + xq0_cryst
            !
          ENDDO
        ENDDO
      ENDDO
      CALL cryst_to_cart(quality3, center_xq, S%bg, 1)
      DO jq = 1, quality3
        CALL freq_phq_safe(center_xq(:,jq), S, fc2, freq(:,2), dummy)
        ! the third vector (why it's minus?)
        xq(:,3) = -(center_xq(:,jq)+xq(:,1))
        CALL freq_phq_safe(xq(:,3), S, fc2, freq(:,3), dummy)
        DO ibnd = 1, S%nat3
          DO jbnd = 1, S%nat3
            index_double = S%nat3*(ibnd-1) + jbnd
            freqs_doubled_X(index_double,jq) = freq(ibnd,2) + freq(jbnd,3)
            freqs_doubled_C(index_double,jq) = freq(jbnd,3) - freq(ibnd,2)
          ENDDO
        ENDDO
      ENDDO
      timer_CALL t_freqd%stop()
      ! PRINT*, "freqs_doubled_X", freqs_doubled_X(25,1)
      X_flag = .false.
      C_flag = .true.
      wgout = 0._dp
      timer_CALL t_thtetra%start()
      DO ibnd = 1, S%nat3
        IF (X_flag .or. ANY(freqs_doubled_X < freq(ibnd,1))) THEN
          X_flag = .true.
          wgout(:,ibnd,X) = SUM(tetra_weights_delta(quality**3, S%nat3**2, freqs_doubled_X, freq(ibnd,1)), dim=2)
          ! PRINT*, "X", SUM(wgout(:,ibnd,X))
        ENDIF

        IF(ANY(C_flag .and. freqs_doubled_C - freq(ibnd,1) < 0)) THEN
          C_flag = .true.
          wgout(:,ibnd,C) = SUM(tetra_weights_delta(quality**3, S%nat3**2, freqs_doubled_C, freq(ibnd,1)), dim=2)
          ! PRINT*, "C", SUM(wgout(:,ibnd,C))
        ENDIF

      END DO
      ! PRINT*, "-----------"
      timer_CALL t_thtetra%stop()

      timer_CALL t_freq%start()
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
      IF(present(lw_UN)) THEN
        xq_aux = refold_bz(xq(:,1), S%bg) &
          +refold_bz(xq(:,2), S%bg) &
          +refold_bz(xq(:,3), S%bg)
        IF(ALL(ABS(xq_aux)<1.d-6)) THEN
          j_un = normal
        ELSE
          j_un = umklapp
        ENDIF
      ENDIF

!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      timer_CALL t_freq%stop()
      !
      ! IF (SUM(weights_C(:,:,iq)) + SUM(weights_C(:,:,iq)) < threshold) THEN
      !    skip = skip + 1
      ! ENDIF
      timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      timer_CALL t_fc3int%stop()
      timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      timer_CALL t_fc3rot%stop()
      timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      timer_CALL t_fc3m2%stop()
      !
      DO it = 1,nconf
        timer_CALL t_bose%start()
        ! Compute bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 2,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        timer_CALL t_bose%stop()
        timer_CALL t_sum%start()
        aux = sum_linewidth_modes_dense(S, freq, bose, V3sq, nu0, wgout)
        lw(:,it) = lw(:,it) + aux
        IF(present(lw_UN)) lw_UN(:,j_un,it) = lw_UN(:,j_un,it) + aux
        timer_CALL t_sum%stop()
      ENDDO ! nconf
      !
    ENDDO ! grid%nq
    ! open(unit=10, file='points.txt', status='replace', action='write')
    ! CALL cryst_to_cart(grid%nq, grid%xq, S%at, -1)
    ! WRITE(10,*) grid%xq
    ! CALL cryst_to_cart(grid%nq, grid%xq, S%bg, 1)
    ! close(10)

    ! ioWRITE(*,*) "numero di skippati", skip
    !
    timer_CALL t_mpicom%start()
    IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf,lw)
    IF(grid%scattered .and. present(lw_UN)) CALL mpi_bsum(S%nat3,2,nconf,lw_UN)

    timer_CALL t_mpicom%stop()
    linewidth_q_dense = lw*pi/2._dp/grid%nqtot
    !
    DEALLOCATE(U, V3sq, D3)
    !
  END FUNCTION linewidth_q_dense

! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION linewidth_q_tetra(xq0, nconf, T, S, grid, fc2, fc3, freq1, U1, lw_UN)
    USE q_grids,          ONLY : q_grid
    USE constants,        ONLY : pi
    USE input_fc,         ONLY : ph_system_info
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE functions,        ONLY : refold_bz
    USE thtetra,          ONLY : tetra_init, tetra_weights_delta
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    !
    REAL(DP),OPTIONAL,INTENT(out)   :: lw_UN(S%nat3,2,nconf) ! Normal/Umklapp contribution
    REAL(DP),OPTIONAL,INTENT(in)    :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    COMPLEX(DP) :: dummy(S%nat3,S%nat3)
    REAL(DP), ALLOCATABLE :: weights_C(:,:,:), weights_X(:,:,:)
    ! 1 is the folded index which hosts freq(ibnd,2) +/- freq(jbnd,3)
    ! 2 is the q point index
    ! 3 is the "fermi energy" of the tetrahedron method, equals to freq1(ibnd) here,
    ! FUNCTION RESULT:
    REAL(DP) :: linewidth_q_tetra(S%nat3,nconf)
    REAL(DP) :: lw(S%nat3,nconf) !aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, it, nu0(3), j_un, ibnd, jbnd, index_double
    INTEGER,PARAMETER :: normal=1, umklapp=2
    ! REAL(DP) ::  xq_cryst(3)
    REAL(DP) :: freqs_doubled_X(S%nat3**2, grid%nq), freqs_doubled_C(S%nat3**2, grid%nq)
    REAL(DP), ALLOCATABLE :: gather_doubled_X(:,:), gather_doubled_C(:,:)
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3), xq_aux(3), aux(S%nat3)
    ! REAL(DP) :: bose_doubled(S%nat3**2, grid%nq,2,nconf)
    !
    CALL tetra_init( grid%n, S%bg, .true., grid%scattered)

    ALLOCATE(weights_X(S%nat3, S%nat3**2, grid%nqtot))
    ALLOCATE(weights_C(S%nat3, S%nat3**2, grid%nqtot))
    ! ioWRITE(*,*) "ntetra", ntetra, "grid%nq", grid%nq, "grid%nqtot", grid%nqtot, "grid%w(1)", grid%w(1)
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    lw = 0._dp
    IF(present(lw_UN)) lw_UN = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    timer_CALL t_freq%start()
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    IF(present(freq1) .and. present(U1)) THEN
      freq(:,1) = freq1
      U(:,:,1)    = U1
    ELSE
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    ENDIF
    timer_CALL t_freq%stop()
    !
    !WRITE(*,*) "Summing over a grid of", grid%nq, " points"

    ! this cycle initializes freqs, as a grid of frequencies
    ! this cycle initializes freqs_doubled, which is freq(ibnd,2) +/- freq(jbnd,3)
    timer_CALL t_freqd%start()
    DO iq = 1, grid%nq
      CALL freq_phq_safe(grid%xq(:,iq), S, fc2, freq(:,2), dummy)
      ! the third vector (why it's minus?)
      xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
      CALL freq_phq_safe(xq(:,3), S, fc2, freq(:,3), dummy)
      DO ibnd = 1, S%nat3
        DO jbnd = 1, S%nat3
          index_double = S%nat3*(ibnd-1) + jbnd
          freqs_doubled_X(index_double,iq) = freq(ibnd,2) + freq(jbnd,3)
          freqs_doubled_C(index_double,iq) = freq(jbnd,3) - freq(ibnd,2)
        ENDDO
      ENDDO
    ENDDO

    IF(grid%scattered) THEN
      CALL allgather_mat(S%nat3**2, grid%nq, freqs_doubled_X, gather_doubled_X)
      CALL allgather_mat(S%nat3**2, grid%nq, freqs_doubled_C, gather_doubled_C)
    ELSE
      gather_doubled_X = freqs_doubled_X
      gather_doubled_C = freqs_doubled_C
    ENDIF

    timer_CALL t_freqd%stop()
    timer_CALL t_thtetra%start()
    DO ibnd = 1, S%nat3
      weights_X(ibnd,:,:) = tetra_weights_delta(grid%nqtot, S%nat3**2, gather_doubled_X, freq(ibnd,1))
      weights_C(ibnd,:,:) = tetra_weights_delta(grid%nqtot, S%nat3**2, gather_doubled_C, freq(ibnd,1))
    ENDDO
    IF(grid%scattered) THEN
      DO iq = 1, grid%nqtot
        CALL mpi_sum_mat(S%nat3, S%nat3**2, weights_X(:,:,iq)) ! MPI_SUM MPI_REDUCE
        CALL mpi_sum_mat(S%nat3, S%nat3**2, weights_C(:,:,iq))
      ENDDO
      CALL scatteri_tns(S%nat3, S%nat3**2, grid%nqtot, weights_X)
      CALL scatteri_tns(S%nat3, S%nat3**2, grid%nqtot, weights_C)
    ENDIF

    timer_CALL t_thtetra%stop()
    ! TODO: qui devo scatterare la griglia e i pesi
    ! threshold = 1d-1*SUM(summ(:))/grid%nq
    ! threshold = (SUM(weights_C) + SUM(weights_X))/(grid%nq)
    ! skip = 0
    DO iq = 1, grid%nq
      !
      timer_CALL t_freq%start()
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
      IF(present(lw_UN)) THEN
        xq_aux = refold_bz(xq(:,1), S%bg) &
          +refold_bz(xq(:,2), S%bg) &
          +refold_bz(xq(:,3), S%bg)
        IF(ALL(ABS(xq_aux)<1.d-6)) THEN
          j_un = normal
        ELSE
          j_un = umklapp
        ENDIF
      ENDIF

!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      timer_CALL t_freq%stop()
      !
      ! IF (SUM(weights_C(:,:,iq)) + SUM(weights_C(:,:,iq)) < threshold) THEN
      !    skip = skip + 1
      ! ENDIF
      timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      timer_CALL t_fc3int%stop()
      timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      timer_CALL t_fc3rot%stop()
      timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      timer_CALL t_fc3m2%stop()
      !
      DO it = 1,nconf
        timer_CALL t_bose%start()
        ! Compute bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 2,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        timer_CALL t_bose%stop()
        timer_CALL t_sum%start()
        aux = sum_linewidth_modes_weights(S, freq, bose, V3sq, nu0, weights_X(:,:,iq), weights_C(:,:,iq))
        lw(:,it) = lw(:,it) + aux
        IF(present(lw_UN)) lw_UN(:,j_un,it) = lw_UN(:,j_un,it) + aux
        timer_CALL t_sum%stop()
      ENDDO
      !
    ENDDO
    ! open(unit=10, file='points.txt', status='replace', action='write')
    ! CALL cryst_to_cart(grid%nq, grid%xq, S%at, -1)
    ! WRITE(10,*) grid%xq
    ! CALL cryst_to_cart(grid%nq, grid%xq, S%bg, 1)
    ! close(10)

    ! ioWRITE(*,*) "numero di skippati", skip
    !
    timer_CALL t_mpicom%start()
    IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf,lw)
    IF(grid%scattered .and. present(lw_UN)) CALL mpi_bsum(S%nat3,2,nconf,lw_UN)

    timer_CALL t_mpicom%stop()
    linewidth_q_tetra = lw*pi/2._dp
    !
    DEALLOCATE(U, V3sq, D3)
    !
  END FUNCTION linewidth_q_tetra

! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION linewidth_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, freq1, U1, lw_UN)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE functions,        ONLY : refold_bz
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    REAL(DP),OPTIONAL,INTENT(out)   :: lw_UN(S%nat3,2,nconf) ! Normal/Umklapp contribution
    REAL(DP),OPTIONAL,INTENT(in)    :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! FUNCTION RESULT:
    REAL(DP) :: linewidth_q(S%nat3,nconf)
    REAL(DP) :: lw(S%nat3,nconf) !aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, it, nu0(3), j_un
    INTEGER,PARAMETER :: normal=1, umklapp=2
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3), xq_aux(3), aux(S%nat3)
    !

    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    lw = 0._dp
    IF(present(lw_UN)) lw_UN = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    timer_CALL t_freq%start()
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    IF(present(freq1) .and. present(U1)) THEN
      freq(:,1) = freq1
      U(:,:,1)    = U1
    ELSE
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    ENDIF
    timer_CALL t_freq%stop()
    !
    !WRITE(*,*) "Summing over a grid of", grid%nq, " points"
    DO iq = 1, grid%nq
      !
      timer_CALL t_freq%start()
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
      IF(present(lw_UN)) THEN
        xq_aux = refold_bz(xq(:,1), S%bg) &
          +refold_bz(xq(:,2), S%bg) &
          +refold_bz(xq(:,3), S%bg)
        IF(ALL(ABS(xq_aux)<1.d-6)) THEN
          j_un = normal
        ELSE
          j_un = umklapp
        ENDIF
      ENDIF

      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
      !$OMP END PARALLEL DO
      timer_CALL t_freq%stop()
      !
      timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      timer_CALL t_fc3int%stop()
      timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      timer_CALL t_fc3rot%stop()
      timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      timer_CALL t_fc3m2%stop()
      !
      DO it = 1,nconf
        timer_CALL t_bose%start()
        ! Compute bose-einstein occupation at q2 and q3
        !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
        !$OMP END PARALLEL DO
        timer_CALL t_bose%stop()
        timer_CALL t_sum%start()
        aux = -0.5_dp * grid%w(iq)*sum_linewidth_modes(S, sigma(it), freq, bose, V3sq, nu0)
        lw(:,it) = lw(:,it) + aux
        IF(present(lw_UN)) lw_UN(:,j_un,it) = lw_UN(:,j_un,it) + aux
        timer_CALL t_sum%stop()
      ENDDO
      !
    ENDDO
    !
    timer_CALL t_mpicom%start()
    IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf,lw)
    IF(grid%scattered .and. present(lw_UN)) CALL mpi_bsum(S%nat3,2,nconf,lw_UN)

    timer_CALL t_mpicom%stop()
    linewidth_q = lw
    !
    DEALLOCATE(U, V3sq, D3)
    !
  END FUNCTION linewidth_q

  ! <<^V^\\=========================================//-//-//========//O\\//
  ! This function returns the complex self energy from the bubble(?) diagram:
  !  its real part is the order contribution to lineshift
  !  its complex part is the intrinsic linewidth
  FUNCTION selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, freq1, U1)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! FUNCTION RESULT:
    COMPLEX(DP) :: selfnrg_q(S%nat3,nconf)
    COMPLEX(DP) :: se(S%nat3,nconf) ! aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, it
    INTEGER :: nu0(3)
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    se = (0._dp, 0._dp)
    !
    timer_CALL t_freq%start()
    ! Compute eigenvalues, eigenmodes at q1
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    IF(present(freq1) .and. present(U1)) THEN
      freq(:,1) = freq1
      U(:,:,1)    = U1
    ELSE
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    ENDIF
    timer_CALL t_freq%stop()

    DO iq = 1, grid%nq
      !
      timer_CALL t_freq%start()
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      timer_CALL t_freq%stop()
      !
      timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      timer_CALL t_fc3int%stop()
      timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      timer_CALL t_fc3rot%stop()
      timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      timer_CALL t_fc3m2%stop()
      !
      DO it = 1,nconf
        ! Compute bose-einstein occupation at q2 and q3
        timer_CALL t_bose%start()
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        timer_CALL t_bose%stop()
        !
        timer_CALL t_sum%start()
        se(:,it) = se(:,it) + grid%w(iq)*sum_selfnrg_modes( S, sigma(it), T(it), freq, bose, V3sq, nu0 )
        timer_CALL t_sum%stop()
        !
      ENDDO
      !
    ENDDO
    !
    timer_CALL t_mpicom%start()
    IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf,se)
    timer_CALL t_mpicom%stop()
    selfnrg_q = -0.5_dp * se
    !
    DEALLOCATE(U, V3sq)
    !
  END FUNCTION selfnrg_q
  ! \/o\________\\\_________________________________________/^>
  ! Simple spectral function, computed as a superposition of Lorentzian functions
  FUNCTION simple_spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1, shift) &
    RESULT(spectralf)
    USE q_grids,            ONLY : q_grid
    USE constants,          ONLY : RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe
    USE fc3_interpolate,    ONLY : forceconst3

    IMPLICIT NONE
    ! ARGUMENTS:
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    LOGICAL,INTENT(in)  :: shift
    !
    REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! FUNCTION RESULT:
    REAL(DP) :: spectralf(ne,S%nat3,nconf)
    !
    COMPLEX(DP) :: selfnrg(S%nat3,nconf)
    !
    INTEGER :: it, i, ie
    REAL(DP) :: gamma(S%nat3,nconf), delta(S%nat3,nconf), omega, denom
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    REAL(DP) :: freq(S%nat3)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    IF(present(freq1) .and. present(U1)) THEN
      freq = freq1
      U    = U1
    ELSE
      !CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      CALL freq_phq_safe(xq0, S, fc2, freq, U)
    ENDIF
    !
    ! use lineshift to compute 3rd order linewidth and lineshift
    IF(shift)THEN
      selfnrg = selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
      gamma(:,:) = -DIMAG(selfnrg(:,:))
      delta(:,:) =   DBLE(selfnrg(:,:))
    ELSE
      gamma = linewidth_q_tetra(xq0, nconf, T, S, grid, fc2, fc3)
      delta = 0._dp
    ENDIF
    !
    DO it = 1,nconf
      WRITE(*,'(30x,6e12.2,5x,6e12.2)') gamma(:,it)*RY_TO_CMM1,delta(:,it)*RY_TO_CMM1
    ENDDO
    !
    delta = 0._dp
    ! Compute and superpose the Lorentzian functions corresponding to each band
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1, ne
          ! Only apply shift (i.e. real part of self nrg) if requested
          omega = freq(i)
          denom =   (ener(ie)**2 -omega**2 -2*omega*delta(i,it))**2 &
            + 4*omega**2 *gamma(i,it)**2
          IF(ABS(denom)/=0._dp)THEN
            spectralf(ie,i,it) = 2*omega*gamma(i,it) / denom
          ELSE
            spectralf(ie,i,it) = 0._dp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
  END FUNCTION simple_spectre_q


  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Self energy for all the phonon bands at q in a range of frequencies
  FUNCTION selfnrg_omega_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
    RESULT(selfnrg_wq)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)           :: grid
    REAL(DP),INTENT(in)               :: sigma(nconf) ! ry
    !
    REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! To interpolate D2 and D3:
    INTEGER :: iq, jq, it
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: nu0(3)
    COMPLEX(DP),ALLOCATABLE :: selfnrg(:,:,:)
    ! FUNCTION RESULT:
    COMPLEX(DP) :: selfnrg_wq(ne,S%nat3,nconf)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    ALLOCATE(selfnrg(ne,S%nat3,nconf))
    !
    selfnrg = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    timer_CALL t_freq%start()
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    IF(present(freq1) .and. present(U1)) THEN
      freq(:,1) = freq1
      U(:,:,1)    = U1
    ELSE
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    ENDIF
    timer_CALL t_freq%stop()
    !
    DO iq = 1, grid%nq
      !CALL print_percent_wall(33.333_dp, 300._dp, iq, grid%nq, (iq==1))
      !
      timer_CALL t_freq%start()
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      timer_CALL t_freq%stop()
      !
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      timer_CALL t_fc3int%stop()
      timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      timer_CALL t_fc3rot%stop()
      timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      timer_CALL t_fc3m2%stop()
      !
      DO it = 1,nconf
        ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
        timer_CALL t_bose%start()
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        timer_CALL t_bose%stop()
        timer_CALL t_sum%start()
        selfnrg(:,:,it) = selfnrg(:,:,it) + grid%w(iq)*sum_selfnrg_spectre( S, sigma(it), freq, bose, V3sq, ne, ener, nu0 )
        timer_CALL t_sum%stop()
        !
      ENDDO
      !
    ENDDO
    !
    timer_CALL t_mpicom%start()
    IF(grid%scattered) CALL mpi_bsum(ne,S%nat3,nconf,selfnrg)
    timer_CALL t_mpicom%stop()
    selfnrg_wq = -0.5_dp * selfnrg
    !
    DEALLOCATE(U, V3sq, D3, selfnrg)
    !
  END FUNCTION selfnrg_omega_q
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Spectral weight function, computed as in eq. 1 of arXiv:1312.7467v1
  FUNCTION spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1, shift) &
    RESULT(spectralf)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    LOGICAL,INTENT(in)  :: shift ! set to false to drop the real part of the self-energy
    !
    REAL(DP),INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie, it
    REAL(DP) :: gamma, delta, omega, denom
    COMPLEX(DP) :: selfnrg(ne,S%nat3,nconf)
    ! FUNCTION RESULT:
    REAL(DP)    :: spectralf(ne,S%nat3,nconf)
    !
    ! Once we have the self-energy, the rest is trivial
    selfnrg = selfnrg_omega_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
    !
    timer_CALL t_mkspf%start()
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1, ne
          gamma =  -DIMAG(selfnrg(ie,i,it))
          IF(shift) THEN
            delta =   DBLE(selfnrg(ie,i,it))
          ELSE
            delta = 0._dp
          ENDIF
          omega = freq1(i)
          denom =   (ener(ie)**2 -omega**2 -2*omega*delta)**2 &
            + 4*omega**2 *gamma**2
          IF(ABS(denom)/=0._dp)THEN
            spectralf(ie,i,it) = 2*omega*gamma / denom
          ELSE
            spectralf(ie,i,it) = 0._dp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    timer_CALL t_mkspf%stop()
    !
  END FUNCTION spectre_q

  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Spectral weight function, computed as in eq. 1 of arXiv:1312.7467v1
  ! test, compute this as imaginary part of epsilon(omega)
  FUNCTION spectre2_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
    RESULT(spectralf)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !LOGICAL,INTENT(in)  :: shift ! set to false to drop the real part of the self-energy
    !
    REAL(DP),INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie, it
    COMPLEX(DP) :: tepsilon(ne,S%nat3,nconf)
    ! FUNCTION RESULT:
    REAL(DP)    :: spectralf(ne,S%nat3,nconf)
    !
    ! Once we have the self-energy, the rest is trivial
    tepsilon = tepsilon_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
    !
    timer_CALL t_mkspf%start()
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1,ne
          !IF(freq1(i)/=0._dp)THEN
          spectralf(ie,i,it) = DIMAG(tepsilon(ie,i,it))
          !ELSE
          !  spectralf(ie,i,it) = 0._dp
          !ENDIF
        ENDDO
      ENDDO
    ENDDO
    timer_CALL t_mkspf%stop()
    !
  END FUNCTION spectre2_q
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  !  Complex \tilde{epsilon} in a range of frequencies
  ! Experimental subroutine, assumes the material cubic and takes
  ! epsilon(1,1) as epsilon_infty
  FUNCTION tepsilon_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
    RESULT(tepsilon)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    USE constants,        ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    REAL(DP),INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie, it
    REAL(DP) :: pref, epsilon_infty, rmass, zeu_i2
    COMPLEX(DP) :: denom
    COMPLEX(DP):: selfnrg(ne,S%nat3,nconf)
    ! FUNCTION RESULT:
    COMPLEX(DP)    :: tepsilon(ne,S%nat3,nconf)
    !
    EXTERNAL errore
    IF(.not.S%lrigid) CALL errore("tepsilon","Cannot compute \tilde{epsilon} without epsilon0", 1)
    ioWRITE(*,*) "BEWARE: tepsilon is only isotropic"
    !
    ! Once we have the self-energy, the rest is trivial
    selfnrg = selfnrg_omega_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
    !
    ! S = 4piZ^2/(volume mass omega0^2)

    rmass = S%amass(1)*S%amass(2)/(S%amass(1)+S%amass(2))
    epsilon_infty = S%epsil(1,1)
    zeu_i2        = S%zeu(1,1,1)**2
    pref =  fpi * zeu_i2 /(S%omega * rmass)
    !pref = 1/freq1(i)**2 ! 4 pi Z^2  / (Volume mu  omega0^2) # mu = reduced mass
!     print*, "rmass", rmass, S%amass(1:2)
!     print*, "epsilon", epsilon_infty
!     print*, "zeu2", zeu_i2
!     print*, "pref", pref
!     print*, "vol", S%omega
!     print*, "denom", rmass*S%omega
!     print*, "fpi * zeu_i2", fpi * zeu_i2
!  denom   1002568.9328649712
!  epsilon   3.1410114605330000
!  fpi * zeu_i2   48.024054760187823
!  pref   4.7901000306236137E-005
!  rmass   8793.6751743770546        22152.652305024902        14582.196429874200
!  vol   114.01023041950073
!  zeu2   3.8216328511998792
!
    timer_CALL t_mkspf%start()
    tepsilon = DCMPLX(0._dp, 0._dp)
    DO it = 1,nconf
      DO i = 1,S%nat3
        IF(freq1(i)==0._dp) CYCLE
        DO ie = 1,ne
          denom = freq1(i)**2 - ener(ie)**2 - 2*freq1(i)*selfnrg(ie,i,it)
          IF(denom==DCMPLX(0._dp,0._dp)) CYCLE
          !tepsilon(ie,i,it) = 1 + pref*freq1(i)**2/denom
          tepsilon(ie,i,it) = epsilon_infty + pref/denom
        ENDDO
      ENDDO
    ENDDO
    timer_CALL t_mkspf%stop()
    !
  END FUNCTION tepsilon_q
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Infra red reflectivity |sqrt(epsilon)+1/sqrt(epsilon)-1|^2
  ! Experimental! Assumes isoropic material and a bunch of other stuff
  FUNCTION ir_reflectivity_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
    RESULT(reflectivity)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    REAL(DP),INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, it
    COMPLEX(DP) :: tepsilon(ne,S%nat3,nconf)
    COMPLEX(DP) :: aux(ne),aux2(ne)
    ! FUNCTION RESULT:
    REAL(DP)    :: reflectivity(ne,S%nat3,nconf)
    !
    ! We use tilde{epsilon}, which itself comes from the self-energy
    tepsilon = tepsilon_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
    !
    timer_CALL t_mkspf%start()
    DO it = 1,nconf
      DO i = 1,S%nat3
        aux  = SQRT(tepsilon(:,i,it))
        aux2 = (aux-1)/(aux+1)
        reflectivity(:,i,it) = DBLE(aux2*CONJG(aux2))
      ENDDO
    ENDDO
    timer_CALL t_mkspf%stop()
    !
  END FUNCTION ir_reflectivity_q
  !
  ! Sum the self energy at the provided ener(ne) input energies
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_selfnrg_spectre(S, sigma, freq, bose, V3sq, ne, ener, nu0)
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    USE functions,          ONLY : sigma_mgo
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma   ! smearing (regularization) (Ry)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
    !
    INTEGER,INTENT(in)  :: ne           ! number of energies...
    REAL(DP),INTENT(in) :: ener(ne)     ! energies for which to compute the spectral function
    INTEGER,INTENT(in)  :: nu0(3)       ! first non-zero phomnon frequency
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations
    REAL(DP) :: factor, sigma_, T=0._dp !freqtotm1
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    COMPLEX(DP) :: ctm_P,ctm_M, reg
    COMPLEX(DP) :: ctm(ne)
    REAL(DP)    :: freqm1(S%nat3,3)  ! phonon energies (Ry)
    !
    INTEGER :: i,j,k, ie
    !
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_selfnrg_spectre(ne,S%nat3)
    COMPLEX(DP),ALLOCATABLE :: spf(:,:)
    !
!     IF(sigma<=0._dp)THEN
!       CALL errore("sum_selfnrg_spectre","spf not implemented in the static limit. "&
!                   //"NEW: To do unshifted spf use 'spf imag'",1)
!     ENDIF
    !
    ALLOCATE(spf(ne,S%nat3))
    spf = (0._dp, 0._dp)
    !
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(ie,i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2, &
!$OMP                     ctm_P,ctm_M,ctm,reg,factor) &
!$OMP             REDUCTION(+: spf) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        bose_P   = 1 + bose(j,2) + bose(k,3)
        omega_P  = freq(j,2)+freq(k,3)
        omega_P2 = omega_P**2
        !
        bose_M   = bose(k,3) - bose(j,2)
        omega_M  = freq(j,2)-freq(k,3)
        omega_M2 = omega_M**2
        !
        ! A little optimization: precompute the parts that depends only on the energy
! ! !         DO ie = 1,ne
! ! !           ! regularization:
! ! !           reg = CMPLX(ener(ie), sigma, kind=DP)**2
! ! !           ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
! ! !           ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
! ! !           ctm(ie) = ctm_P + ctm_M
! ! !         ENDDO
        !
        DO i = 1,S%nat3
          !
          factor = V3sq(i,j,k) *freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
          !
          IF(ABS(sigma-666._dp)< 1._dp)THEN
            sigma_ = sigma_mgo(freq(i,1),T) &
              +sigma_mgo(freq(j,2),T) &
              +sigma_mgo(freq(k,3),T)
            !print*, RY_TO_CMM1*freq(i,1), RY_TO_CMM1*freq(j,2), RY_TO_CMM1*freq(k,3), RY_TO_CMM1*sigma_
          ELSE
            sigma_ = sigma
          ENDIF

          DO ie = 1, ne
            ! regularization:
            reg = CMPLX(ener(ie), sigma_, kind=DP)**2
            !
            ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
            ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
            ctm(ie) = ctm_P + ctm_M
            !
            spf(ie,i) = spf(ie,i) + ctm(ie) * factor !V3sq(i,j,k) * freqtotm1
          ENDDO
! ! !           spf(:,i) = spf(:,i) + ctm(:) * factor !V3sq(i,j,k) * freqtotm1

          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(ne, S%nat3, spf, freq(:,1))
    sum_selfnrg_spectre = spf
    DEALLOCATE(spf)
    !
  END FUNCTION sum_selfnrg_spectre
  !
  ! \/o\________\\\_________________________________________/^>
  ! Sum the self energy in a range of frequencies
  FUNCTION sum_selfnrg_modes(S, sigma, T, freq, bose, V3sq, nu0)
    USE input_fc,           ONLY : ph_system_info
    USE functions,          ONLY : df_bose
    USE merge_degenerate,   ONLY : merge_degen
    USE functions,          ONLY : sigma_mgo
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma, T
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations
    REAL(DP) :: freqtotm1
    REAL(DP) :: omega_P,  omega_M   ! \sigma\omega
    REAL(DP) :: omega_P2, omega_M2  ! \sigma\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg
    REAL(DP) :: freqm1(S%nat3,3), sigma_
    !
    INTEGER :: i,j,k
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_selfnrg_modes(S%nat3)
    COMPLEX(DP) :: se(S%nat3)
    !
    se(:) = (0._dp, 0._dp)
    !
    freqm1 = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_P,bose_M,omega_P,omega_M,&
!$OMP                     omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtotm1) &
!$OMP             REDUCTION(+: se) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        bose_P   = 1 + bose(j,2) + bose(k,3)
        bose_M   = bose(k,3) - bose(j,2)
        omega_P  = freq(j,2)+freq(k,3)
        omega_M  = freq(j,2)-freq(k,3)
        !
        IF(sigma>0._dp)THEN
          omega_P2 = omega_P**2
          omega_M2 = omega_M**2
        ELSE IF(sigma<0._dp)THEN ! (static limit)
          ctm_P = 2 * bose_P *omega_P/(omega_P**2+sigma**2)
          ctm_M = 2 * bose_M *omega_M/(omega_M**2+sigma**2)
        ELSE !IF (sigma==0._dp)THEN
          ! In the static limit with sigma=0 case we have to take the
          ! derivative of (n_3-n2)/(w_2-w_3) when w_2 is close to w_3
          IF(omega_P>0._dp)THEN
            ctm_P = 2 * bose_P /omega_P
          ELSE
            ctm_P = 0._dp
          ENDIF
          !
          IF(ABS(omega_M)>1.e-5_dp)THEN
            ctm_M = 2 * bose_M /omega_M
          ELSE
            IF(T>0._dp)THEN
              ctm_M = -2* df_bose(0.5_dp * omega_P, T)
            ELSE
              ctm_M = 0._dp
            ENDIF
          ENDIF
          !
        ENDIF
        !
        DO i = 1,S%nat3
          !
          !IF(freq(i,1)<1.d-3 .or. freq(j,2)<1.d-3 .or. freq(k,3)<1.d-3) CYCLE
          ! This comes from the definition of u_qj, Ref. 1. (there is an hidden factor 1/8)
          freqtotm1 = freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
          !
          ! regularization:
          IF(sigma>0._dp)THEN
            IF(ABS(sigma-666._dp)<1._dp)THEN
              sigma_ = sigma_mgo(freq(i,1),T) &
                +sigma_mgo(freq(j,2),T) &
                +sigma_mgo(freq(k,3),T)
              !print*, RY_TO_CMM1*freq(i,1), RY_TO_CMM1*freq(j,2), RY_TO_CMM1*freq(k,3), RY_TO_CMM1*sigma_
            ELSE
              sigma_ = sigma
            ENDIF

            reg = CMPLX(freq(i,1), sigma_, kind=DP)**2
            ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
            ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
          ENDIF
          !
          !
          se(i) = se(i) + (ctm_P + ctm_M)*freqtotm1 * V3sq(i,j,k)
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(S%nat3, se, freq(:,1))
    sum_selfnrg_modes = se
    !
  END FUNCTION sum_selfnrg_modes
  !
  ! Sum the Imaginary part of the self energy for the phonon modes, uses gaussian function for energy delta
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_linewidth_modes_weights(S, freq, bose, V3sq, nu0, weights_X, weights_C)
    USE constants, ONLY : pi
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    REAL(DP), INTENT(IN) :: weights_X(S%nat3,S%nat3**2), weights_C(S%nat3,S%nat3**2)
    !
    REAL(DP) :: sum_linewidth_modes_weights(S%nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_C, bose_X ! final/initial state populations
    REAL(DP) :: ctm_C, ctm_X   !
    REAL(DP) :: freqtotm1, freqtotm1_23
    REAL(DP) :: freqm1(S%nat3,3)
    !REAL(DP),SAVE :: leftover_e
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)!, sigma_
    lw(:) = 0._dp
    !
    freqm1 = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: lw) COLLAPSE(2)
    DO j = 1,S%nat3
      DO k = 1,S%nat3
        !
        bose_C = 2* (bose(j,2) - bose(k,3))
        bose_X = bose(j,2) + bose(k,3) + 1
        freqtotm1_23= freqm1(j,2) * freqm1(k,3)
        !
        DO i = 1,S%nat3
          !
          !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
          !
          freqtotm1 = freqm1(i,1) * freqtotm1_23
          !IF(freqtot/=0._dp)THEN
          !
          ! dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
          ctm_C = bose_C * weights_C(i,S%nat3*(j-1) + k) !* f_gauss(dom_C, sigma)
          !
          ! dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
          ctm_X = bose_X * weights_X(i,S%nat3*(j-1) + k)!* f_gauss(dom_X, sigma)
          !
          lw(i) = lw(i) + (ctm_C + ctm_X) * V3sq(i,j,k)*freqtotm1
          !
          !leftover_e = pi*freqtotm1 * (ctm_C*dom_C + ctm_X*dom_X) * V3sq(i,j,k)
          !ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(S%nat3, lw, freq(:,1))
    sum_linewidth_modes_weights = lw
    !
  END FUNCTION sum_linewidth_modes_weights

  FUNCTION sum_linewidth_modes_dense(S, freq, bose, V3sq, nu0, wg)
    USE constants, ONLY : pi
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    REAL(DP), INTENT(IN) :: wg(S%nat3**2,S%nat3,2)
    !
    REAL(DP) :: sum_linewidth_modes_dense(S%nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_C, bose_X ! final/initial state populations
    REAL(DP) :: ctm_C, ctm_X   !
    REAL(DP) :: freqtotm1, freqtotm1_23
    REAL(DP) :: freqm1(S%nat3,3)
    INTEGER :: index_double
    !REAL(DP),SAVE :: leftover_e
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)!, sigma_
    lw(:) = 0._dp
    !
    freqm1 = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: lw) COLLAPSE(2)
    DO j = 1,S%nat3
      DO k = 1,S%nat3
        !
        index_double = S%nat3*(j-1) + k
        bose_C = 2* (bose(j,2) - bose(k,3))
        bose_X = bose(j,2) + bose(k,3) + 1
        freqtotm1_23= freqm1(j,2) * freqm1(k,3)
        !
        DO i = 1,S%nat3

          !
          !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
          !
          freqtotm1 = freqm1(i,1) * freqtotm1_23
          !IF(freqtot/=0._dp)THEN
          !
          ! dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
          ctm_C = bose_C * wg(index_double,i,2) !* f_gauss(dom_C, sigma)
          !
          ! dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
          ctm_X = bose_X * wg(index_double,i,1)!* f_gauss(dom_X, sigma)
          !
          lw(i) = lw(i) + (ctm_C + ctm_X) * V3sq(i,j,k)*freqtotm1
          !
          !leftover_e = pi*freqtotm1 * (ctm_C*dom_C + ctm_X*dom_X) * V3sq(i,j,k)
          !ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(S%nat3, lw, freq(:,1))
    sum_linewidth_modes_dense = lw
    !
  END FUNCTION sum_linewidth_modes_dense

  FUNCTION sum_linewidth_modes(S, sigma, freq, bose, V3sq, nu0)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    REAL(DP) :: sum_linewidth_modes(S%nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_C, bose_X ! final/initial state populations
    REAL(DP) :: dom_C, dom_X   ! \delta\omega
    REAL(DP) :: ctm_C, ctm_X   !
    REAL(DP) :: freqtotm1, freqtotm1_23
    REAL(DP) :: freqm1(S%nat3,3)
    !REAL(DP),SAVE :: leftover_e
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)!, sigma_
    lw(:) = 0._dp
    !
    freqm1 = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: lw) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        bose_C = 2* (bose(j,2) - bose(k,3))
        bose_X = bose(j,2) + bose(k,3) + 1
        freqtotm1_23= freqm1(j,2) * freqm1(k,3)
        !
        DO i = 1,S%nat3
          !
          !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
          !
          freqtotm1 = freqm1(i,1) * freqtotm1_23
          !IF(freqtot/=0._dp)THEN
          !
          dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
          ctm_C = bose_C * f_gauss(dom_C, sigma)
          !
          dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
          ctm_X = bose_X * f_gauss(dom_X, sigma)
          !
          lw(i) = lw(i) - pi * (ctm_C + ctm_X) * V3sq(i,j,k)*freqtotm1
          !
          !leftover_e = pi*freqtotm1 * (ctm_C*dom_C + ctm_X*dom_X) * V3sq(i,j,k)
          !ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(S%nat3, lw, freq(:,1))
    sum_linewidth_modes = lw
    !
  END FUNCTION sum_linewidth_modes

  ! \/o\________\\\_________________________________________/^>
  ! Add the elastic peak of Raman
  SUBROUTINE add_exp_t_factor(nconf, T, ne, nat3, ener, spectralf)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: ne, nat3
    REAL(DP),INTENT(in) :: ener(ne)
    !
    REAL(DP),INTENT(inout) :: spectralf(ne,nat3,nconf)
    !
    REAL(DP) :: factor(ne)
    INTEGER :: it,ie,nu
    !
    DO it = 1,nconf
      factor = (1 + f_bose(ener,T(it))) / ener
      !
      DO nu = 1,nat3
        DO ie = 1,ne
          spectralf(ie,nu,it) = factor(ie)*spectralf(ie,nu,it)
        ENDDO
      ENDDO
      !
    ENDDO
    !
  END SUBROUTINE add_exp_t_factor
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE gauss_convolution(nconf, T, ne, nat3, ener, spectralf)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: ne, nat3
    REAL(DP),INTENT(in) :: ener(ne)
    !
    REAL(DP),INTENT(inout) :: spectralf(ne,nat3,nconf)
    !
    REAL(DP) :: convol(ne)
    !
    INTEGER :: it,ie,nu
    !
    convol = f_bose(ener,T(it))

    DO it = 1,nconf
      DO nu = 1,nat3
        DO ie = 1,ne
          spectralf(ie,nu,it) = convol(ie)*spectralf(ie,nu,it)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE gauss_convolution
  !
  !
END MODULE linewidth
