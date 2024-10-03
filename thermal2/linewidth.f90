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
  USE kinds,            ONLY : DP
  USE mpi_thermal,      ONLY : my_id, num_procs, mpi_bsum, allgather_mat, scatteri_tns, mpi_sum_mat
  USE constants,        ONLY : RY_TO_CMM1
  USE q_grids,          ONLY : q_grid, setup_grid, fc_info
  USE constants,        ONLY : pi
  USE input_fc,         ONLY : ph_system_info
  USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat, d3_mixed
  USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
  USE functions,        ONLY : refold_bz
  USE thtetra,          ONLY : tetra_init, tetra_delta
  USE code_input,       ONLY : code_input_type
  USE timers

  TYPE(code_input_type)           :: input
  TYPE(forceconst2_grid)          :: fc2
  CLASS(forceconst3), POINTER     :: fc3
  TYPE(ph_system_info)            :: S
  TYPE(q_grid)                    :: grid
  INTEGER                         :: nat3
  INTEGER                         :: nconf
  CHARACTER(10)                   :: calc
  REAL(DP)                        :: xq(3,3)
  REAL(DP), ALLOCATABLE           :: V3(:,:,:)
  REAL(DP), ALLOCATABLE           :: weights_C(:,:,:), weights_X(:,:,:), freq(:,:), bose(:,:)
  COMPLEX(DP), ALLOCATABLE        :: U(:,:,:), D3(:,:,:), lw_UN(:,:,:)
  REAL(DP)                        :: energy
  TYPE(d3_mixed)                  :: Dqr

  external :: cryst_to_cart
  EXTERNAL :: errore
CONTAINS
  !
  SUBROUTINE lw_init(input_, fc, calc_)
    TYPE(code_input_type), INTENT(IN) :: input_
    TYPE(fc_info), INTENT(IN) :: fc
    CHARACTER(*), INTENT(IN), OPTIONAL :: calc_
    fc2 = fc%fc2
    fc3 => fc%fc3
    S = fc%S
    grid = fc%grid
    nat3 = fc%S%nat3
    nconf = input%nconf
    if(PRESENT(calc_)) then
      calc = calc_
    else
      calc = input%delta_approx
    endif
    input = input_
    !
    ALLOCATE(weights_C(nat3, nat3**2, grid%nqtot), weights_X(nat3, nat3**2, grid%nqtot))
    ALLOCATE(U(nat3,nat3,3), D3(nat3,nat3,nat3), V3(nat3,nat3,nat3), freq(nat3,3), bose(nat3,3))
    if(input%store_lw) ALLOCATE(lw_UN(nat3,2,input%nconf))
  END SUBROUTINE
  !
  FUNCTION linewidth_q(xq0, input_, fc)
    REAL(DP),INTENT(in) :: xq0(3)
    TYPE(code_input_type), INTENT(IN) :: input_
    TYPE(fc_info), INTENT(IN) :: fc
    !
    REAL(DP) :: linewidth_q(fc%S%nat3,input_%nconf)

    IF(.not.ALLOCATED(D3)) CALL lw_init(input_, fc)
    xq(:,1) = xq0
    SELECT CASE (input%delta_approx)
     CASE ("tetra")
      CALL weights_tetra()
      calc = "lwtetra"
      linewidth_q = REAL(sum_q2(), DP)
     CASE("gauss")
      calc = "lwgauss"
      linewidth_q = REAL(sum_q2(), DP)
     CASE DEFAULT
      CALL errore("linewidth_q", "only delta/gauss as delta_approx are permitted", 1)
    END SELECT

  END FUNCTION linewidth_q

  FUNCTION Pijk(coal)
    LOGICAL, INTENT(IN) :: coal
    !! true if we are computing the coalescence term
    INTEGER :: iq, i,j,k
    !> RETURN VALUE
    REAL(DP) :: Pijk(nat3,nat3,nat3)
    !
    if(coal) THEN
      xq(:,3) = xq(:,1) + xq(:,2)
      CALL fc3%interpolate(xq(:,2), - xq(:,3), S%nat3, D3)
      xq(:,3) = xq(:,1) - xq(:,2)
      CALL fc3%interpolate(- xq(:,2), - xq(:,3), S%nat3, D3)
    ENDIF

    DO iq = 1,3
      CALL freq_phq_safe(xq(:,iq), S, fc2, freq(:,iq), U(:,:,iq))
      CALL bose_phq(300._dp, nat3, freq(:,iq), bose(:,iq))
    ENDDO
    CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
    V3 = REAL( CONJG(D3)*D3 , kind=DP)

    DO i = 1, nat3
      DO j = 1, nat3
        DO k = 1, nat3
          if(coal) then
            Pijk(i,j,k) = weights_C(i,nat3*(j-1) + k, iq) * &
              bose(i,1) * bose(j,2) * (1.0_dp + bose(k,3)) * V3(i,j,k)
          else
            Pijk(i,j,k) = weights_X(i,nat3*(j-1) + k, iq) * &
              bose(i,1) * (1.0_dp + bose(j,2)) * (1.0_dp + bose(k,3)) * V3(i,j,k)
          endif
        ENDDO
      ENDDO
    ENDDO
  END FUNCTION

! <<^V^\\=========================================//-//-//========//O\\//
  SUBROUTINE weights_tetra()
    USE q_grids,          ONLY : q_grid
    USE constants,        ONLY : pi
    USE input_fc,         ONLY : ph_system_info
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE functions,        ONLY : refold_bz
    USE thtetra,          ONLY : tetra_init, tetra_weights_delta
    IMPLICIT NONE
    !
    INTEGER :: iq, ibnd, jbnd, index_double
    REAL(DP) :: freqs_doubled_X(nat3**2, grid%nq), freqs_doubled_C(nat3**2, grid%nq)
    REAL(DP), ALLOCATABLE :: gather_doubled_X(:,:), gather_doubled_C(:,:)
    !
    CALL tetra_init( grid%n, S%bg, .true., grid%scattered)

    ! ALLOCATE(weights_X(nat3, nat3**2, grid%nqtot), weights_C(nat3, nat3**2, grid%nqtot))
    !
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1))
    ! this cycle initializes freqs, as a grid of frequencies
    ! this cycle initializes freqs_doubled, which is freq(ibnd,2) +/- freq(jbnd,3)
    timer_CALL t_freqd%start()
    DO iq = 1, grid%nq
      CALL freq_phq_safe(grid%xq(:,iq), S, fc2, freq(:,2))
      ! the third vector (why it's minus?)
      xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
      CALL freq_phq_safe(xq(:,3), S, fc2, freq(:,3))
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
      gather_doubled_X = freqs_doubled_X
      gather_doubled_C = freqs_doubled_C
    ENDIF

    timer_CALL t_freqd%stop()
    timer_CALL t_thtetra%start()
    DO ibnd = 1, S%nat3
      weights_X(ibnd,:,:) = tetra_weights_delta(grid%nqtot, nat3**2, gather_doubled_X, freq(ibnd,1))
      weights_C(ibnd,:,:) = tetra_weights_delta(grid%nqtot, nat3**2, gather_doubled_C, freq(ibnd,1))
    ENDDO
    IF(grid%scattered) THEN
      DO iq = 1, grid%nqtot
        CALL mpi_sum_mat(nat3, nat3**2, weights_X(:,:,iq)) ! MPI_SUM MPI_REDUCE
        CALL mpi_sum_mat(nat3, nat3**2, weights_C(:,:,iq))
      ENDDO
      CALL scatteri_tns(nat3, nat3**2, grid%nqtot, weights_X)
      CALL scatteri_tns(nat3, nat3**2, grid%nqtot, weights_C)
    ENDIF

    timer_CALL t_thtetra%stop()

  END SUBROUTINE

  PURE function rotate_single_d3(i, j, k, invert) result(D3_s2)
    INTEGER, INTENT(IN) :: i, j, k
    LOGICAL, INTENT(IN) :: invert
    !
    INTEGER :: a, b, c
    COMPLEX(DP) :: D3_s, D3_s2, aux
    D3_s = 0._dp
    DO c = 1, nat3
      DO b = 1, nat3
        if (invert) then
          aux = U(b,i,1)* U(c,k,3)
        else
          aux = U(b,j,2)* U(c,k,3)
        endif
        DO a = 1, nat3
          if (invert) then
            D3_s = D3_s + D3(a,b,c) * CONJG( aux * U(a,j,2) )
          else
            D3_s = D3_s + D3(a,b,c) * CONJG( aux * U(a,i,1) )
          endif
        END DO
      END DO
    END DO
    D3_s2 = REAL( CONJG(D3_s)* D3_s, kind=DP)
  end function

  FUNCTION sum_q2(xq1, coalescence)
    USE functions, ONLY : f_gauss
    USE fc3_interpolate, ONLY : sum_R3, d3_mixed
    USE merge_degenerate,   ONLY : merge_degen

    !! Normal/Umklapp contribution
    REAL(DP), INTENT(IN), OPTIONAL :: xq1(3)
    LOGICAL, INTENT(IN), OPTIONAL :: coalescence
    INTEGER,PARAMETER :: normal=1, umklapp=2
    INTEGER :: iq, jq, j_un, it, nu0(3), i,j,k
    REAL(DP) :: f(3)
    REAL(DP) :: freqm1(nat3,3), freqtotm1_23, freqtotm1, bose_C, bose_X, dom_C, dom_X, sigma
    COMPLEX(DP) :: aux(nat3), sum_q2(nat3,nconf), ctm

    ! ALLOCATE(U(nat3,nat3,3), D3(nat3,nat3,nat3))

    sum_q2 = 0._dp
    IF(ALLOCATED(lw_UN)) lw_UN = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    timer_CALL t_freq%start()
    IF(PRESENT(xq1)) xq(:,1) = xq1
    nu0(1) = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    timer_CALL t_freq%stop()

    CALL fc3%sum_R2(xq(:,1), nat3, Dqr)
    DO iq = 1, grid%nq
      !
      timer_CALL t_freq%start()
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      IF(PRESENT(coalescence)) THEN
        xq(:,2) = grid%xq(:,iq)
        xq(:,3) = -(xq(:,2)+xq(:,1))
      ELSE
        xq(:,2) = - grid%xq(:,iq)
        xq(:,3) = - xq(:,1) + xq(:,2)
      ENDIF
      IF(ALLOCATED(lw_UN)) THEN
        IF(ALL(ABS(refold_bz(xq(:,1), S%bg) &
          +refold_bz(xq(:,2), S%bg) &
          +refold_bz(xq(:,3), S%bg))<1.d-6)) THEN
          j_un = normal
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
      ! CALL fc3%interpolate(xq(:,2), xq(:,3), nat3, D3)
      CALL sum_R3(S, xq(:,3), Dqr, D3)
      print*, D3(1,2,3)
      timer_CALL t_fc3int%stop()
      DO it = 1,nconf
        timer_CALL t_bose%start()
        ! Compute bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(input%T(it),nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        timer_CALL t_bose%stop()
        timer_CALL t_sum%start()
        freqm1 = 0._dp
        aux = 0._dp
        DO i = 1,nat3
          do j = 1,3
            IF(i>=nu0(j)) freqm1(i,j) = 0.5_dp/freq(i,j)
          ENDDO
        ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: aux) COLLAPSE(2)
        DO k = 1,nat3
          DO j = 1,nat3
            !
            bose_C = 2* (bose(j,2) - bose(k,3))
            bose_X = bose(j,2) + bose(k,3) + 1
            freqtotm1_23= freqm1(j,2) * freqm1(k,3)
            !
            DO i = 1,nat3
              !
              !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
              !
              freqtotm1 = freqm1(i,1) * freqtotm1_23
              !IF(freqtot/=0._dp)THEN
              !
              sigma = input%sigma(it)
              SELECT CASE (calc)
               CASE ("lwgauss")
                dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
                dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
                ctm = bose_C * f_gauss(dom_C, sigma) + bose_X * f_gauss(dom_X, sigma)
               CASE ("lwtetra")
                ctm = bose_C * weights_C(i,nat3*(j-1) + k, iq) + bose_X * weights_X(i,nat3*(j-1) + k, iq)
               CASE ("selfnrg")
                f(1) = freq(i,1)
                f(2) = freq(j,2)
                f(3) = freq(k,3)
                ctm = ctm_selfnrg(sigma, input%T(it), f, bose_C, bose_X)
               CASE ("selfnrg_sp")
                f(1) = freq(i,1)
                f(2) = freq(j,2)
                f(3) = freq(k,3)
                ctm = ctm_selfnrg_spectre(sigma, f, energy, bose_C, bose_X)
               CASE("cgp_C")
                !> the prefactor 2 is to annul the 1/2 in the definition of the linewidth, but maybe
                !> there's another prefactor
                ctm = 2 * bose(i,1) * bose(j,2) * (1.0_dp + bose(k,3)) * weights_C(i,nat3*(j-1) + k, iq)
               CASE DEFAULT
                CALL errore("sum_rotate_lw", "you should give sigma/tetra_weights", 1)
              END SELECT
              !
              IF(TRIM(input%calculation) == "cgp" .or. TRIM(input%calculation) == "exact") &
                ctm = ctm * 2 * bose(i,1) * (bose(i,1) + 1)
              ! IF(REAL(ctm, DP) < 1) CYCLE
              ! ci sono tanti negativi che contribuiscono sulla terza cifra, mica da poco
              !
              !> rotation of single D3 only where we know that we have scattering
              timer_CALL t_fc3rot%start()
              aux(i) = aux(i) + ctm * rotate_single_d3(i,j,k, invert=.true.) * freqtotm1
              timer_CALL t_fc3rot%stop()
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        !
        CALL merge_degen(nat3, aux, freq(:,1))

        sum_q2(:,it) = sum_q2(:,it) + aux * grid%w(iq)
        IF(ALLOCATED(lw_UN)) lw_UN(:,j_un,it) = lw_UN(:,j_un,it) + aux
        timer_CALL t_sum%stop()
      ENDDO
      !
    ENDDO

    timer_CALL t_mpicom%start()
    IF(grid%scattered) CALL mpi_bsum(nat3,input%nconf,sum_q2)
    IF(grid%scattered .and. ALLOCATED(lw_UN)) CALL mpi_bsum(nat3,2,input%nconf,lw_UN)
    timer_CALL t_mpicom%stop()

    sum_q2 = sum_q2 * pi/2
  END FUNCTION

  ! \/o\________\\\_________________________________________/^>
  ! Simple spectral function, computed as a superposition of Lorentzian functions
  FUNCTION simple_spectre_q(xq0, input, fc, ne, ener, freq1, U1, shift) &
    RESULT(spectralf)
    USE q_grids,            ONLY : q_grid
    USE constants,          ONLY : RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe
    USE fc3_interpolate,    ONLY : forceconst3

    IMPLICIT NONE
    ! ARGUMENTS:
    REAL(DP),INTENT(in) :: xq0(3)
    TYPE(code_input_type), INTENT(IN) :: input
    TYPE(fc_info), INTENT(IN) :: fc
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    LOGICAL,INTENT(in)  :: shift
    !
    REAL(DP),OPTIONAL,INTENT(in) :: freq1(fc%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(fc%nat3,fc%nat3)
    CHARACTER(10), PARAMETER :: SELFNRG_STRING = "selfnrg"
    !
    ! FUNCTION RESULT:
    REAL(DP) :: spectralf(ne,fc%nat3,input%nconf)
    !
    COMPLEX(DP) :: selfnrg(fc%nat3,input%nconf)
    !
    INTEGER :: it, i, ie
    REAL(DP) :: gamma(fc%nat3,input%nconf), delta(fc%nat3,input%nconf), omega, denom
    COMPLEX(DP) :: U(fc%nat3, fc%nat3)
    REAL(DP) :: freq(fc%nat3)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    IF(present(freq1) .and. present(U1)) THEN
      freq = freq1
      U    = U1
      !CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      CALL freq_phq_safe(xq0, fc%S, fc%fc2, freq, U)
    ENDIF
    !
    ! use lineshift to compute 3rd order linewidth and lineshift
    IF(shift)THEN
      selfnrg = sum_q2()
      gamma(:,:) = -DIMAG(selfnrg(:,:))
      delta(:,:) =   DBLE(selfnrg(:,:))
      gamma = linewidth_q(xq0, input, fc)
      delta = 0._dp
    ENDIF
    !
    DO it = 1,input%nconf
      WRITE(*,'(30x,6e12.2,5x,6e12.2)') gamma(:,it)*RY_TO_CMM1,delta(:,it)*RY_TO_CMM1
    ENDDO
    !
    delta = 0._dp
    ! Compute and superpose the Lorentzian functions corresponding to each band
    DO it = 1,input%nconf
      DO i = 1,fc%nat3
        DO ie = 1, ne
          ! Only apply shift (i.e. real part of self nrg) if requested
          omega = freq(i)
          denom =   (ener(ie)**2 -omega**2 -2*omega*delta(i,it))**2 &
            + 4*omega**2 *gamma(i,it)**2
          IF(ABS(denom)/=0._dp)THEN
            spectralf(ie,i,it) = 2*omega*gamma(i,it) / denom
            spectralf(ie,i,it) = 0._dp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
  END FUNCTION simple_spectre_q


  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Spectral weight function, computed as in eq. 1 of arXiv:1312.7467v1
  FUNCTION spectre_q(xq0, input, fc, ener, shift) &
    RESULT(spectralf)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    TYPE(code_input_type), INTENT(IN) :: input
    TYPE(fc_info), INTENT(IN) :: fc
    !
    REAL(DP),INTENT(in) :: ener(input%ne)
    !
    LOGICAL,INTENT(in)  :: shift ! set to false to drop the real part of the self-energy
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie, it
    REAL(DP) :: gamma, delta, omega, denom, freq1(fc%nat3)
    COMPLEX(DP) :: selfnrg(input%ne,fc%nat3,input%nconf)
    ! FUNCTION RESULT:
    REAL(DP)    :: spectralf(input%ne,fc%nat3,input%nconf)

    calc = "selfnrg_sp"
    ! Once we have the self-energy, the rest is trivial
    do ie = 1, input%ne
      energy = ener(ie)
      selfnrg(ie,:,:) = sum_q2()
    enddo
    !
    timer_CALL t_mkspf%start()
    DO it = 1,input%nconf
      CALL freq_phq_safe(xq0, fc%S, fc%fc2, freq1)
      DO i = 1,fc%nat3
        DO ie = 1, input%ne
          gamma =  -DIMAG(selfnrg(ie,i,it))
          IF(shift) THEN
            delta =   DBLE(selfnrg(ie,i,it))
            delta = 0._dp
          ENDIF
          omega = freq1(i)
          denom =   (ener(ie)**2 -omega**2 -2*omega*delta)**2 &
            + 4*omega**2 *gamma**2
          IF(ABS(denom)/=0._dp)THEN
            spectralf(ie,i,it) = 2*omega*gamma / denom
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
  FUNCTION spectre2_q(xq0, input, fc, ener) &
    RESULT(spectralf)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    TYPE(code_input_type), INTENT(IN) :: input
    TYPE(fc_info), INTENT(IN) :: fc
    REAL(DP),INTENT(in) :: xq0(3)
    REAL(DP),INTENT(in) :: ener(input%ne)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie, it
    COMPLEX(DP) :: tepsilon(input%ne,fc%nat3,input%nconf)
    ! FUNCTION RESULT:
    REAL(DP)    :: spectralf(input%ne,fc%nat3,input%nconf)
    !
    ! Once we have the self-energy, the rest is trivial
    tepsilon = tepsilon_q(xq0, input, fc, ener)
    !
    timer_CALL t_mkspf%start()
    DO it = 1,input%nconf
      DO i = 1,fc%nat3
        DO ie = 1,input%ne
          !IF(freq1(i)/=0._dp)THEN
          spectralf(ie,i,it) = DIMAG(tepsilon(ie,i,it))
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
  FUNCTION tepsilon_q(xq0, input, fc, ener) &
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
    TYPE(code_input_type), INTENT(IN) :: input
    TYPE(fc_info), INTENT(IN) :: fc
    !
    REAL(DP),INTENT(in) :: ener(input%ne)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie, it
    REAL(DP) :: freq1(fc%nat3)
    REAL(DP) :: pref, epsilon_infty, rmass, zeu_i2
    COMPLEX(DP) :: denom
    COMPLEX(DP):: selfnrg(input%ne,fc%nat3,input%nconf)
    ! FUNCTION RESULT:
    COMPLEX(DP)    :: tepsilon(input%ne,fc%nat3,input%nconf)
    CHARACTER(10), PARAMETER :: SELFNRG_SPECTRE_STRING = "selfnrg_sp"
    !
    IF(.not.fc%S%lrigid) CALL errore("tepsilon","Cannot compute \tilde{epsilon} without epsilon0", 1)
    ioWRITE(*,*) "BEWARE: tepsilon is only isotropic"
    !
    ! Once we have the self-energy, the rest is trivial
    !
    ! Once we have the self-energy, the rest is trivial
    do ie = 1, input%ne
      energy = ener(ie)
      selfnrg(ie,:,:) = sum_q2()
    enddo
    !
    ! S = 4piZ^2/(volume mass omega0^2)

    rmass = fc%S%amass(1)*fc%S%amass(2)/(fc%S%amass(1)+fc%S%amass(2))
    epsilon_infty = fc%S%epsil(1,1)
    zeu_i2        = fc%S%zeu(1,1,1)**2
    pref =  fpi * zeu_i2 /(fc%S%omega * rmass)
    !pref = 1/freq1(i)**2 ! 4 pi Z^2  / (Volume mu  omega0^2) # mu = reduced mass
!     print*, "rmass", rmass, fc%S%amass(1:2)
!     print*, "epsilon", epsilon_infty
!     print*, "zeu2", zeu_i2
!     print*, "pref", pref
!     print*, "vol", fc%S%omega
!     print*, "denom", rmass*fc%S%omega
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
    DO it = 1,input%nconf
      CALL freq_phq_safe(xq0, fc%S, fc%fc2, freq1)
      DO i = 1,fc%S%nat3
        IF(freq1(i)==0._dp) CYCLE
        DO ie = 1,input%ne
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
  FUNCTION ir_reflectivity_q(xq0, input, fc, ener) &
    RESULT(reflectivity)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
    USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    TYPE(code_input_type), INTENT(IN) :: input
    TYPE(fc_info), INTENT(IN) :: fc
    REAL(DP),INTENT(in) :: ener(input%ne)
    !
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, it
    COMPLEX(DP) :: tepsilon(input%ne,fc%S%nat3,input%nconf)
    COMPLEX(DP) :: aux(input%ne),aux2(input%ne)
    ! FUNCTION RESULT:
    REAL(DP)    :: reflectivity(input%ne,fc%S%nat3,input%nconf)
    !
    ! We use tilde{epsilon}, which itself comes from the self-energy
    tepsilon = tepsilon_q(xq0, input, fc, ener)
    !
    timer_CALL t_mkspf%start()
    DO it = 1,input%nconf
      DO i = 1,fc%S%nat3
        aux  = SQRT(tepsilon(:,i,it))
        aux2 = (aux-1)/(aux+1)
        reflectivity(:,i,it) = DBLE(aux2*CONJG(aux2))
      ENDDO
    ENDDO
    timer_CALL t_mkspf%stop()
    !
  END FUNCTION ir_reflectivity_q
  !

  FUNCTION ctm_selfnrg_spectre(sigma, freq, ener, bose_C, bose_X)
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    USE functions,          ONLY : sigma_mgo
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: sigma   ! smearing (regularization) (Ry)
    REAL(DP),INTENT(in) :: freq(3)  ! phonon energies (Ry)
    !
    REAL(DP),INTENT(in) :: ener     ! energies for which to compute the spectral function
    REAL(DP), INTENT(IN) :: bose_C, bose_X
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    COMPLEX(DP) :: ctm_P,ctm_M, reg
    !
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: ctm_selfnrg_spectre
    !
!     IF(sigma<=0._dp)THEN
!       CALL errore("sum_selfnrg_spectre","spf not implemented in the static limit. "&
!                   //"NEW: To do unshifted spf use 'spf imag'",1)
!     ENDIF
    !
    !
    omega_P  = freq(2)+freq(3)
    omega_P2 = omega_P**2
    !
    omega_M  = freq(2)-freq(3)
    omega_M2 = omega_M**2
    !
    !
    reg = CMPLX(ener, sigma, kind=DP)**2
    !
    ctm_P = 2 * bose_C *omega_P/(omega_P2-reg)
    ctm_M = 2 * bose_X *omega_M/(omega_M2-reg)
    !
    ctm_selfnrg_spectre = ctm_P + ctm_M
    !
  END FUNCTION

  !
  ! \/o\________\\\_________________________________________/^>
  ! Sum the self energy in a range of frequencies
  FUNCTION ctm_selfnrg(sigma, T, freq, bose_C, bose_X)
    USE input_fc,           ONLY : ph_system_info
    USE functions,          ONLY : df_bose
    USE merge_degenerate,   ONLY : merge_degen
    USE functions,          ONLY : sigma_mgo
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: sigma, T, freq(3), bose_C, bose_X
    REAL(DP) :: omega_X,  omega_C   ! \sigma\omega
    COMPLEX(DP) :: ctm_X, ctm_C, reg, ctm_selfnrg
    !
    ctm_X = 0._dp
    ctm_C = 0._dp
    omega_X  = freq(2)+freq(3)
    omega_C  = freq(2)-freq(3)
    IF(sigma>0._dp)THEN
      reg = CMPLX(freq(1), sigma, kind=DP)**2
      ctm_X = 2 * bose_X *omega_X/(omega_X**2-reg )
      ctm_C = 2 * bose_C *omega_C/(omega_C**2-reg )
      ctm_X = 2 * bose_X *omega_X/(omega_X**2+sigma**2)
      ctm_C = 2 * bose_C *omega_C/(omega_C**2+sigma**2)
      ! In the static limit with sigma=0 case we have to take the
      ! derivative of (n_3-n2)/(w_2-w_3) when w_2 is close to w_3
      IF(omega_X>0._dp)THEN
        ctm_X = 2 * bose_X /omega_X
        ctm_X = 0._dp
      ENDIF
      !
      IF(ABS(omega_C)>1.e-5_dp)THEN
        ctm_C = 2 * bose_C /omega_C
        IF(T>0._dp)THEN
          ctm_C = -2* df_bose(0.5_dp * omega_X, T)
          ctm_C = 0._dp
        ENDIF
      ENDIF
      !
    ENDIF

    ctm_selfnrg = ctm_C + ctm_X
  END FUNCTION
  !
  FUNCTION sum_rotate_lw(S, freq, bose, D3, U, nu0, calc, sigma, T, weights_C, weights_X)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: freq(S%nat3,3), bose(S%nat3,3)
    COMPLEX(DP),INTENT(in) :: D3(S%nat3,S%nat3,S%nat3), U(S%nat3,S%nat3,3)
    INTEGER,INTENT(in)  :: nu0(3)
    CHARACTER(10), INTENT(IN) :: calc
    REAL(DP),INTENT(in), OPTIONAL :: sigma, T
    REAL(DP), INTENT(IN), OPTIONAL :: weights_C(S%nat3,S%nat3**2), weights_X(S%nat3,S%nat3**2)
    !
    COMPLEX(DP) :: sum_rotate_lw(S%nat3)
    COMPLEX(DP) :: D3_s
    REAL(DP) :: D3_s2
    INTEGER :: a,b,c
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_C, bose_X ! final/initial state populations
    REAL(DP) :: dom_C, dom_X   ! \delta\omega
    COMPLEX(DP) :: ctm   !
    REAL(DP) :: freqtotm1, freqtotm1_23
    REAL(DP) :: freqm1(S%nat3,3)
    COMPLEX(DP) ::  aux
    !REAL(DP),SAVE :: leftover_e
    !
    INTEGER :: i,j,k
    !
    freqm1 = 0._dp
    sum_rotate_lw = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: sum_rotate_lw) COLLAPSE(2)
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
          SELECT CASE (calc)
           CASE ("gauss")
            dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
            dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
            ctm = (bose_C * f_gauss(dom_C, sigma) + bose_X * f_gauss(dom_X, sigma))
           CASE ("tetra")
            ctm = bose_C * weights_C(i,S%nat3*(j-1) + k) + bose_X * weights_X(i,S%nat3*(j-1) + k)
           CASE ("selfnrg")
            ctm = ctm_selfnrg(sigma, T, freq, bose_C, bose_X)
           CASE DEFAULT
            CALL errore("sum_rotate_lw", "you should give sigma/tetra_weights", 1)
          END SELECT
          ! IF(REAL(ctm, DP) < 1) CYCLE
          ! ci sono tanti negativi che contribuiscono sulla terza cifra, mica da poco
          !
          timer_CALL t_fc3rot%start()
          !> rotation of single D3 only where we know that we have scattering
          D3_s = 0._dp
          DO c = 1, S%nat3
            DO b = 1, S%nat3
              aux = U(b,i,1)* U(c,k,3)
              DO a = 1, S%nat3
                D3_s = D3_s + D3(a,b,c) * CONJG( aux * U(a,j,2) )
              END DO
            END DO
          END DO
          timer_CALL t_fc3rot%stop()
          D3_s2 = REAL( CONJG(D3_s)* D3_s, kind=DP)
          sum_rotate_lw(i) = sum_rotate_lw(i) + ctm * D3_s2 * freqtotm1
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(S%nat3, sum_rotate_lw, freq(:,1))
    !
  END FUNCTION sum_rotate_lw

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
