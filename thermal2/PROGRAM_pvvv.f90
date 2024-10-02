module mpvvv
  USE kinds, ONLY : DP
  use functions, only: f_gauss
  USE fc2_interpolate,    ONLY : bose_phq, freq_phq_safe
  USE input_fc,         ONLY : forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid
  USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
contains
  REAL(DP) FUNCTION Pijk(xq1, xq2, i, j, k, S, fc2, fc3, coal, verbose)
    implicit none

    REAL(DP), INTENT(IN) :: xq1(3), xq2(3)
    INTEGER, INTENT(IN) :: i,j,k
    TYPE(forceconst2_grid), INTENT(IN) :: fc2
    CLASS(forceconst3), INTENT(IN) :: fc3
    TYPE(ph_system_info), INTENT(IN)   :: S
    LOGICAL, INTENT(IN) :: coal
    LOGICAL, OPTIONAL, INTENT(IN) :: verbose
    COMPLEX(DP) :: D3(S%nat3,S%nat3,S%nat3), U(S%nat3,S%nat3,3)
    REAL(DP) :: freq(S%nat3,3), V3(S%nat3,S%nat3,S%nat3), bose(S%nat3,3)
    REAL(DP) :: xq(3,3)
    REAL(DP) :: sigma = 0.1_dp
    INTEGER :: iq

    xq(:,1) = xq1
    xq(:,2) = xq2

    if(coal) THEN
      xq(:,3) = xq(:,1) + xq(:,2)
      CALL fc3%interpolate(xq(:,2), - xq(:,3), S%nat3, D3)
    ELSE
      xq(:,3) = xq(:,1) - xq(:,2)
      CALL fc3%interpolate(- xq(:,2), - xq(:,3), S%nat3, D3)
    ENDIF

    DO iq = 1,3
      CALL freq_phq_safe(xq, S, fc2, freq(:,i), U(:,:,i))
      CALL bose_phq(300._dp, S%nat3, freq(:,i), bose(:,i))
    ENDDO
    CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
    V3 = REAL( CONJG(D3)*D3 , kind=DP)

    if(coal) then
      Pijk = f_gauss(freq(i,1)+freq(j,2)-freq(k,3), sigma) * &
        bose(i,1) * bose(j,2) * (1.0_dp + bose(k,3)) * V3(i,j,k)
    else
      Pijk = f_gauss(freq(i,1)-freq(j,2)-freq(k,3), sigma) * &
        bose(i,1) * (1.0_dp + bose(j,2)) * (1.0_dp + bose(k,3)) * V3(i,j,k)
    endif

    IF (PRESENT(verbose)) THEN
      PRINT*, "V3", V3(i,j,k)
      IF(coal) then
        print*, "delta", f_gauss(freq(i,1)+freq(j,2)-freq(k,3), sigma)
        print*, "bose", bose(i,1) * bose(j,2) * (1.0_dp + bose(k,3))
      else
        print*, "delta", f_gauss(freq(i,1)-freq(j,2)-freq(k,3), sigma)
        print*, "bose", bose(i,1) * (1.0_dp + bose(j,2)) * (1.0_dp + bose(k,3))
      endif
    endif
  END FUNCTION
END MODULE

PROGRAM pvvv
  use kinds, only : DP
  use linewidth
  use functions, only: f_gauss
  USE input_fc,         ONLY : ph_system_info
  USE fc2_interpolate,  ONLY : forceconst2_grid
  USE q_grids,          ONLY : q_grid, setup_simple_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  use mpvvv,             only: Pijk
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  implicit none
  real(DP) :: xq(3,3), P1_23, P32_1, P12_3
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(code_input_type)     :: input
  TYPE(q_grid)      :: qpoints, qgrid

  integer :: i,j,k, iq
  INTEGER, PARAMETER :: GRID_SIZE = 40
  CALL start_mpi()
  CALL init_nanoclock()

  CALL READ_INPUT("LW", input, qpoints, S, fc2, fc3)

  xq(:,1) = qpoints%xq(:,1)
  i = 3; j = 2; k = 1

  CALL setup_simple_grid(S%bg, GRID_SIZE, GRID_SIZE, GRID_SIZE, qgrid)
  DO iq = 1, qgrid%nq
    xq(:,2) = qgrid%xq(:,iq)
    P1_23 = Pijk(xq(:,1), xq(:,2), i, j, k, S, fc2, fc3, .true.)
    IF (P1_23 > 1e-14) THEN
      print*, "found one big P"
      EXIT
    ENDIF
  ENDDO

  P12_3 = Pijk(xq(:,1), xq(:,2), i, j, k, S, fc2, fc3, .true., .true.)
  ! xq(:,3) = xq(:,1) - xq(:,2)
  P12_3 = Pijk(xq(:,2), xq(:,1), j, i, k, S, fc2, fc3, .true., .true.)

  ! P32_1 = Pijk(xq(:,3), xq(:,2), k, j, i, S, fc2, fc3, .true., .true.)

  print*, P1_23, P32_1

  CALL stop_mpi()
END PROGRAM
