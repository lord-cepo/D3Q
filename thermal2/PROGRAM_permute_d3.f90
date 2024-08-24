PROGRAM linewidth

  USE kinds,            ONLY : DP
  USE linewidth_program
  USE input_fc,         ONLY : forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi, ionode
  USE more_constants,   ONLY : print_citations_linewidth
  USE mc_grids,         ONLY : print_optimized_stats
  USE thtetra,          ONLY : tetra_init
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(code_input_type)     :: lwinput
  TYPE(q_grid)      :: qpath
  COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:), D3_perm(:,:,:)
  REAL(DP) :: xq(3,3)

  
  !   CALL mp_world_start(world_comm)
  !   CALL environment_start('LW')
  
  CALL start_mpi()
  CALL init_nanoclock()
  !
  IF(ionode) CALL print_citations_linewidth()
  
  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("LW", lwinput, qpath, S, fc2, fc3)
  xq(:,1) = 0.2_dp
  xq(:,2) = 0.1_dp
  xq(:,3) = -xq(:,1) -xq(:,2)
  CALL cryst_to_cart(3, xq, S%bg, 1)

  ! if (ionode) CALL tetra_init( lwinput%nk, S%bg, .true.)
  !
  ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
  ALLOCATE(D3_perm(S%nat3, S%nat3, S%nat3))
  CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
  CALL fc3%interpolate(xq(:,1), xq(:,3), S%nat3, D3_perm)
  D3_perm = reshape(D3_perm, shape(D3_perm), order=[2,1,3])

  print*, ABS(d3 - d3_perm) < 1e-13

  CALL stop_mpi()

END PROGRAM