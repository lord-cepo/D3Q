!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Code contributions from Giorgia Fugallo and Michele Lazzeri
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM r2rr
  use mpi_thermal, only: start_mpi, stop_mpi
  USE kinds,           ONLY : DP
  USE input_fc,        ONLY : ph_system_info, aux_system, div_mass_fc2, read_fc2
  USE fc3_interpolate, ONLY : grid
  USE f3_bwfft,        ONLY : d3_list, read_d3_matrices, bwfft_d3_interp, test_fwfft_d3
  USE clib_wrappers,        ONLY : memstat
  use fc2_interpolate, only: forceconst2_grid
  use quter_defect, only : fc_sc2RR, r2q_at_once, allocate_fc2_sc, forceconst2_sc
  use q_grids, only: q_grid, setup_simple_grid
  use d3_basis, only : d3_2idx_2_4idx
  USE cmdline_param_module

  IMPLICIT NONE
  TYPE(d3_list),ALLOCATABLE ::  d3grid(:)
  TYPE(grid)                :: fc3
  type(q_grid) :: q2,q3
  INTEGER :: kb
  integer :: nq_trip, nq_grid, i2, i3, comp
  !
  TYPE(ph_system_info) :: S, Sd
  type(forceconst2_grid) :: fc2d, fc2
  type(forceconst2_sc) :: fc2sc
  COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:)
  complex(dp), allocatable :: D(:,:), D4(:,:)
  !
  INTEGER :: far, ios, icriterium
  CHARACTER(len=512) :: argv, filename, dummy, filein, file2
  CHARACTER(len=:),ALLOCATABLE :: cmdline
  LOGICAL :: write_diff, skip_test

  call start_mpi()
  file2 =    cmdline_param_char("i", "mat2R")
  filein =   cmdline_param_char("x", "matSC")
  filename = cmdline_param_char("o", "matRR")
  far      = cmdline_param_int("f", 2)
  write_diff = cmdline_param_logical("w")
  skip_test = cmdline_param_logical("s")
  icriterium = cmdline_param_int("c", 2)
  !
  IF (cmdline_param_logical('h')) THEN
    WRITE(*,*) "Syntax: ls anh*| d3_qq2rr.x NQX NQY NQZ [-o FILEOUT] [-f NFAR] [-w] [-s] [-c ICRIT]"
    WRITE(*,*) ""
    WRITE(*,*) "Selects a grid of (NQX x NQY x NQZ) points from the anh* files"
    WRITE(*,*) "Apply the inverse Fourier transform, and saves it to FILEOUT (default: mat3R)."
    WRITE(*,*) "Check the minimum distance up to NFAR supercells away (default: 2)."
    WRITE(*,*) "Setting NFAR=0 will produce 'periodic' force constants. "
    WRITE(*,*)
    WRITE(*,*) "-w : when performing the FFT test, if the re-computed D3 matrix differs"
    WRITE(*,*) "     significantly from the initial one it will be printed to a file."
    WRITE(*,*) "     The file name will start with prefix 'anh_cmplx' if the test was"
    WRITE(*,*) "     performed with complex force constant and 'anh_real' if the test"
    WRITE(*,*) "     was performed with just the real part of the FCs (they should be real)"
    WRITE(*,*) ""
    WRITE(*,*) "-s : skip the test"
    WRITE(*,*) ""
    WRITE(*,*) "-c ICRIT : specify the localization criteria (1=perimeter,"
    WRITE(*,*) "           2=squared perimeter, 3=radius of incribing circle,"
    WRITE(*,*) "           4=distance from baricenter, 5=dfb squared)"
    STOP 1
  ENDIF

  ! cmdline = cmdline_residual()
  ! READ(cmdline, *, iostat=ios) nq, dummy
  ! IF(ios==0 .and. len(dummy)>0) CALL errore("qq2rr", "too many argument use command '-h' for help",1)
  ! READ(cmdline, *, iostat=ios) nq
  ! IF(ios/=0) CALL errore("qq2rr", "missing argument use command '-h' for help",1)
  !
  ! WRITE(*,*) "Number of neighbours to check for BZ", far



  !
  ! WRITE(*,*) "Reading D3 matrices..."
  ! CALL read_d3_matrices(nq, nq_trip, S, d3grid)
  ! WRITE(*,*) "Reading D3 matrices done"

  CALL read_fc2(file2, S, fc2)
  CALL aux_system(S)
  call div_mass_fc2(S, fc2)
  !
  CALL read_fc2(filein, Sd, fc2d)
  CALL aux_system(Sd)
  ! call div_mass_fc2(Sd, fc2d)
  !
  call allocate_fc2_sc(fc2sc, S, fc2%nq)
  fc2sc%fc = fc_sc2RR(fc2%nq, S, Sd, fc2d%fc)
  nq_grid = fc2%nq(1)*fc2%nq(2)*fc2%nq(3)
  nq_trip = nq_grid**2
  ALLOCATE(d3grid(nq_trip))
  call setup_simple_grid(S%bg, fc2%nq(1), fc2%nq(2), fc2%nq(3), q2)
  call setup_simple_grid(S%bg, fc2%nq(1), fc2%nq(2), fc2%nq(3), q3)
  allocate(D(S%nat3,S%nat3))

  do i2 = 1, q2%nq
    do i3 = 1, q2%nq
      comp = i3 + (i2-1)*q2%nq
      call r2q_at_once(fc2sc, q2%xq(:,i2), q3%xq(:,i3), S%nat3, D)
      allocate(d3grid(comp)%D(1,3,3,1,S%nat,S%nat))
      call d3_2idx_2_4idx(S%nat, d3grid(comp)%D, D)
      d3grid(comp)%xq2 = q2%xq(:,i2)
      d3grid(comp)%xq3 = q3%xq(:,i3)
    enddo
  enddo

  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  WRITE(*,*) "Doing Backward FFT..."
  CALL bwfft_d3_interp(fc2%nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fc3, far, icriterium, 1)
  WRITE(*,*) "Backward FFT done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  IF(filename /="none")THEN
    WRITE(*,*) "Writing FCs to file..."
    fc3%nq = fc2%nq
    print*, fc3%nq
    CALL fc3%write(filename, S, .true.)
  ENDIF
  !
  ! IF(.not.skip_test)THEN
  !   WRITE(*,*) "Testing Forward FFT, with imaginary part..."
  !   WRITE(*,*) "(you can stop the code with CTRL-C to avoid running tests)  "
  !   CALL test_fwfft_d3(1, S, d3grid, fc3, .true., write_diff, "anh_cmplx")
  !   WRITE(*,*) "Testing Forward FFT, without imaginary part..."
  !   DEALLOCATE(fc3%ifc)
  !   CALL test_fwfft_d3(nq_trip, S, d3grid, fc3, .false., write_diff, "anh_real")
  !   WRITE(*,*) "Testing forward FFT done"
  ! ENDIF
  !
  call stop_mpi()
END PROGRAM


