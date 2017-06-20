!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Code contributions from Giorgia Fugallo and Michele Lazzeri
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM qq2rr
  USE kinds,           ONLY : DP
  USE input_fc,        ONLY : ph_system_info
  USE iso_c_binding,   ONLY : c_int
  USE fc3_interpolate, ONLY : grid
  USE f3_bwfft,        ONLY : d3_list, read_d3_matrices, bwfft_d3_interp, test_fwfft_d3
  USE cmdline_param_module
  
  IMPLICIT NONE
  TYPE(d3_list),ALLOCATABLE ::  d3grid(:)
  TYPE(grid)                :: fc3
  INTEGER :: nq(3)
  INTEGER :: nq_trip, nq_grid
  INTEGER(kind=c_int)    :: kb
  !
  TYPE(ph_system_info) :: S
  COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:)
  !
  INTEGER :: far, ios
  CHARACTER(len=512) :: argv, filename
  CHARACTER(len=:),ALLOCATABLE :: cmdline

  filename = cmdline_param_char("o", "mat3R")
  far      = cmdline_param_int("f", 2)
  IF (cmdline_param_logical('h')) THEN
      WRITE(*,*) "Syntax: ls anh*| d3_qq2rr.x NQX NQY NQZ [-o FILEOUT] [-f NFAR]"
      WRITE(*,*) ""
      WRITE(*,*) "Selects a grid of (NQX x NQY x NQZ) points from the anh* files"
      WRITE(*,*) "Apply the inverse Fourier transform, and saves it to FILEOUT (default: mat3R)."
      WRITE(*,*) "Check for shortes perimeter up to NFAR unit cells away (default: 2)."
      STOP 1
  ENDIF
  
  cmdline = cmdline_residual()
  READ(cmdline, *, iostat=ios) nq
  IF(ios/=0) CALL errore("import3py", "missing argument use command '-h' for help",1)
  !
  WRITE(*,*) "Reading grid", nq
  WRITE(*,*) "Number of neighbours to check for BZ", far
  
  nq_grid = nq(1)*nq(2)*nq(3)
  nq_trip = nq_grid**2
  ALLOCATE(d3grid(nq_trip))
  !
  WRITE(*,*) "Reading D3 matrices..."
  CALL read_d3_matrices(nq, nq_trip, S, d3grid)
  WRITE(*,*) "Reading D3 matrices done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  WRITE(*,*) "Doing Backward FFT..."
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fc3, far)
  fc3%nq = nq
  WRITE(*,*) "Backward FFT done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  IF(filename /="none")THEN
    WRITE(*,*) "Writing FCs to file..."
    CALL fc3%write(filename, S)
  ENDIF
  !
  WRITE(*,*) "Testing Forward FFT, with imaginary part..."
  WRITE(*,*) "(you can stop the code with CTRL-C to avoid running tests)  "
  CALL test_fwfft_d3(nq_trip, S, d3grid, fc3)
  WRITE(*,*) "Testing Forward FFT, without imaginary part..."
  DEALLOCATE(fc3%ifc)
  CALL test_fwfft_d3(nq_trip, S, d3grid, fc3)
  WRITE(*,*) "Testing forward FFT done"
  !
END PROGRAM qq2rr


