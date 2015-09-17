!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE more_constants
  USE kinds, ONLY : DP
  REAL(DP),PARAMETER :: RY_TO_JOULE =  0.5* 4.35974394e-18
  REAL(DP),PARAMETER :: RY_TO_SECOND = 2* 2.418884326505e-17
  REAL(DP),PARAMETER :: RY_TO_METER = 5.2917721092e-11
  CHARACTER(len=3),PARAMETER :: INVALID = '///'
  REAL(DP),PARAMETER :: MASS_DALTON_TO_RY = 0.5_dp*1822.88839_dp
  REAL (DP), PARAMETER :: RY_TO_WATT = RY_TO_JOULE / RY_TO_SECOND
  REAL(DP),PARAMETER :: RY_TO_WATTMM1KM1 = RY_TO_WATT / RY_TO_METER
  
  !
  REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
  REAL(DP),PARAMETER :: eps_freq = 1.e-8_dp

  CONTAINS
  !
!   CHARACTER(len=256) &
!   FUNCTION sigma_file_name(prefix, nq1, nq2, nq3, cT, csigma)
!     IMPLICIT NONE
!     CHARACTER(len=256),INTENT(in) :: prefix
!     INTEGER,INTENT(in) :: nq1, nq2, nq3
!     CHARACTER(len=6),INTENT(in) ::  cT, csigma
!     !
!     CHARACTER (LEN=6), EXTERNAL :: int_to_char
!     !
!     sigma_file_name= TRIM(prefix)//&
!                 "."//TRIM(int_to_char(nq1))// &
!                 "."//TRIM(int_to_char(nq2))// &
!                 "."//TRIM(int_to_char(nq3))// &
!                 "."//TRIM(cT)// &
!                 "."//TRIM(csigma)//".out"
!     !
!   END FUNCTION sigma_file_name
!   !
  ! Write a number from a list using as many digits after the dot as the longest in the list.
  CHARACTER(len=6) &
  FUNCTION write_conf(it,nt,T) RESULT(str)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: it, nt
    REAL(DP),INTENT(in) :: T(nt)
    
    INTEGER :: max_digit_left, max_digit_right, jt
    REAL(DP) :: Tx, Tfrac
    CHARACTER(len=64):: fmt1='', fmt2=''
    
    max_digit_right=0 
    max_digit_left=0 
   
    DO jt = 1,nt      
      max_digit_left = MAX( max_digit_left , CEILING(LOG10(T(jt))) )
      IF(T(jt)>0) THEN
        Tfrac = T(jt)-INT(T(jt))
        IF(Tfrac>0._dp) THEN
          Tx = 1._dp/Tfrac
          max_digit_right = MAX( max_digit_right , CEILING(LOG10(Tx)) )
        !rint*, ">>", jt, T(jt), Tx,  CEILING(LOG10(Tx))
        ENDIF
      ENDIF
    ENDDO
    
    str=""
    WRITE(fmt1,'(i6)') max_digit_right
    fmt2 = "(1f6."//TRIM(ADJUSTL(fmt1))//")"
    !print*, fmt1, fmt2, max_digit_left, max_digit_right
    
    WRITE(str,fmt2) T(it)
    str=ADJUSTL(str)
    
  END FUNCTION write_conf
  !
!   CHARACTER(len=6) &
!   FUNCTION write_sigma(it,nat3,nt,sigma) RESULT(str)
!     USE constants, ONLY : eps12
!     IMPLICIT NONE
!     INTEGER,INTENT(in)  :: it, nt, nat3
!     REAL(DP),INTENT(in) :: sigma(nat3,nt)
!     REAL(DP) :: csigma(nt)
!     INTEGER :: jt, is
!     LOGICAL :: same
!     CHARACTER(len=6),EXTERNAL :: int_to_char
! 
!     same = .true.
!     DO is = 1,nat3
!       same=same.and.(ABS(sigma(1,it)-sigma(is,it))<eps12)
!     ENDDO
!     IF(.not.same) THEN
!         str="X"//TRIM(int_to_char(it))//"."
!         RETURN
!     ENDIF
!     !
!     csigma = 0._dp
!     DO jt = 1,nt
!       same = .true.
!       DO is = 1,nat3
!         same=same.and.(ABS(sigma(1,jt)-sigma(is,jt))<eps12)
!       ENDDO
!       IF(same) csigma(jt) = sigma(1,jt)
!     ENDDO
!     csigma = 2*csigma
!     str = write_temperature(it,nt,csigma)
!     
!   END FUNCTION
  
END MODULE more_constants
!
