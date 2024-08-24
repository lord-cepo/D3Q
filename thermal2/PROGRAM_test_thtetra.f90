

PROGRAM test_thtetra
   USE thtetra, ONLY: tetra_init, tetra_weights_delta, opt_tetra_init, deallocate_tetra, tetra_weights_theta
   USE kinds,   ONLY: DP

   IMPLICIT NONE
   INTEGER :: idx, i,j,k, n2
   INTEGER, PARAMETER :: n = 80
   REAL(DP) :: xq(3, n**3), eq(1,n**3), bg(3,3), wq(1,n**3), radius


   bg = 0._dp
   bg(1,1) = 1._dp
   bg(2,2) = 1._dp
   bg(3,3) = 1._dp

   idx = 0
   n2 = n/2
   DO i = -n2, n2-1
      DO j = -n2, n2-1
         DO k = -n2, n2-1
            !
            idx = idx+1
            xq(1,idx) = REAL(i,kind=DP)/REAL(n,kind=DP)
            xq(2,idx) = REAL(j,kind=DP)/REAL(n,kind=DP)
            xq(3,idx) = REAL(k,kind=DP)/REAL(n,kind=DP)
            eq(1,idx) = NORM2(xq(:,idx) - (/0.0178_dp,0.00416_dp,0.00251_dp/)) 
            !
         ENDDO
      ENDDO
   ENDDO

   radius = 0.413_dp
   ! WRITE(*,*) xq(:,197), eq(1,197)
   WRITE(*,*) "-----------------SURFACE --------------------"
   WRITE(*,*) "expected", 4._dp*3.141592653589793*radius**2
   CALL opt_tetra_init((/n,n,n/), bg)
   wq = tetra_weights_delta(n**3, 1, eq, radius)
   WRITE(*,*) "optimized tetra", SUM(wq), "fraction", SUM(wq)/(4._dp*3.141592653589793*radius**2)

   CALL deallocate_tetra()
   CALL tetra_init((/n,n,n/), bg)
   wq = tetra_weights_delta(n**3, 1, eq, radius)
   WRITE(*,*) "linear tetra", SUM(wq)
   WRITE(*,*) "------- VOLUME ----------"
   WRITE(*,*) "expected", 4._dp/3._dp*3.141592653589793*radius**3
   wq = tetra_weights_theta(n**3, 1, eq, radius)
   WRITE(*,*) "linear tetra", SUM(wq)
   
   
   radius = 0.254_dp
   WRITE(*,*) "-----------------SURFACE --------------------"
   ! WRITE(*,*) xq(:,197), eq(1,197)
   WRITE(*,*) "expected", 4._dp*3.141592653589793*radius**2
   CALL opt_tetra_init((/n,n,n/), bg)
   wq = tetra_weights_delta(n**3, 1, eq, radius)
   WRITE(*,*) "optimized tetra", SUM(wq), "fraction", SUM(wq)/(4._dp*3.141592653589793*radius**2)

   CALL deallocate_tetra()
   CALL tetra_init((/n,n,n/), bg)
   wq = tetra_weights_delta(n**3, 1, eq, radius)
   WRITE(*,*) "linear tetra", SUM(wq)
   WRITE(*,*) "------- VOLUME ----------"
   WRITE(*,*) "expected", 4._dp/3._dp*3.141592653589793*radius**3
   wq = tetra_weights_theta(n**3, 1, eq, radius)
   WRITE(*,*) "linear tetra", SUM(wq)

   WRITE(*,*) "---------------------------------------"
   




END PROGRAM test_thtetra
