PROGRAM test_tetra
   USE new_ktetra,  ONLY: &
   opt_tetra_init, &
   deallocate_tetra, &
   tetra_weights_only, &
   opt_tetra_weights_only, &
   tetra_type, &
   wlsm, &
   tetra, &
   opt_tetra_dos_t
   USE thtetra, ONLY: &
   my_opt_tetra_init => opt_tetra_init, &
   my_tetra_init => tetra_init, &
   tetra_weights_theta, &
   tetra_weights_delta, &
   my_deallocate_tetra => deallocate_tetra, &
   my_wlsm => wlsm, &
   my_tetra => tetra
   USE kinds,   ONLY: DP

   IMPLICIT NONE
   INTEGER, PARAMETER :: n = 30
   INTEGER :: s(3,3,48), idx, i, j, k, nk1, nk2, nk3, isk(n**3), t_rev(48)
   REAL(DP) :: at(3,3), xq(3,n**3), eq(1,n**3), radius, wq(1,n**3), dost(2)

   isk = 0

   nk1 = n
   nk2 = n
   nk3 = n
   DO i = 1, nk1
      DO j = 1, nk2
         DO k = 1, nk3
            !  this is nothing but consecutive ordering
            idx = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
            !  xkg are the components of the complete grid in crystal axis
            xq(1,idx) = DBLE(i-1)/nk1
            xq(2,idx) = DBLE(j-1)/nk2
            xq(3,idx) = DBLE(k-1)/nk3
            eq(1,idx) = NORM2(xq(:,idx) - (/0.5178_dp,0.50416_dp,0.50251_dp/))
         ENDDO
      ENDDO
   ENDDO

   s = 0
   s(1,1,1) = 1
   s(2,2,1) = 1
   s(3,3,1) = 1

   at = REAL(s(:,:,1), DP)

   t_rev = 0

   ! WRITE(*,*) "before init"
   radius = 0.413_dp
   
   CALL my_opt_tetra_init((/n,n,n/), at)
   

   ! wlsm e tetra sono le stesse, il problema non sembra essere in opt_tetra_init
   ! tetra_type = 2
   ! CALL opt_tetra_init(1, s, .FALSE., t_rev, at, at, n**3, 0,0,0, n,n,n, n**3, xq, 1)
   ! WRITE(*,*) ALL(my_wlsm == wlsm)
   ! WRITE(*,*) ALL(my_tetra == tetra)
   
   WRITE(*,*) "------- SURFACE ----------"
   WRITE(*,*) "expected", 4._dp*3.141592653589793*radius**2
   ! tetra_old = tetra
   CALL my_tetra_init((/n,n,n/), at)
   wq = tetra_weights_delta(n**3, 1, eq, radius)
   WRITE(*,*) "my opt", SUM(wq)
   
   ! opt_dos_t funziona, a parte un fattore 2
   ! CALL opt_tetra_dos_t(eq, 1, 1, n**3, radius, dost)
   ! WRITE(*,*) "opt_dos_t", dost(1)

   
   ! WRITE(*,*) "------- VOLUME ----------"
   ! WRITE(*,*) "expected", 4._dp/3._dp*3.141592653589793*radius**3

   ! CALL my_tetra_init((/n,n,n/), at)
   ! wq = tetra_weights_theta(n**3, 1, eq, radius)
   ! WRITE(*,*) "tetra_weights_theta", SUM(wq)
   ! tetra_type = 1
   ! CALL opt_tetra_init(1, s, .FALSE., t_rev, at, at, n**3, 0,0,0, n,n,n, n**3, xq, 1)
   ! CALL tetra_weights_only(n**3, 0, 0, isk, 1, radius, eq, radius, wq)
   ! CALL deallocate_tetra()
   ! WRITE(*,*) "tetra_weights_only", SUM(wq)
   ! tetra_type = 2
   ! CALL opt_tetra_init(1, s, .FALSE., t_rev, at, at, n**3, 0,0,0, n,n,n, n**3, xq, 1)
   ! CALL opt_tetra_weights_only(n**3, 0, 1, eq, radius, wq, 0, isk)
   ! CALL deallocate_tetra()
   ! WRITE(*,*) "opt_tetra_weights_only", SUM(wq)

END PROGRAM test_tetra
