!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! #define __OLD_NONCOLIN_GGA
MODULE gc_d3
  !
  USE kinds, ONLY: DP
  !
  USE gc_ph,            ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  !
  REAL(DP), ALLOCATABLE :: &
       dvxc_rrr(:,:,:),          &! dfftp%nnr, nspin, nspin), 
       dvxc_srr(:,:,:),          &! dfftp%nnr, nspin, nspin),
       dvxc_ssr(:,:,:),          &! dfftp%nnr, nspin, nspin), 
       dvxc_sss(:,:,:)            ! dfftp%nnr, nspin, nspin),
  !
  CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE setup_d3gc
    !-----------------------------------------------------------------------
    ! wat20100930:
    ! Allocate and setup all variable needed in the gradient correction case
    !
    !
    !
    USE constants,            ONLY : e2
    USE gvect,                ONLY : ngm, g, nl
    USE lsda_mod,             ONLY : nspin
    USE spin_orb,             ONLY : domag
    USE scf,                  ONLY : rho, rho_core, rhog_core
    USE noncollin_module,     ONLY : noncolin
    USE kinds,                ONLY : DP
    USE funct,                ONLY : dft_is_gradient, gcxc, dgcxc, d3gcxc
!     USE gc_ph,            ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  !   USE gc_d3,            ONLY : dvxc_rrr, dvxc_srr, &
  !                                dvxc_ssr, dvxc_sss
    USE nlcc_ph,              ONLY : nlcc_any
    USE fft_base,             ONLY : dfftp
    USE fft_interfaces,       ONLY : fwfft
    USE wavefunctions_module, ONLY : psic
    USE io_global,            ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER :: ir, is, nspin0
    REAL(DP) :: grho2 (2), fac, &
        sx, sc, v1x, v2x, v1c,v2c, &
        vrrx, vsrx, vssx, vrrc, vsrc, vssc, &
        dvxcrrr, dvxcsrr, dvxcssr, dvxcsss, &
        vrrrx, vsrrx, vssrx, vsssx, &
        vrrrc, vsrrc, vssrc, vsssc
    REAL(DP), ALLOCATABLE :: rho_tot_r(:,:)
    COMPLEX(DP), ALLOCATABLE :: rho_tot_g(:,:)
    REAL (DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10

    !
    WRITE(stdout, '(5x,a)') "Setting up GGA 2nd derivative"
    grho2 = 0._dp
    !
    IF ( .NOT. dft_is_gradient() ) RETURN
    
    nspin0=nspin
    !
    IF ( .NOT. nspin0 == 1) &
        CALL errore ('setup_d3gc', ' Gradient corrections implemented for nspin = 1 only ! ', 1 )
    !
    ALLOCATE (grho    (  3    , dfftp%nnr   , nspin0))
    ALLOCATE (rho_tot_r  (  dfftp%nnr , nspin0))
    ALLOCATE (rho_tot_g  (  ngm , nspin0))
    ALLOCATE (dvxc_rr (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_sr (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_ss (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_s  (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_rrr(  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_srr(  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_ssr(  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_sss(  dfftp%nnr , nspin0 , nspin0))
    !
    grho    (:,:,:) = 0.d0
    dvxc_s  (:,:,:) = 0.d0
    dvxc_rr (:,:,:) = 0.d0
    dvxc_sr (:,:,:) = 0.d0
    dvxc_ss (:,:,:) = 0.d0
    dvxc_rrr(:,:,:) = 0.d0
    dvxc_srr(:,:,:) = 0.d0
    dvxc_ssr(:,:,:) = 0.d0
    dvxc_sss(:,:,:) = 0.d0
    !
    fac = 1.d0 / DBLE (nspin0)
    IF (noncolin.AND.domag) THEN
      call errore('setup_d3gc',' domag not implemented',1)
    ELSE
      !
      IF (.NOT. nlcc_any) THEN
        !
        DO is = 1, nspin0
          rho_tot_r(:,is) = rho%of_r(:,is)
!           rho_tot_g(:,is) = rho%of_g(:,is)
        ENDDO
        !
      ELSE
        !
        DO is = 1, nspin0
          rho_tot_r(:,is) = fac * rho_core(:)  + rho%of_r(:,is)
!           rho_tot_g(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
        ENDDO
        !
      ENDIF
      !
      !
      DO is = 1, nspin0
        psic(:) = rho_tot_r(:,is)
        CALL fwfft ('Dense', psic, dfftp)
        rho_tot_g(:,is) = psic(nl(:))

!           CALL gradrho (nrx1, nrx2, nrx3, nr1, nr2, nr3, dfftp%nnr, rho_tot_g (1, is), &
!               ngm, g, nl, grho (1, 1, is) )          
          CALL gradrho( dfftp%nnr, rho_tot_g(1,is), ngm, g, nl, grho(:,:,is) )
      ENDDO
      !
    END IF
    
    !WHERE(ABS(grho)>1.d+32) grho = 0._dp
    
    DO ir = 1, dfftp%nnr
      !print*, rho%of_r(ir,1), grho (1,ir,1), grho(2,ir,1), grho(3,ir,1)
      grho2(1) = grho (1,ir,1)**2 + grho(2,ir,1)**2 + grho(3,ir,1)**2
      IF (nspin0 == 1) THEN
          IF (ABS (rho_tot_r (ir, 1) ) > epsr .AND. grho2 (1) > epsg) THEN
            call gcxc  (rho_tot_r (ir, 1), grho2(1), sx, sc, v1x, v2x, v1c, v2c)
            call dgcxc (rho_tot_r (ir, 1), grho2(1), vrrx, vsrx, vssx, vrrc, &
                  vsrc, vssc)
            dvxc_rr (ir, 1, 1) = e2 * (vrrx + vrrc)
            dvxc_sr (ir, 1, 1) = e2 * (vsrx + vsrc)
            dvxc_ss (ir, 1, 1) = e2 * (vssx + vssc)
            dvxc_s  (ir, 1, 1) = e2 * (v2x + v2c)
            call d3gcxc (rho_tot_r (ir, 1), grho2(1), vrrrx, vsrrx, vssrx, vsssx, &
                  vrrrc, vsrrc, vssrc, vsssc )
!             write(10001, '(i7,99f12.6)') ir, rho_tot_r(ir, 1),  rho%of_r(ir,is), rho_core(ir), grho2(1), vrrrx, vsrrx, vssrx, vsssx, vrrrc, vsrrc, vssrc, vsssc 
            !
            dvxc_rrr(ir, 1, 1) = e2 * (vrrrx + vrrrc)
            dvxc_srr(ir, 1, 1) = e2 * (vsrrx + vsrrc)
            dvxc_ssr(ir, 1, 1) = e2 * (vssrx + vssrc)
            dvxc_sss(ir, 1, 1) = e2 * (vsssx + vsssc)
            !
          ENDIF
      ELSE
          CALL errore('setup_d3gc',' nspin>1 not implemented',1)
      ENDIF
      !
    ENDDO
    
!     WRITE(10006, '(2f12.6)') dvxc_rrr
!     WRITE(10007, '(2f12.6)') dvxc_srr
!     WRITE(10008, '(2f12.6)') dvxc_ssr
!     WRITE(10009, '(2f12.6)') dvxc_sss
    
    IF (noncolin.AND.domag) &
      CALL errore('setup_d3gc',' domag not implemented',1)
    !
    DEALLOCATE(rho_tot_r, rho_tot_g)
    !
    RETURN
    !
  END SUBROUTINE setup_d3gc

END MODULE gc_d3