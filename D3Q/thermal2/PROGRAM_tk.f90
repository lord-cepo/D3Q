!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE thermalk_program
  !
#include "mpi_thermal.h"
  USE kinds,       ONLY : DP
  USE mpi_thermal, ONLY : ionode
  USE posix_signal,ONLY : check_graceful_termination
  USE timers
  !
  CONTAINS
  ! 
  SUBROUTINE check_negative_lw(lw, nat3, nconf, name)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: lw(nat3, nconf)
    INTEGER,INTENT(in)  :: nat3, nconf
    CHARACTER(len=*),INTENT(in) :: name
    !
    INTEGER :: i,j
    IF(ANY(lw<0._dp))THEN
      DO i = 1,nconf
        DO j = 1,nat3
          IF(lw(j,i) < 0._dp) THEN
            WRITE(*,*) i, j, lw(j,i)
          ENDIF
        ENDDO
      ENDDO
      CALL errore(name, "negative linewidth", 1)
    ENDIF
    RETURN
  END SUBROUTINE
  !
  ! This subroutine computes the SMA thermal conducivity, it is mainly just a driver
  ! that uses other subroutines to compute th intrinsic, isotopic and casimir linewidths,
  ! than it sums everything up and takes care of input/output.
  !
  ! This subroutine is obsoleted by the first iteration of the variational method, 
  ! but we keep it for didactical purposes
  SUBROUTINE TK_SMA(input, out_grid, S, fc2, fc3)
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1, K_BOLTZMANN_RY
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf
    USE q_grids,            ONLY : q_grid, setup_grid !, setup_bz_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel
    USE input_fc,           ONLY : ph_system_info
    USE code_input,         ONLY : code_input_type
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: out_grid
    !
    TYPE(q_grid) :: in_grid
    REAL(DP) :: sigma_ry(input%nconf)
    REAL(DP) :: lw(S%nat3,input%nconf)
    REAL(DP) :: lw_isotopic(S%nat3,input%nconf)
    REAL(DP) :: lw_phph(S%nat3,input%nconf)
    REAL(DP) :: lw_casimir(S%nat3)
    REAL(DP) :: freq(S%nat3)
    REAL(DP) :: bose(S%nat3,input%nconf)
    COMPLEX(DP) :: U(S%nat3,S%nat3)
    !
    REAL(DP) :: tk(3,3,input%nconf)
    REAL(DP) :: vel(3,S%nat3)
    REAL(DP) :: pref
    INTEGER  :: iq, it, a, b, nu
    !
    REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    ! the inner grid (in_grid) is scatterd over MPI
    CALL setup_grid(input%grid_type, S%bg, input%nk_in(1), input%nk_in(2), input%nk_in(3),&
                    in_grid, scatter=.true.)
    !CALL in_grid%scatter()
    !
!     ioWRITE(stdout,'(1x,a,i10,a)') "Integrating over an inner grid of", in_grid%nq, " points"
!     ioWRITE(stdout,'(1x,a,i10,a)') "Integrating over an outer grid of", out_grid%nq, " points"
    !
    IF(ionode)THEN
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              "lw."//TRIM(input%prefix)//&
                               "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                               "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
      ioWRITE(1000+it, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
      ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
      ioFLUSH(1000+it)
      !
      IF(input%isotopic_disorder) THEN
        OPEN(unit=2000+it, file=TRIM(input%outdir)//"/"//&
                                "lwiso."//TRIM(input%prefix)//&
                                    "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                    "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
        ioWRITE(2000+it, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
        ioWRITE(2000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
        ioFLUSH(2000+it)
      ENDIF
    ENDDO
      IF(input%casimir_scattering) THEN
        OPEN(unit=3000, file=TRIM(input%outdir)//"/"//&
                                "lwcas."//TRIM(input%prefix)//".out")
        ioWRITE(3000, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
!         ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
        ioFLUSH(3000)
      ENDIF
    ENDIF
    !
    tk = 0._dp
    !dq = S%Omega/out_grid%nq
    !
      timer_CALL t_tksma%start()
    QPOINT_LOOP : &
    DO iq = 1,out_grid%nq
      ioWRITE(stdout,'(i6,3f15.8)') iq, out_grid%xq(:,iq)
      !
        timer_CALL t_lwphph%start()
      lw_phph = linewidth_q(out_grid%xq(:,iq), input%nconf, input%T,&
                       sigma_ry, S, in_grid, fc2, fc3)
      CALL check_negative_lw(lw_phph, S%nat3, input%nconf, "SMA:phph")
        timer_CALL t_lwphph%stop()
      !
      ! Compute contribution of isotopic disorder
      IF(input%isotopic_disorder)THEN
          timer_CALL t_lwisot%start()
        lw_isotopic = isotopic_linewidth_q(out_grid%xq(:,iq), input%nconf, input%T, &
                                           sigma_ry, S, out_grid, fc2)
        CALL check_negative_lw(lw_isotopic, S%nat3, input%nconf, "SMA:isotopic")
          timer_CALL t_lwisot%stop()
      ELSE
        lw_isotopic = 0._dp
      ENDIF
      !
      !
        timer_CALL t_velcty%start() 
      ! Velocity
      vel = velocity(S, fc2, out_grid%xq(:,iq))
      CALL  freq_phq_safe(out_grid%xq(:,iq), S, fc2, freq, U)
        timer_CALL t_velcty%stop() 
      !
      ! Compute anisotropic Casimir linewidth
      IF(input%casimir_scattering) THEN
          timer_CALL t_lwcasi%start() 
        lw_casimir = casimir_linewidth_vel(vel, input%casimir_length, input%casimir_dir, S%nat3)
          timer_CALL t_lwcasi%stop() 
      ELSE
        lw_casimir = 0._dp
      ENDIF
      !
        timer_CALL t_lwinout%start()
      DO it = 1, input%nconf
        ioWRITE(1000+it,'(3f12.6,99e20.10)') out_grid%xq(:,iq), lw_phph(:,it)*RY_TO_CMM1
        IF(input%isotopic_disorder) THEN
          ioWRITE(2000+it,'(3f12.6,99e20.10)') out_grid%xq(:,iq), lw_isotopic(:,it)*RY_TO_CMM1
        ENDIF
      ENDDO
      IF(input%casimir_scattering) THEN
        ioWRITE(3000,'(3f12.6,99e20.10)') out_grid%xq(:,iq), lw_casimir(:)*RY_TO_CMM1
      ENDIF
        timer_CALL t_lwinout%stop()
      !
        timer_CALL t_tksum%start()
      ! Casimir linewidth is temperature/smearing-independent, sum it to all configurations
      DO it = 1,input%nconf
        lw(:,it) = lw_phph(:,it) + lw_isotopic(:,it) + lw_casimir
      ENDDO
      !
      CONF_LOOP : &
      DO it = 1, input%nconf
        !
        IF(input%T(it)==0._dp) CYCLE CONF_LOOP
        !
        CALL  bose_phq(input%T(it), S%nat3, freq, bose(:,it))
        !
        MODE_LOOP : &
        DO nu = 1, S%nat3
          ! Check if we have zero linewidth and non-zero velocity it is a problem
          ! lw can be NaN when T=0 and xq=0, check for lw>0 insteand, because NaN/=0 is true
          IF(lw(nu,it)<0._dp)THEN ! true for NaN
            WRITE(stdout,"(3x,a,e12.4,3i6)") "WARNING! Negative lw (idx q, mode, conf):", lw(nu,it), iq, nu, it 
            lw(nu,it) = - lw(nu,it)
          ENDIF
          IF(.not. lw(nu,it)>0._dp)THEN ! false for NaN
            IF(ANY(ABS(vel(:,nu))>eps_vel ))THEN
              WRITE(stdout,'(3i6,1e20.10,5x,3e20.10)') iq, nu, it, lw(nu,it), vel(:,nu)
              CALL errore("TK_SMA", "cannot treat this case", 1)
            ELSE
              !ioWRITE(stdout,"(3x,a,3i6)") "skip (iq,nu,it):", iq, nu, it
              CYCLE MODE_LOOP 
            ENDIF
          ENDIF
          !
          pref = freq(nu)**2 *bose(nu,it)*(1+bose(nu,it))&
                             /(input%T(it)**2 *S%Omega*K_BOLTZMANN_RY)&
                             *out_grid%w(iq) /lw(nu,it) 
          !ioWRITE(stdout,"(3x,a,3i6,4e15.6)") "do:", iq, nu, it, pref, freq(nu), bose(nu,it), lw(nu,it)
          DO a = 1,3
          DO b = 1,3
            tk(a,b,it) = tk(a,b,it) + pref*vel(a,nu)*vel(b,nu)
          ENDDO
          ENDDO
        ENDDO MODE_LOOP 
        !
      ENDDO CONF_LOOP 
        timer_CALL t_tksum%stop()
      !
    ENDDO QPOINT_LOOP
      timer_CALL t_tksma%stop()

    IF(ionode)THEN      
    DO it = 1, input%nconf
      CLOSE(1000+it)
      IF(input%isotopic_disorder) CLOSE(2000+it)
    ENDDO
    IF(input%casimir_scattering) CLOSE(3000)
    ENDIF
    !
    ! Write to disk
    IF(ionode) OPEN(unit=10000, file=TRIM(input%outdir)//"/"//&
                                     TRIM(input%prefix)//"."//"out")
    ioWRITE(10000,'(4a)') "#conf  sigma[cmm1]   T[K]  ",&
                        "    K_x            K_y            K_z             ",&
                        "    K_xy           K_xz           K_yz            ",&
                        "    K_yx           K_zx           K_zy       "
    tk = tk*RY_TO_WATTMM1KM1
    DO it = 1,input%nconf
      ioWRITE(10000,"(i3,2f12.6,3(3e15.6,5x))") it, input%sigma(it), input%T(it), &
      tk(1,1,it),tk(2,2,it),tk(3,3,it), &
      tk(1,2,it),tk(1,3,it),tk(2,3,it), &
      tk(2,1,it),tk(3,1,it),tk(3,2,it)
    ENDDO
    IF(ionode) CLOSE(10000)
    !
    ! Write to screen
    ioWRITE(stdout,"(3x,a,/,3x,a)") "************", "SMA thermal conductivity, stored to file:"
    ioWRITE(stdout,'(5x,a)') TRIM(input%outdir)//"/"//TRIM(input%prefix)//"."//"out"
    ioWRITE(stdout,"(3x,a)") "Diagonal components (conf, sigma, T, K_x, K_y, K_z):"
    DO it = 1,input%nconf
      ioWRITE(stdout,"(i3,2f12.6,3e16.8)")  it, input%sigma(it), input%T(it),&
                                      tk(1,1,it), tk(2,2,it), tk(3,3,it)
    ENDDO
    !
#ifdef timer_CALL
    ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
    CALL print_timers_header()
    CALL t_tksma%print()
    ioWRITE(stdout,'(a)') "*** * Contributions to SMA conductivity:"
    CALL t_tksum%print()
    CALL t_lwisot%print()
    CALL t_lwcasi%print()
    CALL t_lwphph%print()
    CALL t_velcty%print()
    CALL t_lwinout%print()
    CALL t_mpicom%print()
    CALL t_readdt%print()
    ioWRITE(stdout,'(a)') "*** * Contributions to ph-ph linewidth time:"
    CALL t_freq%print()
    CALL t_bose%print()
    CALL t_sum%print()
    CALL t_fc3int%print()
    CALL t_fc3m2%print()
    CALL t_fc3rot%print()
    CALL t_merged%print()
#endif
    !
    !
  END SUBROUTINE TK_SMA
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-! 
  SUBROUTINE TK_CG_prec(input, out_grid, S, fc2, fc3)
    USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_q
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE q_grids,            ONLY : q_grid, q_basis, setup_grid, setup_bz_grid, &
                                   prepare_q_basis, qbasis_dot, qbasis_ax, &
                                   qbasis_a_over_b
    USE variational_tk
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: out_grid
    !
    INTEGER :: ix, nu, iq, it, nu0, iter
    !
    TYPE(q_grid)         :: in_grid ! inner grid is MPI-scattered it is used 
                                    ! for integrating the ph-ph scattering terms and linewidth
    !
    TYPE(q_basis) :: qbasis
    REAL(DP),ALLOCATABLE :: A_out(:,:,:), inv_sqrt_A_out(:,:,:), inv_A_out(:,:,:)
    REAL(DP),ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:), h(:,:,:,:), t(:,:,:,:), j(:,:,:,:)
    REAL(DP),ALLOCATABLE :: Af(:,:,:,:)
    REAL(DP),ALLOCATABLE :: g_dot_h(:,:), h_dot_t(:,:), &
                            g_mod2(:,:), g_mod2_old(:,:), &
                            pref(:,:)
    REAL(DP) :: tk(3,3,input%nconf), delta_tk(3,3,input%nconf), &
                tk_old(3,3,input%nconf)
    LOGICAL :: conv
    CHARACTER(len=6) :: what
    CHARACTER(len=256) :: filename
    !
    ! For code readability:
    LOGICAL :: restart_ok
    INTEGER :: nconf, nat3, nq
    nconf = input%nconf
    nat3  = S%nat3
    nq    = out_grid%nq
    
        timer_CALL t_tkprec%start()
    ! make the inner grid on top of the outer one, they are actually identical
    ! but the inner one is MPI-scattered, so we need to separate onjects to hold them
    CALL setup_grid(input%grid_type, S%bg, out_grid%n(1),out_grid%n(2),out_grid%n(3), &
                    in_grid, scatter=.true.)
    CALL prepare_q_basis(out_grid, qbasis, nconf, input%T, S, fc2)
    
    OPEN_FILES : &
    IF(ionode)THEN
      DO it = 0,nconf
        IF(it==0) THEN
          filename=TRIM(input%outdir)//"/"//&
                  TRIM(input%prefix)//".out"
          what="# conf"
        ELSE
          filename=TRIM(input%outdir)//"/"//&
                  TRIM(input%prefix)//"_T"//TRIM(write_conf(it,nconf,input%T))//&
                    "_s"//TRIM(write_conf(it,nconf,input%sigma))//".out"
          what="# iter"
        ENDIF
        !
        IF(input%restart) THEN
          OPEN(unit=10000+it, file=filename, position="append")
        ELSE
          OPEN(unit=10000+it, file=filename)
          ioWRITE(10000+it, '(a)') "# Thermal conductivity from BTE"
          IF(it>0) THEN
            ioWRITE(10000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, &
                    "     T=",input%T(it), "    sigma=", input%sigma(it)
          ENDIF
          ioWRITE(10000+it,'(5a)') what, " sigma[cmm1]   T[K]  ",&
                              "    K_x              K_y              K_z              ",&
                              "    K_xy             K_xz             K_yz             ",&
                              "    K_yx             K_zx             K_zy"
        ENDIF
        ioFLUSH(10000+it)
      ENDDO
    ENDIF OPEN_FILES
    !
    ALLOCATE(A_out(nconf, nat3, nq))
    ALLOCATE(inv_sqrt_A_out(nconf, nat3, nq))
    ALLOCATE(f(3, nconf, nat3, nq))
    ALLOCATE(g(3, nconf, nat3, nq))
    ALLOCATE(h(3, nconf, nat3, nq))
    ALLOCATE(t(3, nconf, nat3, nq))
    ALLOCATE(   g_dot_h(3, nconf) )
    ALLOCATE(   h_dot_t(3, nconf) )
    ALLOCATE(    g_mod2(3, nconf) )
    ALLOCATE(g_mod2_old(3, nconf) )
    ALLOCATE(      pref(3, nconf) )
      timer_CALL t_tkprec%stop()
    !
    ! Try if restart file can be read (will return false if input%restart is false):
      timer_CALL t_restart%start()
    restart_ok = read_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq)
      timer_CALL t_restart%stop()
    !
    IF(.not. restart_ok) THEN
      ! Compute A_out diagonal matrix
        timer_CALL t_tkaout%start()
      CALL compute_A_out(A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
          timer_CALL t_tkaout%stop()
    ENDIF
    
    ! Compute A_out^(-1/2) and A_out^-1 from A_out
          timer_CALL t_tkprec%start()
    CALL compute_inv_sqrt_A_out(A_out, inv_sqrt_A_out, nconf, nat3, nq)
          timer_CALL t_tkprec%stop()
    ! Go to the reduced variables:
    ! \tilde{b} = A_out^(-1/2) b 
    qbasis%b = A_diag_f(inv_sqrt_A_out, qbasis%b, nconf, nat3, nq)
    !      
    !
    ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") "iter ", 0
    ioWRITE(stdout,'(5a)') "       ", " sigma[cmm1]   T[K]  ",&
                          "    K_x              K_y              K_z              "
    !
    IF(.not. restart_ok) THEN
      ! \tilde{f0} = A_out^(1/2) f0 = A_out^(-1/2) b = \tilde{b}
      f = qbasis%b
      ! "tilde" label dropped from here on
      !
      ! Compute SMA tk = lambda f0.b
          timer_CALL t_tkprec%start()
      g = 0._dp
      tk_old = calc_tk_gf(g, f, qbasis%b, input%T, out_grid%w, S%omega, nconf, nat3, nq)
      CALL print_tk(tk_old, input%sigma, input%T, nconf, "SMA TK", 10000, -1)
          timer_CALL t_tkprec%stop()
      !
      ! Compute g0 = Af0 - b
        timer_CALL t_tkain%start()
      CALL tilde_A_times_f(f, g, inv_sqrt_A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
        timer_CALL t_tkain%stop()
      !
        timer_CALL t_tkprec%start()
      g =  g - qbasis%b
      !
      ! Trivially set h0 = -g0
      h = -g
        timer_CALL t_tkprec%stop()
      !
        timer_CALL t_restart%start()
      CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq)
        timer_CALL t_restart%stop()
    ENDIF
    !
    ! Compute initial variational tk0 =-\lambda (f0.g0-f0.b)
    ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") "iter ", 0
    tk = calc_tk_gf(g, f, qbasis%b, input%T, out_grid%w, S%omega, nconf, nat3, nq)
    CALL print_tk(tk, input%sigma, input%T, nconf, "TK from 1/2(fg-fb)", 10000, 0)
    delta_tk = tk-tk_old
    CALL print_tk(delta_tk, input%sigma, input%T, nconf, "Delta TK - initial")
    !
    ! Check initial gradient
    conv = check_conv_tk(input%thr_tk, nconf, delta_tk)
    g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
    !
    INSTANT_CONVERGENCE : &
    IF(conv) THEN
      ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
                    "Converged at first iteration"
    ELSE INSTANT_CONVERGENCE
      CGP_ITERATION : &
      DO iter = 1,input%niter_max
        CALL check_graceful_termination()
        ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") "iter ", iter
        ! t = Ah = [A_out^(-1/2) (1+A_in) A_out^(-1/2)] h
          timer_CALL t_tkain%start()
        CALL tilde_A_times_f(h, t, inv_sqrt_A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
          timer_CALL t_tkain%stop()
        !
          timer_CALL t_tkcg%start()
        ! Do a CG step:
        ! f_(i+1) = f_i - (g_i.h_i / h_i.t_i) h_i
        ! g_(i+1) = g_i - (g_i.h_i / h_i.t_i) t_i
        g_dot_h = qbasis_dot(g, h,  nconf, nat3, nq )   ! g.h
        h_dot_t = qbasis_dot(h, t, nconf, nat3, nq )    ! h.t
        pref = qbasis_a_over_b(g_dot_h, h_dot_t, nconf) ! g.h/h.t
        f = f - qbasis_ax(pref, h, nconf, nat3, nq)
        g = g - qbasis_ax(pref, t, nconf, nat3, nq)
        ! 
        !In case you want to compute explicitly the gradient (i.e. for testing):
  !       CALL tilde_A_times_f(f, j, inv_sqrt_A_out, input, qbasis, &
  !                            out_grid, in_grid, S, fc2, fc3)
  !       j = j-qbasis%b
        !
        tk_old = tk
        !
        !tk = -\lambda (f.g-f.b)
        tk = calc_tk_gf(g, f, qbasis%b, input%T, out_grid%w, S%omega, nconf, nat3, nq)
        CALL print_tk(tk, input%sigma, input%T, nconf, "TK from 1/2(fg-fb)", 10000, iter)
        ! also compute the variation of tk and check for convergence
        delta_tk = tk-tk_old
        CALL print_tk(delta_tk, input%sigma, input%T, nconf, "Delta TK")
        conv = check_conv_tk(input%thr_tk, nconf, delta_tk)
        !
        IF(conv) THEN
          ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
                        "Convergence achieved"
          CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq)
          EXIT CGP_ITERATION
        ENDIF
        !
        ! h_(i+1) = -g_(i+1) + (g_(i+1).g_(i+1)) / (g_i.g_i) h_i
        g_mod2_old = g_mod2
        g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
        !
        ! Compute the new conjugate gradient: 
        ! h_(i+1) = (g_(i+1).g_(i+1) / g_i.g_i) h_i - g_(i+1)
        pref = qbasis_a_over_b(g_mod2, g_mod2_old, nconf)
        h = -g + qbasis_ax(pref, h, nconf, nat3, nq)
          timer_CALL t_tkcg%stop()
        !
        ! Store to file for restart
          timer_CALL t_restart%start()
        CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq)
          timer_CALL t_restart%stop()
        !
      ENDDO &
      CGP_ITERATION
      !
      IF(iter>input%niter_max) THEN
        ioWRITE (stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
                      "Maximum number of iterations reached."
      ENDIF
    !
    ENDIF &
    INSTANT_CONVERGENCE 

    IF(ionode)THEN
    DO it = 0,nconf
      CLOSE(10000+it)
    ENDDO
    ENDIF
    !
#ifdef timer_CALL
      ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
      CALL print_timers_header()
      CALL t_tkprec%print()
      CALL t_tkaout%print()
      CALL t_tkain%print()
      CALL t_tkcg%print()
      CALL t_restart%print()
      ioWRITE(*,'(a)') "*** * High level components:"
      CALL t_tktld%print()
      CALL t_lwisot%print()
      CALL t_lwcasi%print()
      CALL t_lwphph%print()
      CALL t_lwchk%print()
      ioWRITE(*,'(a)') "*** * Low level subroutines: "
      CALL t_freq%print() 
      CALL t_bose%print() 
      CALL t_sum%print() 
      CALL t_fc3int%print() 
      CALL t_fc3dint%print() 
      CALL t_fc3m2%print() 
      CALL t_fc3rot%print() 
      CALL t_mpicom%print() 
      CALL t_merged%print()
#endif
    END SUBROUTINE TK_CG_prec
  END MODULE thermalk_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!             !               !               !               !               !               !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM thermalk

  USE kinds,            ONLY : DP
  USE thermalk_program
!   USE environment,      ONLY : environment_start, environment_end
!   USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid !, setup_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE nanoclock,        ONLY : init_nanoclock
  USE more_constants,   ONLY : print_citations_linewidth
  !
  USE posix_signal,       ONLY : set_TERMINATE_GRACEFULLY
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid)     :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)       :: S
  TYPE(code_input_type)      :: tkinput
  TYPE(q_grid)               :: out_grid

!   CALL mp_world_start(world_comm)
!   CALL environment_start('TK')
  CALL init_nanoclock()
  CALL start_mpi()
  CALL print_citations_linewidth()
  CALL set_TERMINATE_GRACEFULLY() !print_timers_and_die)

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("TK", tkinput, out_grid, S, fc2, fc3)
  !
  IF(TRIM(tkinput%calculation) == "sma") THEN
    !
    CALL TK_SMA(tkinput, out_grid, S, fc2, fc3)
    !
  ELSEIF(TRIM(tkinput%calculation) == "cgp" .or. TRIM(tkinput%calculation) == "exact") THEN
    !
    CALL TK_CG_prec(tkinput, out_grid, S, fc2, fc3)
    !
  ELSE
    CALL errore("tk", "Unknown calculation type: "//TRIM(tkinput%calculation), 1)
  ENDIF
  !
  IF(ionode) CALL print_citations_linewidth()
  CALL stop_mpi()
 
END PROGRAM thermalk
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

