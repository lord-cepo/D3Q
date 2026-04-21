module symm_q_mat
  use kinds, only : DP
  USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, irt, copy_sym, nrot, inverse_s
  USE cell_base,          ONLY :  at, bg, celldm, ibrav, omega
  USE ions_base,          ONLY : ityp, ntyp => nsp, atm, tau, amass
  use thutils,            only : cryst2cart
  USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq

  ! ====================================================================
  ! START OF STAR AVERAGING BLOCK
  ! ====================================================================
  ! At this point, we have:
  ! - xq: the irreducible q-point
  ! - s(:,:,1:nsym): all crystal symmetries
  ! - invs(:): inverses of the symmetry operations
  ! - irt, rtau: atom rotations and fractional translations
  ! - isq(:): maps each symmetry operation to a point in the star
  ! - sxq(:,iq): the actual q-vectors of the star
  ! - phi_star(:,:,:,:,iq): the un-averaged matrices at each q-point in the star
  !
  ! We want to compute: phi (the perfectly averaged matrix at xq)
  ! ====================================================================
  public :: apply_sym
contains
  subroutine apply_sym(at_, bg_, nat, ityp_, tau_, dyn_in, equiv, grid, average)
    real(dp), intent(in) :: at_(3,3), bg_(3,3)
    integer, intent(in) :: nat
    real(dp), intent(in) :: tau_(3,nat)
    integer, intent(in) :: ityp_(nat)
    integer, intent(in) :: equiv(:)
    real(dp), intent(in) :: grid(:,:)
    complex(dp), allocatable, intent(inout) :: dyn_in(:,:,:)
    logical, intent(in) :: average
    !
    real(dp) :: sxq(3,48), diff(3), xq(3)
    logical :: sym(48)
    logical :: found
    integer :: isym, iq, i, j, nqs, iiq, iq_sym
    integer :: nq, iq_irr, isq(48), iq_isym, nq_irr, invsm, imq
    complex(dp) :: phi_in(3,3,nat,nat, size(dyn_in,3))
    complex(DP) :: phi_tmp(3,3,nat,nat), d2(3*nat, 3*nat)
    complex(DP), allocatable :: phi_avg(:,:,:,:,:), dyn_star(:,:,:)
    real(dp), allocatable :: m_loc(:,:)
    integer, allocatable :: irr_map(:)
    !
    at = at_
    bg = bg_
    ityp = ityp_
    tau = tau_
    !
    if(.not. allocated(rtau)) allocate(rtau(3, 48, nat))
    ALLOCATE(m_loc(3,nat))
    m_loc = 0._dp
    do iq = 1, size(dyn_in, 3)
      call scompact_dyn(nat, dyn_in(:,:,iq), phi_in(:,:,:,:,iq))
      call trntnsc_ats(phi_in(:,:,:,:,iq), at, bg, -1)
    enddo
    nq = size(grid, 2)
    !
    nq_irr = maxval(equiv)
    allocate(irr_map(nq_irr))
    allocate(phi_avg(3,3,nat,nat,nq_irr))
    !
    ! ! ######################### symmetry setup #########################
    ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~
    !
    ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~
    CALL set_sym_bl ( )
    CALL find_sym ( nat, tau, ityp, .false., m_loc )
    !
    CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
    !
    phi_avg = 0._dp
    irr_map = 0
    do iq_irr = 1, nq
      if(irr_map(equiv(iq_irr)) /= 0) then
        cycle
      else
        irr_map(equiv(iq_irr)) = iq_irr
      end if
      !
      if(average) then
        CALL star_q(grid(:,iq_irr), at, bg, nsym, s, invs, nqs, sxq, isq, imq, .false. )
        ! Loop over all symmetry operations of the crystal
        do isym = 1, nsym
          ! Find which q-point in the star this symmetry operation generates
          found = .false.
          do iq = 1, nq
            diff = cryst2cart(sxq(:,isq(isym)) - grid(:,iq), at, -1)
            if (norm2(diff-NINT(diff)) < 1e-5_dp) then
              iq_isym = iq
              found = .true.
              exit
            end if
          enddo
          if(.not. found) then
            print*, cryst2cart(sxq(:,isq(isym)), at, -1)
            print*, "-------------------------------"
            print*, cryst2cart(grid(:,iq_irr), at, -1)
            call errore("apply_sym", "not all the star is in the grid", 1)
          endif
          ! 1. Take the matrix at that q-point
          ! 2. Transform it to crystal coordinates

          ! 3. Rotate it BACK to the irreducible point xq
          ! We use the INVERSE operation (invs(isym)) to map q_eq -> xq.
          ! Note: rotate_and_add_dyn handles the rtau phases internally.
          ! The last argument is the wavevector we are coming FROM (sxq(:, iq)).
          call rotate_and_add_dyn (phi_in(:,:,:,:,iq_isym), phi_avg(:,:,:,:,equiv(iq_irr)), &
            nat, invs(isym), s, invs, irt, rtau, sxq(:,isq(isym)) )
        enddo
      endif
      ! ====================================================================
      ! END OF STAR AVERAGING BLOCK
      ! ====================================================================
    enddo
    phi_avg = phi_avg / real(nsym, dp)
    !
    if (.not. average) then
      do iq = 1, nq_irr
        phi_avg(:,:,:,:,iq) = phi_in(:,:,:,:,iq)
      enddo
      deallocate(dyn_in)
      allocate(dyn_in(3*nat, 3*nat, nq))
    endif
    !
    dyn_in = (1._dp, 0._dp)
    minus_q = .true.
    do iq_irr = 1, nq_irr
      CALL set_sym_bl ( )
      CALL find_sym ( nat, tau, ityp, .false., m_loc )
      !
      xq = grid(:,irr_map(iq_irr))
      sym = .false.
      sym(1:nsym) = .true.
      CALL smallg_q(xq, 0, at, bg, nsym, s, sym, minus_q)
      nsymq = copy_sym(nsym, sym)
      ! recompute the inverses as the order of sym.ops. has changed
      CALL inverse_s ( )
      ! part 2: this computes gi, gimq
      call set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
      !
      ! finally this does some of the above again and also computes rtau...
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      ! ######################### star of q #########################
      CALL symdynph_gq_no_herm(xq, phi_avg(:,:,:,:,iq_irr), s, invs, rtau, irt, nsymq, nat, &
        irotmq, minus_q)
      !
      CALL star_q(xq, at, bg, nsym, s, invs, nqs, sxq, isq, imq, .false. )
      !
      call compact_dyn(nat, d2, phi_avg(:,:,:,:,iq_irr))
      allocate(dyn_star(3*nat, 3*nat, nqs))
      CALL q2qstar_ph_nowrite(d2, at, bg, nat, nsym, s, invs, irt, rtau, &
        nqs, sxq, isq, imq, 1, dyn_star)
      ! dyn_star(:,:,1) = d2

      do iq = 1, nqs
        found = .false.
        do iiq = 1, nq
          diff = cryst2cart(sxq(:,iq) - grid(:,iiq), at, -1)
          if (norm2(diff-NINT(diff)) < 1e-5_dp) then
            iq_sym = iiq
            found = .true.
            exit
          end if
        enddo
        if(.not. found) call errore("apply_sym", "(2) not all the star is in the grid", 1)
        if(any(abs(dyn_in(:,:,iq_sym))-1._dp > 1e-10_dp)) call errore("apply_sym", "dyn_in is not empty", 1)
        dyn_in(:,:,iq_sym) = dyn_star(:,:,iq)
      enddo
      deallocate(dyn_star)
    enddo
  end subroutine
  !
  subroutine trntnsc_ats(phi_, at_, bg_, sign)
    complex(dp), intent(inout) :: phi_(:,:,:,:)
    real(dp), intent(in) :: at_(3,3), bg_(3,3)
    integer, intent(in) :: sign
    !
    integer :: na, nb, nat_
    nat_ = size(phi_, 3)
    !
    do na = 1, nat_
      do nb = 1, nat_
        call trntnsc (phi_(1,1,na,nb), at_, bg_, sign)
      enddo
    enddo
  end subroutine
!
  subroutine q2qstar_ph_nowrite(dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
    nq, sxq, isq, imq, iudyn, dyn_grid)
    !-----------------------------------------------------------------------
    !! Generates the dynamical matrices for the star of q and writes them on
    !! disk for later use.
    !! If there is a symmetry operation such that \(q \rightarrow -q+G \) then
    !! imposes on dynamical matrix those conditions related to time reversal
    !! symmetry.
    !
    USE kinds, only : DP
    USE io_dyn_mat, only : write_dyn_mat
    USE control_ph, only : xmldyn
    implicit none
    !
    integer :: nat
    !! number of atoms in the unit cell
    integer :: nsym
    !! number of symmetry operations
    integer :: s(3,3,48)
    !! the symmetry operations
    integer :: invs(48)
    !! index of the inverse operations
    integer :: irt(48,nat)
    !! index of the rotated atom
    integer :: nq
    !! degeneracy of the star of q
    complex(dp), intent(out) :: dyn_grid(3*nat, 3*nat, nq)
    integer :: isq(48)
    !! symmetry op. giving the rotated q
    integer :: imq
    !! index of -q in the star (0 if non present)
    integer :: iudyn
    !! unit number
    complex(DP), intent(inout) :: dyn(3*nat,3*nat)
    !! the input matrix used to generate the star. If \(\text{imq}\) is
    !! different from 0, the \(-q\) partner is generated by conjugation
    real(DP) :: at (3,3)
    !! direct lattice vectors
    real(DP) :: bg (3,3)
    !! reciprocal lattice vectors
    real(DP) :: rtau (3,48,nat)
    !! for each atom and rotation gives the R vector involved
    real(DP) :: sxq (3,48)
    !! list of q in the star
    !
    ! ... local variables
    !
    integer :: na, nb, iq, nsq, isym, icar, jcar, i, j
    ! counters
    ! nsq: number of sym.op. giving each q in the list

    complex(DP) :: phi (3, 3, nat, nat), phi2 (3, 3, nat, nat)
    ! work space
    complex(dp) :: d2(3*nat, 3*nat)
    !
    ! Sets number of symmetry operations giving each q in the list
    !
    nsq = nsym / nq
    if (nsq * nq /= nsym) call errore ('q2star_ph', 'wrong degeneracy', 1)
    !
    ! Writes dyn.mat. dyn(3*nat,3*nat) on the 4-index array phi(3,3,nat,nat)
    !
    CALL scompact_dyn(nat, dyn, phi)
    !
    ! Go to crystal coordinates
    !
    ! If -q is in the list impose first of all the conditions coming from
    ! time reversal symmetry
    !
    ! if (imq /= 0) then
    !   phi2 (:,:,:,:) = (0.d0, 0.d0)
    !   isym = 1
    !   do while (isq (isym) /= imq)
    !     isym = isym + 1
    !   enddo
    !   call rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
    !     rtau, sxq (1, imq) )
    !   do na = 1, nat
    !     do nb = 1, nat
    !       do i = 1, 3
    !         do j = 1, 3
    !           phi (i, j, na, nb) = 0.5d0 * (phi (i, j, na, nb) + &
    !             CONJG(phi2(i, j, na, nb) ) )
    !         enddo
    !       enddo
    !     enddo
    !   enddo
    !   phi2 (:,:,:,:) = phi (:,:,:,:)
    !   ! !
    !   ! Back to cartesian coordinates
    !   !
    !   do na = 1, nat
    !     do nb = 1, nat
    !       call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
    !     enddo
    !   enddo
    !   !
    !   ! Saves 4-index array phi2(3,3,nat,nat) on the dyn.mat. dyn(3*nat,3*nat)
    !   !
    !   CALL compact_dyn(nat, dyn, phi2)
    ! endif
    !
    ! For each q of the star rotates phi with the appropriate sym.op. -> phi
    !
    do iq = 1, nq
      phi2 (:,:,:,:) = (0.d0, 0.d0)
      do isym = 1, nsym
        if (isq (isym) == iq) then
          call rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
            rtau, sxq (1, iq) )
        endif
      enddo
      phi2 (:,:,:,:) = phi2 (:,:,:,:) / DBLE (nsq)
      !
      ! Back to cartesian coordinates
      !
      do na = 1, nat
        do nb = 1, nat
          call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
        enddo
      enddo
      call compact_dyn(nat, d2, phi2)
      dyn_grid(:,:,iq) = d2
      !
      ! Writes the dynamical matrix in cartesian coordinates on file
      !
      ! IF (xmldyn) THEN
      !   call write_dyn_mat(nat, counter, sxq(1,iq), phi2)
      ! ELSE
      !   call write_dyn_on_file (sxq (1, iq), phi2, nat, iudyn)
      ! ENDIF
      if (imq == 0) then
        !
        print*, "WARNING: no inversion symmetry"
        ! if -q is not in the star recovers its matrix by time reversal
        !
        ! do na = 1, nat
        !   do nb = 1, nat
        !     do i = 1, 3
        !       do j = 1, 3
        !         phi2 (i, j, na, nb) = CONJG(phi2 (i, j, na, nb) )
        !       enddo
        !     enddo
        !   enddo
        ! enddo
        !
        ! and writes it (changing temporarily sign to q)
        !
        ! sxq (:, iq) = - sxq (:, iq)
        ! IF (xmldyn) THEN
        !   call write_dyn_mat(nat, counter, sxq(1,iq), phi2)
        ! ELSE
        !   call write_dyn_on_file (sxq (1, iq), phi2, nat, iudyn)
        ! ENDIF
        ! sxq (:, iq) = - sxq (:, iq)
      endif
    enddo
    !
    return
  end subroutine
!
  subroutine symdynph_gq_no_herm( xq, phi, s, invs, rtau, irt, nsymq, &
    nat, irotmq, minus_q )
    !-----------------------------------------------------------------------
    !! This routine receives as input an unsymmetrized dynamical
    !! matrix expressed on the crystal axes and imposes the symmetry
    !! of the small group of q. Furthermore it imposes also the symmetry
    !! q -> -q+G if present.
    !! February 2020: Update (A. Urru) to include the symmetry operations
    !! that require the time reversal operator (meaning that TS is a
    !! symmetry of the crystal). For more information please see:
    !! Phys. Rev. B 100, 045115 (2019).
    !
    USE kinds, only : DP
    USE constants, ONLY: tpi
    USE symm_base, ONLY : t_rev
    !
    implicit none
    !
    integer :: nat
    !! input: the number of atoms
    integer :: s(3,3,48)
    !! input: the symmetry matrices
    integer :: irt(48,nat)
    !! input: the rotated of each vector
    integer :: invs(48)
    !! input: the inverse of each matrix
    integer :: nsymq
    !! input: the order of the small group
    integer :: irotmq
    !! input: the rotation sending q ->-q+G
    real(DP) :: xq(3)
    !! input: the q point
    real(DP) :: rtau(3,48,nat)
    !! input: the R associated at each t
    logical :: minus_q
    !! input: true if a symmetry q->-q+G
    complex(DP) :: phi(3,3,nat,nat)
    !! inp/out: the matrix to symmetrize
    !
    ! ... local variables
    !
    integer :: isymq, sna, snb, irot, na, nb, ipol, jpol, lpol, kpol, &
      iflb (nat, nat)
    ! counters, indices, work space

    real(DP) :: arg
    ! the argument of the phase

    complex(DP) :: work (3, 3), faseq (48)
    ! work space, phase factors
    !
    !
    if ( (nsymq == 1) .and. (.not.minus_q) ) return
    !
    !    Here we symmetrize with respect to the small group of q
    !
    if (nsymq == 1) return

    iflb (:, :) = 0
    do na = 1, nat
      do nb = 1, nat
        if (iflb (na, nb) == 0) then
          work(:,:) = (0.d0, 0.d0)
          do isymq = 1, nsymq
            irot = isymq
            sna = irt (irot, na)
            snb = irt (irot, nb)
            arg = 0.d0
            do ipol = 1, 3
              arg = arg + (xq (ipol) * (rtau (ipol, irot, na) - &
                rtau (ipol, irot, nb) ) )
            enddo
            arg = arg * tpi
            faseq (isymq) = CMPLX(cos (arg), sin (arg) ,kind=DP)
            do ipol = 1, 3
              do jpol = 1, 3
                do kpol = 1, 3
                  do lpol = 1, 3
                    IF (t_rev(isymq)==1) THEN
                      work (ipol, jpol) = work (ipol, jpol) + &
                        s (ipol, kpol, irot) * s (jpol, lpol, irot) &
                        * CONJG(phi (kpol, lpol, sna, snb) * faseq (isymq))
                    ELSE
                      work (ipol, jpol) = work (ipol, jpol) + &
                        s (ipol, kpol, irot) * s (jpol, lpol, irot) &
                        * phi (kpol, lpol, sna, snb) * faseq (isymq)
                    ENDIF
                  enddo
                enddo
              enddo
            enddo
          enddo
          do isymq = 1, nsymq
            irot = isymq
            sna = irt (irot, na)
            snb = irt (irot, nb)
            do ipol = 1, 3
              do jpol = 1, 3
                phi (ipol, jpol, sna, snb) = (0.d0, 0.d0)
                do kpol = 1, 3
                  do lpol = 1, 3
                    IF (t_rev(isymq)==1) THEN
                      phi(ipol,jpol,sna,snb)=phi(ipol,jpol,sna,snb) &
                        + s(ipol,kpol,invs(irot))*s(jpol,lpol,invs(irot))&
                        * CONJG(work (kpol, lpol)*faseq (isymq))
                    ELSE
                      phi(ipol,jpol,sna,snb)=phi(ipol,jpol,sna,snb) &
                        + s(ipol,kpol,invs(irot))*s(jpol,lpol,invs(irot))&
                        * work (kpol, lpol) * CONJG(faseq (isymq) )
                    ENDIF
                  enddo
                enddo
              enddo
            enddo
            iflb (sna, snb) = 1
          enddo
        endif
      enddo
    enddo
    phi (:, :, :, :) = phi (:, :, :, :) / DBLE(nsymq)
    return
  end subroutine
!
end module
