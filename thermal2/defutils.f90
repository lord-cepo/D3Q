module defutils
    use constants, only : dp
contains
  function flatten_RR_cmplx(mat) result(mat_flat)
    complex(dp), intent(in) :: mat(:,:,:,:)
    integer :: nat3i, nat3j, nqi, nqj
    !
    complex(dp), allocatable :: mat_flat(:,:)
    integer :: i, j, na1, na2, n1, n2
    !
    nat3i = size(mat,1)
    nat3j = size(mat,2)
    nqi = size(mat,3)
    nqj = size(mat,4)

    allocate(mat_flat(nat3i*nqi, nat3j*nqj))
    do n1 = 1, nqi
      do n2 = 1, nqj
        do na2 = 1, nat3j
          j = na2 + (n2-1)*nat3j
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat_flat(i,j) = mat(na1,na2,n1,n2)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function flatten_RR_real(mat) result(mat_flat)
    real(dp), intent(in) :: mat(:,:,:,:)
    integer :: nat3i, nat3j, nqi, nqj
    !
    real(dp), allocatable :: mat_flat(:,:)
    integer :: i, j, na1, na2, n1, n2
    !
    nat3i = size(mat,1)
    nat3j = size(mat,2)
    nqi = size(mat,3)
    nqj = size(mat,4)

    allocate(mat_flat(nat3i*nqi, nat3j*nqj))
    do n1 = 1, nqi
      do n2 = 1, nqj
        do na2 = 1, nat3j
          j = na2 + (n2-1)*nat3j
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat_flat(i,j) = mat(na1,na2,n1,n2)
          enddo
        enddo
      enddo
    enddo
  end function
  !
  function unflatten_RR_cmplx(mat_flat, nqi, nqj) result(mat)
    complex(dp), intent(in) :: mat_flat(:,:)
    integer,  intent(in) :: nqi, nqj
    integer :: nat3i, nat3j
    complex(dp), allocatable :: mat(:,:,:,:)
    integer :: n1, n2, na1, na2, i, j

    nat3i = size(mat_flat,1) / nqi
    nat3j = size(mat_flat,2) / nqj
    allocate(mat(nat3i, nat3j, nqi, nqj))

    do n2 = 1, nqj
      do na2 = 1, nat3j
        j = na2 + (n2-1)*nat3j
        do n1 = 1, nqi
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat(na1,na2,n1,n2) = mat_flat(i,j)
          end do
        end do
      end do
    end do
  end function

  function unflatten_RR_real(mat_flat, nqi, nqj) result(mat)
    real(dp), intent(in) :: mat_flat(:,:)
    integer,  intent(in) :: nqi, nqj
    integer :: nat3i, nat3j
    real(dp), allocatable :: mat(:,:,:,:)
    integer :: n1, n2, na1, na2, i, j

    nat3i = size(mat_flat,1) / nqi
    nat3j = size(mat_flat,2) / nqj
    allocate(mat(nat3i, nat3j, nqi, nqj))

    do n2 = 1, nqj
      do na2 = 1, nat3j
        j = na2 + (n2-1)*nat3j
        do n1 = 1, nqi
          do na1 = 1, nat3i
            i = na1 + (n1-1)*nat3i
            mat(na1,na2,n1,n2) = mat_flat(i,j)
          end do
        end do
      end do
    end do
  end function
  !
  !
  subroutine diag_degen_cmplx(nat3, dim_deg, V1, Udeg, e1)
    use fc2_interpolate, only : mat2_diag
    !
    integer, intent(in) :: nat3, dim_deg
    complex(dp), intent(in) :: V1(nat3,nat3)
    complex(dp), intent(inout) :: Udeg(nat3,dim_deg)
    complex(dp), intent(out), optional :: e1(dim_deg)
    !
    complex(dp) :: Vdeg(dim_deg,dim_deg)
    !
    Vdeg = matmul(TRANSPOSE(CONJG(Udeg)), matmul(V1, Udeg))
    call mat2_diag(dim_deg, Vdeg, e1) !freq(i:i+dim_deg-1))
    Udeg = matmul(Udeg, Vdeg)
  end subroutine
  !
  subroutine diag_degen_real(nat3, dim_deg, V1, Udeg, e1)
    use fc2_interpolate, only : mat2_diag
    !
    integer, intent(in) :: nat3, dim_deg
    complex(dp), intent(in) :: V1(nat3,nat3)
    complex(dp), intent(inout) :: Udeg(nat3,dim_deg)
    real(dp), intent(out), optional :: e1(dim_deg)
    !
    complex(dp) :: Vdeg(dim_deg,dim_deg)
    !
    Vdeg = matmul(TRANSPOSE(CONJG(Udeg)), matmul(V1, Udeg))
    call mat2_diag(dim_deg, Vdeg, e1) !freq(i:i+dim_deg-1))
    Udeg = matmul(Udeg, Vdeg)
  end subroutine
  !
  subroutine lift_degen(nat3, freq, V1, U, e1)
    use thutils, only: near, id_mat, outer_product, near, braket
    use functions, only: quicksort_idx
    use fc2_interpolate, only : mat2_diag
    !
    integer, intent(in) :: nat3
    real(dp), intent(in) :: freq(nat3)
    complex(dp), intent(in) :: V1(nat3,nat3)
    complex(dp), intent(inout) :: U(nat3,nat3)
    real(dp), intent(out), optional :: e1(nat3)
    !
    REAL(DP),PARAMETER :: epsq = 1.e-8_dp
    integer :: i,j, dim_deg
    real(dp) :: e1_(nat3)
    !
    i = 1
    do ! i cycle
      do j = i+1, nat3
        if (.not. near(freq(i), freq(j))) exit
      enddo
      dim_deg = j-i
      if(dim_deg > 1) then
        call diag_degen_real(nat3, dim_deg, V1, U(:,i:i+dim_deg-1), e1_(i:i+dim_deg-1))
      else
        e1_(i) = REAL(braket(U(:,i), V1), DP)
      endif
      i = i + dim_deg
      if(i > nat3) exit
    enddo
    if (present(e1)) e1 = e1_
  end subroutine
  !
  SUBROUTINE freq_phq_degen(xq, S, fc2, freq, V1, U, e1)
    use fc2_interpolate, only : mat2_diag, fftinterp_mat2
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants,          ONLY : RY_TO_CMM1
    use thutils, only: near
    IMPLICIT NONE
    REAL(DP),INTENT(in)               :: xq(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    complex(dp), intent(in)           :: V1(S%nat3,S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    real(dp), intent(out), optional   :: e1(S%nat3)
    REAL(DP),PARAMETER :: epsq = 1.e-8_dp
    REAL(DP) :: cq(3), chk(3)
    LOGICAL :: gamma
    !
    ! RAF
    !U = CONJG(U)
    CALL fftinterp_mat2(xq, S, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    cq = xq
    CALL cryst_to_cart(1,cq,S%at,-1)
    gamma = ALL( ABS(cq-NINT(cq))<epsq)
    IF( gamma )THEN
      freq(1:3) = 0._dp
      U(:,1:3) = (0._dp, 0._dp)
    ENDIF

    !> I lift the degeneracy by diagonalizing the perturbation
    !> we obtain new kets and frequencies, needed for 2nd order (FGR)
    if (present(e1)) then
      call lift_degen(S%nat3, freq, V1, U, e1)
    else
      call lift_degen(S%nat3, freq, V1, U)
    endif

    chk(:) = cq(:)*fc2%nq(:)
    chk=chk-NINT(chk(:))
    IF(fc2%periodic .AND. SUM(ABS(chk)) > 1.d-8) CALL errore("interp","cannot interpolate with periodic matrices",1 )

    IF(ANY(freq<0._dp)) THEN
      WRITE(*,*) gamma
      WRITE(*,"(e12.5,e12.5,e12.5)") cq
      WRITE(*,"(e12.5,e12.5,e12.5)") xq
      WHERE    (freq >  0.0)
        freq = DSQRT(freq)
      ELSEWHERE(freq < 0.0)
        freq = -DSQRT(-freq)
      ENDWHERE
      WRITE(*,"('negative freq = ',12e12.4)")freq*RY_TO_CMM1
      CALL errore("freq_phq_safe", "cannot continue with negative frequencies",1)
    ELSE
      !
      freq = DSQRT(freq)
    END IF
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE quter_cmplx(nq1, nq2, nq3, nat,tau,at,bg, matq, gridq, fc, xR, far)
    USE input_fc,  ONLY : forceconst2_grid, allocate_fc2_grid, write_fc2
    USE constants, ONLY : tpi
    USE functions, ONLY : default_if_not_present
    use quter_module, only : R_list_idx, expand_matR
    IMPLICIT NONE
    ! Dummy arguments
    INTEGER,INTENT(in)     :: nq1, nq2, nq3 ! dimensions of the q-points grid
    INTEGER,INTENT(in)     :: nat ! number of atoms
    REAL(DP),INTENT(in)    :: tau(3,nat) ! atom positions (alat units)
    REAL(DP),INTENT(in)    :: at(3,3), bg(3,3) ! real, reciprocal lattice
    COMPLEX(DP),INTENT(in) :: matq(3,3,nat,nat,nq1*nq2*nq3)
    REAL(DP),INTENT(in)    :: gridq(3,nq1*nq2*nq3)
    complex(dp), allocatable, INTENT(out) :: fc(:,:,:)
    INTEGER,OPTIONAL,INTENT(in) :: far
    REAL(DP), allocatable, INTENT(out)    :: xR(:,:)
    !
    REAL(DP),PARAMETER :: eps = 1.d-8
    INTEGER :: i, na1, na2, j1,j2, jn1,jn2,  iiq, iR, l1,l2,l3, far_
    INTEGER :: nRbig, nqt
    REAL(DP),ALLOCATABLE :: Rbig(:,:)
    REAL(DP) :: arg, dist(3), totalweight, Rx(3)
    complex(dp) :: aus
    !
    ! For the list of used R vectors:
    REAL(DP),POINTER :: Rout(:,:) => null()
    INTEGER :: iout, nRout
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=5000
    INTEGER :: nrws
    REAL(DP) :: atws(3,3) ! supercell of size nq1 x nq2 x nq3
    REAL(DP) :: wg, rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    !
    COMPLEX(DP),POINTER :: matR(:,:,:,:,:) => null()
    !
    far_ = default_if_not_present(2, far)

    nqt = nq1*nq2*nq3
    atws(:,1) = nq1*at(:,1)
    atws(:,2) = nq2*at(:,2)
    atws(:,3) = nq3*at(:,3)
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,atws)
    !
    ! Construct a big enough lattice of Supercell vectors
    IF(far_>0)THEN
      nRbig = (2*far_*nq1+1)*(2*far_*nq2+1)*(2*far_*nq3+1)
      ALLOCATE(Rbig(3,nRbig))
      nRbig=0
      DO l1=-far_*nq1,far_*nq1
        DO l2=-far_*nq2,far_*nq2
          DO l3=-far_*nq3,far_*nq3
            nRbig=nRbig+1
            Rbig(:, nRbig) = at(:,1)*l1 +at(:,2)*l2 +at(:,3)*l3
          END DO
        END DO
      END DO
      IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)
      !       WRITE(*,*) "seeking over ", nRbig," vectors"
    ELSE
      nRbig = nq1*nq2*nq3
      ALLOCATE(Rbig(3,nRbig))
      nRbig=0
      DO l1=0,nq1-1
        DO l2=0,nq2-1
          DO l3=0,nq3-1
            nRbig=nRbig+1
            Rbig(:, nRbig) = at(:,1)*l1 +at(:,2)*l2 +at(:,3)*l3
          END DO
        END DO
      END DO
      IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)
!        WRITE(*,*) "nfar==0 => standard FT over ", nRbig," vectors"
    ENDIF
    !
    ! dyn.mat. FFT
    !
    nRout=0
    DO na1=1,nat
      DO na2=1,nat
        DO j1=1,3
          DO j2=1,3
            totalweight = 0._dp
            DO iR= 1,nRbig
              !
              aus= 0._dp
              IF(far_>0)THEN
                dist(:) = Rbig(:,iR) + tau(:,na1) - tau(:,na2)
                wg = wsweight(dist,rws,nrws)
                wg = wg/DFLOAT(nqt)
              ELSE
                wg = 1._dp/DFLOAT(nqt)
              ENDIF
              IF(wg /= 0) THEN
                !
                DO iiq=1,nqt
                  !
                  arg=tpi*SUM(gridq(:,iiq)*Rbig(:,iR))
                  aus = aus + CMPLX(cos(arg),sin(arg),kind=DP)*matq(j1,j2,na1,na2,iiq)
                  !
                END DO
                !
                Rx = Rbig(:,iR)
                CALL cryst_to_cart(1,Rx,bg,-1)
                ! find if we have already used this R point ...
                iout = R_list_idx(nRout,Rout,Rx)
                ! ... if not, increase storage
                CALL expand_matR(iout,nRout,nat,matR)
                matR(j1,j2,na1,na2,iout) = aus * wg
                !
                totalweight=totalweight+wg
                !
              END IF
            END DO
            IF(ABS(totalweight-1._dp)>eps) THEN
              print*, totalweight, na1, na2
              CALL errore('main','wrong totalweight',1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    if(allocated(fc)) deallocate(fc)
    ALLOCATE( fc(3*nat,3*nat,nRout) )
    DO na1=1,nat
      DO na2=1,nat
        DO j1=1,3
          jn1 = j1 + (na1-1)*3
          DO j2=1,3
            jn2 = j2 + (na2-1)*3
            !
            DO i = 1, nRout
              fc(jn1,jn2,i) = matR(j1,j2,na1,na2,i)
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    if(allocated(xR)) deallocate(xR)
    allocate(xR(3,nRout))
    xR = Rout
    CALL cryst_to_cart(nRout, xR, at, 1)
    !
  END SUBROUTINE
  !
  SUBROUTINE fftinterp_mat2_cmplx(xq, S, fc, xR, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq(:)
    TYPE(ph_system_info),INTENT(in) :: S
    complex(dp), INTENT(in) :: fc(:,:,:)
    real(dp), INTENT(in) :: xR(:,:)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    !
    INTEGER :: i, mu,nu
    REAL(DP) :: varg(size(fc,3)), vcos(size(fc,3)), vsin(size(fc,3))
    COMPLEX(DP) :: vphase(size(fc,3))
    !
    D = (0._dp, 0._dp)
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
    FORALL(i=1:size(fc,3)) varg(i) =  tpi * SUM(xq(:)*xR(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
    CALL vdCos(size(fc,3), varg, vcos)
    CALL vdSin(size(fc,3), varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    DO i = 1, size(fc,3)
      DO mu= 1, S%nat3
        DO nu= 1, S%nat3
          D(nu,mu) = D(nu,mu) + vphase(i) * fc(nu,mu, i)
        ENDDO
      ENDDO
    END DO
    !
  END SUBROUTINE
end module
