 !
 ! Copyright (C) 2016-2023 EPW-Collaboration
 ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
 !
 ! This file is distributed under the terms of the GNU General Public
 ! License. See the file `LICENSE' in the root directory of the
 ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
 !
 !----------------------------------------------------------------------
MODULE EPW_utilities
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines associated with Broyden's method,
  !! Pade' approximants, DOS and Fermi level determination
  !!
  IMPLICIT NONE
  !
CONTAINS
  !
  subroutine mix_broyden_cmplx(ndim, deltaout, deltain, alphamix, iter, n_iter, df, dv)
    !
    ! Wrapper to call the real mix_broyden subroutine for complex arrays
    !
    USE kinds,         ONLY : DP
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ndim
    INTEGER, INTENT(in) :: iter
    INTEGER, INTENT(in) :: n_iter
    REAL(DP), INTENT(in) :: alphamix
    COMPLEX(DP), INTENT(inout) :: deltaout(ndim)
    COMPLEX(DP), INTENT(inout) :: deltain(ndim)
    COMPLEX(DP), INTENT(inout) :: df(ndim, n_iter)
    COMPLEX(DP), INTENT(inout) :: dv(ndim, n_iter)
    !
    REAL(DP), dimension(ndim) :: deltaout_r, deltain_r, deltaout_i, deltain_i
    REAL(DP), dimension(ndim, n_iter) :: df_r, dv_r, df_i, dv_i
    !
    df_r = real(df, dp)
    dv_r = real(dv, dp)
    deltaout_r = real(deltaout, dp)
    deltain_r = real(deltain, dp)
    df_i = aimag(df)
    dv_i = aimag(dv)
    deltaout_i = aimag(deltaout)
    deltain_i = aimag(deltain)
    !
    CALL mix_broyden(ndim, deltaout_r, deltain_r, alphamix, iter, n_iter, df_r, dv_r)
    CALL mix_broyden(ndim, deltaout_i, deltain_i, alphamix, iter, n_iter, df_i, dv_i)
    !
    deltain = cmplx(deltain_r, deltain_i, dp)
    df = cmplx(df_r, df_i, dp)
    dv = cmplx(dv_r, dv_i, dp)
    deltaout = cmplx(deltaout_r, deltaout_i, dp)
    !
  end subroutine
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mix_broyden(ndim, deltaout, deltain, alphamix, iter, n_iter, df, dv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    !
    USE kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(DP), INTENT(inout) :: deltaout(ndim)
    !! output delta at current iteration
    REAL(DP), INTENT(inout) :: deltain(ndim)
    !! delta at previous iteration
    REAL(DP), INTENT(inout) :: df(ndim, n_iter)
    !! arrays containing info from previous iterations
    REAL(DP), INTENT(inout) :: dv(ndim, n_iter)
    !! arrays containing info from previous iterations
    !
    ! Local variables
    INTEGER, PARAMETER :: maxter = 8
    !! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER :: n
    !! Counter on deltain/deltaout dimension (1 to nmin)
    INTEGER :: i, j
    !! Counter on iterations (1 to iter_used)
    INTEGER :: iter_used
    !! number of iterations used in mixing
    INTEGER :: ipos
    !! position at which results from the present iteraction are stored
    INTEGER :: inext
    !! position at which results for the next iteration are stored
    INTEGER :: ierr
    !! Error status
    INTEGER :: info
    !! Exit info DSYTRF and DSYTRI lapack subroutines
    INTEGER :: iwork(maxter)
    !! Workspace array in DSYTRF and DSYTRI lapack subroutines
    !
    REAL(DP) :: gammamix
    !! Mixing parameter
    REAL(DP) :: norm
    !! norm of df
    REAL(DP) :: inv_norm
    !! 1.0/norm. Defined for efficiency reasons
    REAL(DP) :: wg0
    !!
    REAL(DP) :: work(maxter)
    !!
    REAL(DP) :: wg(maxter)
    !!
    REAL(DP), ALLOCATABLE :: deltainsave(:)
    !! Array to store deltain from previous iteration
    REAL(DP) :: beta(maxter, maxter)
    !!
    REAL(DP), EXTERNAL :: DDOT
    !! Inner product of two vectors
    REAL(DP), EXTERNAL :: DNRM2
    !! Norm of a vector
    !!
    !
    ! adjustable parameters as suggested in the original paper
    wg0 = 1e-2_dp
    wg  = maxter
    !
    IF (iter < 1) CALL errore('mix_broyden', 'iter is smaller than 1', 1)
    IF (n_iter > maxter) CALL errore('mix_broyden', 'n_iter is too big', 1)
    IF (ndim <= 0) CALL errore('mix_broyden', 'ndim <= 0', 1)
    !
    ALLOCATE(deltainsave(ndim), STAT = ierr)
    IF (ierr /= 0) CALL errore('mix_broyden', 'Error allocating deltainsave', 1)
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = MIN(iter - 1, n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos = iter - 1 until ipos = n_iter, then back to 1, 2,...
    !
    ipos = iter - 1 - ((iter - 2) / n_iter) * n_iter
    !
    DO n = 1, ndim
      deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
      DO n = 1, ndim
        df(n, ipos) = deltaout(n) - df(n, ipos)
        dv(n, ipos) = deltain(n)  - dv(n, ipos)
      ENDDO
      norm = (DNRM2(ndim, df(1, ipos), 1))**2
      norm = DSQRT(norm)
      inv_norm = 1._dp / norm
      ! DSCAL scales df and dv by inv_norm
      CALL DSCAL(ndim, inv_norm, df(1, ipos), 1)
      CALL DSCAL(ndim, inv_norm, dv(1, ipos), 1)
    ENDIF
    !
    DO i = 1, iter_used
      DO j = i + 1, iter_used
        beta(i, j) = wg(i) * wg(j) * DDOT(ndim, df(1, j), 1, df(1, i), 1)
      ENDDO
      beta(i, i) = wg0**2 + wg(i)**2
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix
    !
    CALL DSYTRF('u', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('mix_broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('u', iter_used, beta, maxter, iwork, work, info)
    CALL errore('mix_broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
      DO j = i + 1, iter_used
        beta(j, i) = beta(i, j)
      ENDDO
    ENDDO
    !
    DO i = 1, iter_used
      work(i) = DDOT(ndim, df(1, i), 1, deltaout, 1)
    ENDDO
    !
    DO n = 1, ndim
      deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
      gammamix = 0._dp
      DO j = 1, iter_used
        gammamix = gammamix + beta(j, i) * wg(j) * work(j)
      ENDDO
      !
      DO n = 1, ndim
        deltain(n) = deltain(n) - wg(i) * gammamix * (alphamix * df(n, i) + dv(n, i))
      ENDDO
    ENDDO
    !
    inext = iter - ((iter - 1) / n_iter) * n_iter
    df(:, inext) = deltaout(:)
    dv(:, inext) = deltainsave(:)
    !
    DEALLOCATE(deltainsave, STAT = ierr)
    IF (ierr /= 0) CALL errore('mix_broyden', 'Error deallocating deltainsave', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE
  !-----------------------------------------------------------------------------

SUBROUTINE mix_broyden_full(ndim, deltaout, deltain, alphamix, iter, n_iter, df, dv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    !
    USE kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    complex(DP), INTENT(inout) :: deltaout(ndim)
    !! output delta at current iteration
    complex(DP), INTENT(inout) :: deltain(ndim)
    !! delta at previous iteration
    complex(DP), INTENT(inout) :: df(ndim, n_iter)
    !! arrays containing info from previous iterations
    complex(DP), INTENT(inout) :: dv(ndim, n_iter)
    !! arrays containing info from previous iterations
    !
    ! Local variables
    INTEGER, PARAMETER :: maxter = 8
    !! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER :: n
    !! Counter on deltain/deltaout dimension (1 to nmin)
    INTEGER :: i, j
    !! Counter on iterations (1 to iter_used)
    INTEGER :: iter_used
    !! number of iterations used in mixing
    INTEGER :: ipos
    !! position at which results from the present iteraction are stored
    INTEGER :: inext
    !! position at which results for the next iteration are stored
    INTEGER :: ierr
    !! Error status
    INTEGER :: info
    !! Exit info DSYTRF and DSYTRI lapack subroutines
    INTEGER :: iwork(maxter)
    !! Workspace array in DSYTRF and DSYTRI lapack subroutines
    !
    complex(DP) :: gammamix
    !! Mixing parameter
    REAL(DP) :: norm
    !! norm of df
    REAL(DP) :: inv_norm
    !! 1.0/norm. Defined for efficiency reasons
    REAL(DP) :: wg0
    !!
    complex(DP) :: work(maxter)
    !!
    REAL(DP) :: wg(maxter)
    !!
    complex(DP), ALLOCATABLE :: deltainsave(:)
    !! Array to store deltain from previous iteration
    complex(DP) :: beta(maxter, maxter)
    real(dp), external :: dznrm2
    complex(dp), external :: ZDOTC
    !!
    !
    ! adjustable parameters as suggested in the original paper
    wg0 = 1e-2_dp
    wg  = maxter
    !
    IF (iter < 1) CALL errore('mix_broyden', 'iter is smaller than 1', 1)
    IF (n_iter > maxter) CALL errore('mix_broyden', 'n_iter is too big', 1)
    IF (ndim <= 0) CALL errore('mix_broyden', 'ndim <= 0', 1)
    !
    ALLOCATE(deltainsave(ndim), STAT = ierr)
    IF (ierr /= 0) CALL errore('mix_broyden', 'Error allocating deltainsave', 1)
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = MIN(iter - 1, n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos = iter - 1 until ipos = n_iter, then back to 1, 2,...
    !
    ipos = iter - 1 - ((iter - 2) / n_iter) * n_iter
    !
    DO n = 1, ndim
      deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
      DO n = 1, ndim
        df(n, ipos) = deltaout(n) - df(n, ipos)
        dv(n, ipos) = deltain(n)  - dv(n, ipos)
      ENDDO
      norm = dznrm2(ndim, df(1, ipos), 1)
      inv_norm = 1._dp / norm
      ! DSCAL scales df and dv by inv_norm
      CALL ZDSCAL(ndim, inv_norm, df(1, ipos), 1)
      CALL ZDSCAL(ndim, inv_norm, dv(1, ipos), 1)
    ENDIF
    !
    DO i = 1, iter_used
      DO j = i + 1, iter_used
        beta(i, j) = wg(i) * wg(j) * ZDOTC(ndim, df(1, j), 1, df(1, i), 1)
      ENDDO
      beta(i, i) = wg0**2 + wg(i)**2
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix
    !
    CALL ZHETRF('u', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('mix_broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL ZHETRI('u', iter_used, beta, maxter, iwork, work, info)
    CALL errore('mix_broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
      DO j = i + 1, iter_used
        beta(j, i) = conjg(beta(i, j))
      ENDDO
    ENDDO
    !
    DO i = 1, iter_used
      work(i) = ZDOTC(ndim, df(1, i), 1, deltaout, 1)
    ENDDO
    !
    DO n = 1, ndim
      deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
      gammamix = 0._dp
      DO j = 1, iter_used
        gammamix = gammamix + beta(j, i) * wg(j) * work(j)
      ENDDO
      !
      DO n = 1, ndim
        deltain(n) = deltain(n) - wg(i) * gammamix * (alphamix * df(n, i) + dv(n, i))
      ENDDO
    ENDDO
    !
    inext = iter - ((iter - 1) / n_iter) * n_iter
    df(:, inext) = deltaout(:)
    dv(:, inext) = deltainsave(:)
    !
    DEALLOCATE(deltainsave, STAT = ierr)
    IF (ierr /= 0) CALL errore('mix_broyden', 'Error deallocating deltainsave', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE
  !-----------------------------------------------------------------------------
END MODULE