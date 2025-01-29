!
! Written by Raja Sen (2024) IMPMC @ UPMC / CNRS UMR7590
!  This subroutine is taken from elphbolt code.
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE tetra_raja
  !
  USE kinds, ONLY : DP
  USE mpi_thermal, ONLY : ionode
  real(dp), parameter :: EPSILON_INEQUALITY = 1.0e-10_dp
  !
  !
CONTAINS
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  !
  SUBROUTINE form_tetrahedra_3d(nk, mesh, tetra, tetracount, tetramap)
    !! Form all the tetrahedra of a 3d FBZ mesh.
    !!
    !! nk Number of points in the list of FBZ wave vectors
    !! mesh Wave vector grid
    !! tetra List of the tetrahedra vertices
    !! tetracount Number of tetrahedra in which a wave vector belongs
    !! tetramap Wave vector to (tetrahedron, vertex) mapping
    !! blocks Is the FBZ wave vector list full or energy restricted?
    !! indexlist List of muxed indices of the FBZ wave vectors

    IMPLICIT NONE

    INTEGER, INTENT(in)    :: nk, mesh(3)
    INTEGER, INTENT(inout) :: tetra(6*nk, 4), tetracount(nk), tetramap(2, nk, 24)

    !Local variables
    INTEGER :: ik, i, j, k, ijk(3), ii, jj, kk, tk, tl, aux, count
    INTEGER :: ip1, jp1, kp1, n1, n2, n3, tmp
    INTEGER :: tetra_vertices_labels(6, 4)
    INTEGER :: scvol_vertices(8, 3) ! subcell volume vertices

    n1 = mesh(1)
    n2 = mesh(2)
    n3 = mesh(3)

    !Label of the vertices of the tetrahedra for a given subcell
    tetra_vertices_labels = reshape([ &
      1, 2, 3, 6, &
      1, 3, 5, 6, &
      3, 5, 6, 7, &
      3, 6, 7, 8, &
      3, 4, 6, 8, &
      2, 3, 4, 6 ], &
      shape(tetra_vertices_labels), order = [2, 1])

    tetra(:,:) = 0
    tetracount(:) = 0
    tetramap(:,:,:) = 0
    count = 1 !tetrahedron counter

    do ik = 1, nk !Run over all wave vectors in FBZ
      call demux_vector(ik, ijk, mesh, 1)
      !call demux_vector(ik, ijk, mesh)
      i = ijk(1)
      j = ijk(2)
      k = ijk(3)

      !Apply periodic boundary condition
      if (i == n1) then
        ip1 = 1
      else
        ip1 = i + 1
      end if

      if (j == n2) then
        jp1 = 1
      else
        jp1 = j + 1
      end if

      if (k == n3) then
        kp1 = 1
      else
        kp1 = k + 1
      end if

      !For each subcell save the vertices
      scvol_vertices = reshape([ &
        i,   j,   k,   &
        ip1, j,   k,   &
        i,   jp1, k,   &
        ip1, jp1, k,   &
        i,   j,   kp1, &
        ip1, j,   kp1, &
        i,   jp1, kp1, &
        ip1, jp1, kp1 ], &
        shape(scvol_vertices), order = [2, 1])
      !scvol_vertices = reshape([ &
      !      i,   j,   k,     &
      !      i,   j,   kp1,   &
      !      i,   jp1, k,     &
      !      i,   jp1, kp1,   &
      !      ip1, j,   k,     &
      !      ip1, j,   kp1,   &
      !      ip1, jp1, k,     &
      !      ip1, jp1, kp1 ], &
      !      shape(scvol_vertices), order = [2, 1])
      !
      do tk = 1, 6 !Run over 6 tetrahedra
        do tl = 1, 4 !Run over the labels of the vertices that
          !make up each tetrahedron
          aux = tetra_vertices_labels(tk, tl)
          ii = scvol_vertices(aux,1)
          jj = scvol_vertices(aux,2)
          kk = scvol_vertices(aux,3)
          aux = mux_vector([ii, jj, kk], mesh, 1)
          tmp = aux !Guaranteed to be > 0
          tetra(count, tl) = tmp ! tmp is wave vector index

          if(tmp > 0) then
            !Save the mapping of a wave vector index to a (tetrahedron, vertex)
            ! Each wave vector is associated with 24 different tetrahedra.
            ! So tetracount(tmp) <= 24
            tetracount(tmp) = tetracount(tmp) + 1
            tetramap(1, tmp, tetracount(tmp)) = count
            tetramap(2, tmp, tetracount(tmp)) = tl
          else
            call errore ('tmp is less than zero',1)
          end if
        end do
        count = count + 1
      end do
    end do
  END SUBROUTINE form_tetrahedra_3d
  !
  subroutine fill_tetrahedra_3d(nk, pbands, tetra, evals, tetra_evals)
    !! Populate the (sorted along the vertices) eigenvalues on all the vertices of the tetrahedra
    !!
    !! tetra List of the tetrahedra vertices
    !! evals List of eigenvalues
    !! tetra_evals Tetrahedra populated with the eigenvalues
    integer, intent(in) :: nk
    integer, intent(in) :: pbands
    integer, intent(in) :: tetra(6*nk, 4)
    real(DP), intent(in) :: evals(pbands, nk)
    real(DP), intent(out) :: tetra_evals(6*nk, pbands, 4)

    !Local variables
    integer :: iv, it, ib, numbands, aux, numtetra

    !numtetra = size(tetra(:, 1))
    !numbands = size(evals(:, 1))
    numtetra = 6*nk
    numbands = pbands

    !allocate(tetra_evals(numtetra, numbands, 4))

    !Note: Eigenvalues outside the transport active window is taken to be zero.
    !      As such, close to the transport window boundary, this method is
    !      inaccurate.
    !      A large enough transport window must be chosen to obtain accurate
    !      transport coefficients.
    tetra_evals(:,:,:) = 0._dp

    do it = 1, numtetra !Run over tetrahedra
      !do ib = 1, numbands !Run over bands
      do iv = 1, 4 !Run over vertices
        aux = tetra(it, iv)
        if(aux > 0) then !Only eigenvalues inside transport active region
          tetra_evals(it, :, iv) = evals(:, aux)
        else
          call errore ('aux is less than zero',1)
        end if
      end do
    end do

    do it = 1, numtetra
      do ib = 1, numbands
        call sort(tetra_evals(it, ib, :))
      end do
    end do
  end subroutine fill_tetrahedra_3d
  !
  function delta_fn_tetra(e, ik, ib, mesh, tetramap, tetracount, tetra_evals)
    !! Calculate delta function using the tetraheron method.
    !!
    !! e Sample energy
    !! ik Wave vector index
    !! ib Band index
    !! mesh Wave vector grid
    !! tetramap Wave vector to (tetrahedron, vertex) mapping
    !! tetracount Number of tetrahedra in which a wave vector belongs
    !! tetra_evals Tetrahedra populated with the eigenvalues

    !$acc routine seq

    real(DP), intent(in) :: e
    integer, intent(in) :: ik, ib
    integer, intent(in) :: mesh(3), tetramap(:,:,:), tetracount(:)
    real(DP), intent(in) :: tetra_evals(:,:,:)

    !Local variables
    integer :: iv, it, itk, num, numtetra
    logical :: c1, c2, c3
    real(DP) :: delta_fn_tetra
    real(DP) :: e1, e2, e3, e4, e1e, e2e, e3e, e4e, &
      e21, e31, e41, e32, e42, e43, tmp ! eji \equiv ej - ei

    tmp = 0._dp
    delta_fn_tetra = 0._dp

    !Total number of tetrahedra in the system
    numtetra = product(mesh)*6

    !Grab number of tetrahedra in which wave vector belongs
    ! Num is always 24.
    num = tetracount(ik)

    do itk = 1, num !Run over tetrahedra
      it = tetramap(1, ik, itk) !Grab tetrahedron
      iv = tetramap(2, ik, itk) !Grab vertex

      !Grab vertex energies
      e1 = tetra_evals(it, ib, 1)
      e2 = tetra_evals(it, ib, 2)
      e3 = tetra_evals(it, ib, 3)
      e4 = tetra_evals(it, ib, 4)

      !Define the energy differences
      e1e = e1 - e
      e2e = e2 - e
      e3e = e3 - e
      e4e = e4 - e
      e21 = e2 - e1
      e31 = e3 - e1
      e41 = e4 - e1
      e32 = e3 - e2
      e42 = e4 - e2
      e43 = e4 - e3

      !Evaluate the three cases
      c1 = e1 <= e .and. e <= e2
      c2 = e2 <= e .and. e <= e3
      c3 = e3 <= e .and. e <= e4

      if(.not. (e < e1 .or. e > e4)) then
        !Evaluate the expressions for the three cases
        select case(iv)
         case(1)
          if(c1) then
            tmp = (e2e/e21 + e3e/e31 + e4e/e41)*(e1e**2)/e41/e31/e21

            if(e1 == e2) then
              tmp = 0._dp
            end if
          else if(c2) then
            tmp = -0.5_dp*(e3e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41) &
              + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))

            if(e2 == e3) then
              tmp = -0.5_dp*(e4e*e1e/e41/e42 + e1e/e41 &
                + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
            end if
          else if(c3) then
            tmp = (e4e**3)/(e41**2)/e42/e43

            if(e3 == e4) then
              tmp = (e4e**2)/(e41**2)/e42
            end if
          end if
         case(2)
          if(c1) then
            tmp = -(e1e**3)/(e21**2)/e31/e41

            if(e1 == e2) then
              tmp = 0.0_dp
            end if
          else if(c2) then
            tmp = -0.5_dp*(e3e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
              + e4e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41))

            if(e2 == e3) then
              tmp = -0.5_dp*(0._dp + e4e/e42/e41 + 0._dp &
                + e4e/(e42**2)*(0._dp + e4e*e1e/e41/e31 + 1._dp))
            end if
          else if(c3) then
            tmp = (e4e**3)/e41/(e42**2)/e43

            if(e3 == e4) then
              tmp = 0.0_dp
            end if
          end if
         case(3)
          if(c1) then
            tmp = -(e1e**3)/e21/(e31**2)/e41

            if(e1 == e2) then
              tmp = 0.0_dp
            end if
          else if(c2) then
            tmp = 0.5_dp*(e2e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
              + e1e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41))

            if(e2 == e3) then
              tmp = 0.5_dp*(0._dp + e4e/e42/e41 + e1e/e31/e41 &
                + e1e/(e31**2)*(0._dp + e4e*e1e/e41/e42 + e1e/e41))
            end if
          else if(c3) then
            tmp = (e4e**3)/e41/e42/(e43**2)

            if(e3 == e4) then
              tmp = 0.0_dp
            end if
          end if
         case(4)
          if(c1) then
            tmp = -(e1e**3)/e21/e31/(e41**2)
            if(e1 == e2) then
              tmp = 0.0_dp
            end if
          else if(c2) then
            tmp = 0.5_dp*(e2e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41) &
              + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))

            if(e2 == e3) then
              tmp = 0.5_dp*(0._dp &
                + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
            end if
          else if(c3) then
            tmp = -(e3e/e43 + e2e/e42 + e1e/e41)*(e4e**2)/e41/e42/e43

            if(e3 == e4) then
              tmp = 0.0_dp
            end if
          end if
        end select

        if ((e1 == e2) .and. (e1 == e3) .and. (e1 == e4) .and. (e == e1)) then
          tmp = 0.25_dp
        end if

        delta_fn_tetra = delta_fn_tetra + tmp
      end if ! .not. (e <= e1 .or. e >= e4)
    end do !itk

    if(delta_fn_tetra < 1.0e-12) delta_fn_tetra = 0.0_dp

    !Normalize with the total number of tetrahedra
    delta_fn_tetra = delta_fn_tetra/numtetra
  end function delta_fn_tetra
  !
  function real_tetra(e, ik, ib, mesh, tetramap, tetracount, tetra_evals)
    !! Calculate the real part of the matrix elements of the resolvent operator
    !! using the analytic tetraheron method.
    !! Lambin and Vigneron Phys. Rev. B 29 6 1984 Eqs. A3-A6
    !! Note that typos in Eqs. A4 and A5 have been corrected.
    !! Here we use the expressions given in
    !! V. Eyert The Augmented Spherical Wave Method DOI
    !10.1007/978-3-642-25864-0.
    !!
    !! e Sample energy
    !! ik Wave vector index
    !! ib Band index
    !! mesh Wave vector grid
    !! tetramap Wave vector to (tetrahedron, vertex) mapping
    !! tetracount Number of tetrahedra in which a wave vector belongs
    !! tetra_evals Tetrahedra populated with the eigenvalues

    real(DP), intent(in) :: e
    integer, intent(in) :: ik, ib
    integer, intent(in) :: mesh(3), tetramap(:,:,:), tetracount(:)
    real(DP), intent(in) :: tetra_evals(:,:,:)

    !Local variables
    integer :: iv, it, itk, num, numtetra
    logical :: c1, c2, c3, c4, c5, c6, c7, l01, l12, l23
    real(DP) :: real_tetra
    real(DP) :: e0, e1, e2, e3, &
      ee0, ee1, ee2, ee3, &
      logabs_ee0, logabs_ee1, logabs_ee2, logabs_ee3, &
      e01, e02, e03, e12, e13, e23, tmp

    tmp = 0.0_dp
    real_tetra = 0.0_dp

    !Total number of tetrahedra in the system
    numtetra = product(mesh)*6

    !Grab number of tetrahedra in which wave vector belongs
    num = tetracount(ik)

    do itk = 1, num !Run over tetrahedra
      it = tetramap(1, ik, itk) !Grab tetrahedron
      iv = tetramap(2, ik, itk) !Grab vertex

      !Grab vertex energies
      e0 = tetra_evals(it, ib, 1)
      e1 = tetra_evals(it, ib, 2)
      e2 = tetra_evals(it, ib, 3)
      e3 = tetra_evals(it, ib, 4)

      !Define the energy differences
      ee0 = e - e0
      ee1 = e - e1
      ee2 = e - e2
      ee3 = e - e3
      e01 = e0 - e1
      e02 = e0 - e2
      e03 = e0 - e3
      e12 = e1 - e2
      e13 = e1 - e3
      e23 = e2 - e3

      !Precalculate all the log(abs(e - e_vertex))
      logabs_ee0 = 0.0_dp
      if(less_than(0.0_dp, abs(ee0/e03))) logabs_ee0 = log(abs(ee0))

      logabs_ee1 = 0.0_dp
      if(less_than(0.0_dp, abs(ee0/e03))) logabs_ee1 = log(abs(ee1))

      logabs_ee2 = 0.0_dp
      if(less_than(0.0_dp, abs(ee0/e03))) logabs_ee2 = log(abs(ee2))

      logabs_ee3 = 0.0_dp
      if(less_than(0.0_dp, abs(ee0/e03))) logabs_ee3 = log(abs(ee3))

      ! logabs_ee0 = 0.0_dp
      ! if(ee0 /= 0.0_dp) logabs_ee0 = log(abs(ee0))

      ! logabs_ee1 = 0.0_dp
      ! if(ee1 /= 0.0_dp) logabs_ee1 = log(abs(ee1))

      ! logabs_ee2 = 0.0_dp
      ! if(ee2 /= 0.0_dp) logabs_ee2 = log(abs(ee2))

      ! logabs_ee3 = 0.0_dp
      ! if(ee3 /= 0.0_dp) logabs_ee3 = log(abs(ee3))


      !Evaluate the seven cases
      l01 = less_than(0.0_dp, ABS(e01/e03))
      l12 = less_than(0.0_dp, ABS(e12/e03))
      l23 = less_than(0.0_dp, ABS(e23/e03))


      c1 = l01 .and. l12 .and. l23
      c2 = .not. l01 .and. l12 .and. l23
      c3 = l01 .and. .not. l12 .and. l23
      c4 = l01 .and. l12 .and. .not. l23
      c5 = .not. l01 .and. .not. l12 .and. l23
      c6 = .not. l01 .and. l12 .and. .not. l23
      c7 = l01 .and. .not. l12 .and. .not. l23

      if(.not. (e < e0 .or. e > e3)) then
        !Evaluate the expressions for the seven cases
        select case(iv) !tetrahedron vertex number
         case(1)
          if(c1) then !Eq. 9.5.124 [x][x]
            tmp = -ee0**2/(e01*e02*e03) &
              *( 1.0_dp + (ee1/e01 + ee2/e02 + ee3/e03)*logabs_ee0 ) &
              + ee1**3/(e01**2*e12*e13)*logabs_ee1 &
              - ee2**3/(e02**2*e12*e23)*logabs_ee2 &
              + ee3**3/(e03**2*e13*e23)*logabs_ee3
          else if(c2) then !Eq. 9.5.130 [x][x]
            tmp = eval_Eq9_5_130()
          else if(c3) then !Eq. 9.5.133 [x][x]
            tmp = -ee0**2/(e01**2*e03)*( 1.0_dp + (2.0_dp*ee1/e01 + ee3/e03)*logabs_ee0 ) &
              - ee1**2/(e01**2*e13)*( 1.0_dp + (-2.0_dp*ee0/e01 + ee3/e13)*logabs_ee1 ) &
              + ee3**3/(e03*e13)**2*logabs_ee3
          else if(c4) then !Eq. 9.5.136 [x][x]
            tmp = -ee0**2/(e02**2*e01)*( 1.0_dp + (2.0_dp*ee2/e02 + ee1/e01)*logabs_ee0 ) &
              + ee2**2/(e02**2*e12)*( 1.0_dp - (2.0_dp*ee0/e02 + ee1/e12)*logabs_ee2 ) &
              + ee1**3/(e01*e12)**2*logabs_ee1
          else if(c5) then !Eq. 9.5.139 [x][x]
            tmp = eval_Eq9_5_139()
          else if(c6) then  !Eq. 9.5.141 [x][x]
            tmp = eval_Eq9_5_141()
          else if(c7) then !Eq. 9.5.143 [x][x]
            tmp = 3.0_dp*ee0**2*ee1/e01**4*(logabs_ee1 - logabs_ee0) &
              - 1.5_dp*ee1*(2.0_dp*ee0 - e01)/e01**3 &
              - 1.0_dp/e01
          end if
         case(2)
          if(c1) then !Eq. 9.5.125 [x][x]
            tmp = ee1**2/(e01*e12*e13) &
              *( 1.0_dp + (-ee0/e01 + ee2/e12 + ee3/e13)*logabs_ee1 ) &
              + ee0**3/(e01**2*e02*e03)*logabs_ee0 &
              - ee2**3/(e02*e12**2*e23)*logabs_ee2 &
              + ee3**3/(e03*e13**2*e23)*logabs_ee3
          else if(c2) then !Eq. 9.5.130 [x][x]
            tmp = eval_Eq9_5_130()
          else if(c3) then !Eq. 9.5.134 [x][x]
            tmp = eval_Eq9_5_134()
          else if(c4) then !Eq. 9.5.137 [x][x]
            tmp = ee1**2/(e12**2*e01)*( 1.0_dp + (2.0_dp*ee2/e12 - ee0/e01)*logabs_ee1 ) &
              + ee2**2/(e12**2*e02)*( 1.0_dp - (2.0_dp*ee1/e12 + ee0/e02)*logabs_ee2 ) &
              + ee0**3/(e01*e02)**2*logabs_ee0
          else if(c5) then !Eq. 9.5.139 [x][x]
            tmp = eval_Eq9_5_139()
          else if(c6) then  !Eq. 9.5.141 [x][x]
            tmp = eval_Eq9_5_141()
          else if(c7) then !Eq. 9.5.144 [x][x]
            tmp = eval_Eq9_5_144()
          end if
         case(3)
          if(c1) then !Eq. 9.5.126 [x][x]
            tmp = -ee2**2/(e02*e12*e23) &
              *( 1.0_dp + (-ee0/e02 - ee1/e12 + ee3/e23)*logabs_ee2 ) &
              + ee0**3/(e01*e02**2*e03)*logabs_ee0 &
              - ee1**3/(e01*e12**2*e13)*logabs_ee1 &
              + ee3**3/(e03*e13*e23**2)*logabs_ee3
          else if(c2) then !Eq. 9.5.131 [x][x]
            tmp = -ee2**2/(e02**2*e23)*( 1.0_dp + (-2.0_dp*ee0/e02 + ee3/e23)*logabs_ee2 ) &
              - ee0**2/(e02**2*e03)*( 1.0_dp + (2.0_dp*ee2/e02 + ee3/e03)*logabs_ee0 ) &
              + ee3**3/(e23*e03)**2*logabs_ee3
          else if(c3) then !Eq. 9.5.134 [x][x]
            tmp = eval_Eq9_5_134()
          else if(c4) then !Eq. 9.5.138 [x][x]
            tmp = eval_Eq9_5_138()
          else if(c5) then !Eq. 9.5.139 [x][x]
            tmp = eval_Eq9_5_139()
          else if(c6) then !Eq. 9.5.142 [x][x]
            tmp = eval_Eq9_5_142()
          else if(c7) then !Eq. 9.5.144 [x][x]
            tmp = eval_Eq9_5_144()
          end if
         case(4)
          if(c1) then !Eq. 9.5.127 [x][x]
            tmp = ee3**2/(e03*e13*e23) &
              *( 1.0_dp + (-ee0/e03 - ee1/e13 - ee2/e23)*logabs_ee3 ) &
              + ee0**3/(e01*e02*e03**2)*logabs_ee0 &
              - ee1**3/(e01*e12*e13**2)*logabs_ee1 &
              + ee2**3/(e02*e12*e23**2)*logabs_ee2
          else if(c2) then !Eq. 9.5.132 [x][x]
            tmp = ee3**2/(e03**2*e23)*( 1.0_dp - (2.0_dp*ee0/e03 + ee2/e23)*logabs_ee3 ) &
              - ee0**2/(e03**2*e02)*( 1.0_dp + (2.0_dp*ee3/e03 + ee2/e02)*logabs_ee0 ) &
              + ee2**3/(e23*e02)**2*logabs_ee2
          else if(c3) then !Eq. 9.5.135 [x][x]
            tmp = ee3**2/(e13**2*e03)*( 1.0_dp - (2.0_dp*ee1/e13 + ee0/e03)*logabs_ee3 ) &
              + ee1**2/(e13**2*e01)*( 1.0_dp + (2.0_dp*ee3/e13 - ee0/e01)*logabs_ee1 ) &
              + ee0**3/(e03*e01)**2*logabs_ee0
          else if(c4) then !Eq. 9.5. 138 [x][x]
            tmp = eval_Eq9_5_138()
          else if(c5) then !Eq. 9.5. 140 [x][x]
            tmp =  3.0_dp*ee0*ee3**2/e03**4*(logabs_ee0 - logabs_ee3) &
              + 1.5_dp*ee0*(2.0_dp*ee3 + e03)/e03**3 &
              + 1.0_dp/e03
          else if(c6) then !Eq. 9.5.142 [x][x]
            tmp = eval_Eq9_5_142()
          else if(c7) then !Eq. 9.5.144 [x][x]
            tmp = eval_Eq9_5_144()
          end if
        end select

        if(e0 == e1 .and. e1 == e2 .and. e2 == e3) tmp = 0.25_dp/ee0

        real_tetra = real_tetra + tmp
      end if
    end do !itk

    !Normalize with the total number of tetrahedra
    real_tetra = real_tetra/numtetra

  contains

    ![x][x]
    pure real(dp) function eval_Eq9_5_130()
      !! Right hand side of Eq. 9.5.130 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_130 = -ee2**3/(e23*e02**3)*logabs_ee2 &
        + ee3**3/(e23*e03**3)*logabs_ee3 &
        + ee0/(e02*e03)*( 0.5_dp + ee2/e02 + ee3/e03 &
        + ((ee2/e02)**2 + (ee3/e03)**2 + ee2*ee3/(e02*e03))*logabs_ee0 )
    end function eval_Eq9_5_130

    ![x][x]
    pure real(dp) function eval_Eq9_5_134()
      !! Right hand side of Eq. 9.5.134 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_134 = ee0**3/(e03*e01**3)*logabs_ee0 &
        + ee3**3/(e03*e13**3)*logabs_ee3 &
        - ee1/(e01*e13)*( 0.5_dp - ee0/e01 + ee3/e13 + &
        ((ee0/e01)**2 + (ee3/e13)**2 - ee0*ee3/(e01*e13))*logabs_ee1 )
    end function eval_Eq9_5_134

    ![x][x]
    pure real(dp) function eval_Eq9_5_138()
      !! Right hand side of Eq. 9.5.138 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_138 = ee0**3/(e01*e02**3)*logabs_ee0 &
        - ee1**3/(e01*e12**3)*logabs_ee1 &
        + ee2/(e02*e12)*( 0.5_dp - ee0/e02 - ee1/e12 + &
        ((ee0/e02)**2 + (ee1/e12)**2 + ee0*ee1/(e02*e12))*logabs_ee2 )
    end function eval_Eq9_5_138

    ![x][x]
    pure real(dp) function eval_Eq9_5_139()
      !! Right hand side of Eq. 9.5.139 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_139 = ee3**3/e03**4*(logabs_ee3 - logabs_ee0) &
        - (ee3**2 + 0.5_dp*ee3*e03 + e03**2/3.0_dp)/e03**3
    end function eval_Eq9_5_139

    ![x][x]
    pure real(dp) function eval_Eq9_5_141()
      !! Right hand side of Eq. 9.5.141 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_141 = 3.0_dp*ee0*ee2**2/e02**4*(logabs_ee0 - logabs_ee2) &
        + 1.5_dp*ee0*(2.0_dp*ee2 + e02)/e02**3 &
        + 1.0_dp/e02
    end function eval_Eq9_5_141

    ![x][x]
    pure real(dp) function eval_Eq9_5_142()
      !! Right hand side of Eq. 9.5.142 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_142 = 3.0_dp*ee0**2*ee2/e02**4*(logabs_ee2 - logabs_ee0) &
        - 1.5_dp*ee2*(2.0_dp*ee0 - e02)/e02**3 &
        - 1.0_dp/e02
    end function eval_Eq9_5_142

    ![x][x]
    pure real(dp) function eval_Eq9_5_144()
      !! Right hand side of Eq. 9.5.144 of
      !! V. Eyert The Augmented Spherical Wave Method DOI
      !10.1007/978-3-642-25864-0.

      eval_Eq9_5_144 = ee0**3/e01**4*(logabs_ee0 - logabs_ee1) &
        + (ee0**2 - 0.5_dp*ee0*e01 + e01**2/3.0_dp)/e01**3
    end function eval_Eq9_5_144
  end function real_tetra
  !
  SUBROUTINE demux_vector(i, v, mesh, base)
    !! Demultiplex index of a single wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base chooses whether v has 0- or 1-based indexing.

    INTEGER, INTENT(in) :: i, mesh(3), base
    INTEGER, INTENT(out) :: v(3)
    INTEGER :: aux

    if(base < 0 .or. base > 1) &
      call errore ('Base has to be either 0 or 1 in tetrahedron_scheme.f90:demux_vector',1)

    call int_div_tetra(i - 1, mesh(1), aux, v(1))
    call int_div_tetra(aux, mesh(2), v(3), v(2))
    if(base == 1) v = v + 1
  END SUBROUTINE demux_vector
  !
  SUBROUTINE int_div_tetra(num, denom, q, r)
    !! Quotient(q) and remainder(r) of the integer division num/denom.

    integer, intent(in) :: num, denom
    integer, intent(out) :: q, r

    q = num/denom
    r = mod(num, denom)
  END SUBROUTINE int_div_tetra
  !
  !  SUBROUTINE demux_vector(i, v, mesh)
  !    !! Demultiplex index of a single wave vector.
  !    !! i is the multiplexed index of a wave vector (always 1-based).
  !    !! v is the demultiplexed triplet of a wave vector.
  !    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
  !    IMPLICIT NONE
  !
  !    INTEGER, INTENT(IN) :: i, mesh(3)
  !    INTEGER, INTENT(OUT) :: v(3)
  !    INTEGER :: total_combinations
  !    INTEGER :: index
  !
  !    ! Number of total combinations
  !    total_combinations = mesh(1) * mesh(2) * mesh(3)
  !
  !    IF (i < 1 .OR. i > total_combinations) THEN
  !      call errore ('Index out of range in tetrahedron_scheme.f90:demux_vector',1)
  !    END IF
  !
  !    ! Convert 1-based index to 0-based for easier calculations
  !    index = i - 1
  !
  !    ! Compute each component based on index
  !    v(1) = (index / (mesh(2) * mesh(3))) + 1
  !    v(2) = (MOD(index / mesh(3), mesh(2))) + 1
  !    v(3) = (MOD(index, mesh(3))) + 1
  !
  !  END SUBROUTINE demux_vector
  !
  FUNCTION mux_vector(v, mesh, base)
    !! Multiplex index of a single wave vector.
    !! Output is always 1-based.
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base states whether v has 0- or 1-based indexing.

    !$acc routine seq

    integer, intent(in) :: v(3), mesh(3), base
    integer :: mux_vector

    !if(base == 0) then
    !   mux_vector = (v(3)*mesh(2) + v(2))*mesh(1) + v(1) + 1
    !else
    !mux_vector = ((v(3) - 1)*mesh(2) + (v(2) - 1))*mesh(1) + v(1)
    !mux_vector = ((v(1) - 1)*mesh(2) + (v(2) - 1))*mesh(3) + v(3)
    mux_vector = ((v(1) - 1) * mesh(2) * mesh(3)) + ((v(2) - 1) * mesh(3)) + (v(3) - 1) + 1

    !end if
  END FUNCTION mux_vector
  !
  subroutine sort(list)
    !! Swap sort list of reals

    real(DP), intent(inout) :: list(:)
    real(DP) :: aux, tmp
    integer :: i, j, n

    n = size(list)

    do i = 1, n
      aux = list(i)
      do j = i + 1, n
        if (aux > list(j)) then
          tmp = list(j)
          list(j) = aux
          list(i) = tmp
          aux = tmp
        end if
      end do
    end do
  end subroutine sort
  !
  pure function less_than(a, b) result(res)
    !! Check if a < b
    !!
    !! a First number
    !! b Second number

    real(dp), intent(in) :: a, b
    logical :: res

    res = a + EPSILON_INEQUALITY < b
  end function less_than
END MODULE

