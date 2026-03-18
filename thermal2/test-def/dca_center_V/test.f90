program test
  use constants, only: dp
  use dca, only: center_V
  use input_fc, only: read_fc2, aux_system, ph_system_info, &
    forceconst2_grid, div_mass_fc2
  use quter_defect, only: forceconst2_sc
  use quter_defect, only: asr3, fc_RR2sc, fc_sc2RR, fc_uc2RR
  use thutils, only: grid_vec_cart
  use test_print, only: allclose
  use q_grids, only: symmetrize_system
  use asr2_module, only: impose_asr2
  use test_print, only: test_log
  !
  type(ph_system_info) :: S, Sd(2)
  type(forceconst2_grid) :: fc2, fc2d(2)
  type(forceconst2_sc) :: fc2_sc(2)
  integer, parameter :: Q_LENGTH = 2
  real(dp) :: DRR(6,6,27,27), summ
  integer :: i, sc_grid(3), iq, jq, jn1, jn2, j, k, jmax, kmax, iR1, iR2
  character(len=20) :: filename(2)
  real(dp) :: xq(3,Q_LENGTH**3), max_diff, diff
  complex(dp), allocatable :: Vqqs1(:,:,:,:,:)
  complex(dp) :: Vqqs2(6,6)
  logical :: test_passed
  !
  sc_grid = [3,3,3]
  CALL read_fc2('mat2R', S, fc2)
  call aux_system(S)
  call impose_asr2('diff', S%nat, fc2)
  call div_mass_fc2(S, fc2)
  call symmetrize_system(S)
  xq = grid_vec_cart([Q_LENGTH,Q_LENGTH,Q_LENGTH], S%bg, divide=.true., natural=.true.)
  filename = ['mat2D-normal ', 'mat2D-rotated']
  !
  do i = 1, 2
    call read_fc2(trim(filename(i)), Sd(i), fc2d(i))
    ! do jn1 = 1, size(fc2d(i)%fc, 1)
    !   do jn2 = 1, jn1
    !     summ = (fc2d(i)%fc(jn1,jn2,1) + fc2d(i)%fc(jn2,jn1,1))/2
    !     fc2d(i)%fc(jn1,jn2,1) = summ
    !     fc2d(i)%fc(jn2,jn1,1) = summ
    !   enddo
    ! enddo
    call aux_system(Sd(i))
    DRR = fc_sc2RR(sc_grid, S, Sd(i), fc2d(i)%fc(:,:,1))
    call asr3(DRR)
    fc2d(i)%fc(:,:,1) = fc_RR2sc(sc_grid, S, Sd(i), DRR)
    call div_mass_fc2(Sd(i), fc2d(i))
    call fc2_sc(i)%allocate(S, Sd(i), sc_grid)
    fc2_sc(i)%fc = fc_sc2RR(sc_grid, S, Sd(i), fc2d(i)%fc) - fc_uc2RR(fc2)
  enddo
  i = 2
  !
  ! print*, fc2_sc(1)%fc(1,1,1,1), fc2_sc(2)%fc(4,4,1,1)
  ! print*, fc2_sc(1)%fc(1,1,1,1), fc2_sc(2)%fc(4,4,1,1)
  !
  call center_V(xq, S, Sd(1), fc2_sc(1), Vqqs1)
  !
  call fc2_sc(i)%center(sc_grid, S)
  !
  test_passed = .true.
  do iq = 1, size(xq,2)
    call fc2_sc(i)%r2q(xq(:,iq))
    do jq = 1, size(xq,2)
      call fc2_sc(i)%r2q(xq(:,jq), Vqqs2)
      if( any(abs(Vqqs1(:,:,iq,jq,2) - Vqqs2) > 1e-10_dp + 1e-1_dp * abs(Vqqs2))) then
        ! print*, maxval(abs(Vqqs1(:,:,iq,jq,2) - Vqqs2)), max(abs(Vqqs1(:,:,iq,jq,2)), abs(Vqqs2))
        test_passed = .false.
        do j=1,6
          do k=1,6
            diff = abs(Vqqs1(j,k,iq,jq,2) - Vqqs2(j,k))
            if (diff > 1e-1_dp*abs(Vqqs2(j,k))) then
              print*, "violates:", j,k,diff,1e-1_dp*abs(Vqqs2(j,k))
            end if
          end do
        end do
      endif
    enddo
  enddo
  !
  call test_log("dca_center_V", test_passed)
  !
end program
