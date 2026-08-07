! tensor_geom_mod.f90
! Module: tensor-product (square / [-1,1]^2) reference-element geometry.
! Parallel to koorn_geom_mod (which covers Koornwinder / triangle).
!
! Public:
!   line3quadr_3dline_T(x_uvs, ordert, nquad, tgl, wgl, Dgl,
!                       sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd)
!     Tensor-product analog of koorn_geom_mod::line3quadr_3dline.
!
module tensor_geom_mod
#ifdef BIESOLVER_R64_ONLY
  use quatapproximation_mod, only: r64, r128, gauss_r64
#else
  use quatapproximation_mod, only: r64, r128, gauss_r64, gauss_r128
#endif
  implicit none
  private
  public :: line3quadr_3dline_T
#ifndef BIESOLVER_R64_ONLY
  public :: line3quadr_3dline_T_r128
#endif

contains

  ! ----------------------------------------------------------------
  ! line3quadr_3dline_T
  ! Generate GL boundary nodes/weights/tangents/speeds for a
  ! tensor-product (square) source patch's boundary.
  !
  ! Inputs:
  !   x_uvs(3, ordert*ordert)  source patch values at GL tensor grid
  !   ordert                   1D tensor order
  !   nquad                    GL order per panel
  !   tgl, wgl, Dgl            nquad GL setup
  !   sbdnp                    number of panels along [0, 2pi]
  !   tpan(sbdnp+1)            panel breakpoints in [0, 2pi]
  !   nbd = sbdnp*nquad
  !
  ! Outputs:
  !   sxbd(3, nbd)             boundary point positions
  !   swbd(nbd)                speed-weighted GL quadrature weights
  !   stangbd(3, nbd)          unit tangent vectors
  !   sspbd(nbd)               speeds |x'(t)|
  !
  ! Algorithm (analog of koorn version):
  !   1. Build 1D Legendre vals -> coeffs matrix Legmat for ordert GL nodes.
  !   2. coefs(:,:,c) = Legmat * x_grid(:,:,c) * Legmat^T   (2D vals -> coeffs)
  !   3. Per panel: GL t-nodes -> (u, v) on box boundary via tparam_to_uv_square
  !      -> P_n(u), P_m(v) via 3-term recurrence -> contract with coefs -> sxbd.
  !   4. xp via Dgl per panel; then sspbd, stangbd, swbd.
  ! ----------------------------------------------------------------
  subroutine line3quadr_3dline_T(x_uvs, ordert, nquad,    &
                                  tgl, wgl, Dgl, sbdnp, tpan, nbd, &
                                  sxbd, swbd, stangbd, sspbd)
    integer(8), intent(in)    :: ordert, nquad, sbdnp, nbd
    real(r64),  intent(in)    :: x_uvs(3, ordert*ordert)
    real(r64),  intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r64),  intent(in)    :: tpan(sbdnp+1)
    real(r64),  intent(inout) :: sxbd(3,nbd), swbd(nbd)
    real(r64),  intent(inout) :: stangbd(3,nbd), sspbd(nbd)

    real(r64)  :: Legmat(ordert, ordert)
    real(r64)  :: tgl_t(ordert), wgl_t(ordert), Dgl_t(ordert, ordert)
    real(r64)  :: x_grid(ordert, ordert, 3)
    real(r64)  :: coefs (ordert, ordert, 3)
    real(r64)  :: tlo(sbdnp), thi(sbdnp), pt(sbdnp)
    real(r64)  :: pu(ordert), pv(ordert)
    real(r64)  :: t(nbd), sxp(3, nbd), sx_ell(3, nquad)
    real(r64)  :: uv(2), tmp_mat(ordert, ordert), accum
    integer(8) :: i, j, k, c, ell, ii, itmp, ii_start, ii_end

    ! ---- 1. 1D Legendre vals -> coeffs matrix at ordert GL nodes
    call gauss_r64(ordert, tgl_t, wgl_t, Dgl_t)
    do j = 1, ordert
      call legepols(tgl_t(j), ordert, pu)        ! pu(k) = P_{k-1}(tgl_t(j))
      do k = 1, ordert
        Legmat(k, j) = 0.5_r64 * (2.0_r64*real(k, r64) - 1.0_r64) * wgl_t(j) * pu(k)
      end do
    end do

    ! ---- 2. Reshape x_uvs and project to 2D Legendre coeffs per xyz
    do c = 1, 3
      do j = 1, ordert
        do i = 1, ordert
          x_grid(i, j, c) = x_uvs(c, (j-1)*ordert + i)
        end do
      end do
    end do
    do c = 1, 3
      tmp_mat        = matmul(Legmat, x_grid(:, :, c))
      coefs(:, :, c) = matmul(tmp_mat, transpose(Legmat))
    end do

    ! ---- 3. Per-panel t grid and bare GL weights (speed factor folded in below)
    tlo = tpan(1:sbdnp)
    thi = tpan(2:sbdnp+1)
    pt  = thi - tlo
    do ell = 1, sbdnp
      itmp = (ell-1)*nquad
      do i = 1, nquad
        ii        = itmp + i
        t(ii)     = tlo(ell) + 0.5_r64*(1.0_r64 + tgl(i)) * pt(ell)
        swbd(ii)  = 0.5_r64 * wgl(i) * pt(ell)
      end do
    end do

    ! ---- 4. Boundary chart + 2D Legendre evaluation per node.
    !       coefs is indexed (n_fast, n_slow): first axis = degree in the fast
    !       direction (matches MATLAB's meshgrid where x2 varies fast), second
    !       axis = degree in the slow direction (x1).  uv(1) is the slow coord,
    !       uv(2) is the fast coord (mirrors MATLAB rboxZ -> [real; imag]).
    sxbd = 0.0_r64
    do k = 1, nbd
      call tparam_to_uv_square(t(k), uv)
      call legepols(uv(1), ordert, pu)   ! pu = P_n at slow coord (uv(1))
      call legepols(uv(2), ordert, pv)   ! pv = P_n at fast coord (uv(2))
      do c = 1, 3
        accum = 0.0_r64
        do j = 1, ordert
          do i = 1, ordert
            accum = accum + coefs(i, j, c) * pv(i) * pu(j)
          end do
        end do
        sxbd(c, k) = accum
      end do
    end do

    ! ---- 5. Derivative via Dgl per panel, then sspbd, stangbd, swbd
    sxp = 0.0_r64
    do ell = 1, sbdnp
      ii_start = (ell-1)*nquad + 1
      ii_end   = ell*nquad
      sx_ell   = sxbd(:, ii_start:ii_end)
      sxp(1, ii_start:ii_end) = (2.0_r64/pt(ell)) * matmul(Dgl, sx_ell(1,:))
      sxp(2, ii_start:ii_end) = (2.0_r64/pt(ell)) * matmul(Dgl, sx_ell(2,:))
      sxp(3, ii_start:ii_end) = (2.0_r64/pt(ell)) * matmul(Dgl, sx_ell(3,:))
    end do
    sspbd = sqrt(sxp(1,:)**2 + sxp(2,:)**2 + sxp(3,:)**2)
    stangbd(1,:) = sxp(1,:) / sspbd
    stangbd(2,:) = sxp(2,:) / sspbd
    stangbd(3,:) = sxp(3,:) / sspbd
    swbd = swbd * sspbd

  end subroutine line3quadr_3dline_T

  ! ----------------------------------------------------------------
  ! legepols (private)
  ! P_0(x)..P_{n-1}(x) via the standard 3-term recurrence.
  ! ----------------------------------------------------------------
  subroutine legepols(x, n, pols)
    real(r64), intent(in)  :: x
    integer(8), intent(in) :: n
    real(r64), intent(out) :: pols(n)
    integer(8) :: k
    if (n >= 1) pols(1) = 1.0_r64
    if (n >= 2) pols(2) = x
    do k = 3, n
      pols(k) = ((2.0_r64*real(k-2,r64) + 1.0_r64) * x * pols(k-1) &
                - real(k-2,r64) * pols(k-2)) / real(k-1, r64)
    end do
  end subroutine legepols

  ! ----------------------------------------------------------------
  ! tparam_to_uv_square (private)
  ! Box boundary chart: t in [0, 2pi] -> (u, v) on the boundary of
  ! [-1, 1]^2, counterclockwise starting at (-1, -1).  Mirrors
  ! qotential's rboxZ.
  ! ----------------------------------------------------------------
  subroutine tparam_to_uv_square(t, uv)
    real(r64), intent(in)  :: t
    real(r64), intent(out) :: uv(2)
    real(r64), parameter :: PI_   = 4.0_r64*atan(1.0_r64)
    real(r64), parameter :: TWOPI = 2.0_r64*PI_
    real(r64) :: tt
    tt = mod(t, TWOPI)
    if (tt < 0.0_r64) tt = tt + TWOPI
    if (tt < 0.5_r64*PI_) then            ! bottom edge: v=-1, u in [-1, 1]
      uv(1) = -1.0_r64 + (4.0_r64/PI_)*tt
      uv(2) = -1.0_r64
    else if (tt < PI_) then               ! right edge: u=1, v in [-1, 1]
      uv(1) =  1.0_r64
      uv(2) = -1.0_r64 + (4.0_r64/PI_)*(tt - 0.5_r64*PI_)
    else if (tt < 1.5_r64*PI_) then       ! top edge: v=1, u in [1, -1]
      uv(1) =  1.0_r64 - (4.0_r64/PI_)*(tt - PI_)
      uv(2) =  1.0_r64
    else                                  ! left edge: u=-1, v in [1, -1]
      uv(1) = -1.0_r64
      uv(2) =  1.0_r64 - (4.0_r64/PI_)*(tt - 1.5_r64*PI_)
    end if
  end subroutine tparam_to_uv_square

#ifndef BIESOLVER_R64_ONLY
  ! ================================================================
  ! r128 siblings -- bit-for-bit mirror of line3quadr_3dline_T and its
  ! two private helpers (legepols, tparam_to_uv_square) at real(16).
  ! No mex wrappers; intended for the future Lap3dDLP_closepanel_r128
  ! orchestration.
  ! ================================================================

  subroutine line3quadr_3dline_T_r128(x_uvs, ordert, nquad,    &
                                       tgl, wgl, Dgl, sbdnp, tpan, nbd, &
                                       sxbd, swbd, stangbd, sspbd)
    integer(8), intent(in)    :: ordert, nquad, sbdnp, nbd
    real(r128), intent(in)    :: x_uvs(3, ordert*ordert)
    real(r128), intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r128), intent(in)    :: tpan(sbdnp+1)
    real(r128), intent(inout) :: sxbd(3,nbd), swbd(nbd)
    real(r128), intent(inout) :: stangbd(3,nbd), sspbd(nbd)

    real(r128) :: Legmat(ordert, ordert)
    real(r128) :: tgl_t(ordert), wgl_t(ordert), Dgl_t(ordert, ordert)
    real(r128) :: x_grid(ordert, ordert, 3)
    real(r128) :: coefs (ordert, ordert, 3)
    real(r128) :: tlo(sbdnp), thi(sbdnp), pt(sbdnp)
    real(r128) :: pu(ordert), pv(ordert)
    real(r128) :: t(nbd), sxp(3, nbd), sx_ell(3, nquad)
    real(r128) :: uv(2), tmp_mat(ordert, ordert), accum
    integer(8) :: i, j, k, c, ell, ii, itmp, ii_start, ii_end

    call gauss_r128(ordert, tgl_t, wgl_t, Dgl_t)
    do j = 1, ordert
      call legepols_r128(tgl_t(j), ordert, pu)
      do k = 1, ordert
        Legmat(k, j) = 0.5_r128 * (2.0_r128*real(k, r128) - 1.0_r128) * wgl_t(j) * pu(k)
      end do
    end do

    do c = 1, 3
      do j = 1, ordert
        do i = 1, ordert
          x_grid(i, j, c) = x_uvs(c, (j-1)*ordert + i)
        end do
      end do
    end do
    do c = 1, 3
      tmp_mat        = matmul(Legmat, x_grid(:, :, c))
      coefs(:, :, c) = matmul(tmp_mat, transpose(Legmat))
    end do

    tlo = tpan(1:sbdnp)
    thi = tpan(2:sbdnp+1)
    pt  = thi - tlo
    do ell = 1, sbdnp
      itmp = (ell-1)*nquad
      do i = 1, nquad
        ii        = itmp + i
        t(ii)     = tlo(ell) + 0.5_r128*(1.0_r128 + tgl(i)) * pt(ell)
        swbd(ii)  = 0.5_r128 * wgl(i) * pt(ell)
      end do
    end do

    sxbd = 0.0_r128
    do k = 1, nbd
      call tparam_to_uv_square_r128(t(k), uv)
      call legepols_r128(uv(1), ordert, pu)
      call legepols_r128(uv(2), ordert, pv)
      do c = 1, 3
        accum = 0.0_r128
        do j = 1, ordert
          do i = 1, ordert
            accum = accum + coefs(i, j, c) * pv(i) * pu(j)
          end do
        end do
        sxbd(c, k) = accum
      end do
    end do

    sxp = 0.0_r128
    do ell = 1, sbdnp
      ii_start = (ell-1)*nquad + 1
      ii_end   = ell*nquad
      sx_ell   = sxbd(:, ii_start:ii_end)
      sxp(1, ii_start:ii_end) = (2.0_r128/pt(ell)) * matmul(Dgl, sx_ell(1,:))
      sxp(2, ii_start:ii_end) = (2.0_r128/pt(ell)) * matmul(Dgl, sx_ell(2,:))
      sxp(3, ii_start:ii_end) = (2.0_r128/pt(ell)) * matmul(Dgl, sx_ell(3,:))
    end do
    sspbd = sqrt(sxp(1,:)**2 + sxp(2,:)**2 + sxp(3,:)**2)
    stangbd(1,:) = sxp(1,:) / sspbd
    stangbd(2,:) = sxp(2,:) / sspbd
    stangbd(3,:) = sxp(3,:) / sspbd
    swbd = swbd * sspbd

  end subroutine line3quadr_3dline_T_r128

  subroutine legepols_r128(x, n, pols)
    real(r128), intent(in)  :: x
    integer(8), intent(in)  :: n
    real(r128), intent(out) :: pols(n)
    integer(8) :: k
    if (n >= 1) pols(1) = 1.0_r128
    if (n >= 2) pols(2) = x
    do k = 3, n
      pols(k) = ((2.0_r128*real(k-2,r128) + 1.0_r128) * x * pols(k-1) &
                - real(k-2,r128) * pols(k-2)) / real(k-1, r128)
    end do
  end subroutine legepols_r128

  subroutine tparam_to_uv_square_r128(t, uv)
    real(r128), intent(in)  :: t
    real(r128), intent(out) :: uv(2)
    real(r128), parameter :: PI_   = 4.0_r128*atan(1.0_r128)
    real(r128), parameter :: TWOPI = 2.0_r128*PI_
    real(r128) :: tt
    tt = mod(t, TWOPI)
    if (tt < 0.0_r128) tt = tt + TWOPI
    if (tt < 0.5_r128*PI_) then
      uv(1) = -1.0_r128 + (4.0_r128/PI_)*tt
      uv(2) = -1.0_r128
    else if (tt < PI_) then
      uv(1) =  1.0_r128
      uv(2) = -1.0_r128 + (4.0_r128/PI_)*(tt - 0.5_r128*PI_)
    else if (tt < 1.5_r128*PI_) then
      uv(1) =  1.0_r128 - (4.0_r128/PI_)*(tt - PI_)
      uv(2) =  1.0_r128
    else
      uv(1) = -1.0_r128
      uv(2) =  1.0_r128 - (4.0_r128/PI_)*(tt - 1.5_r128*PI_)
    end if
  end subroutine tparam_to_uv_square_r128
#endif

end module tensor_geom_mod
