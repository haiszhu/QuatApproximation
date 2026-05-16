! quatapproximation_mex.f90
! Standalone (non-module) MEX-facing wrappers for quatapproximation_mod.
! Universal qaq_ prefix per design spec §9.

! ------------------------------------------------------------------
! qaq_gauss_mex
! Thin r64 wrapper for quatapproximation_mod::gauss_r64.
! ------------------------------------------------------------------
subroutine qaq_gauss_mex(n, tgl, wgl, Dgl)
  use quatapproximation_mod, only: gauss_r64
  implicit none
  integer(8), intent(in)    :: n
  real(8),    intent(inout) :: tgl(n), wgl(n), Dgl(n, n)

  call gauss_r64(n, tgl, wgl, Dgl)
end subroutine qaq_gauss_mex

! ------------------------------------------------------------------
! qaq_bclaginterpweights_mex
! Thin r64 wrapper for quatapproximation_mod::bclaginterpweights_r64.
! ------------------------------------------------------------------
subroutine qaq_bclaginterpweights_mex(n, tgl, w_bclag)
  use quatapproximation_mod, only: bclaginterpweights_r64
  implicit none
  integer(8), intent(in)    :: n
  real(8),    intent(in)    :: tgl(n)
  real(8),    intent(inout) :: w_bclag(n)

  call bclaginterpweights_r64(n, tgl, w_bclag)
end subroutine qaq_bclaginterpweights_mex

! ------------------------------------------------------------------
! qaq_simplex2equil_with_detj_mex
! Same map as qaq_simplex2equil_mex but additionally returns detJ.
! ------------------------------------------------------------------
subroutine qaq_simplex2equil_with_detj_mex(npts, h_side, uvs, uvs_eq, detJ)
  use quatapproximation_mod, only: simplex2equil_with_detj_r64
  implicit none
  integer(8), intent(in)    :: npts
  real(8),    intent(in)    :: h_side
  real(8),    intent(in)    :: uvs(2, npts)
  real(8),    intent(inout) :: uvs_eq(2, npts)
  real(8),    intent(inout) :: detJ

  call simplex2equil_with_detj_r64(npts, h_side, uvs, uvs_eq, detJ)
end subroutine qaq_simplex2equil_with_detj_mex

! ------------------------------------------------------------------
! qaq_paraboloidparam_mex
! Downward-opening paraboloid eval with derivatives + normal + speed.
! ------------------------------------------------------------------
subroutine qaq_paraboloidparam_mex(npts, h, a, t, p, x, nx, sp, &
                                    rts, rps, rtts, rpps, rtps)
  use quatapproximation_mod, only: paraboloidparam_r64
  implicit none
  integer(8), intent(in)    :: npts
  real(8),    intent(in)    :: h, a
  real(8),    intent(in)    :: t(npts), p(npts)
  real(8),    intent(inout) :: x(3, npts), nx(3, npts), sp(npts)
  real(8),    intent(inout) :: rts(3, npts), rps(3, npts)
  real(8),    intent(inout) :: rtts(3, npts), rpps(3, npts), rtps(3, npts)

  call paraboloidparam_r64(npts, h, a, t, p, x, nx, sp, &
                            rts, rps, rtts, rpps, rtps)
end subroutine qaq_paraboloidparam_mex

! ------------------------------------------------------------------
! qaq_subdivide_simplex_mex
! Single uniform-refinement step: parent (Cx, Cy, scale, orientation)
! -> 4 child triangles as a (4, 4) matrix [Cx_i, Cy_i, scale_i, ori_i].
! ------------------------------------------------------------------
subroutine qaq_subdivide_simplex_mex(Cx, Cy, scale, orientation, sub_tris)
  use quatapproximation_mod, only: subdivide_simplex_r64
  implicit none
  real(8), intent(in)    :: Cx, Cy, scale, orientation
  real(8), intent(inout) :: sub_tris(4, 4)

  call subdivide_simplex_r64(Cx, Cy, scale, orientation, sub_tris)
end subroutine qaq_subdivide_simplex_mex

! ------------------------------------------------------------------
! qaq_assemble_subdivided_nodes_mex
! Lay out N reference simplex nodes inside each of M subdivided
! triangles; output is the flat (2, N*M) concatenation.
! ------------------------------------------------------------------
subroutine qaq_assemble_subdivided_nodes_mex(N, M, uvs_simplex, sub_tris, uvs_all)
  use quatapproximation_mod, only: assemble_subdivided_nodes_r64
  implicit none
  integer(8), intent(in)    :: N, M
  real(8),    intent(in)    :: uvs_simplex(2, N)
  real(8),    intent(in)    :: sub_tris(M, 4)
  real(8),    intent(inout) :: uvs_all(2, N*M)

  call assemble_subdivided_nodes_r64(N, M, uvs_simplex, sub_tris, uvs_all)
end subroutine qaq_assemble_subdivided_nodes_mex

! ------------------------------------------------------------------
! qaq_get_subdivided_triangles_mex
! Recursively apply subdivide_simplex L times.  Caller pre-allocates
! all_tris as (4^L, 4).
! ------------------------------------------------------------------
subroutine qaq_get_subdivided_triangles_mex(C0, s0, o0, L, all_tris)
  use quatapproximation_mod, only: get_subdivided_triangles_r64
  implicit none
  real(8),    intent(in)    :: C0(2), s0, o0
  integer(8), intent(in)    :: L
  real(8),    intent(inout) :: all_tris(4**L, 4)

  call get_subdivided_triangles_r64(C0, s0, o0, L, all_tris)
end subroutine qaq_get_subdivided_triangles_mex

! ------------------------------------------------------------------
! qaq_simplex2equil_mex
! Thin wrapper for quatapproximation_mod::simplex2equil_r64.
! ------------------------------------------------------------------
subroutine qaq_simplex2equil_mex(npts, h_side, uvs, uvs_eq)
  use quatapproximation_mod, only: simplex2equil_r64
  implicit none
  integer(8), intent(in)    :: npts
  real(8),    intent(in)    :: h_side
  real(8),    intent(in)    :: uvs(2, npts)
  real(8),    intent(inout) :: uvs_eq(2, npts)

  call simplex2equil_r64(npts, h_side, uvs, uvs_eq)
end subroutine qaq_simplex2equil_mex

! ------------------------------------------------------------------
! qaq_paraboloid_uv_mex
! Thin wrapper for quatapproximation_mod::paraboloid_uv_r64.
! ------------------------------------------------------------------
subroutine qaq_paraboloid_uv_mex(npts, H, u, v, x)
  use quatapproximation_mod, only: paraboloid_uv_r64
  implicit none
  integer(8), intent(in)    :: npts
  real(8),    intent(in)    :: H
  real(8),    intent(in)    :: u(npts), v(npts)
  real(8),    intent(inout) :: x(3, npts)

  call paraboloid_uv_r64(npts, H, u, v, x)
end subroutine qaq_paraboloid_uv_mex
