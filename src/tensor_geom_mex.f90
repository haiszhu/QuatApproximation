! tensor_geom_mex.f90 -- top-level mwrap wrappers around tensor_geom_mod.
! Symbol prefix: qatg_  (QuatApproximation, tensor geometry; parallels qakg_).

subroutine qatg_line3quadr_3dline_T_mex(x_uvs, ordert, nquad,    &
                                         tgl, wgl, Dgl, sbdnp, tpan, nbd, &
                                         sxbd, swbd, stangbd, sspbd)
  use quatapproximation_mod, only: r64
  use tensor_geom_mod, only: line3quadr_3dline_T
  implicit none
  integer(8), intent(in)    :: ordert, nquad, sbdnp, nbd
  real(r64),  intent(in)    :: x_uvs(3, ordert*ordert)
  real(r64),  intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad, nquad)
  real(r64),  intent(in)    :: tpan(sbdnp+1)
  real(r64),  intent(inout) :: sxbd(3, nbd), swbd(nbd)
  real(r64),  intent(inout) :: stangbd(3, nbd), sspbd(nbd)
  call line3quadr_3dline_T(x_uvs, ordert, nquad, tgl, wgl, Dgl, &
                           sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd)
end subroutine qatg_line3quadr_3dline_T_mex
