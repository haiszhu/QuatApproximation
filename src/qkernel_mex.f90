! qkernel_mex.f90 -- top-level mwrap wrappers around qkernel_mod.
! Symbol prefix: qak_ (QuatApproximation, kernel).

subroutine qak_qnm_i_mex(idx, N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
  use quatapproximation_mod, only: r64
  use qkernel_mod, only: qak_qnm_i_r64
  implicit none
  integer(8), intent(in)    :: idx, N, ncoeff, lptype_id
  real(r64),  intent(in)    :: sx(3, N)
  real(r64),  intent(in)    :: F (N, ncoeff)
  real(r64),  intent(in)    :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
  real(r64),  intent(in)    :: gradxyz(3, N)
  real(r64),  intent(inout) :: q_i(N, ncoeff, 5)
  real(r64),  intent(inout) :: q_j(N, ncoeff, 5)
  real(r64),  intent(inout) :: q_k(N, ncoeff, 5)
  call qak_qnm_i_r64(idx, N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
end subroutine qak_qnm_i_mex
