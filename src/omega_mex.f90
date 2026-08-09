! omega_mex.f90 -- top-level mwrap wrapper around omega_mod.
! Symbol prefix: qao_ (QuatApproximation, omega).

subroutine qao_omeganm_i_mex(N, ncoeff, nslots, r, dr, q_i, q_j, q_k, onm)
  use quatapproximation_mod, only: r64
  use omega_mod, only: qao_omeganm_i_r64
  implicit none
  integer(8), intent(in)    :: N, ncoeff, nslots
  real(r64),  intent(in)    :: r (3, N), dr(3, N)
  real(r64),  intent(in)    :: q_i(N, ncoeff, nslots)
  real(r64),  intent(in)    :: q_j(N, ncoeff, nslots)
  real(r64),  intent(in)    :: q_k(N, ncoeff, nslots)
  real(r64),  intent(inout) :: onm(N, ncoeff, nslots)
  call qao_omeganm_i_r64(N, ncoeff, nslots, r, dr, q_i, q_j, q_k, onm)
end subroutine qao_omeganm_i_mex

subroutine qao_omegaall_mex(m, dim1, n, h_dim, morder, r0, M_all, &
                            onm0, onm1, onm2, onm3, ijIdx, omega)
  use quatapproximation_mod, only: r64
  use omega_mod, only: qao_omegaall_r64
  implicit none
  integer(8), intent(in)    :: m, dim1, n, h_dim, morder
  real(r64),  intent(in)    :: r0(3, m)
  real(r64),  intent(in)    :: M_all(dim1, m)
  real(r64),  intent(in)    :: onm0(n*h_dim, 4), onm1(n*h_dim, 4)
  real(r64),  intent(in)    :: onm2(n*h_dim, 4), onm3(n*h_dim, 4)
  integer(8), intent(in)    :: ijIdx(2, h_dim)
  real(r64),  intent(inout) :: omega(m, 4*h_dim)
  call qao_omegaall_r64(m, dim1, n, h_dim, morder, r0, M_all, &
                        onm0, onm1, onm2, onm3, ijIdx, omega)
end subroutine qao_omegaall_mex

subroutine qao_omegasdlp_mex(m, nterms, ncoeff, h_dim, r0, Ichi, &
                             Ialpha, omega_slp, omega)
  use quatapproximation_mod, only: r64
  use omega_mod, only: qao_omegasdlp_r64
  implicit none
  integer(8), intent(in) :: m, nterms, ncoeff, h_dim
  real(r64), intent(in) :: r0(3,m)
  complex(r64), intent(in) :: Ichi(m,ncoeff,4), Ialpha(m,ncoeff,4)
  real(r64), intent(out) :: omega_slp(h_dim,m), omega(h_dim,m,4)
  call qao_omegasdlp_r64(m, nterms, ncoeff, h_dim, r0, Ichi, Ialpha, &
                         omega_slp, omega)
end subroutine qao_omegasdlp_mex
