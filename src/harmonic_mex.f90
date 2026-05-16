! harmonic_mex.f90
! Standalone (non-module) MEX-facing wrappers for harmonic_mod.
!
! gfortran exports module procedures with module-mangled names
! (__harmonic_mod_MOD_*), but mwrap binds to plain external symbols.
! Each wrapper here is a file-scope subroutine that imports the module
! procedure under a local alias and delegates.  All wrappers carry the
! universal qah_<routine>_mex prefix per the design spec §9 (Hard rule:
! zero un-prefixed top-level _mex.f90 wrappers).

! ------------------------------------------------------------------
! qah_evaltensorproductharmonicgrad_mex
! Thin wrapper for harmonic_mod::evaltensorproductharmonicgrad_r64.
! Tensor-product harmonic polynomial basis values f + Cartesian gradient
! (fx, fy, fz) at nt 3D points, plus the (i, j) index pairs in ijidx.
! ------------------------------------------------------------------
subroutine qah_evaltensorproductharmonicgrad_mex(nt, r, order, fx, fy, fz, f, ijidx)
  use harmonic_mod, only: evaltensorproductharmonicgrad_r64
  implicit none
  integer(8), intent(in)    :: nt, order
  real(8),    intent(in)    :: r(3, nt)
  real(8),    intent(inout) :: fx(nt, order*order), fy(nt, order*order), &
                               fz(nt, order*order), f (nt, order*order)
  integer(8), intent(inout) :: ijidx(2, order*order)

  call evaltensorproductharmonicgrad_r64(nt, r, order, fx, fy, fz, f, ijidx)
end subroutine qah_evaltensorproductharmonicgrad_mex

! ------------------------------------------------------------------
! qah_l3dtavecevalmat_mex
! Thin r64 wrapper for harmonic_mod::l3dtavecevalmat_r64.
! ------------------------------------------------------------------
subroutine qah_l3dtavecevalmat_mex(ztargs, nt, nterms, F, Fx, Fy, Fz, ier)
  use harmonic_mod, only: l3dtavecevalmat_r64
  implicit none
  integer(8),  intent(in)    :: nt, nterms
  integer(8),  intent(inout) :: ier
  real(8),     intent(in)    :: ztargs(3, nt)
  complex(8),  intent(inout) :: F (nt, (nterms+1)**2)
  complex(8),  intent(inout) :: Fx(nt, (nterms+1)**2)
  complex(8),  intent(inout) :: Fy(nt, (nterms+1)**2)
  complex(8),  intent(inout) :: Fz(nt, (nterms+1)**2)

  call l3dtavecevalmat_r64(ztargs, nt, nterms, F, Fx, Fy, Fz, ier)
end subroutine qah_l3dtavecevalmat_mex

! ------------------------------------------------------------------
! qah_l3dtavecevalmat_r128_mex
! r128 evaluation with the design-spec §10 trailing-flag contract.
!
!   flag == 0 : cast r64 inputs up to r128, run harmonic_mod::l3dtavecevalmat_r128,
!               cast r128 outputs back down to r64 for return to MATLAB.
!   flag == 1 : reserved for HDF5 round-trip with a downstream r128
!               consumer.  Layout of qah_l3dtavecevalmat_r128.h5 is left
!               undefined here so it can be co-designed with the first
!               consumer (likely qak_q0kl_r128_mex or similar).  Calling
!               with flag=1 currently aborts.
!
! MATLAB cannot carry r128 complex across mwrap, so even on flag=0 the
! per-element output ends up cast back to complex(8); use the r64 path
! for fast iteration and the (future) flag=1 chain when r128 fidelity is
! required through to a downstream consumer.
! ------------------------------------------------------------------
subroutine qah_l3dtavecevalmat_r128_mex(ztargs, nt, nterms, F, Fx, Fy, Fz, ier, flag)
  use harmonic_mod, only: l3dtavecevalmat_r128, r128, c128
  implicit none
  integer(8), intent(in)    :: nt, nterms, flag
  integer(8), intent(inout) :: ier
  real(8),    intent(in)    :: ztargs(3, nt)
  complex(8), intent(inout) :: F (nt, (nterms+1)**2)
  complex(8), intent(inout) :: Fx(nt, (nterms+1)**2)
  complex(8), intent(inout) :: Fy(nt, (nterms+1)**2)
  complex(8), intent(inout) :: Fz(nt, (nterms+1)**2)

  real(r128)    :: ztargs_r128(3, nt)
  complex(c128) :: F_r128 (nt, (nterms+1)**2)
  complex(c128) :: Fx_r128(nt, (nterms+1)**2)
  complex(c128) :: Fy_r128(nt, (nterms+1)**2)
  complex(c128) :: Fz_r128(nt, (nterms+1)**2)

  if (flag == 0_8) then
    ztargs_r128 = real(ztargs, r128)
    F_r128  = cmplx(F,  kind=c128)
    Fx_r128 = cmplx(Fx, kind=c128)
    Fy_r128 = cmplx(Fy, kind=c128)
    Fz_r128 = cmplx(Fz, kind=c128)

    call l3dtavecevalmat_r128(ztargs_r128, nt, nterms, &
                              F_r128, Fx_r128, Fy_r128, Fz_r128, ier)

    F  = cmplx(F_r128,  kind=8)
    Fx = cmplx(Fx_r128, kind=8)
    Fy = cmplx(Fy_r128, kind=8)
    Fz = cmplx(Fz_r128, kind=8)
  else
    error stop 'qah_l3dtavecevalmat_r128_mex: flag=1 path not yet implemented &
                &(HDF5 dataset layout pending downstream r128 consumer design)'
  end if
end subroutine qah_l3dtavecevalmat_r128_mex
