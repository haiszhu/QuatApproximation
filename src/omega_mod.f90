module omega_mod
  ! omega^{(n,m)}_i boundary-integrand assembly. Mirrors
  ! qotential/utils/omeganm_i.m at r64.
  !
  ! Source MATLAB (canonical):
  !   onm =  (dx.*r3)'.*q_j - (dx.*r2)'.*q_k
  !         +(dy.*r1)'.*q_k - (dy.*r3)'.*q_i
  !         +(dz.*r2)'.*q_i - (dz.*r1)'.*q_j
  !
  ! Equivalent grouped form (Fortran body below):
  !   coeff(q_i) = dz*r2 - dy*r3
  !   coeff(q_j) = dx*r3 - dz*r1
  !   coeff(q_k) = dy*r1 - dx*r2
  use quatapproximation_mod, only: r64, r128
  implicit none
  private
  public :: qao_omeganm_i_r64, qao_omegaall_r64
  public :: qao_omeganm_i_r128, qao_omegaall_r128

contains

  subroutine qao_omeganm_i_r64(N, ncoeff, nslots, r, dr, q_i, q_j, q_k, onm)
    integer(8), intent(in)  :: N, ncoeff, nslots
    real(r64),  intent(in)  :: r (3, N), dr(3, N)
    real(r64),  intent(in)  :: q_i(N, ncoeff, nslots)
    real(r64),  intent(in)  :: q_j(N, ncoeff, nslots)
    real(r64),  intent(in)  :: q_k(N, ncoeff, nslots)
    real(r64),  intent(out) :: onm(N, ncoeff, nslots)
    real(r64)  :: a1, a2, a3
    integer(8) :: n_, c, k
    do k = 1, nslots
      do c = 1, ncoeff
        do n_ = 1, N
          ! coefficient of q_i: dz*r2 - dy*r3
          a1 = dr(3, n_)*r(2, n_) - dr(2, n_)*r(3, n_)
          ! coefficient of q_j: dx*r3 - dz*r1
          a2 = dr(1, n_)*r(3, n_) - dr(3, n_)*r(1, n_)
          ! coefficient of q_k: dy*r1 - dx*r2
          a3 = dr(2, n_)*r(1, n_) - dr(1, n_)*r(2, n_)
          onm(n_, c, k) = a1*q_i(n_, c, k) + a2*q_j(n_, c, k) + a3*q_k(n_, c, k)
        end do
      end do
    end do
  end subroutine qao_omeganm_i_r64

  ! ------------------------------------------------------------------
  ! qao_omegaall_r64
  ! Target-dependent assembly of the close-eval Omega tensor across all
  ! four kernel coefficient channels (omega0..omega3).  Mirrors
  ! qotential/utils/f/moments.f::omegaall at r64; one OpenMP-parallel
  ! loop over targets, no algorithm change.
  !
  ! Inputs:
  !   r0   (3, m)             target points
  !   M_all (dim1, m)         flattened (n, morder) moments per target
  !                           with dim1 = n*morder, morder = 2*order+2
  !   onm0..onm3 (n*h_dim, 4) flattened (n, h_dim, 4) omega^nm tensors
  !                           from qao_omeganm_i_r64; the 4 trailing
  !                           slots are the affine-in-target coefficients
  !   ijIdx (2, h_dim)        (i, j) basis indices for each column
  !
  ! Output:
  !   omega (m, 4*h_dim)      Omega per target, concatenated as
  !                           [omega0 | -omega1 | -omega2 | -omega3]
  ! ------------------------------------------------------------------
  subroutine qao_omegaall_r64(m, dim1, n, h_dim, morder, r0, M_all, &
                              onm0, onm1, onm2, onm3, ijIdx, omega)
    integer(8), intent(in)  :: m, dim1, n, h_dim, morder
    real(r64),  intent(in)  :: r0(3, m)
    real(r64),  intent(in)  :: M_all(dim1, m)
    real(r64),  intent(in)  :: onm0(n*h_dim, 4), onm1(n*h_dim, 4)
    real(r64),  intent(in)  :: onm2(n*h_dim, 4), onm3(n*h_dim, 4)
    integer(8), intent(in)  :: ijIdx(2, h_dim)
    real(r64),  intent(out) :: omega(m, 4*h_dim)

    integer(8) :: ijIdxsum(h_dim)
    integer(8) :: j, k
    real(r64)  :: r0_j(3)
    real(r64)  :: mk_j(n, morder)
    real(r64)  :: mkp1(n, h_dim), mkp0(n, h_dim), otmp(n, h_dim)
    real(r64)  :: omega0_j(h_dim), omega1_j(h_dim)
    real(r64)  :: omega2_j(h_dim), omega3_j(h_dim)
    real(r64)  :: onm0_rs(n, h_dim, 4), onm1_rs(n, h_dim, 4)
    real(r64)  :: onm2_rs(n, h_dim, 4), onm3_rs(n, h_dim, 4)
    real(r64), parameter :: PI = 4.0_r64 * atan(1.0_r64)
    real(r64), parameter :: COEFF = -1.0_r64 / (4.0_r64 * PI)

    ijIdxsum = sum(ijIdx, dim=1)

    onm0_rs = reshape(onm0, (/n, h_dim, 4_8/))
    onm1_rs = reshape(onm1, (/n, h_dim, 4_8/))
    onm2_rs = reshape(onm2, (/n, h_dim, 4_8/))
    onm3_rs = reshape(onm3, (/n, h_dim, 4_8/))

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(j, r0_j, mk_j, mkp1, mkp0, otmp) &
!$OMP   PRIVATE(omega0_j, omega1_j, omega2_j, omega3_j)
    do j = 1, m
      r0_j = r0(:, j)
      mk_j = reshape(M_all(:, j), (/n, morder/))

      mkp1 = mk_j(:, ijIdxsum + 3_8)
      mkp0 = mk_j(:, ijIdxsum + 2_8)

      otmp = onm0_rs(:, :, 2)*r0_j(1) + onm0_rs(:, :, 3)*r0_j(2) &
           + onm0_rs(:, :, 4)*r0_j(3)
      omega0_j = COEFF * sum(mkp1*onm0_rs(:, :, 1) + mkp0*otmp, dim=1)

      otmp = onm1_rs(:, :, 2)*r0_j(1) + onm1_rs(:, :, 3)*r0_j(2) &
           + onm1_rs(:, :, 4)*r0_j(3)
      omega1_j = COEFF * sum(mkp1*onm1_rs(:, :, 1) + mkp0*otmp, dim=1)

      otmp = onm2_rs(:, :, 2)*r0_j(1) + onm2_rs(:, :, 3)*r0_j(2) &
           + onm2_rs(:, :, 4)*r0_j(3)
      omega2_j = COEFF * sum(mkp1*onm2_rs(:, :, 1) + mkp0*otmp, dim=1)

      otmp = onm3_rs(:, :, 2)*r0_j(1) + onm3_rs(:, :, 3)*r0_j(2) &
           + onm3_rs(:, :, 4)*r0_j(3)
      omega3_j = COEFF * sum(mkp1*onm3_rs(:, :, 1) + mkp0*otmp, dim=1)

      do k = 1, h_dim
        omega(j, k          ) =  omega0_j(k)
        omega(j, k +   h_dim) = -omega1_j(k)
        omega(j, k + 2*h_dim) = -omega2_j(k)
        omega(j, k + 3*h_dim) = -omega3_j(k)
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine qao_omegaall_r64

  ! ================================================================
  ! r128 siblings -- bit-for-bit mirrors of the r64 routines above.
  ! No mex wrappers; intended for the future Lap3dDLP_closepanel_r128
  ! orchestration.
  ! ================================================================

  subroutine qao_omeganm_i_r128(N, ncoeff, nslots, r, dr, q_i, q_j, q_k, onm)
    integer(8), intent(in)  :: N, ncoeff, nslots
    real(r128), intent(in)  :: r (3, N), dr(3, N)
    real(r128), intent(in)  :: q_i(N, ncoeff, nslots)
    real(r128), intent(in)  :: q_j(N, ncoeff, nslots)
    real(r128), intent(in)  :: q_k(N, ncoeff, nslots)
    real(r128), intent(out) :: onm(N, ncoeff, nslots)
    real(r128) :: a1, a2, a3
    integer(8) :: n_, c, kk
    do kk = 1, nslots
      do c = 1, ncoeff
        do n_ = 1, N
          a1 = dr(3, n_)*r(2, n_) - dr(2, n_)*r(3, n_)
          a2 = dr(1, n_)*r(3, n_) - dr(3, n_)*r(1, n_)
          a3 = dr(2, n_)*r(1, n_) - dr(1, n_)*r(2, n_)
          onm(n_, c, kk) = a1*q_i(n_, c, kk) + a2*q_j(n_, c, kk) + a3*q_k(n_, c, kk)
        end do
      end do
    end do
  end subroutine qao_omeganm_i_r128

  subroutine qao_omegaall_r128(m, dim1, n, h_dim, morder, r0, M_all, &
                               onm0, onm1, onm2, onm3, ijIdx, omega)
    integer(8), intent(in)  :: m, dim1, n, h_dim, morder
    real(r128), intent(in)  :: r0(3, m)
    real(r128), intent(in)  :: M_all(dim1, m)
    real(r128), intent(in)  :: onm0(n*h_dim, 4), onm1(n*h_dim, 4)
    real(r128), intent(in)  :: onm2(n*h_dim, 4), onm3(n*h_dim, 4)
    integer(8), intent(in)  :: ijIdx(2, h_dim)
    real(r128), intent(out) :: omega(m, 4*h_dim)

    integer(8) :: ijIdxsum(h_dim)
    integer(8) :: j, kk
    real(r128) :: r0_j(3)
    real(r128) :: mk_j(n, morder)
    real(r128) :: mkp1(n, h_dim), mkp0(n, h_dim), otmp(n, h_dim)
    real(r128) :: omega0_j(h_dim), omega1_j(h_dim)
    real(r128) :: omega2_j(h_dim), omega3_j(h_dim)
    real(r128) :: onm0_rs(n, h_dim, 4), onm1_rs(n, h_dim, 4)
    real(r128) :: onm2_rs(n, h_dim, 4), onm3_rs(n, h_dim, 4)
    real(r128), parameter :: PI    = 4.0_r128 * atan(1.0_r128)
    real(r128), parameter :: COEFF = -1.0_r128 / (4.0_r128 * PI)

    ijIdxsum = sum(ijIdx, dim=1)

    onm0_rs = reshape(onm0, (/n, h_dim, 4_8/))
    onm1_rs = reshape(onm1, (/n, h_dim, 4_8/))
    onm2_rs = reshape(onm2, (/n, h_dim, 4_8/))
    onm3_rs = reshape(onm3, (/n, h_dim, 4_8/))

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(j, r0_j, mk_j, mkp1, mkp0, otmp) &
!$OMP   PRIVATE(omega0_j, omega1_j, omega2_j, omega3_j)
    do j = 1, m
      r0_j = r0(:, j)
      mk_j = reshape(M_all(:, j), (/n, morder/))

      mkp1 = mk_j(:, ijIdxsum + 3_8)
      mkp0 = mk_j(:, ijIdxsum + 2_8)

      otmp = onm0_rs(:, :, 2)*r0_j(1) + onm0_rs(:, :, 3)*r0_j(2) &
           + onm0_rs(:, :, 4)*r0_j(3)
      omega0_j = COEFF * sum(mkp1*onm0_rs(:, :, 1) + mkp0*otmp, dim=1)

      otmp = onm1_rs(:, :, 2)*r0_j(1) + onm1_rs(:, :, 3)*r0_j(2) &
           + onm1_rs(:, :, 4)*r0_j(3)
      omega1_j = COEFF * sum(mkp1*onm1_rs(:, :, 1) + mkp0*otmp, dim=1)

      otmp = onm2_rs(:, :, 2)*r0_j(1) + onm2_rs(:, :, 3)*r0_j(2) &
           + onm2_rs(:, :, 4)*r0_j(3)
      omega2_j = COEFF * sum(mkp1*onm2_rs(:, :, 1) + mkp0*otmp, dim=1)

      otmp = onm3_rs(:, :, 2)*r0_j(1) + onm3_rs(:, :, 3)*r0_j(2) &
           + onm3_rs(:, :, 4)*r0_j(3)
      omega3_j = COEFF * sum(mkp1*onm3_rs(:, :, 1) + mkp0*otmp, dim=1)

      do kk = 1, h_dim
        omega(j, kk          ) =  omega0_j(kk)
        omega(j, kk +   h_dim) = -omega1_j(kk)
        omega(j, kk + 2*h_dim) = -omega2_j(kk)
        omega(j, kk + 3*h_dim) = -omega3_j(kk)
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine qao_omegaall_r128

end module omega_mod
