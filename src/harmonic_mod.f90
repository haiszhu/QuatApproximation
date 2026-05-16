module harmonic_mod
  ! ------------------------------------------------------------------
  ! 3D Laplace harmonic basis machinery for QuatApproximation-legacy.
  !
  ! Currently hosts the value+gradient evaluator l3dtavecevalmat plus its
  ! three private helpers (cart2polarl, ylgndrini, ylgndru2sf), ported
  ! byte-faithfully from
  !   /Users/hzhu/Documents/OSSD/RCIP/surfquadr/utils/ssq/linequad.f
  ! and converted from F77 fixed-form to F90 module style.
  !
  ! Per the parallel-independent rule (project_parallel_legacy_packages),
  ! r64/r128 kind parameters are defined locally; nothing crosses into
  ! LineQuaaadrature.
  ! ------------------------------------------------------------------
  implicit none

  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16
  integer, parameter :: c64  = 8
  integer, parameter :: c128 = 16

contains

  ! ------------------------------------------------------------------
  ! l3dtavecevalmat_r64
  ! 3D Laplace solid spherical harmonics: values F and Cartesian gradient
  ! (Fx, Fy, Fz) at target points ztargs(3, nt).  Each output is shaped
  ! (nt, (nterms+1)**2), with one column per basis index (l, m) for
  !   0 <= l <= nterms,  -l <= m <= l.
  ! Column index packing (matches the legacy enumeration):
  !   for n = 0..nterms:
  !     tmpidx0 = n*n
  !     col(n, 0)  = tmpidx0 + n + 1
  !     col(n, +m) = tmpidx0 + n + 1 + m     (m = 1..n)
  !     col(n, -m) = tmpidx0 + n + 1 - m     (m = 1..n)
  ! Ported from subroutine l3dtavecevalmat in linequad.f (lines 3716-3908).
  ! ------------------------------------------------------------------
  subroutine l3dtavecevalmat_r64(ztargs, nt, nterms, F, Fx, Fy, Fz, ier)
    implicit none
    integer(8),  intent(in)    :: nt, nterms
    integer(8),  intent(inout) :: ier
    real(r64),   intent(in)    :: ztargs(3, nt)
    complex(c64), intent(inout) :: F (nt, (nterms+1)**2)
    complex(c64), intent(inout) :: Fx(nt, (nterms+1)**2)
    complex(c64), intent(inout) :: Fy(nt, (nterms+1)**2)
    complex(c64), intent(inout) :: Fz(nt, (nterms+1)**2)

    integer(8) :: i, m, n, k, tmpidx0, tmpidx1, tmpidx2
    real(r64)  :: rscale, center(3), zdiff(3), ztarg(3)
    real(r64)  :: fruse, done, r, theta, phi, d
    real(r64)  :: ctheta, stheta, cphi, sphi
    real(r64)  :: phix, phiy, phiz, rx, ry, rz
    real(r64)  :: thetax, thetay, thetaz
    real(r64)  :: sthetax, sthetay, sthetaz
    real(r64)  :: frdernpp, fruseppdn
    complex(c64) :: eye, ephi1, ephipm
    complex(c64) :: ztmp1, ur2, utheta2, uphi2
    complex(c64) :: ur2ephipm, utheta2ephipm, uphi2ephipm

    real(r64),   allocatable :: pp(:,:), ppd(:,:), fr(:), frder(:)
    real(r64),   allocatable :: rat1(:,:), rat2(:,:)
    complex(c64), allocatable :: ephi(:)

    rscale = 1.0_r64
    center = 0.0_r64
    eye    = (0.0_r64, 1.0_r64)
    done   = 1.0_r64
    ier    = 0_8

    allocate(pp  (0:nterms, 0:nterms), ppd (0:nterms, 0:nterms))
    allocate(rat1(0:nterms, 0:nterms), rat2(0:nterms, 0:nterms))
    allocate(fr(0:nterms+1), frder(0:nterms+1))
    allocate(ephi(-nterms-1 : nterms+1))

    call ylgndrini_r64(nterms, rat1, rat2)

    F  = (0.0_r64, 0.0_r64)
    Fx = (0.0_r64, 0.0_r64)
    Fy = (0.0_r64, 0.0_r64)
    Fz = (0.0_r64, 0.0_r64)

    do k = 1, nt
      ztarg    = ztargs(:, k)
      zdiff(1) = ztarg(1) - center(1)
      zdiff(2) = ztarg(2) - center(2)
      zdiff(3) = ztarg(3) - center(3)

      call cart2polarl_r64(zdiff, r, theta, phi)
      d      = rscale * r
      ctheta = cos(theta)
      stheta = sqrt(done - ctheta*ctheta)
      cphi   = cos(phi)
      sphi   = sin(phi)
      ephi1  = cmplx(cphi, sphi, kind=c64)

      ! e^{i m phi} array for m = -nterms-1 .. nterms+1
      ephi( 0) = (1.0_r64, 0.0_r64)
      ephi( 1) = ephi1
      ephi(-1) = conjg(ephi1)
      fr(0) = 1.0_r64
      fr(1) = d
      do i = 2, nterms + 1
        fr(i)    = fr(i-1) * d
        ephi( i) = ephi( i-1) * ephi1
        ephi(-i) = ephi(-i+1) * ephi(-1)
      end do
      frder(0) = 0.0_r64
      do i = 1, nterms + 1
        frder(i) = real(i, r64) * fr(i-1) * rscale
      end do

      ! Spherical-to-Cartesian gradient coefficients (1/sin(theta) and 1/r
      ! factors absorbed by the scaled Ynm and frder; see legacy comment).
      rx       =  stheta * cphi
      thetax   =  ctheta * cphi
      phix     = -sphi
      ry       =  stheta * sphi
      thetay   =  ctheta * sphi
      phiy     =  cphi
      rz       =  ctheta
      thetaz   = -stheta
      phiz     = 0.0_r64
      sthetax  = stheta * thetax
      sthetay  = stheta * thetay
      sthetaz  = stheta * thetaz

      ! Scaled associated Legendre functions and derivatives.
      call ylgndru2sf_r64(nterms, ctheta, pp, ppd, rat1, rat2)
      F(k, 1) = fr(0)

      do n = 1, nterms
        tmpidx0   = n * n
        fruse     = fr(n-1) * rscale
        frdernpp  = frder(n) * pp(n, 0)
        fruseppdn = fruse    * ppd(n, 0)
        tmpidx1   = tmpidx0 + n + 1
        F (k, tmpidx1) = fr(n) * pp(n, 0)
        Fx(k, tmpidx1) = frdernpp * rx - fruseppdn * sthetax
        Fy(k, tmpidx1) = frdernpp * ry - fruseppdn * sthetay
        Fz(k, tmpidx1) = frdernpp * rz - fruseppdn * sthetaz
        do m = 1, n
          ephipm        = ephi(m)
          ztmp1         = fr(n)    * pp(n, m) * stheta
          ur2           = frder(n) * pp(n, m) * stheta
          utheta2       = fruse    * ppd(n, m)
          ur2ephipm     = ur2     * ephipm
          utheta2ephipm = utheta2 * ephipm
          uphi2         = fruse * pp(n, m) * eye * real(m, r64)
          uphi2ephipm   = uphi2 * ephipm
          tmpidx1       = tmpidx0 + n + 1 + m
          tmpidx2       = tmpidx0 + n + 1 - m
          F (k, tmpidx1) = ztmp1 * ephipm
          F (k, tmpidx2) = conjg(F (k, tmpidx1))
          Fx(k, tmpidx1) =  ur2ephipm * rx - utheta2ephipm * thetax &
                          + uphi2ephipm * phix
          Fx(k, tmpidx2) = conjg(Fx(k, tmpidx1))
          Fy(k, tmpidx1) =  ur2ephipm * ry - utheta2ephipm * thetay &
                          + uphi2ephipm * phiy
          Fy(k, tmpidx2) = conjg(Fy(k, tmpidx1))
          Fz(k, tmpidx1) =  ur2ephipm * rz - utheta2ephipm * thetaz &
                          + uphi2ephipm * phiz
          Fz(k, tmpidx2) = conjg(Fz(k, tmpidx1))
        end do
      end do
    end do

    deallocate(pp, ppd, rat1, rat2, fr, frder, ephi)

  end subroutine l3dtavecevalmat_r64

  ! ------------------------------------------------------------------
  ! l3dtavecevalmat_r128
  ! r128 twin of l3dtavecevalmat_r64.  Kind-only substitution; no
  ! algorithmic differences.  Numeric literals all carry the _r128 suffix.
  ! ------------------------------------------------------------------
  subroutine l3dtavecevalmat_r128(ztargs, nt, nterms, F, Fx, Fy, Fz, ier)
    implicit none
    integer(8),    intent(in)    :: nt, nterms
    integer(8),    intent(inout) :: ier
    real(r128),    intent(in)    :: ztargs(3, nt)
    complex(c128), intent(inout) :: F (nt, (nterms+1)**2)
    complex(c128), intent(inout) :: Fx(nt, (nterms+1)**2)
    complex(c128), intent(inout) :: Fy(nt, (nterms+1)**2)
    complex(c128), intent(inout) :: Fz(nt, (nterms+1)**2)

    integer(8) :: i, m, n, k, tmpidx0, tmpidx1, tmpidx2
    real(r128) :: rscale, center(3), zdiff(3), ztarg(3)
    real(r128) :: fruse, done, r, theta, phi, d
    real(r128) :: ctheta, stheta, cphi, sphi
    real(r128) :: phix, phiy, phiz, rx, ry, rz
    real(r128) :: thetax, thetay, thetaz
    real(r128) :: sthetax, sthetay, sthetaz
    real(r128) :: frdernpp, fruseppdn
    complex(c128) :: eye, ephi1, ephipm
    complex(c128) :: ztmp1, ur2, utheta2, uphi2
    complex(c128) :: ur2ephipm, utheta2ephipm, uphi2ephipm

    real(r128),    allocatable :: pp(:,:), ppd(:,:), fr(:), frder(:)
    real(r128),    allocatable :: rat1(:,:), rat2(:,:)
    complex(c128), allocatable :: ephi(:)

    rscale = 1.0_r128
    center = 0.0_r128
    eye    = (0.0_r128, 1.0_r128)
    done   = 1.0_r128
    ier    = 0_8

    allocate(pp  (0:nterms, 0:nterms), ppd (0:nterms, 0:nterms))
    allocate(rat1(0:nterms, 0:nterms), rat2(0:nterms, 0:nterms))
    allocate(fr(0:nterms+1), frder(0:nterms+1))
    allocate(ephi(-nterms-1 : nterms+1))

    call ylgndrini_r128(nterms, rat1, rat2)

    F  = (0.0_r128, 0.0_r128)
    Fx = (0.0_r128, 0.0_r128)
    Fy = (0.0_r128, 0.0_r128)
    Fz = (0.0_r128, 0.0_r128)

    do k = 1, nt
      ztarg    = ztargs(:, k)
      zdiff(1) = ztarg(1) - center(1)
      zdiff(2) = ztarg(2) - center(2)
      zdiff(3) = ztarg(3) - center(3)

      call cart2polarl_r128(zdiff, r, theta, phi)
      d      = rscale * r
      ctheta = cos(theta)
      stheta = sqrt(done - ctheta*ctheta)
      cphi   = cos(phi)
      sphi   = sin(phi)
      ephi1  = cmplx(cphi, sphi, kind=c128)

      ephi( 0) = (1.0_r128, 0.0_r128)
      ephi( 1) = ephi1
      ephi(-1) = conjg(ephi1)
      fr(0) = 1.0_r128
      fr(1) = d
      do i = 2, nterms + 1
        fr(i)    = fr(i-1) * d
        ephi( i) = ephi( i-1) * ephi1
        ephi(-i) = ephi(-i+1) * ephi(-1)
      end do
      frder(0) = 0.0_r128
      do i = 1, nterms + 1
        frder(i) = real(i, r128) * fr(i-1) * rscale
      end do

      rx       =  stheta * cphi
      thetax   =  ctheta * cphi
      phix     = -sphi
      ry       =  stheta * sphi
      thetay   =  ctheta * sphi
      phiy     =  cphi
      rz       =  ctheta
      thetaz   = -stheta
      phiz     = 0.0_r128
      sthetax  = stheta * thetax
      sthetay  = stheta * thetay
      sthetaz  = stheta * thetaz

      call ylgndru2sf_r128(nterms, ctheta, pp, ppd, rat1, rat2)
      F(k, 1) = fr(0)

      do n = 1, nterms
        tmpidx0   = n * n
        fruse     = fr(n-1) * rscale
        frdernpp  = frder(n) * pp(n, 0)
        fruseppdn = fruse    * ppd(n, 0)
        tmpidx1   = tmpidx0 + n + 1
        F (k, tmpidx1) = fr(n) * pp(n, 0)
        Fx(k, tmpidx1) = frdernpp * rx - fruseppdn * sthetax
        Fy(k, tmpidx1) = frdernpp * ry - fruseppdn * sthetay
        Fz(k, tmpidx1) = frdernpp * rz - fruseppdn * sthetaz
        do m = 1, n
          ephipm        = ephi(m)
          ztmp1         = fr(n)    * pp(n, m) * stheta
          ur2           = frder(n) * pp(n, m) * stheta
          utheta2       = fruse    * ppd(n, m)
          ur2ephipm     = ur2     * ephipm
          utheta2ephipm = utheta2 * ephipm
          uphi2         = fruse * pp(n, m) * eye * real(m, r128)
          uphi2ephipm   = uphi2 * ephipm
          tmpidx1       = tmpidx0 + n + 1 + m
          tmpidx2       = tmpidx0 + n + 1 - m
          F (k, tmpidx1) = ztmp1 * ephipm
          F (k, tmpidx2) = conjg(F (k, tmpidx1))
          Fx(k, tmpidx1) =  ur2ephipm * rx - utheta2ephipm * thetax &
                          + uphi2ephipm * phix
          Fx(k, tmpidx2) = conjg(Fx(k, tmpidx1))
          Fy(k, tmpidx1) =  ur2ephipm * ry - utheta2ephipm * thetay &
                          + uphi2ephipm * phiy
          Fy(k, tmpidx2) = conjg(Fy(k, tmpidx1))
          Fz(k, tmpidx1) =  ur2ephipm * rz - utheta2ephipm * thetaz &
                          + uphi2ephipm * phiz
          Fz(k, tmpidx2) = conjg(Fz(k, tmpidx1))
        end do
      end do
    end do

    deallocate(pp, ppd, rat1, rat2, fr, frder, ephi)

  end subroutine l3dtavecevalmat_r128

  ! ------------------------------------------------------------------
  ! cart2polarl_r64
  ! Convert a Cartesian 3-vector to spherical coordinates (r, theta, phi).
  ! Matches the legacy cart2polarl in linequad.f.
  ! ------------------------------------------------------------------
  subroutine cart2polarl_r64(zat, r, theta, phi)
    implicit none
    real(r64), intent(in)  :: zat(3)
    real(r64), intent(out) :: r, theta, phi

    real(r64) :: proj
    r    = sqrt(zat(1)**2 + zat(2)**2 + zat(3)**2)
    proj = sqrt(zat(1)**2 + zat(2)**2)
    theta = atan2(proj, zat(3))
    if (abs(zat(1)) == 0.0_r64 .and. abs(zat(2)) == 0.0_r64) then
      phi = 0.0_r64
    else
      phi = atan2(zat(2), zat(1))
    end if
  end subroutine cart2polarl_r64

  subroutine cart2polarl_r128(zat, r, theta, phi)
    implicit none
    real(r128), intent(in)  :: zat(3)
    real(r128), intent(out) :: r, theta, phi

    real(r128) :: proj
    r    = sqrt(zat(1)**2 + zat(2)**2 + zat(3)**2)
    proj = sqrt(zat(1)**2 + zat(2)**2)
    theta = atan2(proj, zat(3))
    if (abs(zat(1)) == 0.0_r128 .and. abs(zat(2)) == 0.0_r128) then
      phi = 0.0_r128
    else
      phi = atan2(zat(2), zat(1))
    end if
  end subroutine cart2polarl_r128

  ! ------------------------------------------------------------------
  ! ylgndrini_r64
  ! Precompute the recurrence coefficients rat1(n,m), rat2(n,m) used by
  ! ylgndru2sf for the fast associated-Legendre evaluation.
  ! ------------------------------------------------------------------
  subroutine ylgndrini_r64(nmax, rat1, rat2)
    implicit none
    integer(8), intent(in)  :: nmax
    real(r64),  intent(out) :: rat1(0:nmax, 0:nmax), rat2(0:nmax, 0:nmax)

    integer(8) :: m, n
    rat1(0, 0) = 1.0_r64
    rat2(0, 0) = 1.0_r64
    do m = 0, nmax
      if (m > 0)        rat1(m, m)   = sqrt((2.0_r64*m - 1.0_r64) / (2.0_r64*m))
      if (m > 0)        rat2(m, m)   = 1.0_r64
      if (m < nmax)     rat1(m+1, m) = sqrt(2.0_r64*m + 1.0_r64)
      if (m < nmax)     rat2(m+1, m) = 1.0_r64
      do n = m+2, nmax
        rat1(n, m) = real(2*n - 1, r64)
        rat2(n, m) = sqrt(real(n+m-1, r64) * real(n-m-1, r64))
        rat1(n, m) = rat1(n, m) / sqrt(real(n-m, r64) * real(n+m, r64))
        rat2(n, m) = rat2(n, m) / sqrt(real(n-m, r64) * real(n+m, r64))
      end do
    end do
  end subroutine ylgndrini_r64

  subroutine ylgndrini_r128(nmax, rat1, rat2)
    implicit none
    integer(8), intent(in)  :: nmax
    real(r128), intent(out) :: rat1(0:nmax, 0:nmax), rat2(0:nmax, 0:nmax)

    integer(8) :: m, n
    rat1(0, 0) = 1.0_r128
    rat2(0, 0) = 1.0_r128
    do m = 0, nmax
      if (m > 0)        rat1(m, m)   = sqrt((2.0_r128*m - 1.0_r128) / (2.0_r128*m))
      if (m > 0)        rat2(m, m)   = 1.0_r128
      if (m < nmax)     rat1(m+1, m) = sqrt(2.0_r128*m + 1.0_r128)
      if (m < nmax)     rat2(m+1, m) = 1.0_r128
      do n = m+2, nmax
        rat1(n, m) = real(2*n - 1, r128)
        rat2(n, m) = sqrt(real(n+m-1, r128) * real(n-m-1, r128))
        rat1(n, m) = rat1(n, m) / sqrt(real(n-m, r128) * real(n+m, r128))
        rat2(n, m) = rat2(n, m) / sqrt(real(n-m, r128) * real(n+m, r128))
      end do
    end do
  end subroutine ylgndrini_r128

  ! ------------------------------------------------------------------
  ! ylgndru2sf_r64
  ! Evaluate scaled normalized associated Legendre functions y(n,m) and
  ! their derivatives d(n,m) for 0 <= n <= nmax, 0 <= m <= n.
  !
  !   Ynm(x) = sqrt((n-m)!/(n+m)!) Pnm(x).
  !
  ! For m>0, the returned y is scaled by 1/sqrt(1-x**2); the returned d
  ! is scaled by sqrt(1-x**2).  The legacy comment notes this scaling is
  ! compensated by the (rx, thetax, phix)-side construction in
  ! l3dtavecevalmat above.
  ! ------------------------------------------------------------------
  subroutine ylgndru2sf_r64(nmax, x, y, d, rat1, rat2)
    implicit none
    integer(8), intent(in)  :: nmax
    real(r64),  intent(in)  :: x
    real(r64),  intent(out) :: y(0:nmax, 0:nmax), d(0:nmax, 0:nmax)
    real(r64),  intent(in)  :: rat1(0:nmax, 0:nmax), rat2(0:nmax, 0:nmax)

    integer(8) :: n, m
    real(r64)  :: u, u2

    u2 = (1.0_r64 - x) * (1.0_r64 + x)
    u  = -sqrt(u2)
    y(0, 0) = 1.0_r64
    d(0, 0) = 0.0_r64

    ! Standard Legendre polynomials, m = 0.
    m = 0
    if (m < nmax) y(m+1, m) = x * y(m, m) * rat1(m+1, m)
    if (m < nmax) d(m+1, m) = (x * d(m, m) + y(m, m)) * rat1(m+1, m)
    do n = m+2, nmax
      y(n, m) = rat1(n, m) * x * y(n-1, m) - rat2(n, m) * y(n-2, m)
      d(n, m) = rat1(n, m) * (x * d(n-1, m) + y(n-1, m)) - rat2(n, m) * d(n-2, m)
    end do

    ! Scaled associated Legendre functions, m >= 1.
    do m = 1, nmax
      if (m == 1) y(m, m) = y(m-1, m-1) * (-1.0_r64) * rat1(m, m)
      if (m >  1) y(m, m) = y(m-1, m-1) * u           * rat1(m, m)
      if (m >  0) d(m, m) = y(m, m) * real(-m, r64) * x

      if (m < nmax) y(m+1, m) = x * y(m, m) * rat1(m+1, m)
      if (m < nmax) d(m+1, m) = (x * d(m, m) + u2 * y(m, m)) * rat1(m+1, m)
      do n = m+2, nmax
        y(n, m) = rat1(n, m) * x * y(n-1, m) - rat2(n, m) * y(n-2, m)
        d(n, m) = rat1(n, m) * (x * d(n-1, m) + u2 * y(n-1, m)) - rat2(n, m) * d(n-2, m)
      end do
    end do
  end subroutine ylgndru2sf_r64

  subroutine ylgndru2sf_r128(nmax, x, y, d, rat1, rat2)
    implicit none
    integer(8), intent(in)  :: nmax
    real(r128), intent(in)  :: x
    real(r128), intent(out) :: y(0:nmax, 0:nmax), d(0:nmax, 0:nmax)
    real(r128), intent(in)  :: rat1(0:nmax, 0:nmax), rat2(0:nmax, 0:nmax)

    integer(8) :: n, m
    real(r128) :: u, u2

    u2 = (1.0_r128 - x) * (1.0_r128 + x)
    u  = -sqrt(u2)
    y(0, 0) = 1.0_r128
    d(0, 0) = 0.0_r128

    m = 0
    if (m < nmax) y(m+1, m) = x * y(m, m) * rat1(m+1, m)
    if (m < nmax) d(m+1, m) = (x * d(m, m) + y(m, m)) * rat1(m+1, m)
    do n = m+2, nmax
      y(n, m) = rat1(n, m) * x * y(n-1, m) - rat2(n, m) * y(n-2, m)
      d(n, m) = rat1(n, m) * (x * d(n-1, m) + y(n-1, m)) - rat2(n, m) * d(n-2, m)
    end do

    do m = 1, nmax
      if (m == 1) y(m, m) = y(m-1, m-1) * (-1.0_r128) * rat1(m, m)
      if (m >  1) y(m, m) = y(m-1, m-1) * u            * rat1(m, m)
      if (m >  0) d(m, m) = y(m, m) * real(-m, r128) * x

      if (m < nmax) y(m+1, m) = x * y(m, m) * rat1(m+1, m)
      if (m < nmax) d(m+1, m) = (x * d(m, m) + u2 * y(m, m)) * rat1(m+1, m)
      do n = m+2, nmax
        y(n, m) = rat1(n, m) * x * y(n-1, m) - rat2(n, m) * y(n-2, m)
        d(n, m) = rat1(n, m) * (x * d(n-1, m) + u2 * y(n-1, m)) - rat2(n, m) * d(n-2, m)
      end do
    end do
  end subroutine ylgndru2sf_r128

  ! ------------------------------------------------------------------
  ! evaltensorproductharmonicgrad_r64
  ! Tensor-product harmonic polynomial basis H_{i,j} and its Cartesian
  ! gradient (fx, fy, fz) evaluated at nt 3D points r(3, nt).  Outputs:
  !   fx, fy, fz, f       (nt, order**2)   ! grad and value per (i,j)
  !   ijidx               (2, order**2)    ! (i, j) pair for each column
  ! Total-degree ordering of (i, j) pairs (orderx = ordery = order).
  ! Ported from qotential/utils/f/harmonics.f (Hai, 01/01/23): same
  ! algorithm (Laplacian-power expansion of x^i y^j z^{i+j+1}), F77
  ! fixed-form -> F90 module style, intent(out) -> intent(inout) per
  ! mwrap convention, default integer -> integer(8) for -i8 boundary.
  ! ------------------------------------------------------------------
  subroutine evaltensorproductharmonicgrad_r64(nt, r, order, fx, fy, fz, f, ijidx)
    integer(8), intent(in)    :: nt, order
    real(r64),  intent(in)    :: r(3, nt)
    real(r64),  intent(inout) :: fx(nt, order*order), fy(nt, order*order), &
                                 fz(nt, order*order), f (nt, order*order)
    integer(8), intent(inout) :: ijidx(2, order*order)

    integer(8) :: k
    real(r64),  allocatable :: rxpow(:,:), rypow(:,:), rzpow(:,:)
    real(r64),  allocatable :: coeffs_all(:,:)
    integer(8), allocatable :: powx_all(:,:), powy_all(:,:), powz_all(:,:)
    integer(8), allocatable :: lens(:)
    real(r64),  allocatable :: coeffs_ext(:)
    integer(8), allocatable :: powx_ext(:), powy_ext(:), powz_ext(:)
    integer(8) :: lens_k
    real(r64),  allocatable :: fx_k(:), fy_k(:), fz_k(:), f_k(:)

    allocate(rxpow(nt, order+1), rypow(nt, order+1), rzpow(nt, 2*order+1))
    allocate(coeffs_all(order*order, order*order))
    allocate(powx_all(order*order, order*order), &
             powy_all(order*order, order*order), &
             powz_all(order*order, order*order))
    allocate(lens(order*order))
    allocate(coeffs_ext(order*order))
    allocate(powx_ext(order*order), powy_ext(order*order), powz_ext(order*order))
    allocate(fx_k(nt), fy_k(nt), fz_k(nt), f_k(nt))

    rxpow = 1.0_r64
    rypow = 1.0_r64
    rzpow = 1.0_r64
    do k = 1, order
      rxpow(:, k+1) = rxpow(:, k) * r(1, :)
      rypow(:, k+1) = rypow(:, k) * r(2, :)
    end do
    do k = 1, 2*order
      rzpow(:, k+1) = rzpow(:, k) * r(3, :)
    end do

    call hijcoeffsall_r64(order, coeffs_all, powx_all, powy_all, powz_all, lens, ijidx)

    fx = 0.0_r64
    fy = 0.0_r64
    fz = 0.0_r64
    f  = 0.0_r64
    do k = 1, order*order
      coeffs_ext = coeffs_all(k, :)
      powx_ext   = powx_all  (k, :)
      powy_ext   = powy_all  (k, :)
      powz_ext   = powz_all  (k, :)
      lens_k     = lens(k)

      call fxfyfzf_r64(order*order, coeffs_ext, powx_ext, powy_ext, powz_ext, &
                       lens_k, nt, order, rxpow, rypow, rzpow, &
                       fx_k, fy_k, fz_k, f_k)

      fx(:, k) = fx_k
      fy(:, k) = fy_k
      fz(:, k) = fz_k
      f (:, k) = f_k
    end do

    deallocate(rxpow, rypow, rzpow)
    deallocate(coeffs_all, powx_all, powy_all, powz_all, lens)
    deallocate(coeffs_ext, powx_ext, powy_ext, powz_ext)
    deallocate(fx_k, fy_k, fz_k, f_k)

  end subroutine evaltensorproductharmonicgrad_r64

  ! ------------------------------------------------------------------
  ! fxfyfzf_r64
  ! Inner kernel: given coefficient/power tables for a single basis
  ! index, evaluate the polynomial and its gradient at nt points using
  ! the precomputed power tables rxpow, rypow, rzpow.
  ! ------------------------------------------------------------------
  subroutine fxfyfzf_r64(order2, coeffs, powx, powy, powz, lens_k, &
                         nt, order, rxpow, rypow, rzpow, &
                         fx_k, fy_k, fz_k, f_k)
    integer(8), intent(in)    :: order2, lens_k, nt, order
    real(r64),  intent(in)    :: coeffs(order2)
    integer(8), intent(in)    :: powx(order2), powy(order2), powz(order2)
    real(r64),  intent(in)    :: rxpow(nt, order+1), rypow(nt, order+1)
    real(r64),  intent(in)    :: rzpow(nt, 2*order+1)
    real(r64),  intent(inout) :: fx_k(nt), fy_k(nt), fz_k(nt), f_k(nt)

    integer(8), allocatable :: pow_k(:,:)
    integer(8), allocatable :: powgradx(:,:), powgrady(:,:), powgradz(:,:)
    real(r64),  allocatable :: coeffs_k(:)
    real(r64),  allocatable :: coeffsgradx(:), coeffsgrady(:), coeffsgradz(:)
    integer(8) :: j

    allocate(pow_k(3, lens_k))
    allocate(powgradx(3, lens_k), powgrady(3, lens_k), powgradz(3, lens_k))
    allocate(coeffs_k(lens_k))
    allocate(coeffsgradx(lens_k), coeffsgrady(lens_k), coeffsgradz(lens_k))

    pow_k(1, :) = powx(1:lens_k)
    pow_k(2, :) = powy(1:lens_k)
    pow_k(3, :) = powz(1:lens_k)
    coeffs_k    = coeffs(1:lens_k)

    powgradx = pow_k
    powgrady = pow_k
    powgradz = pow_k
    powgradx(1, :) = powgradx(1, :) - 1
    powgrady(2, :) = powgrady(2, :) - 1
    powgradz(3, :) = powgradz(3, :) - 1
    do j = 1, lens_k
      if (powgradx(1, j) == -1_8) powgradx(1, j) = 0_8
      if (powgrady(2, j) == -1_8) powgrady(2, j) = 0_8
      if (powgradz(3, j) == -1_8) powgradz(3, j) = 0_8
    end do

    coeffsgradx = coeffs_k * real(pow_k(1, :), r64)
    coeffsgrady = coeffs_k * real(pow_k(2, :), r64)
    coeffsgradz = coeffs_k * real(pow_k(3, :), r64)

    fx_k = 0.0_r64
    fy_k = 0.0_r64
    fz_k = 0.0_r64
    f_k  = 0.0_r64
    do j = 1, lens_k
      fx_k = fx_k + coeffsgradx(j) * rxpow(:, powgradx(1, j) + 1) &
                                   * rypow(:, powgradx(2, j) + 1) &
                                   * rzpow(:, powgradx(3, j) + 1)
      fy_k = fy_k + coeffsgrady(j) * rxpow(:, powgrady(1, j) + 1) &
                                   * rypow(:, powgrady(2, j) + 1) &
                                   * rzpow(:, powgrady(3, j) + 1)
      fz_k = fz_k + coeffsgradz(j) * rxpow(:, powgradz(1, j) + 1) &
                                   * rypow(:, powgradz(2, j) + 1) &
                                   * rzpow(:, powgradz(3, j) + 1)
      f_k  = f_k  + coeffs_k(j)    * rxpow(:, pow_k(1, j) + 1) &
                                   * rypow(:, pow_k(2, j) + 1) &
                                   * rzpow(:, pow_k(3, j) + 1)
    end do

    deallocate(pow_k, powgradx, powgrady, powgradz)
    deallocate(coeffs_k, coeffsgradx, coeffsgrady, coeffsgradz)

  end subroutine fxfyfzf_r64

  ! ------------------------------------------------------------------
  ! hijcoeffsall_r64
  ! Build the coefficient/power tables for all order**2 (i, j) basis
  ! pairs in total-degree order.
  ! ------------------------------------------------------------------
  subroutine hijcoeffsall_r64(order, coeffs, powx, powy, powz, lens, ijidx)
    integer(8), intent(in)    :: order
    real(r64),  intent(inout) :: coeffs(order*order, order*order)
    integer(8), intent(inout) :: powx(order*order, order*order), &
                                 powy(order*order, order*order), &
                                 powz(order*order, order*order)
    integer(8), intent(inout) :: lens(order*order)
    integer(8), intent(inout) :: ijidx(2, order*order)

    integer(8) :: orderx, ordery, k, i, j, lens_k
    real(r64),  allocatable :: coeffs_k(:)
    integer(8), allocatable :: pow_k(:,:)

    orderx = order
    ordery = order

    allocate(coeffs_k(order*order), pow_k(3, order*order))

    call tdordering_int(orderx, ordery, ijidx)
    coeffs = 0.0_r64
    powx   = 0_8
    powy   = 0_8
    powz   = 0_8
    lens   = 0_8
    do k = 1, order*order
      i = ijidx(1, k)
      j = ijidx(2, k)
      call hijcoeffs0_r64(i, j, order*order, coeffs_k, pow_k, lens_k)
      lens(k)      = lens_k
      coeffs(k, :) = coeffs_k
      powx(k, :)   = pow_k(1, :)
      powy(k, :)   = pow_k(2, :)
      powz(k, :)   = pow_k(3, :)
    end do

    deallocate(coeffs_k, pow_k)

  end subroutine hijcoeffsall_r64

  ! ------------------------------------------------------------------
  ! hijcoeffs0_r64
  ! Pad hijcoeffs_r64 output (length (i+1)*(j+1)) to fixed length
  ! orderxy.  Returns the actual length in `len`.
  ! ------------------------------------------------------------------
  subroutine hijcoeffs0_r64(i, j, orderxy, coeffs_ext, pow_ext, len)
    integer(8), intent(in)    :: i, j, orderxy
    real(r64),  intent(inout) :: coeffs_ext(orderxy)
    integer(8), intent(inout) :: pow_ext(3, orderxy)
    integer(8), intent(out)   :: len

    real(r64),  allocatable :: coeffs(:)
    integer(8), allocatable :: pow(:,:)
    integer(8) :: l

    allocate(coeffs((i+1)*(j+1)), pow(3, (i+1)*(j+1)))

    call hijcoeffs_r64(i, j, coeffs, pow)
    len = size(coeffs, kind=8)
    coeffs_ext = 0.0_r64
    pow_ext    = 0_8
    coeffs_ext( [ (l, l = 1, len) ] ) = coeffs
    pow_ext(1, [ (l, l = 1, len) ] ) = pow(1, :)
    pow_ext(2, [ (l, l = 1, len) ] ) = pow(2, :)
    pow_ext(3, [ (l, l = 1, len) ] ) = pow(3, :)

    deallocate(coeffs, pow)

  end subroutine hijcoeffs0_r64

  ! ------------------------------------------------------------------
  ! tdordering_int
  ! Total-degree ordering: produce (i, j) pairs in increasing i+j.
  ! ------------------------------------------------------------------
  subroutine tdordering_int(orderx, ordery, ijidx)
    integer(8), intent(in)    :: orderx, ordery
    integer(8), intent(inout) :: ijidx(2, orderx*ordery)

    integer(8), allocatable :: iidx(:,:), jidx(:,:)
    integer(8), allocatable :: iidx_rs(:), jidx_rs(:)
    logical,    allocatable :: tmpidx(:)
    integer(8) :: l, j, kstart

    allocate(iidx(orderx, ordery), jidx(orderx, ordery))
    allocate(iidx_rs(orderx*ordery), jidx_rs(orderx*ordery))
    allocate(tmpidx(orderx*ordery))

    call meshgrid_int(ordery, [ (l, l = 0_8, ordery-1_8) ], &
                      orderx, [ (l, l = 0_8, orderx-1_8) ], jidx, iidx)
    ijidx   = 0_8
    iidx_rs = reshape(iidx, [orderx*ordery])
    jidx_rs = reshape(jidx, [orderx*ordery])

    kstart = 1_8
    do l = 0_8, ordery + orderx - 2_8
      tmpidx = (jidx_rs + iidx_rs) == l
      do j = 1, orderx*ordery
        if (tmpidx(j)) then
          ijidx(1, kstart) = iidx_rs(j)
          ijidx(2, kstart) = jidx_rs(j)
          kstart = kstart + 1_8
        end if
      end do
    end do

    deallocate(iidx, jidx, iidx_rs, jidx_rs, tmpidx)

  end subroutine tdordering_int

  ! ------------------------------------------------------------------
  ! hijcoeffs_r64
  ! Coefficients and powers (x, y, z exponents) for the (i, j) harmonic
  ! polynomial built by Laplacian-power expansion of x^i y^j z^{i+j+1}.
  ! Output sizes: coeffs((i+1)*(j+1)), pow(3, (i+1)*(j+1)).
  ! ------------------------------------------------------------------
  subroutine hijcoeffs_r64(i, j, coeffs, pow)
    integer(8), intent(in)    :: i, j
    real(r64),  intent(inout) :: coeffs((i+1)*(j+1))
    integer(8), intent(inout) :: pow(3, (i+1)*(j+1))

    integer(8), allocatable :: xpow(:,:), ypow(:,:)
    integer(8) :: l, k, length_loc
    real(r64),  allocatable :: coeffs_updt(:)

    allocate(xpow(i+1, j+1), ypow(i+1, j+1))
    allocate(coeffs_updt((i+1)*(j+1)))

    call meshgrid_int(j+1, [ (l, l = 0_8, j) ], &
                      i+1, [ (l, l = 0_8, i) ], ypow, xpow)
    pow(1, :) = reshape(xpow, [(i+1)*(j+1)])
    pow(2, :) = reshape(ypow, [(i+1)*(j+1)])
    pow(3, :) = i + j + 1 - (pow(1, :) + pow(2, :))

    coeffs = 0.0_r64
    do k = 0, floor(real(i+j, r64) / 2.0_r64, kind=8)
      call lenofnckcoeffs0_int(i, j, k, length_loc)
      call hijcoeffsk0_r64(i, j, k, pow, length_loc, coeffs_updt)
      coeffs = coeffs + coeffs_updt
    end do

    deallocate(xpow, ypow, coeffs_updt)

  end subroutine hijcoeffs_r64

  ! ------------------------------------------------------------------
  ! meshgrid_int
  ! Integer meshgrid: xx(i, j) = x(j), yy(i, j) = y(i).
  ! ------------------------------------------------------------------
  subroutine meshgrid_int(n, x, m, y, xx, yy)
    integer(8), intent(in)    :: n, m
    integer(8), intent(in)    :: x(n), y(m)
    integer(8), intent(inout) :: xx(m, n), yy(m, n)

    xx = spread(x, 1, size(y))
    yy = spread(y, 2, size(x))

  end subroutine meshgrid_int

  ! ------------------------------------------------------------------
  ! hijcoeffsk0_r64
  ! Update part of hijcoeffs from the k-th Laplacian power.  Accumulates
  ! coefficients into coeffs_updt indexed by the (xpow, ypow) match
  ! with the polynomial's (pow(1,:), pow(2,:)).
  ! ------------------------------------------------------------------
  subroutine hijcoeffsk0_r64(i, j, k, pow, length_loc, coeffs_updt)
    integer(8), intent(in)    :: i, j, k, length_loc
    integer(8), intent(in)    :: pow(3, (i+1)*(j+1))
    real(r64),  intent(inout) :: coeffs_updt((i+1)*(j+1))

    real(r64),  allocatable :: nck_coeffs(:)
    integer(8), allocatable :: nck_pow(:,:), nck_der(:,:), nck_idx(:)
    integer(8) :: nck_der_init(2, k+1), nck_idx_init(k+1)
    logical    :: idx_flag0(k+1)
    integer(8) :: l, l_idx

    allocate(nck_coeffs(length_loc))
    allocate(nck_pow(3, length_loc), nck_der(2, length_loc))
    allocate(nck_idx(length_loc))

    idx_flag0 = (2_8 * [ (l, l = k, 0_8, -1_8) ] <= i) .and. &
                (2_8 * [ (l, l = 0_8, k,  1_8) ] <= j)
    nck_der_init(1, :) = 2_8 * [ (l, l = k, 0_8, -1_8) ]
    nck_der_init(2, :) = 2_8 * [ (l, l = 0_8, k,  1_8) ]
    nck_idx_init       =        [ (l, l = 1_8, k+1_8, 1_8) ]
    nck_der(1, :) = pack(nck_der_init(1, :), idx_flag0)
    nck_der(2, :) = pack(nck_der_init(2, :), idx_flag0)
    nck_idx       = pack(nck_idx_init,       idx_flag0)

    call polylapder2d_r64(i, j, k, length_loc, nck_der, nck_idx, nck_coeffs, nck_pow)
    coeffs_updt = 0.0_r64
    do l = 1, length_loc
      do l_idx = 1, (i+1)*(j+1)
        if (pow(1, l_idx) == nck_pow(1, l) .and. &
            pow(2, l_idx) == nck_pow(2, l)) then
          coeffs_updt(l_idx) = coeffs_updt(l_idx) + nck_coeffs(l)
        end if
      end do
    end do

    deallocate(nck_coeffs, nck_pow, nck_der, nck_idx)

  end subroutine hijcoeffsk0_r64

  ! ------------------------------------------------------------------
  ! lenofnckcoeffs0_int
  ! Count the number of (n, k) coefficients that survive the
  ! 2*l <= i / 2*l <= j filter for ell = ell.
  ! ------------------------------------------------------------------
  subroutine lenofnckcoeffs0_int(i, j, ell, length_loc)
    integer(8), intent(in)  :: i, j, ell
    integer(8), intent(out) :: length_loc

    logical    :: idx_flag0(ell+1)
    integer(8) :: k

    idx_flag0 = (2_8 * [ (k, k = ell, 0_8, -1_8) ] <= i) .and. &
                (2_8 * [ (k, k = 0_8, ell,  1_8) ] <= j)
    length_loc = count(idx_flag0, kind=8)

  end subroutine lenofnckcoeffs0_int

  ! ------------------------------------------------------------------
  ! polylapder2d_r64
  ! Coefficients and powers for the 2D Laplacian derivative applied
  ! to the polynomial x^i y^j z^{i+j+1}.  Author flagged gamma()
  ! factorial usage as inefficient; preserved verbatim per port.
  ! ------------------------------------------------------------------
  subroutine polylapder2d_r64(i, j, ell, length_loc, nck_der, nck_idx, &
                              nck_coeffs, nck_pow)
    integer(8), intent(in)    :: i, j, ell, length_loc
    integer(8), intent(in)    :: nck_der(2, length_loc), nck_idx(length_loc)
    real(r64),  intent(inout) :: nck_coeffs(length_loc)
    integer(8), intent(inout) :: nck_pow(3, length_loc)

    integer(8) :: k
    real(r64)  :: rec_coeffs

    do k = 1, length_loc
      nck_coeffs(k) = gamma(real(ell, r64) + 1.0_r64) &
                    / gamma(real(nck_idx(k), r64))    &
                    / gamma(real(ell, r64) - real(nck_idx(k), r64) + 2.0_r64)
    end do

    rec_coeffs = real((-1_8)**ell, r64) * 1.0_r64 / gamma(2.0_r64*real(ell, r64) + 2.0_r64)
    do k = 1, length_loc
      nck_coeffs(k) = gamma(real(i, r64) + 1.0_r64) &
                    / gamma(real(i, r64) - real(nck_der(1, k), r64) + 1.0_r64) &
                    * gamma(real(j, r64) + 1.0_r64) &
                    / gamma(real(j, r64) - real(nck_der(2, k), r64) + 1.0_r64) &
                    * nck_coeffs(k)
    end do
    nck_coeffs = rec_coeffs * nck_coeffs

    nck_pow(1, :) = i - nck_der(1, :)
    nck_pow(2, :) = j - nck_der(2, :)
    nck_pow(3, :) = 2_8*ell + 1_8

  end subroutine polylapder2d_r64

end module harmonic_mod
