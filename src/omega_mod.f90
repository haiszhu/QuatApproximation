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
  use harmonic_mod, only: cart2polarl_r64
  implicit none
  private
  public :: qao_omeganm_i_r64, qao_omegaall_r64
  public :: qao_omeganm_i_r128, qao_omegaall_r128
  public :: qao_omegasdlp_r64

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

  subroutine qao_omegasdlp_r64(m, nterms, ncoeff, h_dim, r0, Ichi, Ialpha, &
                               omega_slp, omega)
    integer(8),   intent(in)  :: m, nterms, ncoeff, h_dim
    real(r64),    intent(in)  :: r0(3,m)
    complex(r64), intent(in)  :: Ichi(m,ncoeff,4), Ialpha(m,ncoeff,4)
    real(r64),    intent(out) :: omega_slp(h_dim,m), omega(h_dim,m,4)

    complex(r64), parameter :: IMA   = (0.0_r64, 1.0_r64)
    complex(r64), parameter :: CZERO = (0.0_r64, 0.0_r64)

    real(r64)  :: carray(0:4*nterms,0:4*nterms), dc(0:4*nterms,0:4*nterms)
    real(r64)  :: factinv(0:nterms)
    integer(8) :: fact(0:nterms)

    real(r64)    :: c_new(3), rvec(3), d, theta, phi
    complex(r64) :: ephi(-nterms-1:nterms+1), ephi1
    real(r64)    :: rotmatf(-nterms:nterms,-nterms:nterms,0:nterms)
    real(r64)    :: rotmatb(-nterms:nterms,-nterms:nterms,0:nterms)
    real(r64)    :: fr(0:nterms+1)
    complex(r64) :: Icold(ncoeff), Iaold(4,ncoeff)
    complex(r64) :: Icrb(ncoeff), Icsz(ncoeff)
    complex(r64) :: Iarb(4,ncoeff), Iasz(4,ncoeff)
    complex(r64) :: Icnew(h_dim), Ianew(4,h_dim)
    complex(r64) :: erb, erf, erbn, erbp, erfn, erfp
    complex(r64) :: c0j, cc0j, a0j(4), ca0j(4)
    complex(r64) :: szj, cszj, aszj(4), caszj(4)
    real(r64)    :: frdc, frdc1, frdc2, tv, tv2
    integer(8)   :: j, k, l, mm, mp, ij, ll, nmax
    integer(8)   :: col0, col0r, col02, idx, idx2, col, col2, colr

    nmax = 2*nterms
    do l = 0, 2*nmax
      carray(l,0) = 1.0_r64;  dc(l,0) = 1.0_r64
    end do
    do mm = 1, 2*nmax
      carray(mm,mm) = 1.0_r64;  dc(mm,mm) = 1.0_r64
      do l = mm+1, 2*nmax
        carray(l,mm) = carray(l-1,mm) + carray(l-1,mm-1)
        dc(l,mm)     = sqrt(carray(l,mm))
      end do
    end do
    fact(0) = 1_8;  fact(1) = 1_8
    factinv(0) = 1.0_r64;  factinv(1) = 1.0_r64
    do l = 2, nterms
      fact(l)    = fact(l-1)*l
      factinv(l) = 1.0_r64/real(fact(l), r64)
    end do

    do j = 1, m

      c_new(1) = r0(2,j)
      c_new(2) = r0(3,j)
      c_new(3) = r0(1,j)
      rvec = c_new
      call cart2polarl_r64(rvec, d, theta, phi)

      ephi(0)  = (1.0_r64, 0.0_r64)
      ephi1    = exp(IMA*phi)
      ephi(1)  = ephi1
      ephi(-1) = conjg(ephi1)
      do l = 1, nterms
        ephi(l+1)  = ephi(l)*ephi(1)
        ephi(-1-l) = conjg(ephi(l+1))
      end do

      block
        real(r64), parameter :: PRECIS = 1.0e-20_r64
        real(r64)  :: rd1(-nterms:nterms,-nterms:nterms)
        real(r64)  :: rd2(-nterms:nterms,-nterms:nterms)
        real(r64)  :: sqc(0:2*nterms,2)
        real(r64)  :: ctsqc(0:2*nterms), stsqc(0:2*nterms)
        real(r64)  :: cpsqc(0:2*nterms), cnsqc(0:2*nterms)
        integer(8) :: sgn1(0:nterms), sgn2(-nterms:nterms)
        integer(8) :: sgn12(-nterms:nterms,-nterms:nterms)
        real(r64)  :: ww, ct, st, hst, ctp, ctn, dd, sc, ijinv
        integer(8) :: im, imp, q, r

        rd1 = 0.0_r64;  rd2 = 0.0_r64
        ww = sqrt(0.5_r64)
        do q = 0, 2*nterms
          sqc(q,1) = sqrt(real(q, r64))
        end do
        sqc(0,2) = 0.0_r64
        if (nterms > 0_8) sqc(1,2) = 0.0_r64
        do q = 2, 2*nterms
          sqc(q,2) = sqrt(real(q, r64)*real(q-1, r64)*0.5_r64)
        end do

        ct = cos(theta);   if (abs(ct) <= PRECIS) ct = 0.0_r64
        st = sin(-theta);  if (abs(st) <= PRECIS) st = 0.0_r64
        hst = ww*st
        ctp =  2.0_r64*ww*cos(theta*0.5_r64)**2
        ctn = -2.0_r64*ww*sin(theta*0.5_r64)**2
        ctsqc = ct*sqc(:,1);   stsqc = st*sqc(:,1)
        cpsqc = ctp*sqc(:,2);  cnsqc = ctn*sqc(:,2)

        rd1(0,0) = 1.0_r64
        rotmatf(0,0,0) = 1.0_r64
        rotmatb(0,0,0) = 1.0_r64

        sgn1(0) = 1_8;  sgn2(0) = 1_8
        do q = 1, nterms
          sgn1(q)  = -sgn1(q-1)
          sgn2(q)  =  sgn1(q)
          sgn2(-q) =  sgn1(q)
        end do
        do q = -nterms, nterms
          sgn12(:,q) = sgn2(q)*sgn2
        end do

        do r = 1, nterms
          ijinv = 1.0_r64/real(r, r64)

          do im = -r, -1
            rd2(0,im) = -sqc(r-im,2)*rd1(0,im+1)
            if (im > 1-r) rd2(0,im) = rd2(0,im) + sqc(r+im,2)*rd1(0,im-1)
            rd2(0,im) = rd2(0,im)*hst
            if (im > -r) rd2(0,im) = rd2(0,im) &
                                   + rd1(0,im)*ctsqc(r+im)*sqc(r-im,1)
            rd2(0,im) = rd2(0,im)*ijinv
          end do
          rd2(0,0) = rd1(0,0)*ct
          if (r > 1_8) rd2(0,0) = rd2(0,0) &
                                + hst*sqc(r,2)*(2.0_r64*rd1(0,-1))*ijinv
          do im = 1, r
            rd2(0,im) = rd2(0,-im)
            if (mod(im,2_8) == 0_8) then
              rd2(im,0) =  rd2(0,im)
            else
              rd2(im,0) = -rd2(0,im)
            end if
          end do

          do imp = 1, r
            sc = ww/sqc(r+imp,2)
            do im = imp, r
              rd2(imp, im) = rd1(imp-1, im-1)*cpsqc(r+im)
              rd2(imp,-im) = rd1(imp-1,-im+1)*cnsqc(r+im)
              if (im < r-1) then
                rd2(imp, im) = rd2(imp, im) - rd1(imp-1, im+1)*cnsqc(r-im)
                rd2(imp,-im) = rd2(imp,-im) - rd1(imp-1,-im-1)*cpsqc(r-im)
              end if
              if (im < r) then
                dd = stsqc(r+im)*sqc(r-im,1)
                rd2(imp, im) = rd2(imp, im) + rd1(imp-1, im)*dd
                rd2(imp,-im) = rd2(imp,-im) + rd1(imp-1,-im)*dd
              end if
              rd2(imp, im) = rd2(imp, im)*sc
              rd2(imp,-im) = rd2(imp,-im)*sc
              if (im > imp) then
                if (mod(imp+im,2_8) == 0_8) then
                  rd2(im, imp) =  rd2(imp, im)
                  rd2(im,-imp) =  rd2(imp,-im)
                else
                  rd2(im, imp) = -rd2(imp, im)
                  rd2(im,-imp) = -rd2(imp,-im)
                end if
              end if
            end do
          end do

          rd1 = rd2
          do im = -r, r
            do imp = -r, -1
              rd2(imp,im) = rd2(-imp,-im)
            end do
          end do
          rotmatf(:,:,r) = rd2
          rotmatb(:,:,r) = real(sgn12, r64)*rd2
        end do
      end block

      fr(0) = 1.0_r64;  fr(1) = d
      do l = 2, nterms+1
        fr(l) = fr(l-1)*d
      end do

      Icold = Ichi(j,:,1)
      do k = 1, 4
        Iaold(k,:) = Ialpha(j,:,k)
      end do
      Icrb = CZERO;  Icsz = CZERO;  Icnew = CZERO
      Iarb = CZERO;  Iasz = CZERO;  Ianew = CZERO

      col0 = 0;  col0r = 0
      idx = 1;  col = 1;  mm = 0;  mp = 0
      erb = ephi(-mm)*rotmatb(mp,mm,0)
      Icrb(col) = Icrb(col) + erb*Icold(idx)
      Icsz(idx) = Icsz(idx) + Icrb(idx)
      do l = 1, nterms
        ll   = l
        frdc = fr(l)*dc(ll+mm,l)*dc(ll-mm,l)
        col  = (ll+1)*ll/2 + ll + mm + 1
        Icsz(col) = Icsz(col) + frdc*Icrb(idx)
      end do

      ij = 1;  col0 = 1;  col0r = 0
      idx2 = col0
      mm = -1
      idx2 = idx2 + 1;  col2 = col0 + 1
      c0j = Icold(idx2);    cc0j = conjg(c0j)
      a0j = Iaold(:,idx2);  ca0j = conjg(a0j)
      do mp = -1, 0
        erbn = ephi(-mm)*rotmatb(mp, mm,ij)
        erbp = ephi( mm)*rotmatb(mp,-mm,ij)
        Icrb(col2)   = Icrb(col2)   + erbn*c0j + erbp*cc0j
        Iarb(:,col2) = Iarb(:,col2) + erbn*a0j + erbp*ca0j
        col2 = col2 + 1
      end do
      mm = 0
      idx2 = idx2 + 1;  col2 = col0 + 1
      c0j = Icold(idx2);  a0j = Iaold(:,idx2)
      do mp = -ij, 0
        erb = ephi(-mm)*rotmatb(mp,mm,ij)
        Icrb(col2)   = Icrb(col2)   + erb*c0j
        Iarb(:,col2) = Iarb(:,col2) + erb*a0j
        col2 = col2 + 1
      end do

      idx2 = col0
      do mm = -1, 0
        idx2 = idx2 + 1
        Icsz(idx2)   = Icsz(idx2) &
                     + real(fact(ij-1),r64)*(Icrb(idx2)*factinv(ij))
        Iasz(:,idx2) = Iasz(:,idx2) + Iarb(:,idx2)
        do l = 1, nterms-ij
          ll    = l + ij
          frdc  = fr(l)*dc(ll+mm,l)*dc(ll-mm,l)
          frdc1 = real(fact(ij-1)*fact(ll-ij),r64)*(frdc*factinv(ll))
          col   = (ll+1)*ll/2 + ll + mm + 1
          Icsz(col)   = Icsz(col)   + frdc1*Icrb(idx2)
          Iasz(:,col) = Iasz(:,col) + frdc *Iarb(:,idx2)
        end do
        colr = col0r + 1
        szj  = Icsz(idx2);    cszj  = conjg(szj)
        aszj = Iasz(:,idx2);  caszj = conjg(aszj)
        mp = 1
        if (mm < 0_8) then
          erfn = ephi(mp)*rotmatf(mp, mm,ij)
          erfp = ephi(mp)*rotmatf(mp,-mm,ij)
          Icnew(colr)   = Icnew(colr)   + erfn*szj  + erfp*cszj
          Ianew(:,colr) = Ianew(:,colr) + erfn*aszj + erfp*caszj
        else
          erf = ephi(mp)*rotmatf(mp,mm,ij)
          Icnew(colr)   = Icnew(colr)   + erf*szj
          Ianew(:,colr) = Ianew(:,colr) + erf*aszj
        end if
      end do

      do ij = 2, nterms
        col0  = ij**2
        col02 = ij*(ij+1)/2
        col0r = (ij*(ij-1))/2

        idx2 = col02
        do mm = -ij, -1
          idx2 = idx2 + 1
          c0j  = Icold(idx2);    cc0j = conjg(c0j)
          a0j  = Iaold(:,idx2);  ca0j = conjg(a0j)
          col2 = col02 + 1
          do mp = -ij, 0
            erbn = ephi(-mm)*rotmatb(mp, mm,ij)
            erbp = ephi( mm)*rotmatb(mp,-mm,ij)
            Icrb(col2)   = Icrb(col2)   + erbn*c0j + erbp*cc0j
            Iarb(:,col2) = Iarb(:,col2) + erbn*a0j + erbp*ca0j
            col2 = col2 + 1
          end do
        end do
        mm = 0
        idx2 = idx2 + 1;  col2 = col02 + 1
        c0j = Icold(idx2);  a0j = Iaold(:,idx2)
        do mp = -ij, 0
          erb = ephi(-mm)*rotmatb(mp,mm,ij)
          Icrb(col2)   = Icrb(col2)   + erb*c0j
          Iarb(:,col2) = Iarb(:,col2) + erb*a0j
          col2 = col2 + 1
        end do

        idx2 = col02
        tv  = real(fact(ij-2),r64)*factinv(ij-1)
        tv2 = real(fact(ij-1),r64)*factinv(ij)
        do mm = -ij, 0
          idx2 = idx2 + 1
          Icsz(idx2)   = Icsz(idx2)   + tv2*Icrb(idx2)
          Iasz(:,idx2) = Iasz(:,idx2) + tv *Iarb(:,idx2)
          do l = 1, nterms-ij
            ll    = l + ij
            frdc  = fr(l)*dc(ll+mm,l)*dc(ll-mm,l)
            frdc1 = real(fact(ij-1)*fact(ll-ij),r64)*(frdc*factinv(ll))
            frdc2 = real(fact(ij-2)*fact(ll-ij),r64)*(frdc*factinv(ll-1))
            col   = (ll+1)*ll/2 + ll + mm + 1
            Icsz(col)   = Icsz(col)   + frdc1*Icrb(idx2)
            Iasz(:,col) = Iasz(:,col) + frdc2*Iarb(:,idx2)
          end do
          colr = col0r + 1
          szj  = Icsz(idx2);    cszj  = conjg(szj)
          aszj = Iasz(:,idx2);  caszj = conjg(aszj)
          do mp = 1, ij
            if (mm < 0_8) then
              erfn = ephi(mp)*rotmatf(mp, mm,ij)
              erfp = ephi(mp)*rotmatf(mp,-mm,ij)
              Icnew(colr)   = Icnew(colr)   + erfn*szj  + erfp*cszj
              Ianew(:,colr) = Ianew(:,colr) + erfn*aszj + erfp*caszj
            else
              erf = ephi(mp)*rotmatf(mp,mm,ij)
              Icnew(colr)   = Icnew(colr)   + erf*szj
              Ianew(:,colr) = Ianew(:,colr) + erf*aszj
            end if
            colr = colr + 1
          end do
        end do
      end do

      omega_slp(:,j) = aimag(Icnew)
      do k = 1, 4
        omega(:,j,k) = aimag(Ianew(k,:))
      end do

    end do

  end subroutine qao_omegasdlp_r64

end module omega_mod
