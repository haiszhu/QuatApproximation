module quatapproximation_mod
  ! ------------------------------------------------------------------
  ! Primitives module for QuatApproximation-legacy.
  !
  ! Plays the same role as LineQuaaadrature-legacy's linequaaadrature_mod:
  ! owns the package's kind parameters (r64, r128, c64, c128) plus
  ! low-level utilities that are too small or too generic to deserve
  ! their own module.  Other modules in the package (harmonic_mod,
  ! koorn_geom_mod, ...) `use quatapproximation_mod, only: r64, r128, ...`
  ! rather than re-defining the kinds themselves.
  !
  ! Per the parallel-independent rule (project_parallel_legacy_packages),
  ! no cross-package `use` ever appears — this module is the sibling of
  ! linequaaadrature_mod, not a child.
  ! ------------------------------------------------------------------
  implicit none

  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16
  integer, parameter :: c64  = 8
  integer, parameter :: c128 = 16

contains

  ! ------------------------------------------------------------------
  ! simplex2equil_r64
  ! Map npts points uvs(2, npts) on the canonical simplex
  ! {u, v >= 0, u + v <= 1} to an equilateral triangle of side length
  ! h_side centered at the origin.  Mirrors the local MATLAB helper in
  ! test/harmonic_approx/test_paraboloid_cond.m.
  ! ------------------------------------------------------------------
  subroutine simplex2equil_r64(npts, h_side, uvs, uvs_eq)
    integer(8), intent(in)    :: npts
    real(r64),  intent(in)    :: h_side
    real(r64),  intent(in)    :: uvs(2, npts)
    real(r64),  intent(inout) :: uvs_eq(2, npts)

    real(r64) :: h_c, h_b, v1(2), v2(2), v3(2), edge1(2), edge2(2)
    integer(8) :: i

    h_c = sqrt(3.0_r64) / 3.0_r64 * h_side
    h_b = sqrt(3.0_r64) / 6.0_r64 * h_side
    v1  = [ 0.0_r64,         h_c]
    v2  = [-0.5_r64*h_side, -h_b]
    v3  = [ 0.5_r64*h_side, -h_b]
    edge1 = v2 - v1
    edge2 = v3 - v1
    do i = 1, npts
      uvs_eq(:, i) = v1 + edge1 * uvs(1, i) + edge2 * uvs(2, i)
    end do
  end subroutine simplex2equil_r64

  ! ------------------------------------------------------------------
  ! gauss_r64
  ! Gauss-Legendre nodes tgl, weights wgl, and differentiation matrix Dgl
  ! on [-1, 1].  Newton-Raphson on P_n (cosine seed) + Fornberg D-matrix.
  ! Vendored byte-faithfully from LineQuaaadrature-legacy/src/linequaaadrature_mod.f90.
  ! ------------------------------------------------------------------
  subroutine gauss_r64(n, tgl, wgl, Dgl)
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: tgl(n), wgl(n), Dgl(n,n)

    integer(8) :: i, j, k, m
    real(r64)  :: z, z1, p0, p1, p2, pp, pi, tol
    real(r64)  :: a(n)

    pi  = acos(-1.0_r64)
    tol = epsilon(1.0_r64)
    m   = (n + 1_8) / 2_8

    do i = 1, m
      z = cos(pi * (real(i,r64) - 0.25_r64) / (real(n,r64) + 0.5_r64))
      do
        p0 = 1.0_r64;  p1 = z
        do k = 2, n
          p2 = ((2.0_r64*real(k,r64) - 1.0_r64)*z*p1 &
               - (real(k,r64) - 1.0_r64)*p0) / real(k,r64)
          p0 = p1;  p1 = p2
        end do
        pp = real(n,r64) * (z*p1 - p0) / (z*z - 1.0_r64)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= tol * max(1.0_r64, abs(z))) exit
      end do
      tgl(i)       = -z
      tgl(n+1-i)   =  z
      wgl(i)       = 2.0_r64 / ((1.0_r64 - z*z) * pp*pp)
      wgl(n+1-i)   = wgl(i)
    end do

    do k = 1, n
      a(k) = 1.0_r64
      do j = 1, n
        if (j /= k) a(k) = a(k) * (tgl(k) - tgl(j))
      end do
    end do
    do k = 1, n
      do j = 1, n
        if (j /= k) then
          Dgl(j,k) = (a(j)/a(k)) / (tgl(j) - tgl(k))
        end if
      end do
      Dgl(k,k) = 0.0_r64
      do j = 1, n
        if (j /= k) Dgl(k,k) = Dgl(k,k) + 1.0_r64/(tgl(k) - tgl(j))
      end do
    end do

  end subroutine gauss_r64

  ! ------------------------------------------------------------------
  ! gauss_r128
  ! r128 twin of gauss_r64.  Same Newton iteration + barycentric D
  ! construction, just real(16) throughout.  Used by the r128 chain
  ! that mirrors Lap3dDLP_closepanel_demo for future r128 orchestrations.
  ! ------------------------------------------------------------------
  subroutine gauss_r128(n, tgl, wgl, Dgl)
    integer(8), intent(in)    :: n
    real(r128), intent(inout) :: tgl(n), wgl(n), Dgl(n,n)

    integer(8) :: i, j, k, m
    real(r128) :: z, z1, p0, p1, p2, pp, pi, tol
    real(r128) :: a(n)

    pi  = acos(-1.0_r128)
    tol = epsilon(1.0_r128)
    m   = (n + 1_8) / 2_8

    do i = 1, m
      z = cos(pi * (real(i,r128) - 0.25_r128) / (real(n,r128) + 0.5_r128))
      do
        p0 = 1.0_r128;  p1 = z
        do k = 2, n
          p2 = ((2.0_r128*real(k,r128) - 1.0_r128)*z*p1 &
               - (real(k,r128) - 1.0_r128)*p0) / real(k,r128)
          p0 = p1;  p1 = p2
        end do
        pp = real(n,r128) * (z*p1 - p0) / (z*z - 1.0_r128)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= tol * max(1.0_r128, abs(z))) exit
      end do
      tgl(i)       = -z
      tgl(n+1-i)   =  z
      wgl(i)       = 2.0_r128 / ((1.0_r128 - z*z) * pp*pp)
      wgl(n+1-i)   = wgl(i)
    end do

    do k = 1, n
      a(k) = 1.0_r128
      do j = 1, n
        if (j /= k) a(k) = a(k) * (tgl(k) - tgl(j))
      end do
    end do
    do k = 1, n
      do j = 1, n
        if (j /= k) then
          Dgl(j,k) = (a(j)/a(k)) / (tgl(j) - tgl(k))
        end if
      end do
      Dgl(k,k) = 0.0_r128
      do j = 1, n
        if (j /= k) Dgl(k,k) = Dgl(k,k) + 1.0_r128/(tgl(k) - tgl(j))
      end do
    end do

  end subroutine gauss_r128

  ! ------------------------------------------------------------------
  ! bclaginterpweights_r64
  ! Barycentric Lagrange interpolation weights w_k = 1/prod_{j/=k}(t_k - t_j)
  ! for the type-II barycentric formula at the given GL nodes tgl.
  ! Vendored byte-faithfully from LineQuaaadrature-legacy/src/linequaaadrature_mod.f90.
  ! ------------------------------------------------------------------
  subroutine bclaginterpweights_r64(n, tgl, w_bclag)
    integer(8), intent(in)    :: n
    real(r64),  intent(in)    :: tgl(n)
    real(r64),  intent(inout) :: w_bclag(n)

    integer(8) :: j, k
    real(r64)  :: prod

    do k = 1, n
      prod = 1.0_r64
      do j = 1, n
        if (j /= k) prod = prod * (tgl(k) - tgl(j))
      end do
      w_bclag(k) = 1.0_r64 / prod
    end do

  end subroutine bclaginterpweights_r64

  ! ------------------------------------------------------------------
  ! paraboloid_uv_r64
  ! Elliptic paraboloid parametrization r(u, v) = [u; v; (H/2)(u^2+v^2)].
  ! Eq. (5.1) of Zhu-Jiang (lapquad5_cpam_rev1.tex, §5).  Returns
  ! x(3, npts) for npts input (u, v) points.
  ! ------------------------------------------------------------------
  subroutine paraboloid_uv_r64(npts, H, u, v, x)
    integer(8), intent(in)    :: npts
    real(r64),  intent(in)    :: H
    real(r64),  intent(in)    :: u(npts), v(npts)
    real(r64),  intent(inout) :: x(3, npts)

    integer(8) :: i
    do i = 1, npts
      x(1, i) = u(i)
      x(2, i) = v(i)
      x(3, i) = 0.5_r64 * H * (u(i)**2 + v(i)**2)
    end do
  end subroutine paraboloid_uv_r64

  ! ------------------------------------------------------------------
  ! simplex2equil_with_detj_r64
  ! Same map as simplex2equil_r64 but additionally returns detJ, the
  ! Jacobian determinant of the canonical-simplex -> equilateral
  ! transform.  For an equilateral triangle of side h_side this is the
  ! constant sqrt(3)/2 * h_side**2.
  ! ------------------------------------------------------------------
  subroutine simplex2equil_with_detj_r64(npts, h_side, uvs, uvs_eq, detJ)
    integer(8), intent(in)    :: npts
    real(r64),  intent(in)    :: h_side
    real(r64),  intent(in)    :: uvs(2, npts)
    real(r64),  intent(inout) :: uvs_eq(2, npts)
    real(r64),  intent(inout) :: detJ

    call simplex2equil_r64(npts, h_side, uvs, uvs_eq)
    detJ = sqrt(3.0_r64) / 2.0_r64 * h_side * h_side
  end subroutine simplex2equil_with_detj_r64

  ! ------------------------------------------------------------------
  ! paraboloidparam_r64
  ! Downward-opening elliptic paraboloid parametrization used by
  ! test_paraboloid_refinement.m (mirrors the OSSD MATLAB local helper):
  !   r(t, p) = [t; p; h * (1 - (t^2 + p^2)/a^2)],
  !   rt     = [1; 0; -2*h*t/a^2],
  !   rp     = [0; 1; -2*h*p/a^2],
  !   rtt    = [0; 0; -2*h/a^2],
  !   rpp    = [0; 0; -2*h/a^2],
  !   rtp    = [0; 0; 0],
  !   nx     = cross(rp, rt) / |cross(rp, rt)|,    sp = |cross(rp, rt)|.
  ! All outputs are npts-long.  Distinct from paraboloid_uv_r64, which is
  ! the upward-opening (H/2)(u^2+v^2) form from eq. (5.1) of Zhu-Jiang.
  ! ------------------------------------------------------------------
  subroutine paraboloidparam_r64(npts, h, a, t, p, x, nx, sp, &
                                  rts, rps, rtts, rpps, rtps)
    integer(8), intent(in)    :: npts
    real(r64),  intent(in)    :: h, a
    real(r64),  intent(in)    :: t(npts), p(npts)
    real(r64),  intent(inout) :: x(3, npts), nx(3, npts), sp(npts)
    real(r64),  intent(inout) :: rts(3, npts), rps(3, npts)
    real(r64),  intent(inout) :: rtts(3, npts), rpps(3, npts), rtps(3, npts)

    integer(8) :: i
    real(r64)  :: a2_inv, cx(3)

    a2_inv = 1.0_r64 / (a*a)

    do i = 1, npts
      x(1, i) = t(i)
      x(2, i) = p(i)
      x(3, i) = h * (1.0_r64 - (t(i)*t(i) + p(i)*p(i)) * a2_inv)

      rts(1, i) = 1.0_r64
      rts(2, i) = 0.0_r64
      rts(3, i) = -2.0_r64 * h * t(i) * a2_inv

      rps(1, i) = 0.0_r64
      rps(2, i) = 1.0_r64
      rps(3, i) = -2.0_r64 * h * p(i) * a2_inv

      ! cross(rp, rt)
      cx(1) = rps(2,i)*rts(3,i) - rps(3,i)*rts(2,i)
      cx(2) = rps(3,i)*rts(1,i) - rps(1,i)*rts(3,i)
      cx(3) = rps(1,i)*rts(2,i) - rps(2,i)*rts(1,i)
      sp(i) = sqrt(cx(1)*cx(1) + cx(2)*cx(2) + cx(3)*cx(3))
      nx(1, i) = cx(1) / sp(i)
      nx(2, i) = cx(2) / sp(i)
      nx(3, i) = cx(3) / sp(i)

      rtts(1, i) = 0.0_r64
      rtts(2, i) = 0.0_r64
      rtts(3, i) = -2.0_r64 * h * a2_inv

      rpps(1, i) = 0.0_r64
      rpps(2, i) = 0.0_r64
      rpps(3, i) = -2.0_r64 * h * a2_inv

      rtps(1, i) = 0.0_r64
      rtps(2, i) = 0.0_r64
      rtps(3, i) = 0.0_r64
    end do

  end subroutine paraboloidparam_r64

  ! ------------------------------------------------------------------
  ! subdivide_simplex_r64
  ! One uniform-refinement step on a simplex triangle.  A triangle is
  ! encoded as (Cx, Cy, scale, orientation); the four child triangles
  ! tile the parent.  When orientation ~= pi the four children are
  ! three "upward" subtriangles plus one centered "downward" child;
  ! when orientation == pi the offsets are negated and the orientation
  ! flips accordingly.  Mirrors the MATLAB local helper.
  ! ------------------------------------------------------------------
  subroutine subdivide_simplex_r64(Cx, Cy, scale, orientation, sub_tris)
    real(r64), intent(in)    :: Cx, Cy, scale, orientation
    real(r64), intent(inout) :: sub_tris(4, 4)

    real(r64) :: a_sub, offsets(4, 3)
    integer(8) :: i

    a_sub = scale / 2.0_r64
    offsets(1, 1) = -scale / 6.0_r64
    offsets(1, 2) = -scale / 6.0_r64
    offsets(1, 3) = 0.0_r64
    offsets(2, 1) =  scale / 3.0_r64
    offsets(2, 2) = -scale / 6.0_r64
    offsets(2, 3) = 0.0_r64
    offsets(3, 1) = -scale / 6.0_r64
    offsets(3, 2) =  scale / 3.0_r64
    offsets(3, 3) = 0.0_r64
    offsets(4, 1) = 0.0_r64
    offsets(4, 2) = 0.0_r64
    offsets(4, 3) = acos(-1.0_r64)   ! pi

    if (abs(orientation - acos(-1.0_r64)) < 1.0e-10_r64) then
      do i = 1, 4
        offsets(i, 1) = -offsets(i, 1)
        offsets(i, 2) = -offsets(i, 2)
        offsets(i, 3) = modulo(offsets(i, 3) + acos(-1.0_r64), 2.0_r64 * acos(-1.0_r64))
      end do
    end if

    do i = 1, 4
      sub_tris(i, 1) = Cx + offsets(i, 1)
      sub_tris(i, 2) = Cy + offsets(i, 2)
      sub_tris(i, 3) = a_sub
      sub_tris(i, 4) = offsets(i, 3)
    end do
  end subroutine subdivide_simplex_r64

  ! ------------------------------------------------------------------
  ! assemble_subdivided_nodes_r64
  ! Lay out N reference simplex-VR nodes inside each of M subdivided
  ! triangles described by sub_tris(M, 4) = [Cx, Cy, scale, orientation].
  ! For each child triangle, shift the reference nodes by -C0 = -(1/3,1/3),
  ! optionally negate if orientation == pi, then scale and re-center at C.
  ! Output is the flat 2 x (N*M) concatenation in triangle-major order.
  ! ------------------------------------------------------------------
  subroutine assemble_subdivided_nodes_r64(N, M, uvs_simplex, sub_tris, uvs_all)
    integer(8), intent(in)    :: N, M
    real(r64),  intent(in)    :: uvs_simplex(2, N)
    real(r64),  intent(in)    :: sub_tris(M, 4)
    real(r64),  intent(inout) :: uvs_all(2, N*M)

    real(r64) :: C0(2), pts_local(2, N), Ck(2), sk, theta, pi_val
    integer(8) :: k, j, idx_start

    C0(1) = 1.0_r64 / 3.0_r64
    C0(2) = 1.0_r64 / 3.0_r64
    pi_val = acos(-1.0_r64)

    do k = 1, M
      Ck(1) = sub_tris(k, 1)
      Ck(2) = sub_tris(k, 2)
      sk    = sub_tris(k, 3)
      theta = sub_tris(k, 4)
      do j = 1, N
        pts_local(1, j) = uvs_simplex(1, j) - C0(1)
        pts_local(2, j) = uvs_simplex(2, j) - C0(2)
      end do
      if (abs(theta - pi_val) < 1.0e-10_r64) then
        do j = 1, N
          pts_local(1, j) = -pts_local(1, j)
          pts_local(2, j) = -pts_local(2, j)
        end do
      end if
      idx_start = (k-1)*N
      do j = 1, N
        uvs_all(1, idx_start + j) = Ck(1) + sk * pts_local(1, j)
        uvs_all(2, idx_start + j) = Ck(2) + sk * pts_local(2, j)
      end do
    end do
  end subroutine assemble_subdivided_nodes_r64

  ! ------------------------------------------------------------------
  ! get_subdivided_triangles_r64
  ! Recursively apply subdivide_simplex_r64 L times to a root triangle
  ! described by (C0(2), s0, o0).  Output all_tris is (4^L, 4); caller
  ! pre-allocates with that shape.
  ! ------------------------------------------------------------------
  subroutine get_subdivided_triangles_r64(C0, s0, o0, L, all_tris)
    real(r64),  intent(in)    :: C0(2), s0, o0
    integer(8), intent(in)    :: L
    real(r64),  intent(inout) :: all_tris(4**L, 4)

    integer(8) :: ell, num_parents, i
    real(r64), allocatable :: cur(:,:), nxt(:,:)
    real(r64) :: subs(4, 4)

    allocate(cur(1, 4))
    cur(1, 1) = C0(1)
    cur(1, 2) = C0(2)
    cur(1, 3) = s0
    cur(1, 4) = o0

    do ell = 1, L
      num_parents = size(cur, 1, kind=8)
      allocate(nxt(num_parents * 4, 4))
      do i = 1, num_parents
        call subdivide_simplex_r64(cur(i, 1), cur(i, 2), cur(i, 3), cur(i, 4), subs)
        nxt((i-1)*4 + 1 : i*4, :) = subs
      end do
      deallocate(cur)
      allocate(cur(size(nxt, 1, kind=8), 4))
      cur = nxt
      deallocate(nxt)
    end do

    all_tris = cur
    deallocate(cur)
  end subroutine get_subdivided_triangles_r64

end module quatapproximation_mod
