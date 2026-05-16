! test_paraboloid_refinement.f90
!
! Fortran twin of test/harmonic_approx/test_paraboloid_refinement.m:
! quaternion harmonic-approximation error on a paraboloid patch swept
! over polynomial order p in {4,6,8,10,12,14} and refinement level
! L in {0..5} (num_patches = 4^L).  Reproduces the Fig 5.1 RIGHT
! convergence picture from Zhu-Jiang.
!
! For each (p, L):
!   - Subdivide the simplex L times via get_subdivided_triangles
!   - Load VR nodes + weights of order p-1
!   - Lay out node sets across the 4^L sub-triangles
!   - Map to equilateral coords, evaluate paraboloid (downward form)
!     plus 1st/2nd partials, build mean-curvature H per node
!   - Pick the sub-triangle with the highest |H| as the worst patch
!   - Repeat for a uniform reference target raster (11x11 simplex grid
!     filtered to u+v<=1) on the chosen sub-triangle
!   - Build the patch's edge quadrature, compute the circumcircle
!     transform, apply alpha*R*(x-c) inline
!   - Evaluate the 3D Laplace solid spherical harmonics at source nodes
!   - Subset to (l, m), 1 <= m <= l <= p; assemble the 4x4 block A
!   - cond(A) via dgesvd; q_coeffs = A \ [I; 0] * H_src
!   - Evaluate basis at target nodes, reconstruct H_tgt_num, compare
!     to H_tgt_ref, record MyErr(p, L) = max|err| / max|H_tgt_ref|.
!
! Prints MyErr and MyCond tables to stdout.  Vertex permutation
! r_vert(:, [1 3 2]) flips the circumcircle +z choice so the patch bulges
! toward +z.  4.25*sx3min*nc off-surface shift is skipped (matches the
! MATLAB twin).

program test_paraboloid_refinement
  use quatapproximation_mod, only: r64, &
      simplex2equil_with_detj_r64, paraboloidparam_r64, &
      get_subdivided_triangles_r64, assemble_subdivided_nodes_r64, &
      gauss_r64, bclaginterpweights_r64
  use koorn_geom_mod, only: get_vioreanu_nodes, get_vioreanu_wts, &
      koorn_vals2coefs_coefs2vals, line3quadr_3dline, &
      circumcircle_transform_3d
  use harmonic_mod, only: l3dtavecevalmat_r64
  implicit none

  integer(8), parameter :: norders = 6
  integer(8), parameter :: nL      = 6
  integer(8), parameter :: orders(norders) = [4_8, 6_8, 8_8, 10_8, 12_8, 14_8]
  integer(8), parameter :: Ls(nL)         = [0_8, 1_8, 2_8, 3_8, 4_8, 5_8]
  real(r64),  parameter :: h = 1.0_r64
  real(r64),  parameter :: a = 1.0_r64
  logical,    parameter :: if_gauss = .false.   ! mean curvature H by default

  real(r64) :: MyErr (norders, nL)
  real(r64) :: MyCond(norders, nL)

  ! Reference target grid on the canonical simplex (uniform 11x11 raster
  ! filtered to u+v <= 1).
  integer(8), parameter :: Nraster = 11
  real(r64),  allocatable :: uvs_simplex_ref(:,:)
  integer(8) :: nref

  call build_reference_target_grid(Nraster, uvs_simplex_ref, nref)

  call run_sweep()

  call print_tables()

  call write_dat_files()

  call plot_with_gnuplot()

contains

  ! ----------------------------------------------------------------
  subroutine build_reference_target_grid(Nr, uvs, nout)
    integer(8), intent(in)  :: Nr
    real(r64),  allocatable, intent(out) :: uvs(:,:)
    integer(8), intent(out) :: nout
    integer(8) :: i, j, k
    real(r64)  :: u, v
    real(r64), allocatable :: tmp(:,:)
    allocate(tmp(2, Nr*Nr))
    k = 0
    do j = 1, Nr
      v = real(j-1, r64) / real(Nr-1, r64)
      do i = 1, Nr
        u = real(i-1, r64) / real(Nr-1, r64)
        if (u + v <= 1.0_r64 + 1.0e-10_r64) then
          k = k + 1
          tmp(1, k) = u
          tmp(2, k) = v
        end if
      end do
    end do
    nout = k
    allocate(uvs(2, nout))
    uvs = tmp(:, 1:nout)
    deallocate(tmp)
  end subroutine build_reference_target_grid

  ! ----------------------------------------------------------------
  subroutine run_sweep()
    integer(8) :: idx, jL, p, L, n_p_per, num_tris, n_sub, n_ref_sub
    integer(8) :: ier, kpatch, idx_start
    integer(8) :: ij, kk, tmpidx, tmpidx2, nquad, sbdnp, nbd, nbasis_loc, n, kc

    real(r64), allocatable :: uvs_simplex(:,:), wts_simplex(:)
    real(r64), allocatable :: all_tris(:,:)
    real(r64), allocatable :: uvs_simplex_sub(:,:), uvs_eq_sub(:,:)
    real(r64), allocatable :: x_eq_sub(:,:), nx_eq_sub(:,:), sp_eq_sub(:)
    real(r64), allocatable :: rts_s(:,:), rps_s(:,:), rtts_s(:,:), rpps_s(:,:), rtps_s(:,:)
    real(r64), allocatable :: cur_full(:), cur(:)
    real(r64), allocatable :: Es(:), Fs(:), Gs(:), Lss(:), Mss(:), Nss(:)

    real(r64), allocatable :: uvs_simplex_ref_k(:,:), uvs_eq_ref_k(:,:)
    real(r64), allocatable :: tx_eq_ref_k(:,:), nx_eq_ref_k(:,:), sp_ref_k(:)
    real(r64), allocatable :: rts_t(:,:), rps_t(:,:), rtts_t(:,:), rpps_t(:,:), rtps_t(:,:)
    real(r64), allocatable :: cur_ref(:)
    real(r64), allocatable :: Et(:), Ft(:), Gt(:), Ltt(:), Mtt(:), Ntt(:)

    real(r64), allocatable :: sx(:,:), snx(:,:), tx(:,:)
    real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:), w_bclag(:)
    real(r64), allocatable :: umatr(:,:), vmatr(:,:)
    real(r64), allocatable :: tpan(:), sxbd(:,:), swbd(:), stangbd(:,:), sspbd(:)
    real(r64)               :: r_vert(3, 3), r_vert_perm(3, 3)
    real(r64)               :: R_circ(3, 3), c_circ(3), alpha_circ
    real(r64), allocatable :: sxnew(:,:), snxnew(:,:), txnew(:,:)
    real(r64), allocatable :: sxbd2(:,:), swbd2(:), stangbd2(:,:), sspbd2(:)

    integer(8), allocatable :: idxvec(:)
    complex(8), allocatable :: F_c(:,:), Fx_c(:,:), Fy_c(:,:), Fz_c(:,:)
    complex(8), allocatable :: Ft_c(:,:), Ftx_c(:,:), Fty_c(:,:), Ftz_c(:,:)
    real(r64),  allocatable :: F1(:,:), F2(:,:), F3(:,:), F0(:,:)
    real(r64),  allocatable :: Ft1(:,:), Ft2(:,:), Ft3(:,:), Ft0(:,:)
    real(r64),  allocatable :: Amat(:,:), Awork(:,:), rhs(:), q_coeffs(:)
    real(r64),  allocatable :: q_tgt(:), cur_tgt_num(:)
    real(r64)               :: tmpdiff, tmpval, perpatch_max, mxabs

    do idx = 1, norders
      p = orders(idx)
      n_p_per   = p*(p+1)/2

      do jL = 1, nL
        L = Ls(jL)
        num_tris = 4_8**L

        if (allocated(all_tris)) deallocate(all_tris)
        allocate(all_tris(num_tris, 4))
        call get_subdivided_triangles_r64([1.0_r64/3.0_r64, 1.0_r64/3.0_r64], &
                                           1.0_r64, 0.0_r64, L, all_tris)

        if (allocated(uvs_simplex)) deallocate(uvs_simplex)
        if (allocated(wts_simplex)) deallocate(wts_simplex)
        allocate(uvs_simplex(2, n_p_per), wts_simplex(n_p_per))
        call get_vioreanu_nodes(p-1, n_p_per, uvs_simplex)
        call get_vioreanu_wts  (p-1, n_p_per, wts_simplex)

        n_sub = n_p_per * num_tris
        if (allocated(uvs_simplex_sub)) deallocate(uvs_simplex_sub)
        if (allocated(uvs_eq_sub))      deallocate(uvs_eq_sub)
        if (allocated(x_eq_sub))        deallocate(x_eq_sub)
        if (allocated(nx_eq_sub))       deallocate(nx_eq_sub)
        if (allocated(sp_eq_sub))       deallocate(sp_eq_sub)
        if (allocated(rts_s))           deallocate(rts_s, rps_s, rtts_s, rpps_s, rtps_s)
        if (allocated(Es))              deallocate(Es, Fs, Gs, Lss, Mss, Nss, cur_full)
        allocate(uvs_simplex_sub(2, n_sub), uvs_eq_sub(2, n_sub))
        allocate(x_eq_sub(3, n_sub), nx_eq_sub(3, n_sub), sp_eq_sub(n_sub))
        allocate(rts_s (3, n_sub), rps_s (3, n_sub))
        allocate(rtts_s(3, n_sub), rpps_s(3, n_sub), rtps_s(3, n_sub))
        allocate(Es(n_sub), Fs(n_sub), Gs(n_sub), Lss(n_sub), Mss(n_sub), Nss(n_sub))
        allocate(cur_full(n_sub))

        call assemble_subdivided_nodes_r64(n_p_per, num_tris, uvs_simplex, all_tris, uvs_simplex_sub)
        block
          real(r64) :: detJ
          call simplex2equil_with_detj_r64(n_sub, 1.0_r64, uvs_simplex_sub, uvs_eq_sub, detJ)
        end block
        call paraboloidparam_r64(n_sub, h, a, uvs_eq_sub(1,:), uvs_eq_sub(2,:), &
                                  x_eq_sub, nx_eq_sub, sp_eq_sub, &
                                  rts_s, rps_s, rtts_s, rpps_s, rtps_s)

        ! First/second fundamental forms -> K or H per node.
        Es  = sum(rts_s  * rts_s , dim=1)
        Fs  = sum(rts_s  * rps_s , dim=1)
        Gs  = sum(rps_s  * rps_s , dim=1)
        Lss = sum(rtts_s * nx_eq_sub, dim=1)
        Mss = sum(rtps_s * nx_eq_sub, dim=1)
        Nss = sum(rpps_s * nx_eq_sub, dim=1)
        if (if_gauss) then
          cur_full = (Lss*Nss - Mss*Mss) / (Es*Gs - Fs*Fs)
        else
          cur_full = (Es*Nss + Gs*Lss - 2.0_r64*Fs*Mss) / (2.0_r64*(Es*Gs - Fs*Fs))
        end if

        ! Pick sub-triangle with max |cur_full|.
        kpatch = 1
        mxabs  = -1.0_r64
        do kc = 1, num_tris
          idx_start = (kc-1)*n_p_per
          perpatch_max = maxval(abs(cur_full(idx_start+1 : idx_start+n_p_per)))
          if (perpatch_max > mxabs) then
            mxabs  = perpatch_max
            kpatch = kc
          end if
        end do
        idx_start = (kpatch-1)*n_p_per

        if (allocated(cur)) deallocate(cur)
        allocate(cur(n_p_per))
        cur = cur_full(idx_start+1 : idx_start+n_p_per)

        ! Reference targets for the chosen sub-triangle only.
        n_ref_sub = nref
        if (allocated(uvs_simplex_ref_k)) deallocate(uvs_simplex_ref_k, uvs_eq_ref_k)
        if (allocated(tx_eq_ref_k))       deallocate(tx_eq_ref_k, nx_eq_ref_k, sp_ref_k)
        if (allocated(rts_t))             deallocate(rts_t, rps_t, rtts_t, rpps_t, rtps_t)
        if (allocated(Et))                deallocate(Et, Ft, Gt, Ltt, Mtt, Ntt, cur_ref)
        allocate(uvs_simplex_ref_k(2, n_ref_sub), uvs_eq_ref_k(2, n_ref_sub))
        allocate(tx_eq_ref_k(3, n_ref_sub), nx_eq_ref_k(3, n_ref_sub), sp_ref_k(n_ref_sub))
        allocate(rts_t (3, n_ref_sub), rps_t (3, n_ref_sub))
        allocate(rtts_t(3, n_ref_sub), rpps_t(3, n_ref_sub), rtps_t(3, n_ref_sub))
        allocate(Et(n_ref_sub), Ft(n_ref_sub), Gt(n_ref_sub), Ltt(n_ref_sub), Mtt(n_ref_sub), Ntt(n_ref_sub))
        allocate(cur_ref(n_ref_sub))

        call assemble_subdivided_nodes_r64(n_ref_sub, 1_8, uvs_simplex_ref, &
                                            reshape(all_tris(kpatch, :), [1_8, 4_8]), &
                                            uvs_simplex_ref_k)
        block
          real(r64) :: detJ
          call simplex2equil_with_detj_r64(n_ref_sub, 1.0_r64, uvs_simplex_ref_k, uvs_eq_ref_k, detJ)
        end block
        call paraboloidparam_r64(n_ref_sub, h, a, uvs_eq_ref_k(1,:), uvs_eq_ref_k(2,:), &
                                  tx_eq_ref_k, nx_eq_ref_k, sp_ref_k, &
                                  rts_t, rps_t, rtts_t, rpps_t, rtps_t)

        Et  = sum(rts_t  * rts_t , dim=1)
        Ft  = sum(rts_t  * rps_t , dim=1)
        Gt  = sum(rps_t  * rps_t , dim=1)
        Ltt = sum(rtts_t * nx_eq_ref_k, dim=1)
        Mtt = sum(rtps_t * nx_eq_ref_k, dim=1)
        Ntt = sum(rpps_t * nx_eq_ref_k, dim=1)
        if (if_gauss) then
          cur_ref = (Ltt*Ntt - Mtt*Mtt) / (Et*Gt - Ft*Ft)
        else
          cur_ref = (Et*Ntt + Gt*Ltt - 2.0_r64*Ft*Mtt) / (2.0_r64*(Et*Gt - Ft*Ft))
        end if

        ! ---------- GL/bclag/Koornwinder setup ----------
        nquad = p
        if (allocated(tgl)) deallocate(tgl, wgl, Dgl, w_bclag)
        allocate(tgl(nquad), wgl(nquad), Dgl(nquad, nquad), w_bclag(nquad))
        call gauss_r64(nquad, tgl, wgl, Dgl)
        call bclaginterpweights_r64(nquad, tgl, w_bclag)

        n = min(n_p_per, 20_8*21_8/2_8)
        if (allocated(umatr)) deallocate(umatr, vmatr)
        allocate(umatr(n, n), vmatr(n, n))
        umatr = 0.0_r64;  vmatr = 0.0_r64
        call koorn_vals2coefs_coefs2vals(p-1, n, umatr, vmatr)

        if (allocated(sx)) deallocate(sx, snx, tx)
        allocate(sx(3, n_p_per), snx(3, n_p_per))
        sx  = x_eq_sub (:, idx_start+1 : idx_start+n_p_per)
        snx = nx_eq_sub(:, idx_start+1 : idx_start+n_p_per)
        allocate(tx(3, n_ref_sub))
        tx  = tx_eq_ref_k

        ! ---------- patch boundary in raw frame ----------
        sbdnp = 3
        nbd   = sbdnp*nquad
        if (allocated(tpan)) deallocate(tpan, sxbd, swbd, stangbd, sspbd)
        allocate(tpan(sbdnp+1), sxbd(3, nbd), swbd(nbd), stangbd(3, nbd), sspbd(nbd))
        do kc = 0, sbdnp
          tpan(kc+1) = real(kc, r64) * 2.0_r64 * acos(-1.0_r64) / real(sbdnp, r64)
        end do
        call line3quadr_3dline(sx, p-1, n, umatr, nquad, tgl, wgl, Dgl, &
                                sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd, r_vert)

        ! ---------- circumcircle transform, apply ----------
        r_vert_perm(:, 1) = r_vert(:, 1)
        r_vert_perm(:, 2) = r_vert(:, 3)
        r_vert_perm(:, 3) = r_vert(:, 2)
        call circumcircle_transform_3d(r_vert_perm, R_circ, c_circ, alpha_circ)

        if (allocated(sxnew)) deallocate(sxnew, snxnew, txnew)
        allocate(sxnew(3, n_p_per), snxnew(3, n_p_per), txnew(3, n_ref_sub))
        do kc = 1, n_p_per
          sxnew (:, kc) = alpha_circ * matmul(R_circ, sx (:, kc) - c_circ)
          snxnew(:, kc) =              matmul(R_circ, snx(:, kc))
        end do
        do kc = 1, n_ref_sub
          txnew (:, kc) = alpha_circ * matmul(R_circ, tx (:, kc) - c_circ)
        end do

        ! Re-run line quadrature in transformed frame (matches .m).
        if (allocated(sxbd2)) deallocate(sxbd2, swbd2, stangbd2, sspbd2)
        allocate(sxbd2(3, nbd), swbd2(nbd), stangbd2(3, nbd), sspbd2(nbd))
        call line3quadr_3dline(sxnew, p-1, n, umatr, nquad, tgl, wgl, Dgl, &
                                sbdnp, tpan, nbd, sxbd2, swbd2, stangbd2, sspbd2, r_vert)

        ! ---------- (l, m), 1 <= m <= l <= p index subset ----------
        if (allocated(idxvec)) deallocate(idxvec)
        allocate(idxvec(n_p_per))
        idxvec = 0;  tmpidx = 0;  tmpidx2 = 0
        do ij = 0, p
          do kk = -ij, ij
            tmpidx2 = tmpidx2 + 1
            if (ij > 0 .and. kk > 0) then
              tmpidx = tmpidx + 1
              idxvec(tmpidx) = tmpidx2
            end if
          end do
        end do

        ! ---------- harmonic basis values + grad at source nodes ----------
        nbasis_loc = (p+1)*(p+1)
        if (allocated(F_c)) deallocate(F_c, Fx_c, Fy_c, Fz_c)
        allocate(F_c(n_p_per, nbasis_loc), Fx_c(n_p_per, nbasis_loc), Fy_c(n_p_per, nbasis_loc), Fz_c(n_p_per, nbasis_loc))
        F_c = (0.0_r64, 0.0_r64);  Fx_c = (0.0_r64, 0.0_r64)
        Fy_c = (0.0_r64, 0.0_r64); Fz_c = (0.0_r64, 0.0_r64)
        ier = 0
        call l3dtavecevalmat_r64(sxnew, n_p_per, p, F_c, Fx_c, Fy_c, Fz_c, ier)

        if (allocated(F1)) deallocate(F1, F2, F3, F0)
        allocate(F1(n_p_per, n_p_per), F2(n_p_per, n_p_per), F3(n_p_per, n_p_per), F0(n_p_per, n_p_per))
        F0 = 0.0_r64
        do kc = 1, n_p_per
          F1(:, kc) = real(Fx_c(:, idxvec(kc)), r64)
          F2(:, kc) = real(Fy_c(:, idxvec(kc)), r64)
          F3(:, kc) = real(Fz_c(:, idxvec(kc)), r64)
        end do

        ! ---------- 4x4 block A ----------
        if (allocated(Amat)) deallocate(Amat, Awork)
        allocate(Amat(4*n_p_per, 4*n_p_per), Awork(4*n_p_per, 4*n_p_per))
        Amat = 0.0_r64
        Amat(1:n_p_per,           n_p_per+1:2*n_p_per) = -F1
        Amat(1:n_p_per,         2*n_p_per+1:3*n_p_per) = -F2
        Amat(1:n_p_per,         3*n_p_per+1:4*n_p_per) = -F3
        Amat(  n_p_per+1:2*n_p_per,           1:n_p_per) =  F1
        Amat(  n_p_per+1:2*n_p_per, 2*n_p_per+1:3*n_p_per) = -F3
        Amat(  n_p_per+1:2*n_p_per, 3*n_p_per+1:4*n_p_per) =  F2
        Amat(2*n_p_per+1:3*n_p_per,           1:n_p_per) =  F2
        Amat(2*n_p_per+1:3*n_p_per,   n_p_per+1:2*n_p_per) =  F3
        Amat(2*n_p_per+1:3*n_p_per, 3*n_p_per+1:4*n_p_per) = -F1
        Amat(3*n_p_per+1:4*n_p_per,           1:n_p_per) =  F3
        Amat(3*n_p_per+1:4*n_p_per,   n_p_per+1:2*n_p_per) = -F2
        Amat(3*n_p_per+1:4*n_p_per, 2*n_p_per+1:3*n_p_per) =  F1

        Awork = Amat
        MyCond(idx, jL) = svd_cond(4_8*n_p_per, Awork)

        ! ---------- linear solve A * q_coeffs = [I; 0] * cur ----------
        if (allocated(rhs)) deallocate(rhs, q_coeffs)
        allocate(rhs(4*n_p_per), q_coeffs(4*n_p_per))
        rhs              = 0.0_r64
        rhs(1:n_p_per)   = cur
        Awork            = Amat
        q_coeffs         = rhs
        call lu_solve(4_8*n_p_per, Awork, q_coeffs)

        ! ---------- harmonic basis at target nodes ----------
        if (allocated(Ft_c)) deallocate(Ft_c, Ftx_c, Fty_c, Ftz_c)
        allocate(Ft_c(n_ref_sub, nbasis_loc), Ftx_c(n_ref_sub, nbasis_loc), &
                 Fty_c(n_ref_sub, nbasis_loc), Ftz_c(n_ref_sub, nbasis_loc))
        Ft_c = (0.0_r64, 0.0_r64);  Ftx_c = (0.0_r64, 0.0_r64)
        Fty_c = (0.0_r64, 0.0_r64); Ftz_c = (0.0_r64, 0.0_r64)
        ier = 0
        call l3dtavecevalmat_r64(txnew, n_ref_sub, p, Ft_c, Ftx_c, Fty_c, Ftz_c, ier)

        if (allocated(Ft1)) deallocate(Ft1, Ft2, Ft3, Ft0)
        allocate(Ft1(n_ref_sub, n_p_per), Ft2(n_ref_sub, n_p_per), &
                 Ft3(n_ref_sub, n_p_per), Ft0(n_ref_sub, n_p_per))
        Ft0 = 0.0_r64
        do kc = 1, n_p_per
          Ft1(:, kc) = real(Ftx_c(:, idxvec(kc)), r64)
          Ft2(:, kc) = real(Fty_c(:, idxvec(kc)), r64)
          Ft3(:, kc) = real(Ftz_c(:, idxvec(kc)), r64)
        end do

        ! q_tgt = Mt * q_coeffs ; take first n_ref_sub entries
        if (allocated(q_tgt)) deallocate(q_tgt, cur_tgt_num)
        allocate(q_tgt(4*n_ref_sub), cur_tgt_num(n_ref_sub))
        q_tgt(           1:  n_ref_sub) = matmul( Ft0, q_coeffs(           1:  n_p_per)) &
                                          - matmul(Ft1, q_coeffs(  n_p_per+1:2*n_p_per)) &
                                          - matmul(Ft2, q_coeffs(2*n_p_per+1:3*n_p_per)) &
                                          - matmul(Ft3, q_coeffs(3*n_p_per+1:4*n_p_per))
        cur_tgt_num = q_tgt(1:n_ref_sub)

        tmpdiff = maxval(abs(cur_ref - cur_tgt_num))
        tmpval  = maxval(abs(cur_ref))
        MyErr(idx, jL) = tmpdiff / tmpval
      end do
    end do
  end subroutine run_sweep

  ! ----------------------------------------------------------------
  subroutine print_tables()
    integer(8) :: idx, jL
    write(*,'(A)') '=== test_paraboloid_refinement (Fig 5.1 RIGHT reproducer) ==='
    write(*,'(A)') ''
    write(*,'(A)') 'MyErr (relative L_inf error of quaternion harmonic approx)'
    write(*,'(A)', advance='no') '   p \\ L   '
    do jL = 1, nL
      write(*,'(2X,I7)', advance='no') Ls(jL)
    end do
    write(*,'(A)') ''
    do idx = 1, norders
      write(*,'(I7,3X)', advance='no') orders(idx)
      do jL = 1, nL
        write(*,'(2X,ES7.1)', advance='no') MyErr(idx, jL)
      end do
      write(*,'(A)') ''
    end do

    write(*,'(A)') ''
    write(*,'(A)') 'MyCond (cond_2 of the 4x4 block approx matrix A)'
    write(*,'(A)', advance='no') '   p \\ L   '
    do jL = 1, nL
      write(*,'(2X,I7)', advance='no') Ls(jL)
    end do
    write(*,'(A)') ''
    do idx = 1, norders
      write(*,'(I7,3X)', advance='no') orders(idx)
      do jL = 1, nL
        write(*,'(2X,ES7.1)', advance='no') MyCond(idx, jL)
      end do
      write(*,'(A)') ''
    end do
  end subroutine print_tables

  ! ----------------------------------------------------------------
  ! write_dat_files
  ! Dump MyErr and MyCond as ASCII tables for plotting.  One row per L,
  ! columns: L, 4^L (num patches), then MyXxx(:, jL) for each polynomial
  ! order in MyOrder.
  ! ----------------------------------------------------------------
  subroutine write_dat_files()
    integer(8) :: jL, idx, unit_err, unit_cond
    open(newunit=unit_err,  file='paraboloid_refinement_myerr.dat',  status='replace')
    open(newunit=unit_cond, file='paraboloid_refinement_mycond.dat', status='replace')
    write(unit_err,  '(A)') '# L  npatches  MyErr(p=4..14)'
    write(unit_cond, '(A)') '# L  npatches  MyCond(p=4..14)'
    do jL = 1, nL
      write(unit_err,  '(I4, 2X, I8, *(2X, ES14.6))') &
          int(Ls(jL)), 4**int(Ls(jL)), (MyErr (idx, jL), idx = 1, norders)
      write(unit_cond, '(I4, 2X, I8, *(2X, ES14.6))') &
          int(Ls(jL)), 4**int(Ls(jL)), (MyCond(idx, jL), idx = 1, norders)
    end do
    close(unit_err)
    close(unit_cond)
  end subroutine write_dat_files

  ! ----------------------------------------------------------------
  ! plot_with_gnuplot
  ! Spawn gnuplot to render MyErr (Fig 5.1 RIGHT) and MyCond (companion)
  ! into two PDF files alongside the .dat files.  Uses execute_command_line
  ! after writing a small .gp script.
  ! ----------------------------------------------------------------
  subroutine plot_with_gnuplot()
    integer(8) :: u, idx
    character(len=8) :: title
    open(newunit=u, file='paraboloid_refinement.gp', status='replace')
    write(u, '(A)') 'set terminal pdf size 6,4'
    write(u, '(A)') ''
    write(u, '(A)') '# --- MyErr: rel error vs num patches (loglog) ---'
    write(u, '(A)') "set output 'paraboloid_approx_error_vs_refinement.pdf'"
    write(u, '(A)') 'set logscale xy'
    write(u, '(A)') "set xlabel 'Num of Patches'"
    write(u, '(A)') "set ylabel 'Relative Error'"
    write(u, '(A)') 'set yrange [1e-15:1e1]'
    write(u, '(A)') 'set grid'
    write(u, '(A)') 'set key bottom left'
    write(u, '(A)', advance='no') "plot"
    do idx = 1, norders
      write(title, '(A,I0)') 'p=', int(orders(idx))
      write(u, '(A)', advance='no') " 'paraboloid_refinement_myerr.dat' using 2:"
      write(u, '(I0,A,A,A)', advance='no') idx+2, " with linespoints title '", trim(title), "'"
      if (idx < norders) write(u, '(A)', advance='no') ","
    end do
    write(u, '(A)') ''
    write(u, '(A)') ''
    write(u, '(A)') '# --- MyCond: condition number vs L (semilogy) ---'
    write(u, '(A)') "set output 'paraboloid_cond_num_vs_refinement.pdf'"
    write(u, '(A)') 'unset logscale'
    write(u, '(A)') 'set logscale y'
    write(u, '(A)') "set xlabel 'Refinement Level'"
    write(u, '(A)') "set ylabel 'Condition Number'"
    write(u, '(A)') 'set yrange [1e0:1e18]'
    write(u, '(A)') 'set xrange [-0.5:5.5]'
    write(u, '(A)') 'set grid'
    write(u, '(A)') 'set key top right'
    write(u, '(A)', advance='no') "plot"
    do idx = 1, norders
      write(title, '(A,I0)') 'p=', int(orders(idx))
      write(u, '(A)', advance='no') " 'paraboloid_refinement_mycond.dat' using 1:"
      write(u, '(I0,A,A,A)', advance='no') idx+2, " with linespoints title '", trim(title), "'"
      if (idx < norders) write(u, '(A)', advance='no') ","
    end do
    write(u, '(A)') ''
    close(u)
    call execute_command_line('gnuplot paraboloid_refinement.gp', wait=.true.)
    write(*,'(A)') ''
    write(*,'(A)') 'Wrote: paraboloid_approx_error_vs_refinement.pdf'
    write(*,'(A)') '       paraboloid_cond_num_vs_refinement.pdf'
  end subroutine plot_with_gnuplot

  ! ----------------------------------------------------------------
  ! svd_cond
  ! cond_2(A) = sigma_max(A) / sigma_min(A) via LAPACK dgesvd
  ! (Accelerate; integer(4) sigs for LP64 ABI).  A destroyed on exit.
  ! ----------------------------------------------------------------
  function svd_cond(n, A) result(condA)
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: A(n, n)
    real(r64)                 :: condA
    real(r64), allocatable :: S(:), U(:,:), VT(:,:), work(:)
    integer(4)             :: n4, lwork4, info4, ldu4, ldvt4

    n4     = int(n, 4)
    ldu4   = 1_4
    ldvt4  = 1_4
    lwork4 = max(10*n4, 1)
    allocate(S(n), U(1,1), VT(1,1), work(lwork4))
    call dgesvd('N', 'N', n4, n4, A, n4, S, U, ldu4, VT, ldvt4, work, lwork4, info4)
    if (info4 /= 0) then
      write(*,*) 'svd_cond: dgesvd failed info=', info4
      error stop
    end if
    condA = S(1) / S(n)
    deallocate(S, U, VT, work)
  end function svd_cond

  ! ----------------------------------------------------------------
  ! lu_solve
  ! In-place solve of M * x = rhs via LAPACK dgetrf + dgetrs (LP64).
  ! M and rhs are overwritten; rhs becomes the solution on exit.
  ! Fails loudly on any non-zero LAPACK info (no fallback).
  ! ----------------------------------------------------------------
  subroutine lu_solve(n, M, rhs)
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: M(n, n)
    real(r64),  intent(inout) :: rhs(n)
    integer(4)              :: n4, info4
    integer(4), allocatable :: ipiv(:)

    n4 = int(n, 4)
    allocate(ipiv(n4))
    call dgetrf(n4, n4, M, n4, ipiv, info4)
    if (info4 /= 0) then
      write(*,*) 'lu_solve: dgetrf failed info=', info4
      error stop
    end if
    call dgetrs('N', n4, 1_4, M, n4, ipiv, rhs, n4, info4)
    if (info4 /= 0) then
      write(*,*) 'lu_solve: dgetrs failed info=', info4
      error stop
    end if
    deallocate(ipiv)
  end subroutine lu_solve

end program test_paraboloid_refinement
