! test_paraboloid_cond.f90
!
! Fortran twin of test/harmonic_approx/test_paraboloid_cond.m: reproduces
! the cond(A) sweep from Figure 5.1 LEFT of Zhu & Jiang (lapquad5_cpam_rev1.tex,
! §5 "Study of the quaternion approximation").
!
! For each (p, H) in MyOrder x MyH:
!   - Vioreanu nodes on the simplex of order p-1 -> equilateral triangle (side 1) -> paraboloid
!   - Translate so the three vertices lie on the xy-plane
!   - Evaluate 3D Laplace solid spherical harmonics + gradient
!   - Take real parts of the (l, m), 1 <= m <= l <= p subset of Fx, Fy, Fz
!   - Assemble the 4x4 block matrix A from eq. (5.1) of the paper
!   - Compute cond(A) via LAPACK dgesvd (ratio of largest to smallest singular value)
!
! Prints a table of cond(A) values to stdout.
!
! No plotting (that stays in the .m twin). The Fortran path here is the
! gold reference: no MATLAB <-> mex boundary cast, full r64 throughout.

program test_paraboloid_cond
  use quatapproximation_mod, only: r64, simplex2equil_r64, paraboloid_uv_r64
  use harmonic_mod,          only: l3dtavecevalmat_r64
  use koorn_geom_mod,        only: get_vioreanu_nodes
  implicit none

  integer(8), parameter :: norders = 5
  integer(8), parameter :: nH      = 11
  integer(8), parameter :: orders(norders) = [4_8, 8_8, 12_8, 16_8, 20_8]
  real(r64),  parameter :: h_side = 1.0_r64

  real(r64) :: Hs(nH)
  real(r64) :: cond_table(norders, nH)

  integer(8) :: idx, j, p, n_p, nbasis_full, ier
  integer(8) :: ij, kk, tmpidx, tmpidx2, k
  real(r64)  :: H

  real(r64),    allocatable :: uvs_simplex(:,:), uvs_eq(:,:), sx(:,:)
  real(r64)                 :: uvs_vert_simplex(2,3), uvs_vert_eq(2,3), x_vert(3,3)
  complex(8),   allocatable :: F(:,:), Fx(:,:), Fy(:,:), Fz(:,:)
  real(r64),    allocatable :: F1(:,:), F2(:,:), F3(:,:), A(:,:)
  integer(8),   allocatable :: idxvec(:)

  ! Build the H sweep: 0, 1, 2, ..., 10.
  do j = 1, nH
    Hs(j) = real(j-1, r64)
  end do

  do idx = 1, norders
    p           = orders(idx)
    n_p         = p*(p+1)/2
    nbasis_full = (p+1)*(p+1)

    allocate(uvs_simplex(2, n_p), uvs_eq(2, n_p), sx(3, n_p))
    allocate(F (n_p, nbasis_full), Fx(n_p, nbasis_full), &
             Fy(n_p, nbasis_full), Fz(n_p, nbasis_full))
    allocate(F1(n_p, n_p), F2(n_p, n_p), F3(n_p, n_p))
    allocate(A(4*n_p, 4*n_p))
    allocate(idxvec(n_p))

    ! Vioreanu nodes on the simplex (depend only on p, not H).
    call get_vioreanu_nodes(p-1, n_p, uvs_simplex)
    call simplex2equil_r64(n_p, h_side, uvs_simplex, uvs_eq)

    ! Vertex triangle (depends only on p indirectly via h_side; same for all H).
    uvs_vert_simplex = reshape([0.0_r64, 0.0_r64, &
                                 1.0_r64, 0.0_r64, &
                                 0.0_r64, 1.0_r64], [2, 3])
    call simplex2equil_r64(3_8, h_side, uvs_vert_simplex, uvs_vert_eq)

    ! Index subset of the tensor-product basis: keep (l, m) with l > 0 and m > 0.
    idxvec = 0
    tmpidx = 0
    tmpidx2 = 0
    do ij = 0, p
      do kk = -ij, ij
        tmpidx2 = tmpidx2 + 1
        if (ij > 0 .and. kk > 0) then
          tmpidx = tmpidx + 1
          idxvec(tmpidx) = tmpidx2
        end if
      end do
    end do

    do j = 1, nH
      H = Hs(j)

      ! Patch points on the paraboloid, then z-shift so vertices on xy-plane.
      call paraboloid_uv_r64(n_p, H, uvs_eq(1,:), uvs_eq(2,:), sx)
      call paraboloid_uv_r64(3_8, H, uvs_vert_eq(1,:), uvs_vert_eq(2,:), x_vert)
      sx(3,:) = sx(3,:) - sum(x_vert(3,:)) / 3.0_r64

      ! Evaluate harmonic basis values + Cartesian gradient.
      F  = (0.0_r64, 0.0_r64)
      Fx = (0.0_r64, 0.0_r64)
      Fy = (0.0_r64, 0.0_r64)
      Fz = (0.0_r64, 0.0_r64)
      ier = 0
      call l3dtavecevalmat_r64(sx, n_p, p, F, Fx, Fy, Fz, ier)

      ! Take real parts of the (l > 0, m > 0) subset.
      do k = 1, n_p
        F1(:, k) = real(Fx(:, idxvec(k)), r64)
        F2(:, k) = real(Fy(:, idxvec(k)), r64)
        F3(:, k) = real(Fz(:, idxvec(k)), r64)
      end do

      ! Assemble the 4x4 block matrix A (eq. 5.1 of the paper):
      !   A = [ 0   -F1  -F2  -F3
      !        F1   0   -F3   F2
      !        F2   F3   0   -F1
      !        F3  -F2   F1   0  ]
      A = 0.0_r64
      A(1:n_p,         n_p+1:2*n_p)   = -F1
      A(1:n_p,         2*n_p+1:3*n_p) = -F2
      A(1:n_p,         3*n_p+1:4*n_p) = -F3
      A(n_p+1:2*n_p,   1:n_p)         =  F1
      A(n_p+1:2*n_p,   2*n_p+1:3*n_p) = -F3
      A(n_p+1:2*n_p,   3*n_p+1:4*n_p) =  F2
      A(2*n_p+1:3*n_p, 1:n_p)         =  F2
      A(2*n_p+1:3*n_p, n_p+1:2*n_p)   =  F3
      A(2*n_p+1:3*n_p, 3*n_p+1:4*n_p) = -F1
      A(3*n_p+1:4*n_p, 1:n_p)         =  F3
      A(3*n_p+1:4*n_p, n_p+1:2*n_p)   = -F2
      A(3*n_p+1:4*n_p, 2*n_p+1:3*n_p) =  F1

      cond_table(idx, j) = svd_cond(4*n_p, A)
    end do

    deallocate(uvs_simplex, uvs_eq, sx)
    deallocate(F, Fx, Fy, Fz)
    deallocate(F1, F2, F3, A)
    deallocate(idxvec)
  end do

  ! Print the table: rows = polynomial orders, columns = H values.
  write(*,'(A)') '=== test_paraboloid_cond (Fig 5.1 LEFT reproducer) ==='
  write(*,'(A)') 'h_side = 1.0; cond(A) vs H, one row per polynomial order p:'
  write(*,'(A)', advance='no') '   p \\ H'
  do j = 1, nH
    write(*,'(2X,F8.1)', advance='no') Hs(j)
  end do
  write(*,'(A)') ''
  do idx = 1, norders
    write(*,'(I7)', advance='no') orders(idx)
    do j = 1, nH
      write(*,'(2X,ES8.2)', advance='no') cond_table(idx, j)
    end do
    write(*,'(A)') ''
  end do

contains

  ! ----------------------------------------------------------------
  ! svd_cond
  ! cond_2(A) = sigma_max(A) / sigma_min(A) via LAPACK dgesvd.
  ! A is destroyed on exit (caller passes a workspace copy if needed).
  !
  ! Linked against Apple's Accelerate framework on macOS (LP64 LAPACK
  ! interface: 32-bit integer args).  Even though the rest of the
  ! package compiles with -fdefault-integer-8 to match the mwrap -i8
  ! boundary, the LAPACK call here uses explicit integer(4) so it
  ! matches Accelerate's LP64 ABI; condA itself is r64.
  ! ----------------------------------------------------------------
  function svd_cond(n, A) result(condA)
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: A(n, n)
    real(r64)                 :: condA

    real(r64),  allocatable :: S(:), U(:,:), VT(:,:), work(:)
    integer(4)              :: n4, lwork4, info4, ldu4, ldvt4

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

end program test_paraboloid_cond
