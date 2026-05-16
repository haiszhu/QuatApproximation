subroutine qakg_circumcircle_transform_3d_mex(r_vert, R, c, alpha)
  use koorn_geom_mod, only: circumcircle_transform_3d
  implicit none
  real(8), intent(in)    :: r_vert(3,3)
  real(8), intent(inout) :: R(3,3), c(3), alpha

  call circumcircle_transform_3d(r_vert, R, c, alpha)
end subroutine qakg_circumcircle_transform_3d_mex

! ------------------------------------------------------------------
! qakg_circumcircle_transform_3d_r128_mex
! Per-triangle r128 circumcircle transform. Trailing args (ktri, flag):
!   flag = 0 : cast r_vert from r64 args.
!   flag = 1 : load r_vert from lqem_mesh_r128.h5:/tri_vert(:,:,ktri).
! Writes qakg_circumcircle_r128.h5 with /R, /c, /alpha.
! ------------------------------------------------------------------
subroutine qakg_circumcircle_transform_3d_r128_mex(r_vert, R, c, alpha, ktri, flag)
  use koorn_geom_mod, only: qakg_circ_r128 => circumcircle_transform_3d_r128
  use quatapproximation_mod, only: r128
  use iso_c_binding, only: c_char, c_float128, c_int, c_int64_t, c_null_char
  implicit none
  interface
    function hdf5_write_real128_array(file, name, rank, dims, vals, append_flag) bind(C)
      import c_char, c_float128, c_int, c_int64_t
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_int),   value            :: rank
      integer(c_int64_t), intent(in)     :: dims(*)
      real(c_float128), intent(in)       :: vals(*)
      integer(c_int),   value            :: append_flag
      integer(c_int)                     :: hdf5_write_real128_array
    end function hdf5_write_real128_array
    function hdf5_read_real128_array(file, name, rank, dims, vals) bind(C)
      import c_char, c_float128, c_int, c_int64_t
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_int),   value            :: rank
      integer(c_int64_t), intent(in)     :: dims(*)
      real(c_float128), intent(inout)    :: vals(*)
      integer(c_int)                     :: hdf5_read_real128_array
    end function hdf5_read_real128_array
    function hdf5_get_dims_r128(file, name, rank, dims) bind(C)
      import c_char, c_int, c_int64_t
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_int),     value          :: rank
      integer(c_int64_t), intent(inout)  :: dims(*)
      integer(c_int)                     :: hdf5_get_dims_r128
    end function hdf5_get_dims_r128
  end interface
  real(8),    intent(in)    :: r_vert(3,3)
  real(8),    intent(inout) :: R(3,3), c(3), alpha
  integer(8), intent(in)    :: ktri, flag

  real(r128) :: r_vert_r128(3,3), R_r128(3,3), c_r128(3), alpha_r128
  real(r128), allocatable :: tri_vert_full(:,:,:)
  integer(c_int64_t) :: dims3(3), dims2(2), dims1(1), qdims(3)
  integer(c_int)     :: h5_ok
  integer(8)         :: ntri_full

  if (flag == 0_8) then
    r_vert_r128 = real(r_vert, r128)
  else
    h5_ok = hdf5_get_dims_r128('lqem_mesh_r128.h5'//c_null_char, &
        '/tri_vert'//c_null_char, 3_c_int, qdims)
    if (h5_ok == 0) error stop 'qakg_circumcircle_transform_3d_r128_mex: failed to query /tri_vert dims'
    if (qdims(1) /= 3_c_int64_t .or. qdims(2) /= 3_c_int64_t) then
      error stop 'qakg_circumcircle_transform_3d_r128_mex: /tri_vert dims (3,3,*) mismatch'
    end if
    ntri_full = int(qdims(3), 8)
    if (ktri < 1_8 .or. ktri > ntri_full) then
      error stop 'qakg_circumcircle_transform_3d_r128_mex: ktri out of range'
    end if
    allocate(tri_vert_full(3, 3, ntri_full))
    dims3(1) = 3_c_int64_t
    dims3(2) = 3_c_int64_t
    dims3(3) = int(ntri_full, c_int64_t)
    h5_ok = hdf5_read_real128_array('lqem_mesh_r128.h5'//c_null_char, &
        '/tri_vert'//c_null_char, 3_c_int, dims3, tri_vert_full)
    if (h5_ok == 0) error stop 'qakg_circumcircle_transform_3d_r128_mex: failed to read /tri_vert'
    r_vert_r128 = tri_vert_full(:,:,ktri)
    deallocate(tri_vert_full)
  end if

  R_r128     = 0.0_r128
  c_r128     = 0.0_r128
  alpha_r128 = 0.0_r128
  call qakg_circ_r128(r_vert_r128, R_r128, c_r128, alpha_r128)

  dims2(1) = 3_c_int64_t
  dims2(2) = 3_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_circumcircle_r128.h5'//c_null_char, &
      '/R'//c_null_char, 2_c_int, dims2, R_r128, 0_c_int)
  dims1(1) = 3_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_circumcircle_r128.h5'//c_null_char, &
      '/c'//c_null_char, 1_c_int, dims1, c_r128, 1_c_int)
  dims1(1) = 1_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_circumcircle_r128.h5'//c_null_char, &
      '/alpha'//c_null_char, 1_c_int, dims1, [alpha_r128], 1_c_int)

  R     = real(R_r128,     8)
  c     = real(c_r128,     8)
  alpha = real(alpha_r128, 8)
end subroutine qakg_circumcircle_transform_3d_r128_mex

! ------------------------------------------------------------------
! qakg_setup_target_r128_mex
! Per-triangle r128 setup bundle. Produces every quantity needed by
! the kernel-eval and compress mex (sxbd, sxpbd, stangbd, sspbd,
! txnew, snxnew, qhat, kdata) plus circumcircle outputs (R, c, alpha).
! Trailing args (ktri, flag):
!   flag = 0 : cast r_vert, tx, snx, sxbd_in, Dgl from r64 args.
!   flag = 1 : load r_vert, sxbd_in from lqem_mesh_r128.h5 at slice ktri,
!              load Dgl from lqgauss_r128.h5. tx and snx always from MATLAB.
! Writes one HDF5 file: qakg_setup_target_r128.h5 with datasets
!   /R, /c, /alpha, /sxbd, /sxpbd, /stangbd, /sspbd,
!   /txnew, /snxnew, /qhat, /kdata.
! ------------------------------------------------------------------
subroutine qakg_setup_target_r128_mex(m, n, nbd, sbdnp, nquad,                  &
                                      r_vert, tx, snx, sxbd_in, Dgl,            &
                                      R, c, alpha,                              &
                                      sxbd, sxpbd, stangbd, sspbd,              &
                                      txnew, snxnew, qhat, kdata,               &
                                      ktri, flag)
  use koorn_geom_mod, only: qakg_setup_r128 => lqkg_setup_target_r128
  use quatapproximation_mod, only: r128
  use iso_c_binding, only: c_char, c_float128, c_int, c_int64_t, c_null_char
  implicit none
  interface
    function hdf5_write_real128_array(file, name, rank, dims, vals, append_flag) bind(C)
      import c_char, c_float128, c_int, c_int64_t
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_int),   value            :: rank
      integer(c_int64_t), intent(in)     :: dims(*)
      real(c_float128), intent(in)       :: vals(*)
      integer(c_int),   value            :: append_flag
      integer(c_int)                     :: hdf5_write_real128_array
    end function hdf5_write_real128_array
    function hdf5_read_real128_array(file, name, rank, dims, vals) bind(C)
      import c_char, c_float128, c_int, c_int64_t
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_int),   value            :: rank
      integer(c_int64_t), intent(in)     :: dims(*)
      real(c_float128), intent(inout)    :: vals(*)
      integer(c_int)                     :: hdf5_read_real128_array
    end function hdf5_read_real128_array
    function hdf5_get_dims_r128(file, name, rank, dims) bind(C)
      import c_char, c_int, c_int64_t
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_int),     value          :: rank
      integer(c_int64_t), intent(inout)  :: dims(*)
      integer(c_int)                     :: hdf5_get_dims_r128
    end function hdf5_get_dims_r128
  end interface

  integer(8), intent(in)    :: m, n, nbd, sbdnp, nquad
  real(8),    intent(in)    :: r_vert(3,3)
  real(8),    intent(in)    :: tx(3,m), snx(3,n)
  real(8),    intent(in)    :: sxbd_in(3,nbd)
  real(8),    intent(in)    :: Dgl(nquad,nquad)
  real(8),    intent(inout) :: R(3,3), c(3), alpha
  real(8),    intent(inout) :: sxbd(3,nbd), sxpbd(3,nbd)
  real(8),    intent(inout) :: stangbd(3,nbd), sspbd(nbd)
  real(8),    intent(inout) :: txnew(3,m), snxnew(3,n)
  real(8),    intent(inout) :: qhat(3), kdata(3,m)
  integer(8), intent(in)    :: ktri, flag

  real(r128) :: r_vert_r128(3,3), tx_r128(3,m), snx_r128(3,n)
  real(r128) :: sxbd_in_r128(3,nbd), Dgl_r128(nquad,nquad)
  real(r128) :: R_r128(3,3), c_r128(3), alpha_r128
  real(r128) :: sxbd_r128(3,nbd), sxpbd_r128(3,nbd)
  real(r128) :: stangbd_r128(3,nbd), sspbd_r128(nbd)
  real(r128) :: txnew_r128(3,m), snxnew_r128(3,n)
  real(r128) :: qhat_r128(3), kdata_r128(3,m)

  real(r128), allocatable :: xbd_full(:,:,:), tri_vert_full(:,:,:)
  integer(c_int64_t) :: dims3(3), dims2(2), dims1(1), qdims(3)
  integer(c_int)     :: h5_ok
  integer(8)         :: ntri_full

  ! tx, snx always from MATLAB args.
  tx_r128  = real(tx,  r128)
  snx_r128 = real(snx, r128)

  if (flag == 0_8) then
    r_vert_r128  = real(r_vert,  r128)
    sxbd_in_r128 = real(sxbd_in, r128)
    Dgl_r128     = real(Dgl,     r128)
  else
    ! Mesh slice from lqem_mesh_r128.h5
    h5_ok = hdf5_get_dims_r128('lqem_mesh_r128.h5'//c_null_char, &
        '/tri_vert'//c_null_char, 3_c_int, qdims)
    if (h5_ok == 0) error stop 'qakg_setup_target_r128_mex: failed to query /tri_vert dims'
    if (qdims(1) /= 3_c_int64_t .or. qdims(2) /= 3_c_int64_t) then
      error stop 'qakg_setup_target_r128_mex: /tri_vert dims (3,3,*) mismatch'
    end if
    ntri_full = int(qdims(3), 8)
    if (ktri < 1_8 .or. ktri > ntri_full) then
      error stop 'qakg_setup_target_r128_mex: ktri out of range'
    end if
    allocate(tri_vert_full(3, 3, ntri_full))
    dims3(1) = 3_c_int64_t;  dims3(2) = 3_c_int64_t;  dims3(3) = int(ntri_full, c_int64_t)
    h5_ok = hdf5_read_real128_array('lqem_mesh_r128.h5'//c_null_char, &
        '/tri_vert'//c_null_char, 3_c_int, dims3, tri_vert_full)
    if (h5_ok == 0) error stop 'qakg_setup_target_r128_mex: failed to read /tri_vert'
    r_vert_r128 = tri_vert_full(:,:,ktri)
    deallocate(tri_vert_full)

    allocate(xbd_full(3, nbd, ntri_full))
    dims3(1) = 3_c_int64_t;  dims3(2) = int(nbd, c_int64_t);  dims3(3) = int(ntri_full, c_int64_t)
    h5_ok = hdf5_read_real128_array('lqem_mesh_r128.h5'//c_null_char, &
        '/xbd'//c_null_char, 3_c_int, dims3, xbd_full)
    if (h5_ok == 0) error stop 'qakg_setup_target_r128_mex: failed to read /xbd'
    sxbd_in_r128 = xbd_full(:,:,ktri)
    deallocate(xbd_full)

    ! Dgl from lqgauss_r128.h5
    dims2(1) = int(nquad, c_int64_t);  dims2(2) = int(nquad, c_int64_t)
    h5_ok = hdf5_read_real128_array('lqgauss_r128.h5'//c_null_char, &
        '/Dgl'//c_null_char, 2_c_int, dims2, Dgl_r128)
    if (h5_ok == 0) error stop 'qakg_setup_target_r128_mex: failed to read /Dgl'
  end if

  call qakg_setup_r128(m, n, nbd, sbdnp, nquad,                  &
                       r_vert_r128, tx_r128, snx_r128,           &
                       sxbd_in_r128, Dgl_r128,                   &
                       R_r128, c_r128, alpha_r128,               &
                       sxbd_r128, sxpbd_r128, stangbd_r128, sspbd_r128,  &
                       txnew_r128, snxnew_r128, qhat_r128, kdata_r128)

  ! Persist all r128 outputs to a single HDF5 file.
  dims2(1) = 3_c_int64_t;  dims2(2) = 3_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/R'//c_null_char, 2_c_int, dims2, R_r128, 0_c_int)

  dims1(1) = 3_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/c'//c_null_char, 1_c_int, dims1, c_r128, 1_c_int)

  dims1(1) = 1_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/alpha'//c_null_char, 1_c_int, dims1, [alpha_r128], 1_c_int)

  dims2(1) = 3_c_int64_t;  dims2(2) = int(nbd, c_int64_t)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/sxbd'//c_null_char,    2_c_int, dims2, sxbd_r128,    1_c_int)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/sxpbd'//c_null_char,   2_c_int, dims2, sxpbd_r128,   1_c_int)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/stangbd'//c_null_char, 2_c_int, dims2, stangbd_r128, 1_c_int)

  dims1(1) = int(nbd, c_int64_t)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/sspbd'//c_null_char, 1_c_int, dims1, sspbd_r128, 1_c_int)

  dims2(1) = 3_c_int64_t;  dims2(2) = int(m, c_int64_t)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/txnew'//c_null_char, 2_c_int, dims2, txnew_r128, 1_c_int)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/kdata'//c_null_char, 2_c_int, dims2, kdata_r128, 1_c_int)

  dims2(1) = 3_c_int64_t;  dims2(2) = int(n, c_int64_t)
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/snxnew'//c_null_char, 2_c_int, dims2, snxnew_r128, 1_c_int)

  dims1(1) = 3_c_int64_t
  h5_ok = hdf5_write_real128_array('qakg_setup_target_r128.h5'//c_null_char, &
      '/qhat'//c_null_char, 1_c_int, dims1, qhat_r128, 1_c_int)

  ! r64 cast back for MATLAB visibility.
  R       = real(R_r128,       8)
  c       = real(c_r128,       8)
  alpha   = real(alpha_r128,   8)
  sxbd    = real(sxbd_r128,    8)
  sxpbd   = real(sxpbd_r128,   8)
  stangbd = real(stangbd_r128, 8)
  sspbd   = real(sspbd_r128,   8)
  txnew   = real(txnew_r128,   8)
  snxnew  = real(snxnew_r128,  8)
  qhat    = real(qhat_r128,    8)
  kdata   = real(kdata_r128,   8)

end subroutine qakg_setup_target_r128_mex


! ------------------------------------------------------------------
! qakg_get_vioreanu_nodes_mex
! Thin r64 wrapper for koorn_geom_mod::get_vioreanu_nodes. Loads
! precomputed Vioreanu-Rokhlin nodes uvs(2, npols) of order norder,
! where npols = (norder+1)*(norder+2)/2.
! ------------------------------------------------------------------
subroutine qakg_get_vioreanu_nodes_mex(norder, npols, uvs)
  use koorn_geom_mod, only: get_vioreanu_nodes
  implicit none
  integer(8), intent(in)    :: norder, npols
  real(8),    intent(inout) :: uvs(2, npols)
  call get_vioreanu_nodes(norder, npols, uvs)
end subroutine qakg_get_vioreanu_nodes_mex

! ------------------------------------------------------------------
! qakg_get_vioreanu_wts_mex
! Sibling of qakg_get_vioreanu_nodes_mex returning the matching
! Vioreanu-Rokhlin quadrature weights wts(npols).
! ------------------------------------------------------------------
subroutine qakg_get_vioreanu_wts_mex(norder, npols, wts)
  use koorn_geom_mod, only: get_vioreanu_wts
  implicit none
  integer(8), intent(in)    :: norder, npols
  real(8),    intent(inout) :: wts(npols)
  call get_vioreanu_wts(norder, npols, wts)
end subroutine qakg_get_vioreanu_wts_mex

! ------------------------------------------------------------------
! qakg_koorn_vals2coefs_coefs2vals_mex
! Thin r64 wrapper for koorn_geom_mod::koorn_vals2coefs_coefs2vals.
! Returns the Koornwinder vals->coefs matrix (umatr) and the
! coefs->vals Vandermonde (vmatr) at the VR nodes of order korder.
! ------------------------------------------------------------------
subroutine qakg_koorn_vals2coefs_coefs2vals_mex(korder, kpols, umatr, vmatr)
  use koorn_geom_mod, only: koorn_vals2coefs_coefs2vals
  implicit none
  integer(8), intent(in)    :: korder, kpols
  real(8),    intent(inout) :: umatr(kpols, kpols), vmatr(kpols, kpols)

  call koorn_vals2coefs_coefs2vals(korder, kpols, umatr, vmatr)
end subroutine qakg_koorn_vals2coefs_coefs2vals_mex

! ------------------------------------------------------------------
! qakg_line3quadr_3dline_mex
! Thin r64 wrapper for koorn_geom_mod::line3quadr_3dline.  Builds GL
! quadrature nodes/weights/tangents/speeds along the three edges of a
! curved triangle parameterised by Koornwinder coeffs (via umatr).
! ------------------------------------------------------------------
subroutine qakg_line3quadr_3dline_mex(x_uvs, korder, kpols, umatr, nquad,    &
                                       tgl, wgl, Dgl, sbdnp, tpan, nbd,      &
                                       sxbd, swbd, stangbd, sspbd, r_vert)
  use koorn_geom_mod, only: line3quadr_3dline
  implicit none
  integer(8), intent(in)    :: korder, kpols, nquad, sbdnp, nbd
  real(8),    intent(in)    :: x_uvs(3, kpols), umatr(kpols, kpols)
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad, nquad)
  real(8),    intent(in)    :: tpan(sbdnp + 1)
  real(8),    intent(inout) :: sxbd(3, nbd), swbd(nbd)
  real(8),    intent(inout) :: stangbd(3, nbd), sspbd(nbd)
  real(8),    intent(inout) :: r_vert(3, 3)

  call line3quadr_3dline(x_uvs, korder, kpols, umatr, nquad,    &
                          tgl, wgl, Dgl, sbdnp, tpan, nbd,      &
                          sxbd, swbd, stangbd, sspbd, r_vert)
end subroutine qakg_line3quadr_3dline_mex
