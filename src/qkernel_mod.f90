module qkernel_mod
  ! q^{(n,m)}_i kernel-coefficient routines for the Laplace DLP/SLP close-eval
  ! pipeline. Mirrors qotential/utils/q[0-3]kl.m at r64.
  !
  ! lptype_id: 0='d' (DLP/SLPn, 4 slots), 1='s' (SLP, 5 slots),
  !            2='T' (DLPn, 11 slots) -- NOT YET IMPLEMENTED; stubbed.
  use quatapproximation_mod, only: r64, r128
  implicit none
  private
  public :: qak_q0kl_r64, qak_q1kl_r64, qak_q2kl_r64, qak_q3kl_r64
  public :: qak_q0kl_r128, qak_q1kl_r128, qak_q2kl_r128, qak_q3kl_r128
  public :: qak_qnm_i_r64, qak_qnm_i_r128, qak_qnm_i_nslots

  integer(8), parameter, public :: QAK_LPTYPE_D = 0_8
  integer(8), parameter, public :: QAK_LPTYPE_S = 1_8
  integer(8), parameter, public :: QAK_LPTYPE_T = 2_8

contains

  pure function qak_qnm_i_nslots(lptype_id) result(nslots)
    integer(8), intent(in) :: lptype_id
    integer(8)             :: nslots
    select case (lptype_id)
      case (QAK_LPTYPE_D); nslots = 4_8
      case (QAK_LPTYPE_S); nslots = 5_8
      case (QAK_LPTYPE_T); nslots = 11_8
      case default;        nslots = 0_8
    end select
  end function qak_qnm_i_nslots

  subroutine qak_qnm_i_r64(idx, N, ncoeff, sx, lptype_id, &
                           F, F1, F2, F3, gradxyz, &
                           q_i, q_j, q_k)
    integer(8), intent(in)  :: idx, N, ncoeff, lptype_id
    real(r64),  intent(in)  :: sx(3, N)
    real(r64),  intent(in)  :: F (N, ncoeff)
    real(r64),  intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r64),  intent(in)  :: gradxyz(3, N)
    real(r64),  intent(out) :: q_i(N, ncoeff, 5)
    real(r64),  intent(out) :: q_j(N, ncoeff, 5)
    real(r64),  intent(out) :: q_k(N, ncoeff, 5)
    select case (idx)
      case (0_8); call qak_q0kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case (1_8); call qak_q1kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case (2_8); call qak_q2kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case (3_8); call qak_q3kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case default
        write(*,*) 'qak_qnm_i_r64: idx must be 0..3; got ', idx
        error stop
    end select
  end subroutine qak_qnm_i_r64

  subroutine qak_q0kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r64),  intent(in)  :: sx(3, N)
    real(r64),  intent(in)  :: F (N, ncoeff)
    real(r64),  intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r64),  intent(in)  :: gradxyz(3, N)
    real(r64),  intent(out) :: q_i(N, ncoeff, 5)
    real(r64),  intent(out) :: q_j(N, ncoeff, 5)
    real(r64),  intent(out) :: q_k(N, ncoeff, 5)
    integer(8) :: n_, c
    q_i = 0.0_r64;  q_j = 0.0_r64;  q_k = 0.0_r64
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) = -sx(2, n_)*F3(n_, c) + sx(3, n_)*F2(n_, c)
          q_i(n_, c, 3) =  F3(n_, c)
          q_i(n_, c, 4) = -F2(n_, c)
          q_j(n_, c, 1) = -sx(3, n_)*F1(n_, c) + sx(1, n_)*F3(n_, c)
          q_j(n_, c, 2) = -F3(n_, c)
          q_j(n_, c, 4) =  F1(n_, c)
          q_k(n_, c, 1) = -sx(1, n_)*F2(n_, c) + sx(2, n_)*F1(n_, c)
          q_k(n_, c, 2) =  F2(n_, c)
          q_k(n_, c, 3) = -F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) =  sx(1, n_)*F(n_, c)
          q_i(n_, c, 2) = -F(n_, c)
          q_i(n_, c, 5) =  F1(n_, c)
          q_j(n_, c, 1) =  sx(2, n_)*F(n_, c)
          q_j(n_, c, 3) = -F(n_, c)
          q_j(n_, c, 5) =  F2(n_, c)
          q_k(n_, c, 1) =  sx(3, n_)*F(n_, c)
          q_k(n_, c, 4) = -F(n_, c)
          q_k(n_, c, 5) =  F3(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q0kl_r64: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q0kl_r64: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q0kl_r64

  subroutine qak_q1kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r64),  intent(in)  :: sx(3, N)
    real(r64),  intent(in)  :: F (N, ncoeff)
    real(r64),  intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r64),  intent(in)  :: gradxyz(3, N)
    real(r64),  intent(out) :: q_i(N, ncoeff, 5)
    real(r64),  intent(out) :: q_j(N, ncoeff, 5)
    real(r64),  intent(out) :: q_k(N, ncoeff, 5)
    real(r64)  :: mr_dot_f
    integer(8) :: n_, c
    q_i = 0.0_r64;  q_j = 0.0_r64;  q_k = 0.0_r64
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          mr_dot_f = -sx(1, n_)*F1(n_, c) - sx(2, n_)*F2(n_, c) - sx(3, n_)*F3(n_, c)
          q_i(n_, c, 1) =  mr_dot_f + 2.0_r64*sx(1, n_)*F1(n_, c)
          q_i(n_, c, 2) = -F1(n_, c)
          q_i(n_, c, 3) =  F2(n_, c)
          q_i(n_, c, 4) =  F3(n_, c)
          q_j(n_, c, 1) =  sx(1, n_)*F2(n_, c) + sx(2, n_)*F1(n_, c)
          q_j(n_, c, 2) = -F2(n_, c)
          q_j(n_, c, 3) = -F1(n_, c)
          q_k(n_, c, 1) =  sx(1, n_)*F3(n_, c) + sx(3, n_)*F1(n_, c)
          q_k(n_, c, 2) = -F3(n_, c)
          q_k(n_, c, 4) = -F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          ! q_i is identically zero
          q_j(n_, c, 1) =  sx(3, n_)*F(n_, c)
          q_j(n_, c, 4) = -F(n_, c)
          q_j(n_, c, 5) = -F3(n_, c)
          q_k(n_, c, 1) = -sx(2, n_)*F(n_, c)
          q_k(n_, c, 3) =  F(n_, c)
          q_k(n_, c, 5) =  F2(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q1kl_r64: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q1kl_r64: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q1kl_r64

  subroutine qak_q2kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r64),  intent(in)  :: sx(3, N)
    real(r64),  intent(in)  :: F (N, ncoeff)
    real(r64),  intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r64),  intent(in)  :: gradxyz(3, N)
    real(r64),  intent(out) :: q_i(N, ncoeff, 5)
    real(r64),  intent(out) :: q_j(N, ncoeff, 5)
    real(r64),  intent(out) :: q_k(N, ncoeff, 5)
    real(r64)  :: mr_dot_f
    integer(8) :: n_, c
    q_i = 0.0_r64;  q_j = 0.0_r64;  q_k = 0.0_r64
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          mr_dot_f = -sx(1, n_)*F1(n_, c) - sx(2, n_)*F2(n_, c) - sx(3, n_)*F3(n_, c)
          q_i(n_, c, 1) =  sx(2, n_)*F1(n_, c) + sx(1, n_)*F2(n_, c)
          q_i(n_, c, 2) = -F2(n_, c)
          q_i(n_, c, 3) = -F1(n_, c)
          q_j(n_, c, 1) =  mr_dot_f + 2.0_r64*sx(2, n_)*F2(n_, c)
          q_j(n_, c, 2) =  F1(n_, c)
          q_j(n_, c, 3) = -F2(n_, c)
          q_j(n_, c, 4) =  F3(n_, c)
          q_k(n_, c, 1) =  sx(2, n_)*F3(n_, c) + sx(3, n_)*F2(n_, c)
          q_k(n_, c, 3) = -F3(n_, c)
          q_k(n_, c, 4) = -F2(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) = -sx(3, n_)*F(n_, c)
          q_i(n_, c, 4) =  F(n_, c)
          q_i(n_, c, 5) =  F3(n_, c)
          ! q_j is identically zero
          q_k(n_, c, 1) =  sx(1, n_)*F(n_, c)
          q_k(n_, c, 2) = -F(n_, c)
          q_k(n_, c, 5) = -F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q2kl_r64: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q2kl_r64: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q2kl_r64

  subroutine qak_q3kl_r64(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r64),  intent(in)  :: sx(3, N)
    real(r64),  intent(in)  :: F (N, ncoeff)
    real(r64),  intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r64),  intent(in)  :: gradxyz(3, N)
    real(r64),  intent(out) :: q_i(N, ncoeff, 5)
    real(r64),  intent(out) :: q_j(N, ncoeff, 5)
    real(r64),  intent(out) :: q_k(N, ncoeff, 5)
    real(r64)  :: mr_dot_f
    integer(8) :: n_, c
    q_i = 0.0_r64;  q_j = 0.0_r64;  q_k = 0.0_r64
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          mr_dot_f = -sx(1, n_)*F1(n_, c) - sx(2, n_)*F2(n_, c) - sx(3, n_)*F3(n_, c)
          q_i(n_, c, 1) =  sx(3, n_)*F1(n_, c) + sx(1, n_)*F3(n_, c)
          q_i(n_, c, 2) = -F3(n_, c)
          q_i(n_, c, 4) = -F1(n_, c)
          q_j(n_, c, 1) =  sx(3, n_)*F2(n_, c) + sx(2, n_)*F3(n_, c)
          q_j(n_, c, 3) = -F3(n_, c)
          q_j(n_, c, 4) = -F2(n_, c)
          q_k(n_, c, 1) =  mr_dot_f + 2.0_r64*sx(3, n_)*F3(n_, c)
          q_k(n_, c, 2) =  F1(n_, c)
          q_k(n_, c, 3) =  F2(n_, c)
          q_k(n_, c, 4) = -F3(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) =  sx(2, n_)*F(n_, c)
          q_i(n_, c, 3) = -F(n_, c)
          q_i(n_, c, 5) = -F2(n_, c)
          q_j(n_, c, 1) = -sx(1, n_)*F(n_, c)
          q_j(n_, c, 2) =  F(n_, c)
          q_j(n_, c, 5) =  F1(n_, c)
          ! q_k is identically zero
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q3kl_r64: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q3kl_r64: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q3kl_r64

  ! ================================================================
  ! r128 siblings -- bit-for-bit mirrors of the r64 routines above with
  ! the kind switched.  No mex wrappers (deferred); intended for the
  ! future Lap3dDLP_closepanel_r128 orchestration.
  ! 'T' (DLPn) is stubbed in both precisions; not yet implemented.
  ! ================================================================

  subroutine qak_qnm_i_r128(idx, N, ncoeff, sx, lptype_id, &
                            F, F1, F2, F3, gradxyz, &
                            q_i, q_j, q_k)
    integer(8), intent(in)  :: idx, N, ncoeff, lptype_id
    real(r128), intent(in)  :: sx(3, N)
    real(r128), intent(in)  :: F (N, ncoeff)
    real(r128), intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r128), intent(in)  :: gradxyz(3, N)
    real(r128), intent(out) :: q_i(N, ncoeff, 5)
    real(r128), intent(out) :: q_j(N, ncoeff, 5)
    real(r128), intent(out) :: q_k(N, ncoeff, 5)
    select case (idx)
      case (0_8); call qak_q0kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case (1_8); call qak_q1kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case (2_8); call qak_q2kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case (3_8); call qak_q3kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
      case default
        write(*,*) 'qak_qnm_i_r128: idx must be 0..3; got ', idx
        error stop
    end select
  end subroutine qak_qnm_i_r128

  subroutine qak_q0kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r128), intent(in)  :: sx(3, N)
    real(r128), intent(in)  :: F (N, ncoeff)
    real(r128), intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r128), intent(in)  :: gradxyz(3, N)
    real(r128), intent(out) :: q_i(N, ncoeff, 5)
    real(r128), intent(out) :: q_j(N, ncoeff, 5)
    real(r128), intent(out) :: q_k(N, ncoeff, 5)
    integer(8) :: n_, c
    q_i = 0.0_r128;  q_j = 0.0_r128;  q_k = 0.0_r128
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) = -sx(2, n_)*F3(n_, c) + sx(3, n_)*F2(n_, c)
          q_i(n_, c, 3) =  F3(n_, c)
          q_i(n_, c, 4) = -F2(n_, c)
          q_j(n_, c, 1) = -sx(3, n_)*F1(n_, c) + sx(1, n_)*F3(n_, c)
          q_j(n_, c, 2) = -F3(n_, c)
          q_j(n_, c, 4) =  F1(n_, c)
          q_k(n_, c, 1) = -sx(1, n_)*F2(n_, c) + sx(2, n_)*F1(n_, c)
          q_k(n_, c, 2) =  F2(n_, c)
          q_k(n_, c, 3) = -F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) =  sx(1, n_)*F(n_, c)
          q_i(n_, c, 2) = -F(n_, c)
          q_i(n_, c, 5) =  F1(n_, c)
          q_j(n_, c, 1) =  sx(2, n_)*F(n_, c)
          q_j(n_, c, 3) = -F(n_, c)
          q_j(n_, c, 5) =  F2(n_, c)
          q_k(n_, c, 1) =  sx(3, n_)*F(n_, c)
          q_k(n_, c, 4) = -F(n_, c)
          q_k(n_, c, 5) =  F3(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q0kl_r128: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q0kl_r128: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q0kl_r128

  subroutine qak_q1kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r128), intent(in)  :: sx(3, N)
    real(r128), intent(in)  :: F (N, ncoeff)
    real(r128), intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r128), intent(in)  :: gradxyz(3, N)
    real(r128), intent(out) :: q_i(N, ncoeff, 5)
    real(r128), intent(out) :: q_j(N, ncoeff, 5)
    real(r128), intent(out) :: q_k(N, ncoeff, 5)
    real(r128) :: mr_dot_f
    integer(8) :: n_, c
    q_i = 0.0_r128;  q_j = 0.0_r128;  q_k = 0.0_r128
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          mr_dot_f = -sx(1, n_)*F1(n_, c) - sx(2, n_)*F2(n_, c) - sx(3, n_)*F3(n_, c)
          q_i(n_, c, 1) =  mr_dot_f + 2.0_r128*sx(1, n_)*F1(n_, c)
          q_i(n_, c, 2) = -F1(n_, c)
          q_i(n_, c, 3) =  F2(n_, c)
          q_i(n_, c, 4) =  F3(n_, c)
          q_j(n_, c, 1) =  sx(1, n_)*F2(n_, c) + sx(2, n_)*F1(n_, c)
          q_j(n_, c, 2) = -F2(n_, c)
          q_j(n_, c, 3) = -F1(n_, c)
          q_k(n_, c, 1) =  sx(1, n_)*F3(n_, c) + sx(3, n_)*F1(n_, c)
          q_k(n_, c, 2) = -F3(n_, c)
          q_k(n_, c, 4) = -F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_j(n_, c, 1) =  sx(3, n_)*F(n_, c)
          q_j(n_, c, 4) = -F(n_, c)
          q_j(n_, c, 5) = -F3(n_, c)
          q_k(n_, c, 1) = -sx(2, n_)*F(n_, c)
          q_k(n_, c, 3) =  F(n_, c)
          q_k(n_, c, 5) =  F2(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q1kl_r128: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q1kl_r128: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q1kl_r128

  subroutine qak_q2kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r128), intent(in)  :: sx(3, N)
    real(r128), intent(in)  :: F (N, ncoeff)
    real(r128), intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r128), intent(in)  :: gradxyz(3, N)
    real(r128), intent(out) :: q_i(N, ncoeff, 5)
    real(r128), intent(out) :: q_j(N, ncoeff, 5)
    real(r128), intent(out) :: q_k(N, ncoeff, 5)
    real(r128) :: mr_dot_f
    integer(8) :: n_, c
    q_i = 0.0_r128;  q_j = 0.0_r128;  q_k = 0.0_r128
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          mr_dot_f = -sx(1, n_)*F1(n_, c) - sx(2, n_)*F2(n_, c) - sx(3, n_)*F3(n_, c)
          q_i(n_, c, 1) =  sx(2, n_)*F1(n_, c) + sx(1, n_)*F2(n_, c)
          q_i(n_, c, 2) = -F2(n_, c)
          q_i(n_, c, 3) = -F1(n_, c)
          q_j(n_, c, 1) =  mr_dot_f + 2.0_r128*sx(2, n_)*F2(n_, c)
          q_j(n_, c, 2) =  F1(n_, c)
          q_j(n_, c, 3) = -F2(n_, c)
          q_j(n_, c, 4) =  F3(n_, c)
          q_k(n_, c, 1) =  sx(2, n_)*F3(n_, c) + sx(3, n_)*F2(n_, c)
          q_k(n_, c, 3) = -F3(n_, c)
          q_k(n_, c, 4) = -F2(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) = -sx(3, n_)*F(n_, c)
          q_i(n_, c, 4) =  F(n_, c)
          q_i(n_, c, 5) =  F3(n_, c)
          q_k(n_, c, 1) =  sx(1, n_)*F(n_, c)
          q_k(n_, c, 2) = -F(n_, c)
          q_k(n_, c, 5) = -F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q2kl_r128: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q2kl_r128: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q2kl_r128

  subroutine qak_q3kl_r128(N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k)
    integer(8), intent(in)  :: N, ncoeff, lptype_id
    real(r128), intent(in)  :: sx(3, N)
    real(r128), intent(in)  :: F (N, ncoeff)
    real(r128), intent(in)  :: F1(N, ncoeff), F2(N, ncoeff), F3(N, ncoeff)
    real(r128), intent(in)  :: gradxyz(3, N)
    real(r128), intent(out) :: q_i(N, ncoeff, 5)
    real(r128), intent(out) :: q_j(N, ncoeff, 5)
    real(r128), intent(out) :: q_k(N, ncoeff, 5)
    real(r128) :: mr_dot_f
    integer(8) :: n_, c
    q_i = 0.0_r128;  q_j = 0.0_r128;  q_k = 0.0_r128
    select case (lptype_id)
    case (QAK_LPTYPE_D)
      do c = 1, ncoeff
        do n_ = 1, N
          mr_dot_f = -sx(1, n_)*F1(n_, c) - sx(2, n_)*F2(n_, c) - sx(3, n_)*F3(n_, c)
          q_i(n_, c, 1) =  sx(3, n_)*F1(n_, c) + sx(1, n_)*F3(n_, c)
          q_i(n_, c, 2) = -F3(n_, c)
          q_i(n_, c, 4) = -F1(n_, c)
          q_j(n_, c, 1) =  sx(3, n_)*F2(n_, c) + sx(2, n_)*F3(n_, c)
          q_j(n_, c, 3) = -F3(n_, c)
          q_j(n_, c, 4) = -F2(n_, c)
          q_k(n_, c, 1) =  mr_dot_f + 2.0_r128*sx(3, n_)*F3(n_, c)
          q_k(n_, c, 2) =  F1(n_, c)
          q_k(n_, c, 3) =  F2(n_, c)
          q_k(n_, c, 4) = -F3(n_, c)
        end do
      end do
    case (QAK_LPTYPE_S)
      do c = 1, ncoeff
        do n_ = 1, N
          q_i(n_, c, 1) =  sx(2, n_)*F(n_, c)
          q_i(n_, c, 3) = -F(n_, c)
          q_i(n_, c, 5) = -F2(n_, c)
          q_j(n_, c, 1) = -sx(1, n_)*F(n_, c)
          q_j(n_, c, 2) =  F(n_, c)
          q_j(n_, c, 5) =  F1(n_, c)
        end do
      end do
    case (QAK_LPTYPE_T)
      write(*,*) "qak_q3kl_r128: lptype_id=2 (DLPn 'T') not yet implemented"
      error stop
    case default
      write(*,*) 'qak_q3kl_r128: unknown lptype_id=', lptype_id
      error stop
    end select
  end subroutine qak_q3kl_r128

end module qkernel_mod
