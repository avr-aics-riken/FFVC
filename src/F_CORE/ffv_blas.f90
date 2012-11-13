!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan.
!
!********************************************************************

!> @file   ffv_blas.f90
!! @brief  BLAS routine
!! @author onishi
!<

!> ********************************************************************
!! @brief 要素のゼロクリア
!! @param [in,out] x  スカラー
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_clear(x, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

  do k=1-g,kx+g
  do j=1-g,jx+g
  do i=1-g,ix+g
		x(i, j, k) = 0.0
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_clear


!> ********************************************************************
!! @brief コピー
!! @param [out]    y  コピー先
!! @param [in]     x  ソース
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_copy(y, x, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  y, x

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

  do k=1-g,kx+g
  do j=1-g,jx+g
  do i=1-g,ix+g
		y(i, j, k) = x(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_copy


!> ********************************************************************
!! @brief XPAY
!! @param [in,out] y  ベクトル
!! @param [in]     x  ベクトル
!! @param [in]     a  係数
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_xpay(y, x, a, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  y, x
	real																											::  a

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
		y(i, j, k) = x(i, j, k) + a*y(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_xpay


!> ********************************************************************
!! @brief AXPY
!! @param [in,out] y  ベクトル
!! @param [in]     x  ベクトル
!! @param [in]     a  係数
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_axpy(y, x, a, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  y, x
	real																											::  a

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
		y(i, j, k) = a*x(i, j, k) + y(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_axpy


!> ********************************************************************
!! @brief AXPYZ
!! @param [out]    z  ベクトル
!! @param [in]     y  ベクトル
!! @param [in]     x  ベクトル
!! @param [in]     a  係数
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_axpyz(z, x, y, a, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x, y, z
	real																											::  a

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
		z(i, j, k) = a*x(i, j, k) + y(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_axpyz


!> ********************************************************************
!! @brief AXPBYPZ
!! @param [in,out] z  ベクトル
!! @param [in]     y  ベクトル
!! @param [in]     x  ベクトル
!! @param [in]     a  係数
!! @param [in]     b  係数
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_axpbypz(z, x, y, a, b, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x, y, z
	real																											::  a, b

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a, b)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
		z(i, j, k) = a*x(i, j, k) + b*y(i, j, k) + z(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_axpbypz


!> ********************************************************************
!! @brief DOT
!! @param [out]    pq 内積
!! @param [in]     p  ベクトル
!! @param [in]     q  ベクトル
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
subroutine blas_dot(pq, p, q, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, q
	real																											::  pq

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	pq = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:pq) &
!$OMP FIRSTPRIVATE(ix, jx, kx)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
		pq = pq + p(i, j, k)*q(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_dot


!> ********************************************************************
!! @brief AX
!! @param [out]    ap AX
!! @param [in]     p  ベクトル
!! @param [in]     bp BCindexP
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_calcax(ap, p, bp, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ap, p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
    ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
    ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b
    dd = real(ibits(idx, bc_diag, 3))  ! diagonal
    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)
		ap(i, j, k) = (ss - p(i, j, k)*dd) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calcax


!> ********************************************************************
!! @brief 定数項
!! @param [out]    b   定数項
!! @param [in]     s_0 ソース項0
!! @param [in]     s_1 ソース項1
!! @param [in]     bp BCindexP
!! @param [in]     dh 格子幅
!! @param [in]     dt 時間積分幅
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_calcb(b, s_0, s_1, bp, dh, dt, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  b, s_0, s_1
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
	real																											::  dh, dt, c1

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	c1 = dh/dt

!$OMP PARALLEL &
!$OMP PRIVATE(idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx, c1, dh)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
		b(i, j, k) = (c1*s_0(i, j, k) + dh*s_1(i, j, k)) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calcb


!> ********************************************************************
!! @brief 残差
!! @param [out]    r  残差
!! @param [in]     p  圧力
!! @param [in]     b  定数項
!! @param [in]     bp BCindexP
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_calcr(r, p, b, bp, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  r, p, b
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
    ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
    ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b
    dd = real(ibits(idx, bc_diag, 3))  ! diagonal
    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)
		r(i, j, k) = (b(i, j, k) - (ss - p(i, j, k)*dd)) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calcr


!> ********************************************************************
!! @brief 残差の自乗和
!! @param [out]    rr 残差の自乗和
!! @param [in]     p  圧力
!! @param [in]     b  定数項
!! @param [in]     bp BCindexP
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_calcr2(rr, p, b, bp, sz, g)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss
	real																											::  rr, r
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	rr = 0.0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP REDUCTION(+:rr) &
!$OMP PRIVATE(r) &
!$OMP FIRSTPRIVATE(ix, jx, kx)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
    ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
    ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b
    dd = real(ibits(idx, bc_diag, 3))  ! diagonal
    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)
		r = (b(i, j, k) - (ss - p(i, j, k)*dd)) * real(ibits(idx, Active, 1))
		rr = rr + r*r
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL
  return
  end subroutine blas_calcr2


!> ********************************************************************
!! @brief RBSOR によるスムーザー
!! @param [in,out] p     圧力
!! @param [in]     b     定数項
!! @param [in]     bp    BCindexP
!! @param [in]     ip    開始点インデクス
!! @param [in]     color 色（グループ）
!! @param [in]     omg   緩和係数
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル
!<
  subroutine blas_smoother_core (p, b, bp, ip, color, omg, sz, g)
	implicit none
	include 'ffv_f_params.h'
	integer                                                   ::  i, j, k, ix, jx, kx, g, idx
	integer, dimension(3)                                     ::  sz
	real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
	real                                                      ::  omg, dd, ss, dp, pp
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
	integer                                                   ::  ip, color

	ix = sz(1)
	jx = sz(2)
	kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, pp, ss, dp, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx, color, ip, omg)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

	do k=1,kx
	do j=1,jx
	do i=1+mod(k+j+color+ip,2), ix, 2
		idx = bp(i,j,k)
		ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
		ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w 
		ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
		ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
		ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
		ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b
		dd = real(ibits(idx, bc_diag, 3))  ! diagonal
		pp = p(i,j,k)
		ss = ndag_e * p(i+1,j  ,k  ) &
			 + ndag_w * p(i-1,j  ,k  ) &
			 + ndag_n * p(i  ,j+1,k  ) &
			 + ndag_s * p(i  ,j-1,k  ) &
			 + ndag_t * p(i  ,j  ,k+1) &
			 + ndag_b * p(i  ,j  ,k-1)
		dp = -(b(i, j, k) - ss)/dd
		p(i,j,k) = (1.0-omg)*pp + omg*dp
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_smoother_core
