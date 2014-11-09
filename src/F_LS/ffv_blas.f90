!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_blas.f90
!! @brief  BLAS routine
!! @author aics, iis
!<

! onishi version

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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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
!! @param [in,out] y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_xpay(y, x, a, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  y, x
double precision                                          ::  flop, a

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = dble(ix) * dble(jx) * dble(kx) * 2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
y(i, j, k) = x(i, j, k) + a * y(i, j, k)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_xpay


!> ********************************************************************
!! @brief AXPY
!! @param [in,out] y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_axpy(y, x, a, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  y, x
double precision                                          ::  flop, a

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = dble(ix) * dble(jx) * dble(kx) * 2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
y(i, j, k) = a * x(i, j, k) + y(i, j, k)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_axpy


!> ********************************************************************
!! @brief AXPYZ
!! @param [out]    z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_axpyz(z, x, y, a, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x, y, z
double precision                                          ::  flop, a

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = dble(ix) * dble(jx) * dble(kx) * 2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
z(i, j, k) = a * x(i, j, k) + y(i, j, k)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_axpyz


!> ********************************************************************
!! @brief AXPBYPZ
!! @param [in,out] z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     b    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_axpbypz(z, x, y, a, b, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x, y, z
double precision                                          ::  flop, a, b

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = dble(ix) * dble(jx) * dble(kx) * 4.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a, b)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
z(i, j, k) = a * x(i, j, k) + b * y(i, j, k) + z(i, j, k)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_axpbypz



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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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


! keno version


!> ********************************************************************
!! @brief 圧力Poissonの定数項bの計算
!! @param [out] rhs  右辺ベクトルbの自乗和
!! @param [out] b    RHS vector b
!! @param [in]  s_u  \sum {u^*}
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  dt   時間積分幅
!! @param [out] flop flop count
!<
subroutine blas_calc_b (rhs, b, s_u, bp, sz, g, dh, dt, flop)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                       ::  sz
double precision                                            ::  flop, rhs
real                                                        ::  dh, dt, c1, dv
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s_u, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
rhs = 0.0
c1 = dh / dt

flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0 + 8.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*20.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:rhs) &
!$OMP PRIVATE(dv) &
!$OMP FIRSTPRIVATE(ix, jx, kx, c1)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
  dv = c1 * s_u(i,j,k) * real(ibits(bp(i,j,k), Active, 1))
  b(i,j,k) = dv ! \frac{h^2}{\Delta t} \nabla u^*
  rhs = rhs + dble(dv*dv)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_calc_b


!> ********************************************************************
!! @brief DOT1
!! @param [out] r    内積
!! @param [in]  p    ベクトル
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
subroutine blas_dot1(r, p, bp, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
double precision                                          ::  flop, q, r

ix = sz(1)
jx = sz(2)
kx = sz(3)
r  = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*3.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:r) &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(q)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
   q = dble(p(i, j, k))
   r = r + q*q * real(ibits(bp(i,j,k), Active, 1))
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_dot1


!> ********************************************************************
!! @brief DOT2
!! @param [out] r    内積
!! @param [in]  p    ベクトル
!! @param [in]  q    ベクトル
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
subroutine blas_dot2(r, p, q, bp, sz, g, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, q
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  double precision                                          ::  flop, r

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	r  = 0.0

  flop = flop + dble(ix)*dble(jx)*dble(kx)*3.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:r) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
  do j=1,jx
  do i=1,ix
		r = r + dble(p(i, j, k)) * dble(q(i, j, k)) * real(ibits(bp(i,j,k), Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_dot2


!> ********************************************************************
!! @brief AX
!! @param [out] ap   AX
!! @param [in]  p    解ベクトル
!! @param [in]  bp   BCindexP
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
  subroutine blas_calc_ax(ap, p, bp, sz, g, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss, d0, d1, d2
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ap, p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  double precision                                          ::  flop

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*18.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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
    d0 = real(ibits(idx, bc_diag + 0, 1))
    d1 = real(ibits(idx, bc_diag + 1, 1))
    d2 = real(ibits(idx, bc_diag + 2, 1))
    dd = d2*4.0 + d1*2.0 + d0  ! diagonal
    !dd = real(ibits(idx, bc_diag, 3))  ! diagonal
    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)
		ap(i, j, k) = (ss - dd * p(i, j, k)) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calc_ax


!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [out]    r  残差ベクトル
!! @param [in]     p  解ベクトル
!! @param [in]     b  定数項
!! @param [in]     bp BCindexP
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!! @param [out] flop flop count
!<
  subroutine blas_calc_rk(r, p, b, bp, sz, g, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss, d0, d1, d2
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  r, p, b
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  double precision                                          ::  flop

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*19.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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
    !dd = real(ibits(idx, bc_diag, 3))  ! diagonal
    d0 = real(ibits(idx, bc_diag + 0, 1))
    d1 = real(ibits(idx, bc_diag + 1, 1))
    d2 = real(ibits(idx, bc_diag + 2, 1))
    dd = d2*4.0 + d1*2.0 + d0  ! diagonal
    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)
		r(i, j, k) = (b(i, j, k) - (ss - dd * p(i, j, k))) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calc_rk


!> ********************************************************************
!! @brief 残差の自乗和のみ
!! @param [out] res  残差の自乗和
!! @param [in]  p    圧力
!! @param [in]  b    RHS vector
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [out] flop flop count
!<
subroutine blas_calc_r2 (res, p, b, bp, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  dd, ss, dp, d0, d1, d2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
res = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*21.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*39.0d0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, dp, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

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

d0 = real(ibits(idx, bc_diag + 0, 1))
d1 = real(ibits(idx, bc_diag + 1, 1))
d2 = real(ibits(idx, bc_diag + 2, 1))
dd = d2*4.0 + d1*2.0 + d0  ! diagonal
! dd = real(ibits(idx, bc_diag, 3))  iff, K compiler is improved

ss = ndag_e * p(i+1,j  ,k  ) &
   + ndag_w * p(i-1,j  ,k  ) &
   + ndag_n * p(i  ,j+1,k  ) &
   + ndag_s * p(i  ,j-1,k  ) &
   + ndag_t * p(i  ,j  ,k+1) &
   + ndag_b * p(i  ,j  ,k-1)
dp = ( b(i,j,k) - (ss - dd * p(i,j,k)) ) * real(ibits(idx, Active, 1))
res = res + dble(dp*dp)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_calc_r2


!> ********************************************************************
!! @brief 残差
!! @param [out] r    残差
!! @param [in]  p    圧力
!! @param [in]  b    定数項
!! @param [in]  bp   BCindexP
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [in]  pn   係数行列
!! @param [out] flop flop count
!<
subroutine blas_calc_r_naive(r, p, b, bp, sz, g, pn, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  dd, ss
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  r, p, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 7) ::  pn
double precision                                          ::  flop

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix)*dble(jx)*dble(kx)*15.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
idx = bp(i,j,k)
ndag_w = pn(i,j,k,1)  ! w
ndag_e = pn(i,j,k,2)  ! e, non-diagonal
ndag_s = pn(i,j,k,3)  ! s
ndag_n = pn(i,j,k,4)  ! n
ndag_b = pn(i,j,k,5)  ! b
ndag_t = pn(i,j,k,6)  ! t
dd     = pn(i,j,k,7)
ss =  ndag_e * p(i+1,j  ,k  ) &
    + ndag_w * p(i-1,j  ,k  ) &
    + ndag_n * p(i  ,j+1,k  ) &
    + ndag_s * p(i  ,j-1,k  ) &
    + ndag_t * p(i  ,j  ,k+1) &
    + ndag_b * p(i  ,j  ,k-1)
r(i, j, k) = (b(i, j, k) - (ss - dd * p(i, j, k))) * real(ibits(idx, Active, 1))
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_calc_r_naive


!> ********************************************************************
!! @brief AX
!! @param [out] ap   AX
!! @param [in]  p    ベクトル
!! @param [in]  bp   BCindexP
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [in]  pn   係数行列
!! @param [out] flop flop count
!<
subroutine blas_calc_ax_naive(ap, p, bp, sz, g, pn, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  dd, ss
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ap, p
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 7) ::  pn
double precision                                          ::  flop

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix)*dble(jx)*dble(kx)*14.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
idx = bp(i,j,k)
ndag_w = pn(i,j,k,1)  ! w
ndag_e = pn(i,j,k,2)  ! e, non-diagonal
ndag_s = pn(i,j,k,3)  ! s
ndag_n = pn(i,j,k,4)  ! n
ndag_b = pn(i,j,k,5)  ! b
ndag_t = pn(i,j,k,6)  ! t
dd     = pn(i,j,k,7)
ss = ndag_e * p(i+1,j  ,k  ) &
   + ndag_w * p(i-1,j  ,k  ) &
   + ndag_n * p(i  ,j+1,k  ) &
   + ndag_s * p(i  ,j-1,k  ) &
   + ndag_t * p(i  ,j  ,k+1) &
   + ndag_b * p(i  ,j  ,k-1)
ap(i, j, k) = (ss - dd * p(i, j, k)) * real(ibits(idx, Active, 1))
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_calc_ax_naive
