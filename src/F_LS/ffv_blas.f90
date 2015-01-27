!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_blas.f90
!! @brief  BLAS routine
!! @author aics, iis
!<

! based on Onishi version

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
!! @brief AXPYZ
!! @param [out]    z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_triad(z, x, y, a, sz, g, flop)
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
end subroutine blas_triad


!> ********************************************************************
!! @brief BiCGstab 2
!! @param [in,out] z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     b    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_bicg_2(z, x, y, a, b, sz, g, flop)
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
end subroutine blas_bicg_2




!> ********************************************************************
!! @brief 圧力Poissonの定数項bの計算
!! @param [out] rhs  右辺ベクトルbの自乗和
!! @param [out] b    RHS vector b
!! @param [in]  div  div {u^*}
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  dt   時間積分幅
!! @param [out] flop flop count
!<
subroutine blas_calc_b (rhs, b, div, bp, sz, g, dh, dt, flop)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                       ::  sz
double precision                                            ::  flop, rhs
real                                                        ::  dt, c1, d, dx
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp
real, dimension(3)                                          ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)
rhs = 0.0
dx = dh(1)
c1 = dx * dx / dt

flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0 + 10.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:rhs) &
!$OMP PRIVATE(d) &
!$OMP FIRSTPRIVATE(ix, jx, kx, c1)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
  d = c1 * div(i,j,k) * real(ibits(bp(i,j,k), Active, 1))
  b(i,j,k) = d ! \frac{dx^2}{\Delta t} \nabla u^*
  rhs = rhs + dble(d*d)
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
!! @param [in]  dh   格子幅
!! @param [out] flop flop count
!<
  subroutine blas_calc_ax(ap, p, bp, sz, g, dh, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
  real                                                      ::  d_w, d_e, d_s, d_n, d_b, d_t
  real                                                      ::  dd, ss
  real                                                      ::  r_xx, r_xy, r_xz, r_x2, r_y2, r_z2
  real, dimension(3)                                        ::  dh
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ap, p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  double precision                                          ::  flop

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*35.0d0 + 19.0d0

  r_xx = 1.0
  r_xy = dh(1) / dh(2)
  r_xz = dh(1) / dh(3)
  r_x2 = r_xx * r_xx
  r_y2 = r_xy * r_xy
  r_z2 = r_xz * r_xz

!$OMP PARALLEL &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t, dd, ss, idx) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t) &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP FIRSTPRIVATE(r_x2, r_y2, r_z2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    c_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    c_e = real(ibits(idx, bc_ndag_E, 1))  ! e
    c_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    c_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    c_b = real(ibits(idx, bc_ndag_B, 1))  ! b
    c_t = real(ibits(idx, bc_ndag_T, 1))  ! t

    d_w = real(ibits(idx, bc_dn_W, 1))
    d_e = real(ibits(idx, bc_dn_E, 1))
    d_s = real(ibits(idx, bc_dn_S, 1))
    d_n = real(ibits(idx, bc_dn_N, 1))
    d_b = real(ibits(idx, bc_dn_B, 1))
    d_t = real(ibits(idx, bc_dn_T, 1))

    dd = r_x2 * (c_w + c_e) &
       + r_y2 * (c_s + c_n) &
       + r_z2 * (c_b + c_t) &
       + 2.0                &
       *(r_x2 * (d_w + d_e) &
       + r_y2 * (d_s + d_n) &
       + r_z2 * (d_b + d_t) )

    ss = r_x2 * ( c_e * p(i+1,j  ,k  ) + c_w * p(i-1,j  ,k  ) ) &
       + r_y2 * ( c_n * p(i  ,j+1,k  ) + c_s * p(i  ,j-1,k  ) ) &
       + r_z2 * ( c_t * p(i  ,j  ,k+1) + c_b * p(i  ,j  ,k-1) )

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
!! @param [out] r    残差ベクトル
!! @param [in]  p    解ベクトル
!! @param [in]  b    定数項
!! @param [in]  bp   BCindexP
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [in]  dh   格子幅
!! @param [out] flop flop count
!<
  subroutine blas_calc_rk(r, p, b, bp, sz, g, dh, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
  real                                                      ::  d_w, d_e, d_s, d_n, d_b, d_t
  real                                                      ::  dd, ss
  real                                                      ::  r_xx, r_xy, r_xz, r_x2, r_y2, r_z2
  real, dimension(3)                                        ::  dh
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  r, p, b
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  double precision                                          ::  flop

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  r_xx = 1.0
  r_xy = dh(1) / dh(2)
  r_xz = dh(1) / dh(3)
  r_x2 = r_xx * r_xx
  r_y2 = r_xy * r_xy
  r_z2 = r_xz * r_xz

  flop = flop + dble(ix)*dble(jx)*dble(kx)*36.0d0 + 19.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t, dd, ss, idx) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t) &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP FIRSTPRIVATE(r_x2, r_y2, r_z2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    c_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    c_e = real(ibits(idx, bc_ndag_E, 1))  ! e
    c_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    c_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    c_b = real(ibits(idx, bc_ndag_B, 1))  ! b
    c_t = real(ibits(idx, bc_ndag_T, 1))  ! t

    d_w = real(ibits(idx, bc_dn_W, 1))
    d_e = real(ibits(idx, bc_dn_E, 1))
    d_s = real(ibits(idx, bc_dn_S, 1))
    d_n = real(ibits(idx, bc_dn_N, 1))
    d_b = real(ibits(idx, bc_dn_B, 1))
    d_t = real(ibits(idx, bc_dn_T, 1))

    dd = r_x2 * (c_w + c_e) &
       + r_y2 * (c_s + c_n) &
       + r_z2 * (c_b + c_t) &
       + 2.0                &
       *(r_x2 * (d_w + d_e) &
       + r_y2 * (d_s + d_n) &
       + r_z2 * (d_b + d_t) )

    ss = r_x2 * ( c_e * p(i+1,j  ,k  ) + c_w * p(i-1,j  ,k  ) ) &
       + r_y2 * ( c_n * p(i  ,j+1,k  ) + c_s * p(i  ,j-1,k  ) ) &
       + r_z2 * ( c_t * p(i  ,j  ,k+1) + c_b * p(i  ,j  ,k-1) )

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
!! @param [in]  dh   格子幅
!! @param [out] flop flop count
!<
subroutine blas_calc_r2 (res, p, b, bp, sz, g, dh, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res
real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
real                                                      ::  d_w, d_e, d_s, d_n, d_b, d_t
real                                                      ::  dd, ss, dp
real                                                      ::  r_xx, r_xy, r_xz, r_x2, r_y2, r_z2
real, dimension(3)                                        ::  dh
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
res = 0.0

r_xx = 1.0
r_xy = dh(1) / dh(2)
r_xz = dh(1) / dh(3)
r_x2 = r_xx * r_xx
r_y2 = r_xy * r_xy
r_z2 = r_xz * r_xz

flop = flop + dble(ix)*dble(jx)*dble(kx)*38.0d0 + 19.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t, dd, ss, dp, idx) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t) &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP FIRSTPRIVATE(r_x2, r_y2, r_z2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
idx = bp(i,j,k)
c_w = real(ibits(idx, bc_ndag_W, 1))  ! w
c_e = real(ibits(idx, bc_ndag_E, 1))  ! e
c_s = real(ibits(idx, bc_ndag_S, 1))  ! s
c_n = real(ibits(idx, bc_ndag_N, 1))  ! n
c_b = real(ibits(idx, bc_ndag_B, 1))  ! b
c_t = real(ibits(idx, bc_ndag_T, 1))  ! t

d_w = real(ibits(idx, bc_dn_W, 1))
d_e = real(ibits(idx, bc_dn_E, 1))
d_s = real(ibits(idx, bc_dn_S, 1))
d_n = real(ibits(idx, bc_dn_N, 1))
d_b = real(ibits(idx, bc_dn_B, 1))
d_t = real(ibits(idx, bc_dn_T, 1))

dd = r_x2 * (c_w + c_e) &
   + r_y2 * (c_s + c_n) &
   + r_z2 * (c_b + c_t) &
   + 2.0                &
   *(r_x2 * (d_w + d_e) &
   + r_y2 * (d_s + d_n) &
   + r_z2 * (d_b + d_t) )

ss = r_x2 * ( c_e * p(i+1,j  ,k  ) + c_w * p(i-1,j  ,k  ) ) &
   + r_y2 * ( c_n * p(i  ,j+1,k  ) + c_s * p(i  ,j-1,k  ) ) &
   + r_z2 * ( c_t * p(i  ,j  ,k+1) + c_b * p(i  ,j  ,k-1) )

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
!! @brief BiCGstabの部分演算1
!! @param [in,out] p    ベクトル
!! @param [in]     r    ベクトル
!! @param [in]     q    ベクトル
!! @param [in]     beta 係数
!! @param [in]     omg  係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_bicg_1(p, r, q, beta, omg, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, r, q
double precision                                          ::  flop, beta, omg

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = dble(ix) * dble(jx) * dble(kx) * 4.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, beta, omg)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
  p(i,j,k) = r(i,j,k) + beta * ( p(i,j,k) - omg * q(i,j,k) )
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_bicg_1

