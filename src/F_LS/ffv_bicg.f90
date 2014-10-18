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

!> @file   ffv_bicg.f90
!! @brief  BiCG algorithms
!! @author aics, iis
!<


!> ********************************************************************
!! @brief BiCG update x = x + alpha * p + omega * s
!! @param [in,out] x     ベクトル
!! @param [in]     alpha 係数
!! @param [in]     p     ベクトル
!! @param [in]     omega 係数
!! @param [in]     s     ベクトル
!! @param [in]     bp    BCindex P
!! @param [out]    x_l2  解ベクトルの自乗和
!! @param [out]    err   修正ベクトルの自乗和
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル
!! @param [out]    flop  flop count
!<
subroutine bicg_update_x(x, alpha, p, omega, s, bp, x_l2, err, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, aa
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x, p, s
real																											::  alp, omg, pp, dx, xx
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
double precision                                          ::  flop, x_l2, err, alpha, omega

ix = sz(1)
jx = sz(2)
kx = sz(3)
x_l2 = 0.0
err = 0.0
alp = real(alpha)
omg = real(omega)

flop = flop + dble(ix)*dble(jx)*dble(kx)*11.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:x_l2) &
!$OMP REDUCTION(+:err) &
!$OMP FIRSTPRIVATE(ix, jx, kx, alp, omg) &
!$OMP PRIVATE(aa, pp, dx, xx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
  aa = ibits(bp(i,j,k), Active, 1)
  xx = x(i, j, k)
  pp = xx + alp * p(i, j, k) + omg * s(i, j, k)
  dx = pp - xx
  x(i, j, k) = pp
  x_l2 = x_l2 + dble(pp*pp) * dble(aa)
  err  = err  + dble(dx*dx) * dble(aa)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine bicg_update_x


!> ********************************************************************
!! @brief BiCG update p =  r + beta * ( p - omega * q )
!! @param [in,out] p     ベクトル
!! @param [in]     r     ベクトル
!! @param [in]     beta  係数
!! @param [in]     omega 係数
!! @param [in]     q     ベクトル
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル
!! @param [out]    flop  flop count
!<
subroutine bicg_update_p(p, r, beta, omega, q, sz, g, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, q, r
real                                                      ::  omg, bt
double precision                                          ::  flop, beta, omega

ix = sz(1)
jx = sz(2)
kx = sz(3)

bt  = beta
omg = omega

flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, bt, omg)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
  p(i, j, k) = r(i, j, k) + bt * ( p(i, j, k) - omg * q(i, j, k) )
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine bicg_update_p