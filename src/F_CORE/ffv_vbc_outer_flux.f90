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

!> @file   ffv_vbc_outer_flux.f90
!! @brief  流束形式による外部境界条件
!! @author aics
!<


!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv     疑似ベクトルの空間項の評価値
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界処理のときの面番号
!! @param [in]  dh     格子幅
!! @param [in]  rei    Reynolds数の逆数
!! @param [in]  v0     速度ベクトル（n-step）
!! @param [in]  bv     BCindex C
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @param [out] flop   浮動小数点演算数
!! @note vecには，流入条件のとき指定速度
!!  mskで部分的な速度を与える
!<
subroutine vobc_pv_specv (wv, sz, g, m_face, dh, rei, v0, bv, vec, nID, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, g, face, m_face
integer                                                   ::  ix, jx, kx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  Up, Vp, Wp, Ur, Vr, Wr
real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
real                                                      ::  fu, fv, fw, c, ac, msk
real                                                      ::  u_bc, v_bc, w_bc, m
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
real, dimension(3)                                        ::  vec
integer, dimension(0:5)                                   ::  nID
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
face = m_face

dh1= 1.0/dh
dh2= rei*dh1*dh1

! u_bcは境界速度
u_bc = vec(1)
v_bc = vec(2)
w_bc = vec(3)

flop = flop + 13.0d0 ! DP 18 flop

m = 0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP FIRSTPRIVATE(dh1, dh2) &
!$OMP PRIVATE(Up, Vp, Wp, Ur, Vr, Wr) &
!$OMP PRIVATE(fu, fv, fw, EX, EY, EZ, c, ac, msk) &
!$OMP PRIVATE(i, j, k)

FACES : select case (face)

case (X_minus)

i = 1
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
Up = v0(i,j,k,1)
Vp = v0(i,j,k,2)
Wp = v0(i,j,k,3)

Ur = u_bc
Vr = v_bc
Wr = w_bc
c  = u_bc
ac = abs(c)

EX = Up - Ur
EY = Vp - Vr
EZ = Wp - Wr

fu = 0.5*(c*(Up+Ur) - ac*EX)
fv = 0.5*(c*(Vp+Vr) - ac*EY)
fw = 0.5*(c*(Wp+Wr) - ac*EZ)

msk = real(ibits(bv(0,j,k), State, 1))

wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 - EX*dh2 ) * msk
wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 - EY*dh2 ) * msk
wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 - EZ*dh2 ) * msk
m = m + 1.0
endif
end do
end do
!$OMP END DO


case (X_plus)

i = ix
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
Up = v0(i,j,k,1)
Vp = v0(i,j,k,2)
Wp = v0(i,j,k,3)

Ur = u_bc
Vr = v_bc
Wr = w_bc
c  = u_bc
ac = abs(c)

EX = Ur - Up
EY = Vr - Vp
EZ = Wr - Wp

fu = 0.5*(c*(Ur+Up) - ac*EX)
fv = 0.5*(c*(Vr+Vp) - ac*EY)
fw = 0.5*(c*(Wr+Wp) - ac*EZ)

msk = real(ibits(bv(ix+1,j,k), State, 1))

wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 ) * msk
wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 ) * msk
wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 ) * msk
m = m + 1.0
endif
end do
end do
!$OMP END DO


case (Y_minus)

j = 1
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
Up = v0(i,j,k,1)
Vp = v0(i,j,k,2)
Wp = v0(i,j,k,3)

Ur = u_bc
Vr = v_bc
Wr = w_bc
c  = v_bc
ac = abs(c)

EX = Up - Ur
EY = Vp - Vr
EZ = Wp - Wr

fu = 0.5*(c*(Up+Ur) - ac*EX)
fv = 0.5*(c*(Vp+Vr) - ac*EY)
fw = 0.5*(c*(Wp+Wr) - ac*EZ)

msk = real(ibits(bv(i,0,k), State, 1))

wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 - EX*dh2 ) * msk
wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 - EY*dh2 ) * msk
wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 - EZ*dh2 ) * msk
m = m + 1.0
endif
end do
end do
!$OMP END DO


case (Y_plus)

j = jx
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
Up = v0(i,j,k,1)
Vp = v0(i,j,k,2)
Wp = v0(i,j,k,3)

Ur = u_bc
Vr = v_bc
Wr = w_bc
c  = v_bc
ac = abs(c)

EX = Ur - Up
EY = Vr - Vp
EZ = Wr - Wp

fu = 0.5*(c*(Ur+Up) - ac*EX)
fv = 0.5*(c*(Vr+Vp) - ac*EY)
fw = 0.5*(c*(Wr+Wp) - ac*EZ)

msk = real(ibits(bv(i,jx+1,k), State, 1))

wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 ) * msk
wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 ) * msk
wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 ) * msk
m = m + 1.0
endif
end do
end do
!$OMP END DO


case (Z_minus)

k = 1
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
Up = v0(i,j,k,1)
Vp = v0(i,j,k,2)
Wp = v0(i,j,k,3)

Ur = u_bc
Vr = v_bc
Wr = w_bc
c  = w_bc
ac = abs(c)

EX = Up - Ur
EY = Vp - Vr
EZ = Wp - Wr

fu = 0.5*(c*(Up+Ur) - ac*EX)
fv = 0.5*(c*(Vp+Vr) - ac*EY)
fw = 0.5*(c*(Wp+Wr) - ac*EZ)

msk = real(ibits(bv(i,j,0), State, 1))

wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 - EX*dh2 ) * msk
wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 - EY*dh2 ) * msk
wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 - EZ*dh2 ) * msk
m = m + 1.0
endif
end do
end do
!$OMP END DO


case (Z_plus)

k = kx
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
Up = v0(i,j,k,1)
Vp = v0(i,j,k,2)
Wp = v0(i,j,k,3)

Ur = u_bc
Vr = v_bc
Wr = w_bc
c  = w_bc
ac = abs(c)

EX = Ur - Up
EY = Vr - Vp
EZ = Wr - Wp

fu = 0.5*(c*(Ur+Up) - ac*EX)
fv = 0.5*(c*(Vr+Vp) - ac*EY)
fw = 0.5*(c*(Wr+Wp) - ac*EZ)

msk = real(ibits(bv(i,j,kx+1), State, 1))

wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 ) * msk
wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 ) * msk
wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 ) * msk
m = m + 1.0
endif
end do
end do
!$OMP END DO


case default
end select FACES
!$OMP END PARALLEL

flop = flop + dble(m)*28.0d0

return
end subroutine vobc_pv_specv


!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv     疑似ベクトルの空間項の評価値
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界処理のときの面番号
!! @param [in]  dh     格子幅
!! @param [in]  rei    Reynolds数の逆数
!! @param [in]  v0     速度ベクトル（n-step）
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @param [out] flop   浮動小数点演算数
!! @note 対流流束はゼロ，壁面法線方向の1階微分もゼロ
!<
subroutine vobc_pv_wall (wv, sz, g, m_face, dh, rei, v0, vec, nID, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, g, face, m_face
integer                                                   ::  ix, jx, kx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, rix, rjx, rkx
real                                                      ::  u_bc, v_bc, w_bc
real                                                      ::  dh, dh1, dh2, rei
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
real, dimension(3)                                        ::  vec
integer, dimension(0:5)                                   ::  nID

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
face = m_face

dh1= 1.0/dh
dh2= rei*dh1*dh1*2.0

u_bc = vec(1)
v_bc = vec(2)
w_bc = vec(3)


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, dh2, face) &
!$OMP PRIVATE(i, j, k)

FACES : select case (face)

case (X_minus)

i = 1
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
! wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dx = 0
wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dx
wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dx
end do
end do
!$OMP END DO


case (X_plus)

i = ix
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
! wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dx = 0
wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dx
wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dx
end do
end do
!$OMP END DO


case (Y_minus)

j = 1
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dy
! wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dy = 0
wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dy
end do
end do
!$OMP END DO


case (Y_plus)

j = jx
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dy
! wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dy = 0
wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dy
end do
end do
!$OMP END DO


case (Z_minus)

k = 1
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dz
wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dz
! wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dz = 0
end do
end do
!$OMP END DO


case (Z_plus)

k = kx
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dz
wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dz
! wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dz = 0
end do
end do
!$OMP END DO


case default
end select FACES
!$OMP END PARALLEL


rix = dble(jx)*dble(kx)
rjx = dble(ix)*dble(kx)
rkx = dble(ix)*dble(jx)

FACES2 : select case (face)

case (X_minus)
flop = flop + rix*6.0d0

case (X_plus)
flop = flop + rix*6.0d0

case (Y_minus)
flop = flop + rjx*6.0d0

case (Y_plus)
flop = flop + rjx*6.0d0

case (Z_minus)
flop = flop + rkx*6.0d0

case (Z_plus)
flop = flop + rkx*6.0d0

case default
end select FACES2

return
end subroutine vobc_pv_wall

