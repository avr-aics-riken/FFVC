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

!> @file   ffv_vbc_outer_face.f90
!! @brief  外部境界条件
!! @author aics
!<


!> ********************************************************************
!! @brief セルフェイスディリクレ値をセットする
!! @param [out] vf     セルフェイス速度
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界の面番号
!! @param [in]  bv     BCindex C
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 部分的な境界条件の実装のため、ガイドセル部のマスク情報を利用
!<
subroutine vobc_face_drchlt (vf, sz, g, m_face, bv, vec, nID)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, g, face, ix, jx, kx, m_face
integer, dimension(3)                                       ::  sz
real                                                        ::  u_bc, v_bc, w_bc, a
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
real, dimension(3)                                          ::  vec
integer, dimension(0:5)                                     ::  nID

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
face = m_face

! u_bcは境界速度
u_bc = vec(1)
v_bc = vec(2)
w_bc = vec(3)


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP PRIVATE(a)

FACES : select case (face)

case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=1,jx
if ( ibits(bv(1, j, k), bc_face_W, bitw_5) == obc_mask ) then
  a = u_bc * real(ibits(bv(0,j,k), State, 1)) * real(ibits(bv(1,j,k), State, 1))
  vf(0, j, k, 1) = a
endif
end do
end do
!$OMP END DO


case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=1,jx
if ( ibits(bv(ix, j, k), bc_face_E, bitw_5) == obc_mask ) then
  a = u_bc * real(ibits(bv(ix,j,k), State, 1)) * real(ibits(bv(ix+1,j,k), State, 1))
  vf(ix, j, k, 1) = a
endif
end do
end do
!$OMP END DO


case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=1,ix
if ( ibits(bv(i, 1, k), bc_face_S, bitw_5) == obc_mask ) then
  a = v_bc * real(ibits(bv(i,0,k), State, 1)) * real(ibits(bv(i,1,k), State, 1))
  vf(i, 0, k, 2) = a
endif
end do
end do
!$OMP END DO


case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=1,ix
if ( ibits(bv(i, jx, k), bc_face_N, bitw_5) == obc_mask ) then
  a = v_bc * real(ibits(bv(i,jx,k), State, 1)) * real(ibits(bv(i,jx+1,k), State, 1))
  vf(i, jx, k, 2) = a
endif
end do
end do
!$OMP END DO


case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=1,ix
if ( ibits(bv(i, j ,1), bc_face_B, bitw_5) == obc_mask ) then
  a = w_bc * real(ibits(bv(i,j,0), State, 1)) * real(ibits(bv(i,j,1), State, 1))
  vf(i, j, 0, 3) = a
endif
end do
end do
!$OMP END DO


case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=1,ix
if ( ibits(bv(i, j, kx), bc_face_T, bitw_5) == obc_mask ) then
  a = w_bc * real(ibits(bv(i,j,kx), State, 1)) * real(ibits(bv(i,j,kx+1), State, 1))
  vf(i, j, kx, 3) = a
endif
end do
end do
!$OMP END DO

case default
end select FACES

!$OMP END PARALLEL


return
end subroutine vobc_face_drchlt


!> ********************************************************************
!! @brief 外部境界面の流入出量を求める
!! @param [out]    sum    領域境界の流速の積算値 \sum{vf}
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 面番号
!! @param [in]     vf     セルフェイス速度
!! @param [in]     bv     BCindex C
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 有効セルのマスクを掛けて、流量を積算
!<
subroutine vobc_face_massflow (sum, sz, g, m_face, vf, bv, nID)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, g, ix, jx, kx, face, m_face
integer, dimension(3)                                     ::  sz
real                                                      ::  sum, a, s
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vf
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
integer, dimension(0:5)                                   ::  nID

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
face = m_face

a = 0.0   ! sum


!$OMP PARALLEL &
!$OMP REDUCTION(+:a) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face) &
!$OMP PRIVATE(s)

FACES : select case (face)

case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=1,jx
  s = real(ibits(bv(0,j,k), State, 1)) * real(ibits(bv(1,j,k), State, 1))
  a = a + vf(0,j,k,1) * s
end do
end do
!$OMP END DO


case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=1,jx
  s = real(ibits(bv(ix,j,k), State, 1)) * real(ibits(bv(ix+1,j,k), State, 1))
  a = a + vf(ix,j,k,1) * s
end do
end do
!$OMP END DO


case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=1,ix
  s = real(ibits(bv(i,0,k), State, 1)) * real(ibits(bv(i,1,k), State, 1))
  a = a + vf(i,0,k,2) * s
end do
end do
!$OMP END DO


case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=1,ix
  s = real(ibits(bv(i,jx,k), State, 1)) * real(ibits(bv(i,jx+1,k), State, 1))
  a = a + vf(i,jx,k,2) * s
end do
end do
!$OMP END DO


case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=1,ix
  s = real(ibits(bv(i,j,0), State, 1)) * real(ibits(bv(i,j,1), State, 1))
  a = a + vf(i,j,0,3) * s
end do
end do
!$OMP END DO


case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=1,ix
  s = real(ibits(bv(i,j,kx), State, 1)) * real(ibits(bv(i,j,kx+1), State, 1))
  a = a + vf(i,j,kx,3) * s
end do
end do
!$OMP END DO

case default
end select FACES
!$OMP END PARALLEL

sum = a

return
end subroutine vobc_face_massflow
