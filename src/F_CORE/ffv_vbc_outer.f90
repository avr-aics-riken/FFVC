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
! Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
! Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_vbc_outer.f90
!! @brief  外部境界条件
!! @author aics
!<
    

!> ********************************************************************
!! @brief ガイドセルの速度指定境界条件を設定するために必要な参照値をセットする
!! @param [out] v      セルセンタ速度
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界の面番号
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 流束型の境界条件を用いるので，内点の計算に使う参照点に値があればよい（1層）
!<
    subroutine vobc_cc_drchlt (v, sz, g, m_face, vec, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx, gc, m_face
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc, v_bc, w_bc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(3)                                          ::  vec
    integer, dimension(0:5)                                     ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    gc = g
    face = m_face
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face, gc) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)
    
    case (X_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
      do i=1-gc, 0
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (X_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
      do i=ix+1, ix+gc
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
      do j=1-gc, 0
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
      do j=jx+1, jx+gc
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
      do k=1-gc, 0
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
      do k=kx+1, kx+gc
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL
    
    return
    end subroutine vobc_cc_drchlt


!> ********************************************************************
!! @brief ノイマン条件 Cell center
!! @param [in,out] v      速度ベクトル
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!<
    subroutine vobc_cc_neumann (v, sz, g, m_face, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, gc, m_face
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    gc = g
    face = m_face

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, gc, face) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
      do i=1-gc, 0
        v(i, j, k, 1) = v(1-i, j, k, 1)
        v(i, j, k, 2) = v(1-i, j, k, 2)
        v(i, j, k, 3) = v(1-i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO


    case (X_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
      do i=ix+1, ix+gc
        v(i, j, k, 1) = v(2*ix-i+1, j, k, 1)
        v(i, j, k, 2) = v(2*ix-i+1, j, k, 2)
        v(i, j, k, 3) = v(2*ix-i+1, j, k, 3)
      end do
      end do
      end do
!$OMP END DO


    case (Y_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
      do j=1-gc, 0
        v(i, j, k, 1) = v(i, 1-j, k, 1)
        v(i, j, k, 2) = v(i, 1-j, k, 2)
        v(i, j, k, 3) = v(i, 1-j, k, 3)
      end do
      end do
      end do
!$OMP END DO


    case (Y_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
      do j=jx+1, jx+gc
        v(i, j, k, 1) = v(i, 2*jx-j+1, k, 1)
        v(i, j, k, 2) = v(i, 2*jx-j+1, k, 2)
        v(i, j, k, 3) = v(i, 2*jx-j+1, k, 3)
      end do
      end do
      end do
!$OMP END DO


    case (Z_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
      do k=1-gc, 0
        v(i, j, k, 1) = v(i, j, 1-k, 1)
        v(i, j, k, 2) = v(i, j, 1-k, 2)
        v(i, j, k, 3) = v(i, j, 1-k, 3)
      end do
      end do
      end do
!$OMP END DO


    case (Z_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
      do k=kx+1, kx+gc
        v(i, j, k, 1) = v(i, j, 2*kx-k+1, 1)
        v(i, j, k, 2) = v(i, j, 2*kx-k+1, 2)
        v(i, j, k, 3) = v(i, j, 2*kx-k+1, 3)
      end do
      end do
      end do
!$OMP END DO


    case default
    end select FACES
!$OMP END PARALLEL

    return 
    end subroutine vobc_cc_neumann



!> ********************************************************************
!! @brief 外部境界面の速度をコピーする
!! @param [out] v      速度ベクトル（セルセンタ）
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 面番号
!! @param [in]  vc     セルセンタ疑似速度 u^*
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!<
subroutine vobc_cc_copy (v, sz, g, m_face, vc, nID)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, g, ix, jx, kx, face, gc, m_face
integer, dimension(3)                                       ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc
integer, dimension(0:5)                                     ::  nID

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
gc = g
face = m_face

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, face, gc) PRIVATE(i, j, k)

FACES : select case (face)
case (X_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1-gc, 0
v(i, j, k, 1) = vc(i, j, k, 1)
v(i, j, k, 2) = vc(i, j, k, 2)
v(i, j, k, 3) = vc(i, j, k, 3)
end do
end do
end do
!$OMP END DO

case (X_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=ix+1, ix+gc
v(i, j, k, 1) = vc(i, j, k, 1)
v(i, j, k, 2) = vc(i, j, k, 2)
v(i, j, k, 3) = vc(i, j, k, 3)
end do
end do
end do
!$OMP END DO

case (Y_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
do j=1-gc, 0
v(i, j, k, 1) = vc(i, j, k, 1)
v(i, j, k, 2) = vc(i, j, k, 2)
v(i, j, k, 3) = vc(i, j, k, 3)
end do
end do
end do
!$OMP END DO

case (Y_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
do j=jx+1, jx+gc
v(i, j, k, 1) = vc(i, j, k, 1)
v(i, j, k, 2) = vc(i, j, k, 2)
v(i, j, k, 3) = vc(i, j, k, 3)
end do
end do
end do
!$OMP END DO

case (Z_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
do k=1-gc, 0
v(i, j, k, 1) = vc(i, j, k, 1)
v(i, j, k, 2) = vc(i, j, k, 2)
v(i, j, k, 3) = vc(i, j, k, 3)
end do
end do
end do
!$OMP END DO

case (Z_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
do k=kx+1, kx+gc
v(i, j, k, 1) = vc(i, j, k, 1)
v(i, j, k, 2) = vc(i, j, k, 2)
v(i, j, k, 3) = vc(i, j, k, 3)
end do
end do
end do
!$OMP END DO

case default
end select FACES

!$OMP END PARALLEL

return
end subroutine vobc_cc_copy



!> ********************************************************************
!! @brief 速度の外部境界：　トラクションフリー（セルセンタ速度）
!! @param [in,out] v      セルセンタ速度ベクトル
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in,out] vf     セルフェイス速度
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 今のところ、トラクションフリー面は全て流体
!<
subroutine vobc_cc_tfree (v, sz, g, m_face, vf, nID)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, face, g, m_face
integer, dimension(3)                                     ::  sz
real                                                      ::  ut, vt, wt
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vf
integer, dimension(0:5)                                   ::  nID

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
face = m_face


FACES : select case (face)
case (X_minus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(ut, vt, wt) &
!$OMP PRIVATE(j, k)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
ut = 0.5 * (vf(-1, j, k, 1) + vf(0, j  , k  , 1))
vt = 0.5 * (vf( 0, j, k, 2) + vf(0, j-1, k  , 2))
wt = 0.5 * (vf( 0, j, k, 3) + vf(0, j  , k-1, 3))

v(-1, j, k, 1) = ut
v(-1, j, k, 2) = vt
v(-1, j, k, 3) = wt

v(0, j, k, 1) = ut
v(0, j, k, 2) = vt
v(0, j, k, 3) = wt
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (X_plus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(ut, vt, wt) &
!$OMP PRIVATE(j, k)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
ut = 0.5 * (vf(ix  , j, k, 1) + vf(ix+1, j  , k  , 1))
vt = 0.5 * (vf(ix+1, j, k, 2) + vf(ix+1, j-1, k  , 2))
wt = 0.5 * (vf(ix+1, j, k, 3) + vf(ix+1, j  , k-1, 3))

v(ix+1, j, k, 1) = ut
v(ix+1, j, k, 2) = vt
v(ix+1, j, k, 3) = wt

v(ix+2, j, k, 1) = ut
v(ix+2, j, k, 2) = vt
v(ix+2, j, k, 3) = wt
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Y_minus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(ut, vt, wt) &
!$OMP PRIVATE(i, k)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
ut = 0.5 * (vf(i, 0, k, 1) + vf(i-1, 0, k  , 1))
vt = 0.5 * (vf(i,-1, k, 2) + vf(i  , 0, k  , 2))
wt = 0.5 * (vf(i, 0, k, 3) + vf(i  , 0, k-1, 3))

v(i, -1, k, 1) = ut
v(i, -1, k, 2) = vt
v(i, -1, k, 3) = wt

v(i, 0, k, 1) = ut
v(i, 0, k, 2) = vt
v(i, 0, k, 3) = wt
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Y_plus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(ut, vt, wt) &
!$OMP PRIVATE(i, k)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
ut = 0.5 * (vf(i, jx+1, k, 1) + vf(i-1, jx+1, k  , 1))
vt = 0.5 * (vf(i, jx  , k, 2) + vf(i  , jx+1, k  , 2))
wt = 0.5 * (vf(i, jx+1, k, 3) + vf(i  , jx+1, k-1, 3))

v(i, jx+1, k, 1) = ut
v(i, jx+1, k, 2) = vt
v(i, jx+1, k, 3) = wt

v(i, jx+2, k, 1) = ut
v(i, jx+2, k, 2) = vt
v(i, jx+2, k, 3) = wt
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Z_minus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(ut, vt, wt) &
!$OMP PRIVATE(i, j)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
ut = 0.5 * (vf(i, j, 0, 1) + vf(i-1, j  , 0, 1))
vt = 0.5 * (vf(i, j, 0, 2) + vf(i  , j-1, 0, 2))
wt = 0.5 * (vf(i, j,-1, 3) + vf(i  , j  , 0, 3))

v(i, j, -1, 1) = ut
v(i, j, -1, 2) = vt
v(i, j, -1, 3) = wt

v(i, j, 0, 1) = ut
v(i, j, 0, 2) = vt
v(i, j, 0, 3) = wt
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Z_plus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(ut, vt, wt) &
!$OMP PRIVATE(i, j)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
ut = 0.5 * (vf(i, j, kx+1, 1) + vf(i-1, j  , kx+1, 1))
vt = 0.5 * (vf(i, j, kx+1, 2) + vf(i  , j-1, kx+1, 2))
wt = 0.5 * (vf(i, j, kx  , 3) + vf(i  , j  , kx+1, 3))

v(i, j, kx+1, 1) = ut
v(i, j, kx+1, 2) = vt
v(i, j, kx+1, 3) = wt

v(i, j, kx+2, 1) = ut
v(i, j, kx+2, 2) = vt
v(i, j, kx+2, 3) = wt
end do
end do
!$OMP END DO
!$OMP END PARALLEL

case default
end select FACES


return
end subroutine vobc_cc_tfree


!> ********************************************************************
!! @brief 速度の外部境界：　トラクションフリー（セルフェイス速度）
!! @param [in,out] vf     セルフェイス速度
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!! ループは0~ixで回す．0から始めるのは，後半処理でスタガード位置の変数からセンターへ内挿するときに
!! 参照する端点を含めるため．
!<
subroutine vobc_cf_tfree (vf, sz, g, m_face, nID)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, face, g, m_face
integer, dimension(3)                                     ::  sz
integer, dimension(0:5)                                   ::  nID
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vf

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
face = m_face


FACES : select case (face)
case (X_minus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=1,jx
vf(-1, j, k, 1) = vf(0, j, k, 1)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=0,jx
vf(0, j, k, 2) = vf(1, j, k, 2) + vf(0, j+1, k  , 1) - vf(0, j, k, 1)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=0,kx
do j=1,jx
vf(0, j, k, 3) = vf(1, j, k, 3) + vf(0, j  , k+1, 1) - vf(0, j, k, 1)
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (X_plus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=1,jx
vf(ix+1, j, k, 1) = vf(ix, j, k, 1)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=1,kx
do j=0,jx
vf(ix+1, j, k, 2) = vf(ix, j, k, 2) - vf(ix, j+1, k  , 1) + vf(ix, j, k, 1)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
do k=0,kx
do j=1,jx
vf(ix+1, j, k, 3) = vf(ix, j, k, 3) - vf(ix, j  , k+1, 1) + vf(ix, j, k, 1)
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Y_minus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=1,ix
vf(i, -1, k, 2) = vf(i, 0, k, 2)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=0,ix
vf(i, 0, k, 1) = vf(i, 1, k, 1) + vf(i+1, 0, k  , 2) - vf(i, 0, k, 2)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=0,kx
do i=1,ix
vf(i, 0, k, 3) = vf(i, 1, k, 3) + vf(i  , 0, k+1, 2) - vf(i, 0, k, 2)
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Y_plus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=1,ix
vf(i, jx+1, k, 2) = vf(i, jx, k, 2)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=1,kx
do i=0,ix
vf(i, jx+1, k, 1) = vf(i, jx, k, 1) - vf(i+1, jx, k  , 2) + vf(i, jx, k, 2)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
do k=0,kx
do i=1,ix
vf(i, jx+1, k, 3) = vf(i, jx, k, 3) - vf(i  , jx, k+1, 2) + vf(i, jx, k, 2)
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Z_minus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=1,ix
vf(i, j, -1, 3) = vf(i, j, 0, 3)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=0,ix
vf(i, j, 0, 1) = vf(i, j, 1, 1) + vf(i+1, j  , 0, 3) - vf(i, j, 0, 3)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=0,jx
do i=1,ix
vf(i, j, 0, 2) = vf(i, j, 1, 2) + vf(i  , j+1, 0, 3) - vf(i, j, 0, 3)
end do
end do
!$OMP END DO
!$OMP END PARALLEL


case (Z_plus)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=1,ix
vf(i, j, kx+1, 3) = vf(i, j, kx, 3)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=1,jx
do i=0,ix
vf(i, j, kx+1, 1) = vf(i, j, kx, 1) - vf(i+1, j  , kx, 3) + vf(i, j, kx, 3)
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
do j=0,jx
do i=1,ix
vf(i, j, kx+1, 2) = vf(i, j, kx, 2) - vf(i  , j+1, kx, 3) + vf(i, j, kx, 3)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

case default
end select FACES


return
end subroutine vobc_cf_tfree



!> ********************************************************************
!! @brief 外部指定境界条件による速度の発散の修正
!! @param [in,out] div     速度の発散
!! @param [in]     sz      配列長
!! @param [in]     g       ガイドセル長
!! @param [in]     dh      格子幅
!! @param [in]     m_face  面番号
!! @param [in]     bv      BCindex C
!! @param [in]     vec     指定する速度ベクトル
!! @param [in]     nID     隣接ランク番号（nID[]<0の時外部境界面）
!! @note 固体部分は対象外とするのでループ中に判定あり
!!       部分的な境界条件の実装のため、ガイドセル部のマスク情報を利用
!<
    subroutine vobc_div_drchlt (div, sz, g, dh, m_face, bv, vec, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face, bvx, m_face
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u_bc, v_bc, w_bc, a
    real, dimension(3)                                        ::  vec, dh
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face
    
    u_bc = vec(1) / dh(1)
    v_bc = vec(2) / dh(2)
    w_bc = vec(3) / dh(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP PRIVATE(a)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
        bvx = bv(1, j, k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          a = u_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(0, j, k), State, 1))
          div(1, j, k) = div(1, j, k) - a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
        bvx = bv(ix, j, k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          a = u_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(ix+1, j, k), State, 1))
          div(ix, j, k) = div(ix, j, k) + a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
        bvx = bv(i, 1, k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          a = v_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, 0, k), State, 1))
          div(i, 1, k) = div(i, 1, k) - a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
        bvx = bv(i, jx, k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          a = v_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, jx+1, k), State, 1))
          div(i, jx, k) = div(i, jx, k) + a
        endif
      end do
      end do
!$OMP END DO
    
    
    case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
        bvx = bv(i, j, 1)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          a = w_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, j, 0), State, 1))
          div(i, j, 1) = div(i, j, 1) - a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
        bvx = bv(i, j, kx)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          a = w_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, j, kx+1), State, 1))
          div(i, j, kx) = div(i, j, kx) + a
        endif
      end do
      end do
!$OMP END DO
    
    case default
    end select FACES

!$OMP END PARALLEL

    return
    end subroutine vobc_div_drchlt



!> **************************************************************************
!! @brief 対流流出条件 Cell center
!! @param [out] vc      疑似速度ベクトル（セルセンタ）
!! @param [in]  v0      セルセンタ速度 u^n
!! @param [in]  sz      配列長
!! @param [in]  g       ガイドセル長
!! @param [in]  dh      格子幅
!! @param [in]  dt      時間積分幅
!! @param [in]  bv      BCindex C
!! @param [in]  v_cnv   対流流出速度
!! @param [in]  m_face  外部境界の面番号
!! @param [in]  nID     隣接ランク番号（nID[]<0の時外部境界面）
!! @note 壁の場合には、OBC_MASKが設定されていないので、内部壁として認識される
!!       u^nのガイドセルには，u^{n+1}の予測値を入れておく > {}^c \hat{u}_{ix+1}を {}^f \hat{u}_{ix}の計算に使う
!<
subroutine vobc_cc_outflow (vc, v0, sz, g, dh, dt, bv, v_cnv, m_face, nID)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, g, face, ix, jx, kx, bvx, gc, m_face
integer, dimension(3)                                       ::  sz
real                                                        ::  Ue, Uw, Un, Us, Ut, Ub
real                                                        ::  Ve, Vw, Vn, Vs, Vt, Vb
real                                                        ::  We, Ww, Wn, Ws, Wt, Wb
real                                                        ::  dt, v_cnv, up, vp, wp, cx, cy, cz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vc, v0
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
integer, dimension(0:5)                                     ::  nID
real, dimension(3)                                          ::  dh

if ( nID(m_face) >= 0 ) return

ix = sz(1)
jx = sz(2)
kx = sz(3)
gc = g
face = m_face

cx = v_cnv * dt / dh(1)
cy = v_cnv * dt / dh(2)
cz = v_cnv * dt / dh(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, gc, cx, cy, cz, face) &
!$OMP PRIVATE(bvx, up, vp, wp, i, j, k)

FACES : select case (face)

case (X_minus)

if ( cx>0.0 ) cx=0.0


!$OMP DO SCHEDULE(static) PRIVATE(Uw, Vw, Ww) COLLAPSE(2)
do k=1,kx
do j=1,jx
bvx = bv(1, j, k)
if ( (ibits(bvx, State, 1) == id_fluid) .and. (ibits(bvx, bc_face_W, bitw_5) /= obc_mask) ) then

  Uw = v0(0, j  ,k  ,1)
  Vw = v0(0, j  ,k  ,2)
  Ww = v0(0, j  ,k  ,3)

  up = Uw - cx*(v0(1  ,j  ,k  ,1)-Uw)
  vp = Vw - cx*(v0(1  ,j  ,k  ,2)-Vw)
  wp = Ww - cx*(v0(1  ,j  ,k  ,3)-Ww)

  do i=1-gc, 0
    vc(i, j, k, 1) = up
    vc(i, j, k, 2) = vp
    vc(i, j, k, 3) = wp
    v0(i, j, k, 1) = up
    v0(i, j, k, 2) = vp
    v0(i, j, k, 3) = wp
  end do

endif
end do
end do
!$OMP END DO


case (X_plus)

if ( cx<0.0 ) cx=0.0


!$OMP DO SCHEDULE(static) PRIVATE(Ue, Ve, We) COLLAPSE(2)
do k=1,kx
do j=1,jx
bvx = bv(ix, j, k)
if ( (ibits(bvx, State, 1) == id_fluid) .and. (ibits(bvx, bc_face_E, bitw_5) /= obc_mask) ) then

  Ue = v0(ix+1,j  ,k  ,1)
  Ve = v0(ix+1,j  ,k  ,2)
  We = v0(ix+1,j  ,k  ,3)

  up = Ue - cx*(Ue-v0(ix  ,j  ,k  ,1))
  vp = Ve - cx*(Ve-v0(ix  ,j  ,k  ,2))
  wp = We - cx*(We-v0(ix  ,j  ,k  ,3))

  do i=ix+1, ix+gc
    vc(i, j, k, 1) = up
    vc(i, j, k, 2) = vp
    vc(i, j, k, 3) = wp
    v0(i, j, k, 1) = up
    v0(i, j, k, 2) = vp
    v0(i, j, k, 3) = wp
  end do

endif
end do
end do
!$OMP END DO


case (Y_minus)

if ( cy>0.0 ) cy=0.0


!$OMP DO SCHEDULE(static) PRIVATE(Us, Vs, Ws) COLLAPSE(2)
do k=1,kx
do i=1,ix
bvx = bv(i, 1, k)
if ( (ibits(bvx, State, 1) == id_fluid) .and. (ibits(bvx, bc_face_S, bitw_5) /= obc_mask) ) then

  Us = v0(i, 0, k, 1)
  Vs = v0(i, 0, k, 2)
  Ws = v0(i, 0, k, 3)

  up = Us - cy*(v0(i, 1, k, 1)-Us)
  vp = Vs - cy*(v0(i, 1, k, 2)-Vs)
  wp = Ws - cy*(v0(i, 1, k, 3)-Ws)

  do j=1-gc, 0
    vc(i, j, k, 1) = up
    vc(i, j, k, 2) = vp
    vc(i, j, k, 3) = wp
    v0(i, j, k, 1) = up
    v0(i, j, k, 2) = vp
    v0(i, j, k, 3) = wp
  end do

endif
end do
end do
!$OMP END DO


case (Y_plus)

if ( cy<0.0 ) cy=0.0


!$OMP DO SCHEDULE(static) PRIVATE(Un, Vn, Wn) COLLAPSE(2)
do k=1,kx
do i=1,ix
bvx = bv(i, jx, k)
if ( (ibits(bvx, State, 1) == id_fluid) .and. (ibits(bvx, bc_face_N, bitw_5) /= obc_mask) ) then

  Un = v0(i  ,jx+1,k  ,1)
  Vn = v0(i  ,jx+1,k  ,2)
  Wn = v0(i  ,jx+1,k  ,3)

  up = Un - cy*(Un-v0(i  ,jx  ,k  ,1))
  vp = Vn - cy*(Vn-v0(i  ,jx  ,k  ,2))
  wp = Wn - cy*(Wn-v0(i  ,jx  ,k  ,3))

  do j=jx+1, jx+gc
    vc(i, j, k, 1) = up
    vc(i, j, k, 2) = vp
    vc(i, j, k, 3) = wp
    v0(i, j, k, 1) = up
    v0(i, j, k, 2) = vp
    v0(i, j, k, 3) = wp
  end do

endif
end do
end do
!$OMP END DO


case (Z_minus)

if ( cz>0.0 ) cz=0.0


!$OMP DO SCHEDULE(static) PRIVATE(Ub, Vb, Wb) COLLAPSE(2)
do j=1,jx
do i=1,ix
bvx = bv(i, j, 1)
if ( (ibits(bvx, State, 1) == id_fluid) .and. (ibits(bvx, bc_face_B, bitw_5) /= obc_mask) ) then
  Ub = v0(i  ,j  , 0,1)
  Vb = v0(i  ,j  , 0,2)
  Wb = v0(i  ,j  , 0,3)

  up = Ub - cz*(v0(i  ,j  ,1  ,1)-Ub)
  vp = Vb - cz*(v0(i  ,j  ,1  ,2)-Vb)
  wp = Wb - cz*(v0(i  ,j  ,1  ,3)-Wb)

  do k=1-gc,0
    vc(i, j, k, 1) = up
    vc(i, j, k, 2) = vp
    vc(i, j, k, 3) = wp
    v0(i, j, k, 1) = up
    v0(i, j, k, 2) = vp
    v0(i, j, k, 3) = wp
  end do

endif
end do
end do
!$OMP END DO


case (Z_plus)

if ( cz<0.0 ) cz=0.0


!$OMP DO SCHEDULE(static) PRIVATE(Ut, Vt, Wt) COLLAPSE(2)
do j=1,jx
do i=1,ix
bvx = bv(i, j, kx)
if ( (ibits(bvx, State, 1) == id_fluid) .and. (ibits(bvx, bc_face_T, bitw_5) /= obc_mask) ) then
  Ut = v0(i  ,j  ,kx+1,1)
  Vt = v0(i  ,j  ,kx+1,2)
  Wt = v0(i  ,j  ,kx+1,3)

  up = Ut - cz*(Ut-v0(i  ,j  ,kx  ,1))
  vp = Vt - cz*(Vt-v0(i  ,j  ,kx  ,2))
  wp = Wt - cz*(Wt-v0(i  ,j  ,kx  ,3))

  do k=kx+1,kx+gc
    vc(i, j, k, 1) = up
    vc(i, j, k, 2) = vp
    vc(i, j, k, 3) = wp
    v0(i, j, k, 1) = up
    v0(i, j, k, 2) = vp
    v0(i, j, k, 3) = wp
  end do

endif
end do
end do
!$OMP END DO

case default
end select FACES

!$OMP END PARALLEL


return
end subroutine vobc_cc_outflow
