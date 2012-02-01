!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************
!
!> @file cds_poisson.f90
!> @brief Poisson routines for CDS
!> @author keno, FSI Team, VCAD, RIKEN
    
!  ***********************************************************************
!> @subroutine cds_psor (p, sz, g, omg, res, div, bnd, cut, e, para_key)
!! @brief Poissonから導かれた連立一次方程式をSOR法で解く
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param res 残差
!! @param div ソース項
!! @param bnd マスク
!! @param cut カット情報
!! @param e マシンイプシロン
!! @param parakey parallel managerの識別ID
!! @todo
!!    - 2階微分の距離の考慮
!<
    subroutine cds_psor (p, sz, g, omg, res, div, bnd, cut, e, para_key)
    implicit none
    include 'sklparaf.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, dtype
    integer, dimension(3)                                     ::  sz
    real                                                      ::  dp, r1, e, ee, omg, cpd, res, s0, tmp
    real                                                      ::  m1, m2, m3, m4, m5, m6
    real                                                      ::  d1, d2, d3, d4, d5, d6
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, div, bnd
    real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    integer                                                   ::  ierr, iret, para_key

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ee = 1.0+e
    iret = 0
    ierr = 0
    
    r1  = 0.0
    res = 0.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      d1 = cut(1,i,j,k) ! Cut(i)   X+
      d2 = cut(2,i,j,k) ! Cut(i)   X-
      d3 = cut(3,i,j,k) ! Cut(j)   Y+
      d4 = cut(4,i,j,k) ! Cut(j)   Y-
      d5 = cut(5,i,j,k) ! Cut(k)   Z+
      d6 = cut(6,i,j,k) ! Cut(k)   Z-
      
      ! セル内にカットがあればゼロ，それ以外は1
      m1 = aint(d1 + 0.5)
      m2 = aint(d2 + 0.5)
      m3 = aint(d3 + 0.5)
      m4 = aint(d4 + 0.5)
      m5 = aint(d5 + 0.5)
      m6 = aint(d6 + 0.5)
      
      cpd=ee/(e+m1+m2+m3+m4+m5+m6)
      s0 = m1*p(i+1,j  ,k  )+m2*p(i-1,j  ,k  )   &
         + m3*p(i  ,j+1,k  )+m4*p(i  ,j-1,k  )   &
         + m5*p(i  ,j  ,k+1)+m6*p(i  ,j  ,k-1) + div(i,j,k)
      dp = (cpd*s0-p(i,j,k))*bnd(i,j,k)
      p(i,j,k)=p(i,j,k) + omg*dp
      r1 = r1 + dp*dp
    end do
    end do
    end do
    
    res = r1
    call SklIsParallel(iret)
    if ( iret == 1 ) then
      tmp = r1
      call SklAllreduce(tmp, res, 1, SKL_REAL, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
    end if

    return
    end subroutine cds_psor

!  **************************************************************************
!> @subroutine cds_div_cf (div, sz, g, dh, dt, b2, vf, bnd, cut, v00, mode)
!! @brief セルフェイスの値を使って，速度の発散を計算する
!! @param div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param b2 定数ベクトルのノルム
!! @param vf 速度ベクトル
!! @param bnd マスク
!! @param cut カット情報
!! @param v00 参照速度
!! @param mode ノルムを計算するかどうか (0-しない　1-する)
!! @date 2009/03/23
!! @note
!!    - Poissonのソース項で計算できるように-dh/dtの係数をかけてある点に注意
!!    - セルフェイスの値は，評価セルの流体側からの情報を用いて内挿
!<
    subroutine cds_div_cf (div, sz, g, dh, dt, b2, vf, bnd, cut, v00, mode)
    implicit none
    integer                                                     ::  i, j, k, ix, jx, kx, g, mode
    integer, dimension(3)                                       ::  sz
    real                                                        ::  dh, dt, dth, b2, dv
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0
    real                                                        ::  d1, d2, d3, d4, d5, d6
    real                                                        ::  h1, h2, h3, h4, h5, h6
    real                                                        ::  u_ref, v_ref, w_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  bnd, div
		real, dimension(0:3)                                        ::  v00
    real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  cut

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    dth= -dh/dt
    b2 = 0.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      d1 = cut(1,i,j,k) ! Cut(i)   X+
      d2 = cut(2,i,j,k) ! Cut(i)   X-
      d3 = cut(3,i,j,k) ! Cut(j)   Y+
      d4 = cut(4,i,j,k) ! Cut(j)   Y-
      d5 = cut(5,i,j,k) ! Cut(k)   Z+
      d6 = cut(6,i,j,k) ! Cut(k)   Z-
      
      h1 = d1 + 0.5
      h2 = d2 + 0.5
      h3 = d3 + 0.5
      h4 = d4 + 0.5
      h5 = d5 + 0.5
      h6 = d6 + 0.5
      
      Ue = vf(i  , j  , k  , 1)
      Uw = vf(i-1, j  , k  , 1)
      Vn = vf(i  , j  , k  , 2)
      Vs = vf(i  , j-1, k  , 2)
      Wt = vf(i  , j  , k  , 3)
      Wb = vf(i  , j  , k-1, 3)
      
      Ue0 = Ue
      Uw0 = Uw
      Vn0 = Vn
      Vs0 = Vs
      Wt0 = Wt
      Wb0 = Wb
      
      if ( d1 < 1.0 ) Ue = (u_ref + (h1-1.0)*Uw0 )/h1
      if ( d2 < 1.0 ) Uw = (u_ref + (h2-1.0)*Ue0 )/h2
      if ( d3 < 1.0 ) Vn = (v_ref + (h3-1.0)*Vs0 )/h3
      if ( d4 < 1.0 ) Vs = (v_ref + (h4-1.0)*Vn0 )/h4
      if ( d5 < 1.0 ) Wt = (w_ref + (h5-1.0)*Wb0 )/h5
      if ( d6 < 1.0 ) Wb = (w_ref + (h6-1.0)*Wt0 )/h6
      
      dv =( Ue - Uw + Vn - Vs + Wt - Wb ) * dth * bnd(i,j,k)
      div(i,j,k) = dv
      b2 = b2 + dv*dv
    end do
    end do
    end do

    return
    end subroutine cds_div_cf
