!   *****************************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
!
!   *****************************************************************
!
!> @file FB_util.f90
!! @brief FlowBase utilities
!! @author keno, FSI Team, VCAD, RIKEN
!<

!  ******************************************************
!> @subroutine fb_totalp (tp, sz, g, v, p, v00, gs, flop)
!! @brief 全圧を計算する
!! @param tp 全圧
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v 速度ベクトル
!! @param p 圧力
!! @param v00 参照速度
!! @param gs 格子配置 (0-staggered, 1-collocate)
!! @param flop 浮動小数演算数
!<
    subroutine fb_totalp (tp, sz, g, v, p, v00, gs, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g, gs
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u1, u2, u3, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, tp
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    scheme : select case (gs)
      case(0)  ! Staggered
        do k=1,kx
        do j=1,jx
        do i=1,ix
          u1 = 0.5*(v(i,j,k,1) + v(i-1,j  ,k  ,1))-v00(1)
          u2 = 0.5*(v(i,j,k,2) + v(i  ,j-1,k  ,2))-v00(2)
          u3 = 0.5*(v(i,j,k,3) + v(i  ,j  ,k-1,3))-v00(3)
          tp(i,j,k) = 0.5*(u1*u1 + u2*u2 + u3*u3) + p(i,j,k)
        end do
        end do
        end do
        flop = flop + real(ix*jx*kx*13)
        
      case(1)  ! Collocated
        do k=1-g,kx+g
        do j=1-g,jx+g
        do i=1-g,ix+g
          u1 = v(i,j,k,1) - v00(1)
          u2 = v(i,j,k,2) - v00(2)
          u3 = v(i,j,k,3) - v00(3)
          tp(i,j,k) = 0.5*(u1*u1 + u2*u2 + u3*u3) + p(i,j,k)
        end do
        end do
        end do
        flop = flop + real(ix*jx*kx*10)
        
      case default
    end select scheme
    
    return
    end subroutine fb_totalp

!  *****************************************************
!> @subroutine fb_minmax_v (v_min, v_max, sz, g, v00, v)
!! @brief 速度の最小値と最大値を計算する
!! @param v_min 最小値
!! @param v_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param v 速度ベクトル
!<
    subroutine fb_minmax_v (v_min, v_max, sz, g, v00, v)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  v_min, v_max, u1, u2, u3, uu
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    v_min =  1.0e6
    v_max = -1.0e6

    do k=1,kx
    do j=1,jx
    do i=1,ix
			u1 = v(i,j,k,1)-v00(1)
			u2 = v(i,j,k,2)-v00(2)
			u3 = v(i,j,k,3)-v00(3)
			uu = sqrt( u1*u1 + u2*u2 + u3*u3 )
      v_min = min(v_min, uu)
      v_max = max(v_max, uu)
    end do
    end do
    end do

    return
    end subroutine fb_minmax_v

!  ************************************************
!> @subroutine fb_minmax_s (f_min, f_max, sz, g, s)
!! @brief スカラ値の最小値と最大値を計算する
!! @param f_min 最小値
!! @param f_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param s スカラ値
!<
    subroutine fb_minmax_s (f_min, f_max, sz, g, s)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  f_min, f_max
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    f_min =  1.0e6
    f_max = -1.0e6

    do k=1,kx
    do j=1,jx
    do i=1,ix
      f_min = min(f_min, s(i,j,k))
      f_max = max(f_max, s(i,j,k))
    end do
    end do
    end do

    return
    end subroutine fb_minmax_s

!  *******************************************
!> @subroutine fb_set_vector (var, sz, g, val)
!! @brief ベクトル値を設定する
!! @param var ベクトル配列
!! @param sz 配列長
!! @param g ガイドセル長
!! @param val ベクトル値
!<
    subroutine fb_set_vector (var, sz, g, val)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real, dimension(3)                                        ::  val
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  var

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    do k=1-g, kx+g
    do j=1-g, jx+g
    do i=1-g, ix+g
        var(i,j,k,1) = val(1)
        var(i,j,k,2) = val(2)
        var(i,j,k,3) = val(3)
    end do
    end do
    end do

    return
    end subroutine fb_set_vector

!  **************************************
!> @subroutine fb_limit_scalar (t, sz, g)
!! @brief スカラ値の値を[0, 1]に制限する
!! @param t スカラ値
!! @param sz 配列長
!! @param g ガイドセル長
!<
    subroutine fb_limit_scalar (t, sz, g)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  t

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      if (t(i,j,k)<0.0) t(i,j,k)=0.0
      if (t(i,j,k)>1.0) t(i,j,k)=1.0
    end do
    end do
    end do

    return
    end subroutine fb_limit_scalar

!  **************************************************************************************
!> @subroutine subroutine setdt_cp (q, gd, dx, dy, dz, cfl, dt, dxmin, sz, g, jm, km, lm)
!! @brief 時間積分幅を計算する
!! @param q 保存変数ベクトル
!! @param gd ?
!! @param dx x方向の格子幅
!! @param dy y方向の格子幅
!! @param dz z方向の格子幅
!! @param cfl CFL数
!! @param dt 時間積分幅
!! @param dxmin ?
!! @param sz 配列長
!! @param g ガイドセル長
!! @param jm ?
!! @param km ?
!! @param lm ?
!! @note
!!    - 寺島コードのルーチン
!<
    subroutine setdt_cp (q, gd, dx, dy, dz, cfl, dt, dxmin, sz, g, jm, km, lm)
    implicit none
    integer                                                   ::  jt, kt, lt, g, j, k, l, jm, km, lm
    integer, dimension(3)                                     ::  sz
    real                                                      ::  sigmax, ri, u, v, w, er, sndsp2, sndsp
    real                                                      ::  siga, sigb, sigc, sigabc, gd
    real                                                      ::  dxmin, dymin, dzmin, dtx, dty, dtz, dt, cfl
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 5) ::  q
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  dx, dy, dz
    
    jt     = 1
    kt     = 1
    lt     = 1
    sigmax = 0.d0

    do 100 l = 2, lm
    do 100 k = 2, km
    do 100 j = 2, jm
      ri = 1.d0/q(j,k,l,1)
      u  = q(j,k,l,2)*ri
      v  = q(j,k,l,3)*ri
      w  = q(j,k,l,4)*ri
      er = q(j,k,l,5)*ri

      sndsp2 = gd*(er-0.5d0*(u*u+v*v+w*w))
      sndsp  = sqrt(sndsp2)

      siga   = abs(u) + sndsp
      sigb   = abs(v) + sndsp
      sigc   = abs(w) + sndsp
      sigabc = max(siga, sigb, sigc)

      if (sigabc.le.sigmax) go to 100

      jt     = j
      kt     = k
      lt     = l
      sigmax = sigabc
100 continue
 
    dxmin = 1.d+8
    dymin = 1.d+8
    dzmin = 1.d+8
    
    do l = 2, lm
    do k = 2, km
    do j = 2, jm
        if (dx(j,k,l).lt.dxmin) dxmin = dx(j,k,l)
        if (dy(j,k,l).lt.dymin) dymin = dy(j,k,l)
        if (dz(j,k,l).lt.dzmin) dzmin = dz(j,k,l)
    end do
    end do
    end do

    dtx = dxmin*cfl/sigmax
    dty = dymin*cfl/sigmax
    dtz = dzmin*cfl/sigmax

    dt  = min(dtx, dty, dtz)

    return
    end subroutine setdt_cp
