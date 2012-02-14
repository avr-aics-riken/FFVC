!   *****************************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *****************************************************************
!
!> @file FB_util.f90
!! @brief FlowBase utilities
!! @author keno, FSI Team, VCAD, RIKEN
!<

!  ***********************************************************
!> @subroutine fb_delta_v (v_min, v_max, sz, g, v00, v, flop)
!! @brief 速度の最小値と最大値を計算する
!! @param v_min 最小値
!! @param v_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param v 速度ベクトル
!! @param flop 浮動小数演算数
!<
  subroutine fb_delta_v (v_min, v_max, sz, g, v00, v)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  v_min, v_max, u1, u2, u3, uu, flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
  real, dimension(0:3)                                      ::  v00

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  v_min =  1.0e6
  v_max = -1.0e6

  ! 10 + sqrt*1 = 20 ! DP 30
  flop = flop + real(ix)*real(jx)*real(kx)*20.0
  ! flop = flop + real(ix)*real(jx)*real(kx)*30.0 ! DP

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
  end subroutine fb_delta_v

!  **********************************************
!> @subroutine fb_average_v (avr, sz, g, v, flop)
!! @brief スカラ値を加算する
!! @param avr 平均値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v ベクトル値
!! @param flop 浮動小数演算数
!<
  subroutine fb_average_v (avr, sz, g, v, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  flop
  real, dimension(3, sz(1), sz(2), sz(3))                   ::  avr
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + real(ix)*real(jx)*real(kx)*3.0

  do k=1,kx
  do j=1,jx
  do i=1,ix
    avr(1,i,j,k) = avr(1,i,j,k) + v(1,i,j,k)
    avr(2,i,j,k) = avr(2,i,j,k) + v(2,i,j,k)
    avr(3,i,j,k) = avr(3,i,j,k) + v(3,i,j,k)
  end do
  end do
  end do

  return
  end subroutine fb_average_v

!  **********************************************
!> @subroutine fb_average_s (avr, sz, g, s, flop)
!! @brief スカラ値を加算する
!! @param avr 平均値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param s スカラ値
!! @param flop 浮動小数演算数
!<
  subroutine fb_average_s (avr, sz, g, s, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  flop
  real, dimension(sz(1), sz(2), sz(3))                      ::  avr
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + real(ix)*real(jx)*real(kx)*1.0

  do k=1,kx
  do j=1,jx
  do i=1,ix
    avr(i,j,k) = avr(i,j,k) + s(i,j,k)
  end do
  end do
  end do

  return
  end subroutine fb_average_s

!  ***************************************
!> @subroutine fb_copy_real (dst, src, sz)
!! @brief ベクトル値を設定する
!! @param dst コピー先
!! @param src コピー元
!! @param sz 配列長（一次元）
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_copy_real (dst, src, sz)
    implicit none
    integer                                                   ::  i, ix, sz
    real, dimension(sz)                                       ::  dst, src

    ix = sz

    do i=1, ix
      dst(i) = src(i)
    end do

    return
    end subroutine fb_copy_real
    
!  ********************************************
!> @subroutine fb_set_value_int (var, sz, init)
!! @brief ベクトル値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長（一次元）
!! @param val 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_value_int (var, sz, init)
    implicit none
    integer                                                   ::  i, ix, sz
    integer                                                   ::  init
    integer, dimension(sz)                                    ::  var

    ix = sz

    do i=1, ix
      var(i) = init
    end do

    return
    end subroutine fb_set_value_int
    
!  *********************************************
!> @subroutine fb_set_value_real (var, sz, init)
!! @brief ベクトル値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長（一次元）
!! @param val 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_value_real (var, sz, init)
    implicit none
    integer                                                   ::  i, ix, sz
    real                                                      ::  init
    real, dimension(sz)                                       ::  var

    ix = sz

    do i=1, ix
      var(i) = init
    end do

    return
    end subroutine fb_set_value_real
    
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
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
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
          u1 = 0.5*(v(1,i,j,k) + v(1,i-1,j  ,k  ))-v00(1)
          u2 = 0.5*(v(2,i,j,k) + v(2,i  ,j-1,k  ))-v00(2)
          u3 = 0.5*(v(3,i,j,k) + v(3,i  ,j  ,k-1))-v00(3)
          tp(i,j,k) = 0.5*(u1*u1 + u2*u2 + u3*u3) + p(i,j,k)
        end do
        end do
        end do
        
        flop = flop + real(ix)*real(jx)*real(kx)*16.0
        
      case(1)  ! Collocated
        do k=1,kx
        do j=1,jx
        do i=1,ix
          u1 = v(1,i,j,k) - v00(1)
          u2 = v(2,i,j,k) - v00(2)
          u3 = v(3,i,j,k) - v00(3)
          tp(i,j,k) = 0.5*(u1*u1 + u2*u2 + u3*u3) + p(i,j,k)
        end do
        end do
        end do
        
        flop = flop + real(ix)*real(jx)*real(kx)*10.0
        
      case default
    end select scheme
    
    return
    end subroutine fb_totalp

!  ***********************************************************
!> @subroutine fb_minmax_v (v_min, v_max, sz, g, v00, v, flop)
!! @brief 速度の最小値と最大値を計算する
!! @param v_min 最小値
!! @param v_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param v 速度ベクトル
!! @param flop 浮動小数演算数
!<
    subroutine fb_minmax_v (v_min, v_max, sz, g, v00, v)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  v_min, v_max, u1, u2, u3, uu, flop
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    v_min =  1.0e6
    v_max = -1.0e6
    
    ! 10 + sqrt*1 = 20 ! DP 30
    flop = flop + real(ix)*real(jx)*real(kx)*20.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*30.0 ! DP

    do k=1,kx
    do j=1,jx
    do i=1,ix
			u1 = v(1,i,j,k)-v00(1)
			u2 = v(2,i,j,k)-v00(2)
			u3 = v(3,i,j,k)-v00(3)
			uu = sqrt( u1*u1 + u2*u2 + u3*u3 )
      v_min = min(v_min, uu)
      v_max = max(v_max, uu)
    end do
    end do
    end do

    return
    end subroutine fb_minmax_v

!  ******************************************************
!> @subroutine fb_minmax_s (f_min, f_max, sz, g, s, flop)
!! @brief スカラ値の最小値と最大値を計算する
!! @param f_min 最小値
!! @param f_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param s スカラ値
!! @param flop 浮動小数演算数
!<
    subroutine fb_minmax_s (f_min, f_max, sz, g, s, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  f_min, f_max, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    f_min =  1.0e6
    f_max = -1.0e6
    
    flop = flop + real(ix)*real(jx)*real(kx)*20.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*2.0 ! DP

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
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  var

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    do k=1-g, kx+g
    do j=1-g, jx+g
    do i=1-g, ix+g
        var(1,i,j,k) = val(1)
        var(2,i,j,k) = val(2)
        var(3,i,j,k) = val(3)
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
