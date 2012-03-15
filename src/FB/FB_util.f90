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

!  ***************************************************
!> @subroutine fb_delta_v (d, sz, g, vn, vo, bx, flop)
!! @brief 有効セルに対する，1タイムステップ進行時の変化量の2乗和と平均値(RootMean)
!! @param d 戻り値（変化量の2乗和と平均値）
!! @param sz 配列長
!! @param g ガイドセル長
!! @param vn ベクトル値 n+1 step
!! @param vo ベクトル値 n step
!! @param bx BCindex
!! @param flop 浮動小数演算数
!<
  subroutine fb_delta_v (d, sz, g, vn, vo, bx, flop)
  implicit none
  include 'cbc_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  flop, actv
  real                                                      ::  u, v, w, av, rm, x, y, z
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  vn, vo
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
  real, dimension(2)                                        ::  d

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + real(ix)*real(jx)*real(kx)*18.0

  av = 0.0
  rm = 0.0

!$OMP PARALLEL &
!$OMP PRIVATE(actv, u, v, w, x, y, z) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:av) &
!$OMP REDUCTION(+:rm)
  do k=1,kx
  do j=1,jx
  do i=1,ix 
    actv = real(ibits(bx(i,j,k), State, 1))
    
    u = vn(1,i,j,k)
    v = vn(2,i,j,k)
    w = vn(3,i,j,k)
    av = av + (u*u + v*v + w*w)*actv
    
    x = u - vo(1,i,j,k)
    y = v - vo(2,i,j,k)
    z = w - vo(3,i,j,k)
    rm = rm + (x*x + y*y + z*z)*actv
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  d(1) = rm
  d(2) = av

  return
  end subroutine fb_delta_v

!  ***************************************************
!> @subroutine fb_delta_s (d, sz, g, sn, so, bx, flop)
!! @brief 有効セルに対する，1タイムステップ進行時の変化量の2乗和と平均値(RootMean)
!! @param d 戻り値（変化量の2乗和と平均値）
!! @param sz 配列長
!! @param g ガイドセル長
!! @param sn スカラー値 n+1 step
!! @param so スカラー値 n step
!! @param bx BCindex
!! @param flop 浮動小数演算数
!<
  subroutine fb_delta_s (d, sz, g, sn, so, bx, flop)
  implicit none
  include 'cbc_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  flop, actv
  real                                                      ::  a, s, av, rm
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  sn, so
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
  real, dimension(2)                                        ::  d

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + real(ix)*real(jx)*real(kx)*7.0

  av = 0.0
  rm = 0.0

!$OMP PARALLEL &
!$OMP PRIVATE(actv, s, a) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:av) &
!$OMP REDUCTION(+:rm)
  do k=1,kx
  do j=1,jx
  do i=1,ix
    actv = real(ibits(bx(i,j,k), State,  1))
    
    s = sn(i,j,k)
    av = av + s*actv
    
    a = (s - so(i,j,k))*actv
    rm = rm + a*a
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL
  
  d(1) = rm
  d(2) = av

  return
  end subroutine fb_delta_s

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

!$OMP PARALLEL &
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
    avr(1,i,j,k) = avr(1,i,j,k) + v(1,i,j,k)
    avr(2,i,j,k) = avr(2,i,j,k) + v(2,i,j,k)
    avr(3,i,j,k) = avr(3,i,j,k) + v(3,i,j,k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

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

!$OMP PARALLEL &
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
    avr(i,j,k) = avr(i,j,k) + s(i,j,k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_average_s

!  ********************************************
!> @subroutine fb_copy_real_s (dst, src, sz, g)
!! @brief スカラー値をコピーする
!! @param dst コピー先
!! @param src コピー元
!! @param sz 配列長
!! @param g ガイドセル長
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_copy_real_s (dst, src, sz, g)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  dst, src

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
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
      dst(i,j,k) = src(i,j,k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_copy_real_s

!  ********************************************
!> @subroutine fb_copy_real_v (dst, src, sz, g)
!! @brief ベクトル値をコピーする
!! @param dst コピー先
!! @param src コピー元
!! @param sz 配列長
!! @param g ガイドセル長
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_copy_real_v (dst, src, sz, g)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  dst, src

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
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
      dst(1,i,j,k) = src(1,i,j,k)
      dst(2,i,j,k) = src(2,i,j,k)
      dst(3,i,j,k) = src(3,i,j,k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_copy_real_v
    
!  *******************************************
!> @subroutine fb_set_int_s (var, sz, g, init)
!! @brief ベクトル値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長
!! @param g ガイドセル長
!! @param val 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_int_s(var, sz, g, init)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    integer                                                   ::  init
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  var

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, init)

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
      var(i,j,k) = init
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_int_s
    
!  ********************************************
!> @subroutine fb_set_real_s(var, sz, g, init)
!! @brief スカラー値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長
!! @param g guide cell
!! @param val 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_real_s(var, sz, g, init)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  init
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  var

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, init)

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
      var(i,j,k) = init
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_real_s

!  ********************************************
!> @subroutine fb_set_real_v(var, sz, g, init)
!! @brief ベクトル値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長
!! @param g ガイドセル長
!! @param val 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_real_v(var, sz, g, init)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  init
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  var

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, init)

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
      var(1,i,j,k) = init
      var(2,i,j,k) = init
      var(3,i,j,k) = init
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_real_v
    
!  **************************************************
!> @subroutine fb_totalp (tp, sz, g, v, p, v00, flop)
!! @brief 全圧を計算する
!! @param tp 全圧
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v 速度ベクトル
!! @param p 圧力
!! @param v00 参照速度
!! @param flop 浮動小数演算数
!<
    subroutine fb_totalp (tp, sz, g, v, p, v00, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u1, u2, u3, flop, vx, vy, vz
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, tp
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vx = v00(1)
    vy = v00(2)
    vz = v00(3)

!$OMP PARALLEL &
!$OMP PRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

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
      u1 = v(1,i,j,k) - vx
      u2 = v(2,i,j,k) - vy
      u3 = v(3,i,j,k) - vz
      tp(i,j,k) = 0.5*(u1*u1 + u2*u2 + u3*u3) + p(i,j,k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
        
    flop = flop + real(ix)*real(jx)*real(kx)*10.0
    
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
    real                                                      ::  v_min, v_max, u1, u2, u3, uu, flop, vx, vy, vz
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vx = v00(1)
    vy = v00(2)
    vz = v00(3)
    v_min =  1.0e6
    v_max = -1.0e6
    
    ! 10 + sqrt*1 = 20 ! DP 30
    flop = flop + real(ix)*real(jx)*real(kx)*20.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*30.0 ! DP

!$OMP PARALLEL &
!$OMP PRIVATE(u1, u2, u3, uu) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(min:v_min) &
!$OMP REDUCTION(max:v_max)
    do k=1,kx
    do j=1,jx
    do i=1,ix
			u1 = v(1,i,j,k)-vx
			u2 = v(2,i,j,k)-vy
			u3 = v(3,i,j,k)-vz
			uu = sqrt( u1*u1 + u2*u2 + u3*u3 )
      v_min = min(v_min, uu)
      v_max = max(v_max, uu)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

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

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(min:f_min) &
!$OMP REDUCTION(max:f_max)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      f_min = min(f_min, s(i,j,k))
      f_max = max(f_max, s(i,j,k))
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

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
    real                                                      ::  u1, u2, u3
    real, dimension(3)                                        ::  val
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  var

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u1 = val(1)
    u2 = val(2)
    u3 = val(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=1-g, kx+g
    do j=1-g, jx+g
    do i=1-g, ix+g
        var(1,i,j,k) = u1
        var(2,i,j,k) = u2
        var(3,i,j,k) = u3
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

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
    real                                                      ::  tmp
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  t

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(tmp) &
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
      tmp = t(i,j,k)
      if (tmp<0.0) t(i,j,k)=0.0
      if (tmp>1.0) t(i,j,k)=1.0
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_limit_scalar
