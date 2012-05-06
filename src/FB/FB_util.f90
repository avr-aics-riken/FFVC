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

!  **********************************************
!> @subroutine fb_interp_rough_s(dst, sz, g, src)
!! @brief 粗い格子から密な格子への補間
!! @param dst 密な格子系
!! @param sz 配列長
!! @param g ガイドセル長
!! @param src 粗い格子系
!<
  subroutine fb_interp_rough_s(dst, sz, g, src)
  implicit none
  integer                                                      ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                        ::  sz
  real                                                         ::  g_x, g_y, g_z
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)       ::  dst
  real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g) ::  src


  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, scale, rr)

!$OMP DO SCHEDULE(static)

  do k=1, kx/2
    kk = k*2
  do j=1, jx/2
    jj = j*2
  do i=1, ix/2
    ii = i*2

    sp = src(i  , j  , k  )
    sw = src(i-1, j  , k  )
    se = src(i+1, j  , k  )
    ss = src(i  , j-1, k  )
    sn = src(i  , j+1, k  )
    sb = src(i  , j  , k-1)
    st = src(i  , j  , k+1)

    g_x = ( se - sw ) * 0.5
    g_y = ( sn - ss ) * 0.5
    g_z = ( st - sb ) * 0.5

    g_xx = se - 2.0 * sp + sw
    g_yy = sn - 2.0 * sp + ss
    g_zz = st - 2.0 * sp + sb

    g_xy = ( src(i+1, j+1, k  ) - src(i+1, j-1, k  )   &
           - src(i-1, j+1, k  ) + src(i-1, j-1, k  ) ) * 0.25
    g_yz = ( src(i  , j+1, k+1) - src(i  , j+1, k-1)   &
           - src(i  , j-1, k+1) + src(i  , j-1, k-1) ) * 0.25
    g_zx = ( src(i+1, j  , k+1) - src(i-1, j  , k+1)   &
           - src(i+1, j  , k-1) + src(i-1, j  , k-1) ) * 0.25

    dst(ii-1, jj-1, kk-1) = sp 
    dst(ii  , jj  , kk  ) = sp 

  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_interp_rough_s

!  **********************************************
!> @subroutine fb_interp_rough_v(dst, sz, g, src)
!! @brief 粗い格子から密な格子への補間
!! @param dst 密な格子系
!! @param sz 配列長
!! @param g ガイドセル長
!! @param src 粗い格子系
!<
  subroutine fb_interp_rough_s(dst, sz, g, src)
  implicit none
  integer                                                         ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                           ::  sz
  real                                                            ::  gu, gv, gw
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)       ::  dst
  real, dimension(3, 1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g) ::  src


  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, scale, rr)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    gu = 0.5 * ( src(1,i+1,j,k) - src(1,i-1,j,k) )
    dst(1,i,j,k) = ( src(1,i,j,k) * rr + u_ref ) * scale
    dst(2,i,j,k) = ( src(2,i,j,k) * rr + v_ref ) * scale
    dst(3,i,j,k) = ( src(3,i,j,k) * rr + w_ref ) * scale
  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_interp_rough_v

!  ****************************************************************************
!> @subroutine fb_tmp_nd2d (dst, src, sz, Base_tmp, Diff_tmp, klv, scale, flop)
!! @brief 温度値を無次元から有次元へ変換し，scale倍して出力
!! @param dst 有次元
!! @param src 無次元
!! @param sz 配列長（一次元）
!! @param Base_tmp 基準温度(K or C)
!! @param Diff_tmp 代表温度差(K or C)
!! @param klv 絶対温度への変換（入力Cのときklv=273.15, Kのときklv=0）
!! @param scale 倍数（瞬時値のとき1.0）
!! @param flop 浮動小数演算数
!<
  subroutine fb_tmp_nd2d (dst, src, sz, Base_tmp, Diff_tmp, klv, scale, flop)
  implicit none
  integer                                                   ::  i, sz
  real                                                      ::  flop, dp, scale
  real                                                      ::  Base_tmp, Diff_tmp, klv
  real, dimension(sz)                                       ::  dst, src

  dp = scale * abs(Diff_tmp)
  flop = flop + real(sz) * 3.0 + 2.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz, dp, Base_tmp, klv)
   
!$OMP DO SCHEDULE(static)

  do i=1,sz
    dst(i) = src(i) * dp - klv + Base_tmp
  end do
!$OMP END DO

!$OMP END PARALLEL

  return
  end subroutine fb_tmp_nd2d
  
!  ****************************************************************************
!> @subroutine fb_tmp_d2nd (dst, src, sz, Base_tmp, Diff_tmp, klv, scale, flop)
!! @brief 温度値を有次元から無次元へ変換し，scale倍して出力
!! @param dst 無次元
!! @param src 有次元
!! @param sz 配列長（一次元）
!! @param Base_tmp 基準温度(K or C)
!! @param Diff_tmp 代表温度差(K or C)
!! @param klv 絶対温度への変換（入力Cのときklv=273.15, Kのときklv=0）
!! @param scale 倍数（瞬時値のとき1.0）
!! @param flop 浮動小数演算数
!<
  subroutine fb_tmp_d2nd (dst, src, sz, Base_tmp, Diff_tmp, klv, scale, flop)
  implicit none
  integer                                                   ::  i, sz
  real                                                      ::  flop, dp, scale
  real                                                      ::  Base_tmp, Diff_tmp, klv
  real, dimension(sz)                                       ::  dst, src

  dp = scale / abs(Diff_tmp)
  flop = flop + real(sz) * 3.0 + 10.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz, dp, Base_tmp, klv)
   
!$OMP DO SCHEDULE(static)

  do i=1,sz
    dst(i) = ( src(i) + klv - Base_tmp ) * dp
  end do
!$OMP END DO

!$OMP END PARALLEL

  return
  end subroutine fb_tmp_d2nd
  
!  *****************************************************************************
!> @subroutine fb_prs_d2nd (dst, src, sz, Base_prs, Ref_rho, Ref_v, scale, flop)
!! @brief 圧力値を有次元から無次元へ変換し，scale倍して出力
!! @param dst 無次元
!! @param src 有次元
!! @param sz 配列長（一次元）
!! @param Base_prs 基準圧力(Pa) 基準圧がゼロのとき，ゲージ圧
!! @param Ref_rho 代表密度(kg/m^3)
!! @param Ref_v 代表速度(m/s)
!! @param scale 倍数（瞬時値のとき1.0）
!! @param flop 浮動小数演算数
!<
  subroutine fb_prs_d2nd (dst, src, sz, Base_prs, Ref_rho, Ref_v, scale, flop)
  implicit none
  integer                                                   ::  i, sz
  real                                                      ::  flop, dp, scale
  real                                                      ::  Base_prs, Ref_rho, Ref_v
  real, dimension(sz)                                       ::  dst, src

  dp = scale / (Ref_rho * Ref_v * Ref_v)
  flop = flop + real(sz) * 3.0 + 10.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz, dp, Base_prs)
   
!$OMP DO SCHEDULE(static)

  do i=1,sz
    dst(i) = ( src(i) - Base_prs ) * dp
  end do
!$OMP END DO

!$OMP END PARALLEL

  return
  end subroutine fb_prs_d2nd
  
!  *****************************************************************************
!> @subroutine fb_prs_nd2d (dst, src, sz, Base_prs, Ref_rho, Ref_v, scale, flop)
!! @brief 圧力値を無次元から有次元へ変換し，scale倍して出力
!! @param dst 有次元
!! @param src 無次元
!! @param sz 配列長（一次元）
!! @param Base_prs 基準圧力(Pa)
!! @param Ref_rho 代表密度(kg/m^3)
!! @param Ref_v 代表速度(m/s)
!! @param mode 圧力の変換モード（１−絶対圧，2-ゲージ圧）
!! @param scale 倍数（瞬時値のとき1.0）
!! @param flop 浮動小数演算数
!<
  subroutine fb_prs_nd2d (dst, src, sz, Base_prs, Ref_rho, Ref_v, scale, flop)
  implicit none
  integer                                                   ::  i, sz
  real                                                      ::  flop, dp, scale
  real                                                      ::  Base_prs, Ref_rho, Ref_v
  real, dimension(sz)                                       ::  dst, src

  dp = Ref_rho * Ref_v * Ref_v * scale
  flop = flop + real(sz) * 3.0 + 2.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz, dp, Base_prs)

!$OMP DO SCHEDULE(static)

  do i=1,sz
    dst(i) = ( src(i) * dp + Base_prs ) 
  end do
!$OMP END DO

!$OMP END PARALLEL

  return
  end subroutine fb_prs_nd2d
  
!  ************************************************
!> @subroutine fb_xcopy (dst, src, sz, scale, flop)
!! @brief 値をscale倍してコピーする
!! @param dst 出力
!! @param src 入力
!! @param sz 配列長（一次元）
!! @param scale 倍数
!! @param flop 浮動小数演算数
!<
  subroutine fb_xcopy (dst, src, sz, scale, flop)
  implicit none
  integer                                                   ::  i, sz
  real                                                      ::  flop, scale
  real, dimension(sz)                                       ::  dst, src

  flop = flop + real(sz)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz)

!$OMP DO SCHEDULE(static)

  do i=1,sz
    dst(i) = src(i) * scale
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_xcopy
  
!  **********************************************************************
!> @subroutine fb_shift_refv_in (dst, sz, g, src, v00, scale, refv, flop)
!! @brief 速度ベクトルの格子速度変換
!! @param dst 移動系へ変換されたベクトル（平均場の場合は積算値）
!! @param sz 配列長
!! @param g ガイドセル長
!! @param src 静止系のベクトル
!! @param v00 参照速度
!! @param sacle 倍数　（瞬時値の場合には1）
!! @param refv 代表速度
!! @param flop 浮動小数演算数
!! @note dst[] = ( src[]/refv + v00 ) * scale, 有次元のときrefvは次元速度，無次元のとき1.0
!<
  subroutine fb_shift_refv_in (dst, sz, g, src, v00, scale, refv, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  flop, scale, u_ref, v_ref, w_ref, refv, rr
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  dst, src
  real, dimension(0:3)                                      ::  v00

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  
  u_ref = v00(1)
  v_ref = v00(2)
  w_ref = v00(3)
  
  rr = 1.0/refv
  
  flop = flop + real(ix)*real(jx)*real(kx)*9.0 + 8.0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, scale, rr)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix 
    dst(1,i,j,k) = ( src(1,i,j,k) * rr + u_ref ) * scale
    dst(2,i,j,k) = ( src(2,i,j,k) * rr + v_ref ) * scale
    dst(3,i,j,k) = ( src(3,i,j,k) * rr + w_ref ) * scale
  end do
  end do
  end do
  
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_shift_refv_in
  
!  *************************************************************************
!> @subroutine fb_shift_refv_out (dst, sz, g, src, v00, scale, unit_v, flop)
!! @brief 速度ベクトルの格子速度変換をして，scale倍する
!! @param dst 静止系へ変換されたベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param src 移動格子系のベクトル（平均場の場合は積算値）
!! @param v00 参照速度
!! @param scale 倍数（瞬時値の場合には1）
!! @param unit_v 無次元のとき1.0，有次元のとき代表速度(m/s)
!! @param flop 浮動小数演算数
!! @note dst[] = ( src[] * stepAvr ) - v00
!<
  subroutine fb_shift_refv_out (dst, sz, g, src, v00, scale, unit_v, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  flop, u_ref, v_ref, w_ref, unit_v, scale
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  dst, src
  real, dimension(0:3)                                      ::  v00

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  
  u_ref = v00(1)
  v_ref = v00(2)
  w_ref = v00(3)
  
  flop = flop + real(ix)*real(jx)*real(kx)*9.0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, scale, unit_v)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix 
    dst(1,i,j,k) = ( src(1,i,j,k) * scale - u_ref ) * unit_v
    dst(2,i,j,k) = ( src(2,i,j,k) * scale - v_ref ) * unit_v
    dst(3,i,j,k) = ( src(3,i,j,k) * scale - w_ref ) * unit_v
  end do
  end do
  end do
  
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_shift_refv_out
  
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

!$OMP DO SCHEDULE(static) &
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
  !if ((i.ge.150).and.(i.le.202).and.(j.eq.26).and.(k.eq.37)) write(*,*) i,j,k, u, v, w
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

!$OMP DO SCHEDULE(static) &

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
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v, avr

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + real(ix)*real(jx)*real(kx)*3.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

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
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s, avr

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + real(ix)*real(jx)*real(kx)*1.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

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

!  *******************************************
!> @subroutine fb_average (avr, src, sz, flop)
!! @brief 値を加算する
!! @param avr 加算値
!! @param src 元の値
!! @param sz 配列長（一次元）
!! @param flop 浮動小数演算数
!<
  subroutine fb_average (avr, src, sz, flop)
  implicit none
  integer                                                   ::  i, sz
  real                                                      ::  flop
  real, dimension(sz)                                       ::  src, avr

  flop = flop + real(sz)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz)

!$OMP DO SCHEDULE(static)

  do i=1,sz
    avr(i) = avr(i) + src(i)
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_average
  
!  ***************************************
!> @subroutine fb_copy_real (dst, src, sz)
!! @brief 値を一次元インデクスでコピーする
!! @param dst コピー先
!! @param src コピー元
!! @param sz 配列長（一次元）
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_copy_real (dst, src, sz)
    implicit none
    integer                                                   ::  i, sz
    real, dimension(sz)                                       ::  dst, src

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz)

!$OMP DO SCHEDULE(static)

    do i=1,sz
      dst(i) = src(i)
    end do

!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_copy_real

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

!$OMP DO SCHEDULE(static)

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

!$OMP DO SCHEDULE(static)

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

!  **************************************
!> @subroutine fb_set_int (var, sz, init)
!! @brief ベクトル値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長（一次元）
!! @param val 初期値
!<
    subroutine fb_set_int (var, sz, init)
    implicit none
    integer                                                   ::  i, sz, init
    integer, dimension(sz)                                    ::  var

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz, init)

!$OMP DO SCHEDULE(static)

    do i=1,sz
      var(i) = init
    end do

!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_int
    
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

!$OMP DO SCHEDULE(static)

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

!  ***************************************
!> @subroutine fb_set_real (var, sz, init)
!! @brief 値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長（一次元）
!! @param val 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_real (var, sz, init)
    implicit none
    integer                                                   ::  i, sz
    real                                                      ::  init
    real, dimension(sz)                                       ::  var

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(sz, init)

!$OMP DO SCHEDULE(static)

    do i=1,sz
      var(i) = init
    end do

!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_real
    
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

!$OMP DO SCHEDULE(static)

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

!  ******************************************
!> @subroutine fb_set_real_v(var, sz, g, vec)
!! @brief ベクトル値を設定する
!! @param var 配列の先頭ポインタ
!! @param sz 配列長
!! @param g ガイドセル長
!! @param vec 初期値
!! @note realはコンパイラオプションでdouble precision にもなる
!<
    subroutine fb_set_real_v(var, sz, g, vec)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u_ref, v_ref, w_ref
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  var
    real, dimension(3)                                        ::  vec

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    u_ref = vec(1)
    v_ref = vec(2)
    w_ref = vec(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      var(1,i,j,k) = u_ref
      var(2,i,j,k) = v_ref
      var(3,i,j,k) = w_ref
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

!$OMP DO SCHEDULE(static)

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

!$OMP DO SCHEDULE(static) &

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

!$OMP DO SCHEDULE(static) &
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

!$OMP DO SCHEDULE(static)

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

!$OMP DO SCHEDULE(static)

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
