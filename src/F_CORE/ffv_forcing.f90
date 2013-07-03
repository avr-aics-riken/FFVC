!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_forcing.f90
!! @brief  外力タイプの境界条件
!! @author kero
!<

!> ********************************************************************
!! @brief 擬似速度ベクトルの方向の修正、および外力項の付加
!! @param v 速度ベクトル タイムレベルn
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param nv 法線ベクトル
!! @param[out] flop flop count
!<
    subroutine hex_dir (v, sz, g, st, ed, bd, vf, odr, v00, nv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    double precision                                            ::  flop
    real                                                        ::  u_ref, v_ref, w_ref, es, bes
    real                                                        ::  nx, ny, nz, b0, r_bt, uu, u1, u2, u3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  vf
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  nv

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    nx = nv(1)
    ny = nv(2)
    nz = nv(3)

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr) &
!$OMP FIRSTPRIVATE(nx, ny, nz) &
!$OMP PRIVATE(idx, b0, es, bes, u1, u2, u3, r_bt, uu)
!$OMP DO SCHEDULE(static)

    do k=ks,ke
    do j=js,je
    do i=is,ie
      idx = bd(i,j,k)
      b0 = vf(i,j,k)
      es = 0.0

      if ( ibits(idx, 0, bitw_6) == odr ) es=1.0 ! コンポーネントは流体であることが保証されている
      bes = b0*es ! 体積率とコンポーネントの両方でチェック

      ! 擬似速度の方向の強制
      u1 = v(i,j,k,1) - u_ref
      u2 = v(i,j,k,2) - v_ref
      u3 = v(i,j,k,3) - w_ref
      r_bt = 1.0 - bes
      uu = sqrt(u1*u1 + u2*u2 + u3*u3) * bes

      v(i,j,k,1) = (r_bt*u1 + uu*nx) + u_ref
      v(i,j,k,2) = (r_bt*u2 + uu*ny) + v_ref
      v(i,j,k,3) = (r_bt*u3 + uu*nz) + w_ref

    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + dble(ie-is+1)*(je-js+1)*(ke-ks+1)*24.0d0

    return
    end subroutine hex_dir
    
!> ********************************************************************
!! @brief コンポーネントのワーク配列に速度ベクトルを保持
!! @param [out] wk テンポラリのワークベクトル
!! @param [in]  cz コンポーネントの配列長
!! @param [in]  st ループの開始インデクス
!! @param [in]  ed ループの終了インデクス
!! @param [in]  v  速度ベクトル (n+1,k)
!! @param [in]  sz 配列長
!! @param [in]  g  ガイドセル長
!<
    subroutine force_keep_vec (wk, cz, st, ed, v, sz, g)
    implicit none
    integer                                                     ::  i, j, k, g, ii, jj, kk
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, cz, st, ed
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(-1:cz(1)+2, -1:cz(2)+2, -1:cz(3)+2, 3)      ::  wk

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke) &
!$OMP PRIVATE(ii, jj, kk)

! ガイドセルも含めてコピーする
!$OMP DO SCHEDULE(static)

    do k=ks-2,ke+2
      kk = k - ks + 1

    do j=js-2,je+2
      jj = j - js + 1

    do i=is-2,ie+2
      ii = i - is + 1

      wk(ii,jj,kk,1) = v(i,j,k,1)
      wk(ii,jj,kk,2) = v(i,j,k,2)
      wk(ii,jj,kk,3) = v(i,j,k,3)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine force_keep_vec

!> ********************************************************************
!! @brief 圧力損失部におけるPoissonの反復ソース項を計算する
!! @param [in,out] src  反復ソース項 \sum {\beta F}
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     st   ループの開始インデクス
!! @param [in]     ed   ループの終了インデクス
!! @param [in]     bd   BCindex ID
!! @param [in]     vf   コンポーネントの体積率
!! @param [in]     wk   テンポラリのワークベクトル 速度ベクトル
!! @param [in]     cz   コンポーネントの配列長
!! @param [in]     odr  速度境界条件のエントリ
!! @param [in]     v00  参照速度
!! @param [in]     nv   法線ベクトル
!! @param [in]     c    圧力損失部の係数
!! @param [out]    flop flop count
!<
    subroutine hex_psrc (src, sz, g, st, ed, bd, vf, wk, cz, odr, v00, nv, c, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, ii, jj, kk, idx, odr
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed, cz
    double precision                                            ::  flop
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  u_w, u_e, u_s, u_n, u_b, u_t, u_p
    real                                                        ::  v_w, v_e, v_s, v_n, v_b, v_t, v_p
    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t, w_p
    real                                                        ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                        ::  d_w, d_e, d_s, d_n, d_b, d_t, d_p
    real                                                        ::  Fx_w, Fx_e, Fx_p
    real                                                        ::  Fy_s, Fy_n, Fy_p
    real                                                        ::  Fz_b, Fz_t, Fz_p
    real                                                        ::  be, bw, bn, bs, bt, bb, b0, es
    real                                                        ::  re, rw, rn, rs, rt, rb
    real                                                        ::  nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep, pick, qq
    real                                                        ::  q_w, q_e, q_s, q_n, q_b, q_t, q_p
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  vf
    real, dimension(-1:cz(1)+2, -1:cz(2)+2, -1:cz(3)+2, 3)      ::  wk
    real, dimension(6)                                          ::  c
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  nv

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    nx = nv(1)
    ny = nv(2)
    nz = nv(3)
    c1 = c(1)
    c2 = c(2)
    c3 = c(3)
    c4 = c(4)
    ep = c(5)         ! threshold

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, pick, es, ii, jj, kk) &
!$OMP PRIVATE(be, bw, bn, bs, bt, bb, b0) &
!$OMP PRIVATE(u_w, u_e, u_s, u_n, u_b, u_t, u_p) &
!$OMP PRIVATE(v_w, v_e, v_s, v_n, v_b, v_t, v_p) &
!$OMP PRIVATE(w_w, w_e, w_s, w_n, w_b, w_t, w_p) &
!$OMP PRIVATE(g_w, g_e, g_s, g_n, g_b, g_t, g_p) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t, d_p) &
!$OMP PRIVATE(Fx_w, Fx_e, Fx_p) &
!$OMP PRIVATE(Fy_s, Fy_n, Fy_p) &
!$OMP PRIVATE(Fz_b, Fz_t, Fz_p) &
!$OMP PRIVATE(re, rw, rn, rs, rt, rb) &
!$OMP PRIVATE(q_w, q_e, q_s, q_n, q_b, q_t, q_p, qq)

    ! 周囲の1セルを含めてサーチ
!$OMP DO SCHEDULE(static)

    do k=ks-1,ke+1
      kk = k - ks + 1

    do j=js-1,je+1
      jj = j - js + 1

    do i=is-1,ie+1
      ii = i - is + 1

      ! 自セルが流体で、かつ周囲1セルの範囲にコンポーネントが存在する場合のみ有効
      pick = 0.0
      qq = q_p + q_w + q_e + q_s + q_n + q_b + q_t
      if ( qq > 0.0 ) pick = 1.0
      es = pick * real(ibits(idx, State, 1)) ! Fluid=1, Solid=0

      include 'force.h' ! 179 flop

      src(i,j,k) = ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * es ! esはマスク
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    flop = flop + dble(cz(1)+2)*(cz(2)+2)*(cz(3)+2)*200.0d0

    return
    end subroutine hex_psrc

!> ********************************************************************
!! @brief 擬似速度ベクトルの方向の修正、および外力項の付加
!! @param[in,out] vc 擬似速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param v 速度ベクトル タイムレベルn
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param dt 時間積分幅
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param[out] flop flop count
!<
    subroutine hex_force_pvec (vc, sz, g, st, ed, bd, vf, v, odr, v00, dt, nv, c, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    double precision                                            ::  flop
    real                                                        ::  dt, cf
    real                                                        ::  u_ref, v_ref, w_ref, es
    real                                                        ::  u_p, v_p, w_p, g_p, d_p, Fx_p, Fy_p, Fz_p
    real                                                        ::  nx, ny, nz, b0, d_b, r_bt, uu, u1, u2, u3
    real                                                        ::  c1, c2, c3, c4, ep, bes
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  vf
    real, dimension(6)                                          ::  c
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  nv

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    nx = nv(1)
    ny = nv(2)
    nz = nv(3)
    c1 = c(1)
    c2 = c(2)
    c3 = c(3)
    c4 = c(4)
    ep = c(5)         ! threshold
    cf = 0.5*dt

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, cf) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, b0, d_b, u_p, v_p, w_p, g_p, d_p, Fx_p, Fy_p, Fz_p) &
!$OMP PRIVATE(r_bt, uu, u1, u2, u3, es, bes)

!$OMP DO SCHEDULE(static)
    do k=ks,ke
    do j=js,je
    do i=is,ie
      idx = bd(i,j,k)
      b0 = vf(i,j,k)
      es = 0.0

      if ( ibits(idx, 0, bitw_6) == odr ) es=1.0 ! コンポーネントは流体であることが保証されている
      bes = b0*es ! 体積率とコンポーネントの両方でチェック

      ! 圧損による外力の計算(v^n)
      u_p = v(i,j,k,1) - u_ref
      v_p = v(i,j,k,2) - v_ref
      w_p = v(i,j,k,3) - w_ref
      g_p = sqrt(u_p*u_p + v_p*v_p + w_p*w_p)
      d_p = c1*g_p*g_p + c2*g_p + c3
      if (g_p < ep) d_p = c4*g_p*g_p
      Fx_p = -sign(1.0, u_p)*d_p*nx
      Fy_p = -sign(1.0, v_p)*d_p*ny
      Fz_p = -sign(1.0, w_p)*d_p*nz

      ! 擬似速度の方向の強制
      u1 = vc(i,j,k,1) - u_ref
      u2 = vc(i,j,k,2) - v_ref
      u3 = vc(i,j,k,3) - w_ref
      r_bt = 1.0 - bes
      uu = sqrt(u1*u1 + u2*u2 + u3*u3) * bes
      u1 = (r_bt*u1 + uu*nx) + u_ref
      u2 = (r_bt*u2 + uu*ny) + v_ref
      u3 = (r_bt*u3 + uu*nz) + w_ref

      d_b = cf * bes
      vc(i,j,k,1) = u1 + d_b * Fx_p
      vc(i,j,k,2) = u2 + d_b * Fy_p
      vc(i,j,k,3) = u3 + d_b * Fz_p

    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + dble(ie-is+1)*(je-js+1)*(ke-ks+1)*55.0d0

    return
    end subroutine hex_force_pvec

!> ********************************************************************
!! @brief 圧力損失部における速度の修正と発散値の修正
!! @param [in,out] v    速度ベクトル(n+1,k+1)
!! @param [in,out] div  速度の発散
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     st   ループの開始インデクス
!! @param [in]     ed   ループの終了インデクス
!! @param [in]     bd   BCindex ID
!! @param [in]     vf   コンポーネントの体積率
!! @param [in]     wk   テンポラリのワークベクトル 速度ベクトル (n+1,k)
!! @param [in]     cz   コンポーネントの配列長
!! @param [in]     odr  速度境界条件のエントリ
!! @param [in]     v00  参照速度
!! @param [in]     dt   時間積分幅
!! @param [in]     dh   格子幅
!! @param [in]     nv   法線ベクトル
!! @param [in]     c    圧力損失部の係数
!! @param [in]     am   平均速度と圧損量
!! @param [out]    flop flop count
!<
    subroutine hex_force_vec (v, div, sz, g, st, ed, bd, vf, wk, cz, odr, v00, dt, dh, nv, c, am, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, ii, jj, kk, idx, odr
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed, cz
    double precision                                            ::  flop
    real                                                        ::  dt, dh, beta
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  u_w, u_e, u_s, u_n, u_b, u_t, u_p
    real                                                        ::  v_w, v_e, v_s, v_n, v_b, v_t, v_p
    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t, w_p
    real                                                        ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                        ::  d_w, d_e, d_s, d_n, d_b, d_t, d_p
    real                                                        ::  Fx_w, Fx_e, Fx_p
    real                                                        ::  Fy_s, Fy_n, Fy_p
    real                                                        ::  Fz_b, Fz_t, Fz_p
    real                                                        ::  be, bw, bn, bs, bt, bb, b0, es, qq
    real                                                        ::  re, rw, rn, rs, rt, rb
    real                                                        ::  nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep, am1, am2, pick
    real                                                        ::  q_w, q_e, q_s, q_n, q_b, q_t, q_p
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  vf
    real, dimension(-1:cz(1)+2, -1:cz(2)+2, -1:cz(3)+2, 3)      ::  wk
    real, dimension(6)                                          ::  c
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  nv
    real, dimension(2)                                          ::  am

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    nx = nv(1)
    ny = nv(2)
    nz = nv(3)
    c1 = c(1)
    c2 = c(2)
    c3 = c(3)
    c4 = c(4)
    ep = c(5)         ! threshold

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    am1 = 0.0 ! 通過速度の積算値
    am2 = 0.0 ! 圧損量の積算値

!$OMP PARALLEL &
!$OMP REDUCTION(+:am1) &
!$OMP REDUCTION(+:am2) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, dt) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, es, beta, pick, ii, jj, kk) &
!$OMP PRIVATE(be, bw, bn, bs, bt, bb, b0) &
!$OMP PRIVATE(u_w, u_e, u_s, u_n, u_b, u_t, u_p) &
!$OMP PRIVATE(v_w, v_e, v_s, v_n, v_b, v_t, v_p) &
!$OMP PRIVATE(w_w, w_e, w_s, w_n, w_b, w_t, w_p) &
!$OMP PRIVATE(g_w, g_e, g_s, g_n, g_b, g_t, g_p) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t, d_p) &
!$OMP PRIVATE(Fx_w, Fx_e, Fx_p) &
!$OMP PRIVATE(Fy_s, Fy_n, Fy_p) &
!$OMP PRIVATE(Fz_b, Fz_t, Fz_p) &
!$OMP PRIVATE(re, rw, rn, rs, rt, rb) &
!$OMP PRIVATE(q_w, q_e, q_s, q_n, q_b, q_t, q_p)

! 周囲の1セルを含めてサーチ
!$OMP DO SCHEDULE(static)

    do k=ks-1,ke+1
      kk = k - ks + 1

    do j=js-1,je+1
      jj = j - js + 1

    do i=is-1,ie+1
      ii = i - is + 1

      ! 自セルが流体で、かつ周囲1セルの範囲にコンポーネントが存在する場合のみ有効
      pick = 0.0
      qq = q_p + q_w + q_e + q_s + q_n + q_b + q_t
      if ( qq > 0.0 ) pick = 1.0
      es = pick * real(ibits(idx, State, 1)) ! Fluid=1, Solid=0

      include 'force.h' ! 179 flop

      ! 発散値の修正 セルフェイスのフラックスの和
      div(i,j,k) = div(i,j,k) + dt * ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * es ! esはマスク

      beta = dt * b0 * q_p
      v(i,j,k,1) = v(i,j,k,1) + beta * Fx_p
      v(i,j,k,2) = v(i,j,k,2) + beta * Fy_p
      v(i,j,k,3) = v(i,j,k,3) + beta * Fz_p
      am1 = am1 + g_p * q_p
      am2 = am2 + d_p * q_p
      
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + dble(ie-is+3)*(je-js+3)*(ke-ks+3)*213.0d0

    am(1) = am1 ! accumulated velocity
    am(2) = am2 ! accumulated pressure loss

    return
    end subroutine hex_force_vec
