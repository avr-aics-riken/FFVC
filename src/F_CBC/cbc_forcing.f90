!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************

!> @file cbc_forcing.f90
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN

!  ****************************************************************************
!> @subroutine cbc_hex_dir_pvec (v, sz, g, st, ed, bd, vf, odr, v00, vec, flop)
!! @brief 圧力損失部における疑似速度ベクトルの方向修正
!! @param[in/out] v 疑似速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param vec 法線ベクトル
!! @param[out] flop flop count
!! @note テンポラリに if ( (b0 > 0.0) .and. (ibits(idx, 0, bitw_cmp) == odr) ) a[i,j,k]=1.0 else 0.0 ループ分割
!<
    subroutine cbc_hex_dir_pvec (v, sz, g, st, ed, bd, vf, odr, v00, vec, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, u_ref, v_ref, w_ref, r_bt, b0
    real                                                        ::  nx, ny, nz, u1, u2, u3, uu
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  vf
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    nx = vec(1)
    ny = vec(2)
    nz = vec(3)

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, nx, ny, nz) &
!$OMP PRIVATE(idx, b0, r_bt, u1, u2, u3, uu)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=ks,ke
    do j=js,je
    do i=is,ie
      idx = bd(i,j,k)
      b0  = vf(i,j,k)

      if ( (b0 > 0.0) .and. (ibits(idx, 0, bitw_cmp) == odr) ) then ! 体積率があるセルで、コンポーネントのエントリがodrの場合
        u1 = v(1,i,j,k) - u_ref
        u2 = v(2,i,j,k) - v_ref
        u3 = v(3,i,j,k) - w_ref
        r_bt = 1.0-b0

        uu = sqrt(u1*u1 + u2*u2 + u3*u3) * b0
        
        ! 参照座標系上での表現
        v(1,i,j,k) = (r_bt*u1 + uu*nx) + u_ref 
        v(2,i,j,k) = (r_bt*u2 + uu*ny) + v_ref
        v(3,i,j,k) = (r_bt*u3 + uu*nz) + w_ref
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + (ie-is+1)*(je-js+1)*(ke-ks+1)*24.0

    return
    end subroutine cbc_hex_dir_pvec

!  *************************************************************************************
!> @subroutine cbc_hex_psrc (src, sz, g, st, ed, bd, vf, v, odr, v00, coef, nv, c, flop)
!! @brief 圧力損失部におけるPoissonのソース項を計算する
!! @param[out] src ソース項
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param v セルセンターの速度ベクトル n+1
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param coef 係数 dh/dt
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param[out] flop flop count
!! @note 隣接セルが固体の場合，速度はゼロで，外力もゼロ
!<
    subroutine cbc_hex_psrc (src, sz, g, st, ed, bd, vf, v, odr, v00, coef, nv, c, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, pick
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef, cf
    real                                                        ::  u_ref, v_ref, w_ref, b_c
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
    real                                                        ::  c1, c2, c3, c4, ep
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
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
    cf = -0.5*coef

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, cf) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, pick, b_c, es) &
!$OMP PRIVATE(be, bw, bn, bs, bt, bb, b0) &
!$OMP PRIVATE(u_w, u_e, u_s, u_n, u_b, u_t, u_p) &
!$OMP PRIVATE(v_w, v_e, v_s, v_n, v_b, v_t, v_p) &
!$OMP PRIVATE(w_w, w_e, w_s, w_n, w_b, w_t, w_p) &
!$OMP PRIVATE(g_w, g_e, g_s, g_n, g_b, g_t, g_p) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t, d_p) &
!$OMP PRIVATE(Fx_w, Fx_e, Fx_p) &
!$OMP PRIVATE(Fy_s, Fy_n, Fy_p) &
!$OMP PRIVATE(Fz_b, Fz_t, Fz_p) &
!$OMP PRIVATE(re, rw, rn, rs, rt, rb)

    ! 周囲の1セルを含めてサーチ
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=ks-1,ke+1
    do j=js-1,je+1
    do i=is-1,ie+1
      idx = bd(i,j,k)
      pick = 0.0
      
      if ( ibits(idx,             0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i-1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i+1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j-1,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j+1,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j  ,k-1), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j  ,k+1), 0, bitw_cmp) == odr ) pick=1.0
      
      b_c = real(ibits(idx, State, 1)) ! Fluid=1, Solid=0
      es = pick * b_c ! 自セルが流体で、かつ周囲1セルの範囲にコンポーネントが存在する場合のみ有効

      include 'force.h' ! 170 flop

      src(i,j,k) = cf * ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * es ! esはマスク
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    flop = flop + (ie-is+3)*(je-js+3)*(ke-ks+3)*185.0

    return
    end subroutine cbc_hex_psrc

!  ****************************************************************************************
!> @subroutine cbc_hex_force_pvec (vc, sz, g, st, ed, bd, vf, v, odr, v00, dt, nv, c, flop)
!! @brief 圧力損失部における速度ベクトルの射影の修正、および発散値の修正
!! @param[in/out] vc 擬似速度ベクトル
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
    subroutine cbc_hex_force_pvec (vc, sz, g, st, ed, bd, vf, v, odr, v00, dt, nv, c, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, dt, cf
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  u_p, v_p, w_p, g_p, d_p, Fx_p, Fy_p, Fz_p
    real                                                        ::  nx, ny, nz, b0, d_b
    real                                                        ::  c1, c2, c3, c4, ep
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, vc
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
!$OMP PRIVATE(idx, b0, d_b, u_p, v_p, w_p, g_p, d_p, Fx_p, Fy_p, Fz_p)

! ループ範囲は体積率のbbox
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=ks,ke
    do j=js,je
    do i=is,ie
      idx = bd(i,j,k)
      b0 = vf(i,j,k)

      if ( (b0 > 0.0) .and. (ibits(idx, 0, bitw_cmp) == odr) ) then ! 体積率があるセルで、コンポーネントのエントリがodrの場合
        u_p = v(1, i  ,j  ,k  ) - u_ref
        v_p = v(2, i  ,j  ,k  ) - v_ref
        w_p = v(3, i  ,j  ,k  ) - w_ref
        g_p = sqrt(u_p*u_p + v_p*v_p + w_p*w_p)
        d_p = c1*g_p*g_p + c2*g_p + c3
        if (g_p < ep) d_p = c4*g_p*g_p
        Fx_p = -sign(1.0, u_p)*d_p*nx
        Fy_p = -sign(1.0, v_p)*d_p*ny
        Fz_p = -sign(1.0, w_p)*d_p*nz

        d_b = cf * b0
        vc(1,i,j,k) = vc(1,i,j,k) + d_b * Fx_p
        vc(2,i,j,k) = vc(2,i,j,k) + d_b * Fy_p
        vc(3,i,j,k) = vc(3,i,j,k) + d_b * Fz_p
      end if

    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + (ie-is+1)*(je-js+1)*(ke-ks+1)*34.0

    return
    end subroutine cbc_hex_force_pvec

!  ************************************************************************************************
!> @subroutine cbc_hex_force_vec (v, div, sz, g, st, ed, bd, vf, odr, v00, dt, dh, nv, c, am, flop)
!! @brief 圧力損失部における速度の修正と発散値の修正
!! @param[in/out] v 速度ベクトルn+1
!! @param[in/out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param dt 時間積分幅
!! @param dh 格子幅
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param am 平均速度と圧損量
!! @param[out] flop flop count
!<
    subroutine cbc_hex_force_vec (v, div, sz, g, st, ed, bd, vf, odr, v00, dt, dh, nv, c, am, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, pick
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, dt, dh, cf1, cf2
    real                                                        ::  u_ref, v_ref, w_ref, b_c
    real                                                        ::  u_w, u_e, u_s, u_n, u_b, u_t, u_p
    real                                                        ::  v_w, v_e, v_s, v_n, v_b, v_t, v_p
    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t, w_p
    real                                                        ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                        ::  d_w, d_e, d_s, d_n, d_b, d_t, d_p
    real                                                        ::  Fx_w, Fx_e, Fx_p
    real                                                        ::  Fy_s, Fy_n, Fy_p
    real                                                        ::  Fz_b, Fz_t, Fz_p
    real                                                        ::  be, bw, bn, bs, bt, bb, b0, es, beta
    real                                                        ::  re, rw, rn, rs, rt, rb
    real                                                        ::  nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep, am1, am2
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  vf
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
    cf1= 0.5*dt
    cf2= 0.5*dh

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    am1 = 0.0 ! 通過速度の積算値
    am2 = 0.0 ! 圧損量の積算値

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, cf1, cf2) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, pick, b_c, es, beta) &
!$OMP PRIVATE(be, bw, bn, bs, bt, bb, b0) &
!$OMP PRIVATE(u_w, u_e, u_s, u_n, u_b, u_t, u_p) &
!$OMP PRIVATE(v_w, v_e, v_s, v_n, v_b, v_t, v_p) &
!$OMP PRIVATE(w_w, w_e, w_s, w_n, w_b, w_t, w_p) &
!$OMP PRIVATE(g_w, g_e, g_s, g_n, g_b, g_t, g_p) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t, d_p) &
!$OMP PRIVATE(Fx_w, Fx_e, Fx_p) &
!$OMP PRIVATE(Fy_s, Fy_n, Fy_p) &
!$OMP PRIVATE(Fz_b, Fz_t, Fz_p) &
!$OMP PRIVATE(re, rw, rn, rs, rt, rb)

! 周囲の1セルを含めてサーチ
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge) &
#endif
!$OMP REDUCTION(+:am1) &
!$OMP REDUCTION(+:am2)
    do k=ks-1,ke+1
    do j=js-1,je+1
    do i=is-1,ie+1
      idx = bd(i,j,k)
      pick = 0.0

      if ( ibits(idx,             0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i-1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i+1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j-1,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j+1,k  ), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j  ,k-1), 0, bitw_cmp) == odr ) pick=1.0
      if ( ibits(bd(i  ,j  ,k+1), 0, bitw_cmp) == odr ) pick=1.0

      b_c = real(ibits(idx, State, 1)) ! Fluid=1, Solid=0
      es = pick * b_c ! 自セルが流体で、かつ周囲1セルの範囲にコンポーネントが存在する場合のみ有効

      include 'force.h' ! 170 flop

      if ( b0 > 0.0 ) then ! 体積率があるセル
        beta = cf1 * b0
        v(1,i,j,k) = v(1,i,j,k) + beta * Fx_p
        v(2,i,j,k) = v(2,i,j,k) + beta * Fy_p
        v(3,i,j,k) = v(3,i,j,k) + beta * Fz_p
        am1 = am1 + g_p
        am2 = am2 + d_p
      end if

      div(i,j,k) = div(i,j,k) + cf2 * ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * es ! esはマスク
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + (ie-is+3)*(je-js+3)*(ke-ks+3)*199.0

    am(1) = am1
    am(2) = am2

    return
    end subroutine cbc_hex_force_vec
