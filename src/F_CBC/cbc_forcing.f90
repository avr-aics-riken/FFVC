!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
!
!   *********************************************************

!> @file cbc_forcing.f90
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN

!  ************************************************************************
!> @subroutine cbc_pvec_hex (v, sz, g, st, ed, bd, odr, vf, v00, vec, flop)
!! @brief 圧力損失部における疑似速度ベクトルの方向修正
!! @param[in/out] v セルセンターの疑似速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param odr 速度境界条件のエントリ
!! @param vf コンポーネントの体積率
!! @param v00 参照速度
!! @param vec 法線ベクトル
!! @param[out] flop flop count
!<
    subroutine cbc_pvec_hex (v, sz, g, st, ed, bd, odr, vf, v00, vec, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, m, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, u_ref, v_ref, w_ref, bt, r_bt
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

    m = (ie-is+1)*(je-js+1)*(ke-ks+1)
    flop = flop + real(m)*24.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, nx, ny, nz) &
!$OMP PRIVATE(idx) &
!$OMP PRIVATE(u1, u2, u3, uu, bt, r_bt)

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
      if ( ibits(idx, 0, bitw_cmp) == odr ) then ! コンポーネントのエントリがHEXの場合
        u1 = v(1,i,j,k) - u_ref
        u2 = v(2,i,j,k) - v_ref
        u3 = v(3,i,j,k) - w_ref
        bt = vf(i,j,k)
        r_bt = 1.0-bt

        uu = sqrt(u1*u1 + u2*u2 + u3*u3) * bt
        
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

    return
    end subroutine cbc_pvec_hex

!  ***********************************************************************************
!> @subroutine cbc_psrc_hex (src, sz, g, st, ed, dh, bd, vf, odr, v00, nv, c, v, flop)
!! @brief 圧力損失部におけるPoissonのソース項を計算する
!! @param[out] src ソース項
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param v セルセンターの速度ベクトル n+1
!! @param[out] flop flop count
!! @note 隣接セルが固体の場合，速度はゼロで，外力もゼロ
!<
    subroutine cbc_psrc_hex (src, sz, g, st, ed, dh, bd, vf, odr, v00, nv, c, v, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, pick
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, dh
    real                                                        ::  u_ref, v_ref, w_ref, b_c, m
    real                                                        ::  u_w, u_e, u_s, u_n, u_b, u_t, u_p
    real                                                        ::  v_w, v_e, v_s, v_n, v_b, v_t, v_p
    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t, w_p
    real                                                        ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                        ::  d_w, d_e, d_s, d_n, d_b, d_t, d_p
    real                                                        ::  Fx_w, Fx_e, Fx_p
    real                                                        ::  Fy_s, Fy_n, Fy_p
    real                                                        ::  Fz_b, Fz_t, Fz_p
    real                                                        ::  be, bw, bn, bs, bt, bb, b0
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

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    m = 0.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, dh) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, pick, b_c) &
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
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:m)
    do k=ks-1,ke+1
    do j=js-1,je+1
    do i=is-1,ie+1
      idx = bd(i,j,k)
      pick=0
      
      if ( ibits(idx,             0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i-1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i+1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j-1,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j+1,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j  ,k-1), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j  ,k+1), 0, bitw_cmp) == odr ) pick=1
      
      b_c = real(ibits(idx, State, 1)) ! Fluid=1, Solid=0
      
      if ( (pick == 1) .and. (b_c == 1) ) then ! 自セルが流体で、かつ周囲1セルの範囲にコンポーネントが存在する場合力の計算が必要

        include 'force.h'

        src(i,j,k) = -dh * ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * b_c

        m = m+1.0
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    flop = flop + m*185.0

    return
    end subroutine cbc_psrc_hex

!  ***********************************************************************************************
!> @subroutine cbc_update_hex (div, sz, g, st, ed, bd, vf, odr, v00, coef, dt, nv, c, v, flop)
!! @brief 圧力損失部における速度ベクトルの射影の修正、および発散値の修正
!! @param[in/out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param vf コンポーネントの体積率
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param coef 係数
!! @param dt 時間積分幅
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param v セルセンターの速度ベクトル n+1
!! @param[out] flop flop count
!! @note 隣接セルが固体の場合，速度はゼロで，外力もゼロ
!<
    subroutine cbc_update_hex (div, sz, g, st, ed, bd, vf, odr, v00, coef, dt, nv, c, v, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, pick
    integer                                                     ::  is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef
    real                                                        ::  u_ref, v_ref, w_ref, b_c, m
    real                                                        ::  u_w, u_e, u_s, u_n, u_b, u_t, u_p
    real                                                        ::  v_w, v_e, v_s, v_n, v_b, v_t, v_p
    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t, w_p
    real                                                        ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                        ::  d_w, d_e, d_s, d_n, d_b, d_t, d_p
    real                                                        ::  Fx_w, Fx_e, Fx_p
    real                                                        ::  Fy_s, Fy_n, Fy_p
    real                                                        ::  Fz_b, Fz_t, Fz_p
    real                                                        ::  be, bw, bn, bs, bt, bb, b0, db
    real                                                        ::  re, rw, rn, rs, rt, rb
    real                                                        ::  nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
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

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    m = 0.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, coef) &
!$OMP FIRSTPRIVATE(nx, ny, nz, c1, c2, c3, c4, ep) &
!$OMP PRIVATE(idx, pick, b_c) &
!$OMP PRIVATE(be, bw, bn, bs, bt, bb, b0, db) &
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
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:m)
    do k=ks-1,ke+1
    do j=js-1,je+1
    do i=is-1,ie+1
      idx = bd(i,j,k)
      pick=0

      if ( ibits(idx,             0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i-1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i+1,j  ,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j-1,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j+1,k  ), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j  ,k-1), 0, bitw_cmp) == odr ) pick=1
      if ( ibits(bd(i  ,j  ,k+1), 0, bitw_cmp) == odr ) pick=1

      b_c = real(ibits(idx, State, 1)) ! Fluid=1, Solid=0

      if ( (pick == 1) .and. (b_c == 1) ) then ! 自セルが流体で、かつ周囲1セルの範囲にコンポーネントが存在する場合力の計算が必要

        include 'force.h'

        div(i,j,k) = div(i,j,k) + dt * ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * coef * b_c

        db = dt * b0
        v(1,i,j,k) = v(1,i,j,k) + db * Fx_p
        v(2,i,j,k) = v(2,i,j,k) + db * Fy_p
        v(3,i,j,k) = v(3,i,j,k) + db * Fz_p

        m = m+1.0
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + m*193.0

    return
    end subroutine cbc_update_hex
