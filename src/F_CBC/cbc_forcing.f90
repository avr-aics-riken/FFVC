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

!  ********************************************************************
!> @subroutine cbc_pvec_hex (v, sz, g, st, ed, bd, odr, v00, vec, flop)
!! @brief 圧力損失部における疑似速度ベクトルの方向修正
!! @param[in/out] v セルセンターの疑似速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bd BCindex ID
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param vec 法線ベクトル
!! @param[out] flop flop count
!<
    subroutine cbc_pvec_hex (v, sz, g, st, ed, bd, odr, v00, vec, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, m, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  bt, flop, qtz
    real                                                        ::  u_ref, v_ref, w_ref, b_c
    real                                                        ::  nx, ny, nz, u1, u2, u3, uu
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    nx = vec(1)
    ny = vec(2)
    nz = vec(3)
    qtz = 1.0/255.0
    m = 0

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      idx = bd(i,j,k)
      if ( ibits(idx, 0, bitw_cmp) == odr ) then ! コンポーネントのエントリがHEXの場合
        b_c = real(ibits(idx, State, 1))
        
        u1 = v(i,j,k,1) - u_ref
        u2 = v(i,j,k,2) - v_ref
        u3 = v(i,j,k,3) - w_ref
        
        bt = real(ibits( idx, top_vf, bitw_vf )) * qtz
        uu = sqrt(u1*u1 + u2*u2 + u3*u3) * bt
        
        ! 参照座標系上での表現
        v(i,j,k,1) = ((1.0-bt)*u1 + uu*nx) + u_ref 
        v(i,j,k,2) = ((1.0-bt)*u2 + uu*ny) + v_ref
        v(i,j,k,3) = ((1.0-bt)*u3 + uu*nz) + w_ref

        m = m+1
      end if
    end do
    end do
    end do
    flop = flop + real(m*27)

    return
    end subroutine cbc_pvec_hex

!  *******************************************************************************
!> @subroutine cbc_psrc_hex (src, sz, g, st, ed, dh, bd, odr, v00, nv, c, v, flop)
!! @brief 圧力損失部におけるPoissonのソース項を計算する
!! @param[out] src ソース項
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param bd BCindex ID
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param v セルセンターの速度ベクトル n+1
!! @param[out] flop flop count
!! @note 隣接セルが固体の場合，速度はゼロで，外力もゼロ
!<
    subroutine cbc_psrc_hex (src, sz, g, st, ed, dh, bd, odr, v00, nv, c, v, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, m, odr, pick
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, dh
    real                                                        ::  u_ref, v_ref, w_ref, b_c
    real                                                        ::  u_w, u_e, u_s, u_n, u_b, u_t, u_p
    real                                                        ::  v_w, v_e, v_s, v_n, v_b, v_t, v_p
    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t, w_p
    real                                                        ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                        ::  d_w, d_e, d_s, d_n, d_b, d_t, d_p
    real                                                        ::  Fx_w, Fx_e, Fx_p
    real                                                        ::  Fy_s, Fy_n, Fy_p
    real                                                        ::  Fz_b, Fz_t, Fz_p
    real                                                        ::  be, bw, bn, bs, bt, bb, b0, qtz
    real                                                        ::  re, rw, rn, rs, rt, rb
    real                                                        ::  nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
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
    qtz = 1.0/255.0
    m = 0

    ! 周囲の1セルを含めてサーチ
    do k=st(3)-1,ed(3)+1
    do j=st(2)-1,ed(2)+1
    do i=st(1)-1,ed(1)+1
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
      
      if ( (pick == 1) .and. (b_c == 1) ) then ! コンポーネントの周囲1セルの範囲で流体
        b0 = real(ibits( idx,            top_vf, bitw_vf )) * qtz
        
        u_p = v(i  ,j  ,k  ,1) - u_ref
        v_p = v(i  ,j  ,k  ,2) - v_ref
        w_p = v(i  ,j  ,k  ,3) - w_ref
        g_p = sqrt(u_p*u_p + v_p*v_p + w_p*w_p)
        d_p = c1*g_p*g_p + c2*g_p + c3
        if (g_p < ep) d_p = c4*g_p*g_p
        Fx_p = -sign(1.0, u_p)*d_p*nx
        Fy_p = -sign(1.0, v_p)*d_p*ny
        Fz_p = -sign(1.0, w_p)*d_p*nz
        
        u_w = v(i-1,j  ,k  ,1) - u_ref
        v_w = v(i-1,j  ,k  ,2) - v_ref
        w_w = v(i-1,j  ,k  ,3) - w_ref
        g_w = sqrt(u_w*u_w + v_w*v_w + w_w*w_w)
        d_w = c1*g_w*g_w + c2*g_w + c3
        if (g_w < ep) d_w = c4*g_w*g_w
        Fx_w = -sign(1.0, u_w)*d_w*nx
        
        u_e = v(i+1,j  ,k  ,1) - u_ref
        v_e = v(i+1,j  ,k  ,2) - v_ref
        w_e = v(i+1,j  ,k  ,3) - w_ref
        g_e = sqrt(u_e*u_e + v_e*v_e + w_e*w_e)
        d_e = c1*g_e*g_e + c2*g_e + c3
        if (g_e < ep) d_e = c4*g_e*g_e
        Fx_e = -sign(1.0, u_e)*d_e*nx
        
        u_s = v(i  ,j-1,k  ,1) - u_ref
        v_s = v(i  ,j-1,k  ,2) - v_ref
        w_s = v(i  ,j-1,k  ,3) - w_ref
        g_s = sqrt(u_s*u_s + v_s*v_s + w_s*w_s)
        d_s = c1*g_s*g_s + c2*g_s + c3
        if (g_s < ep) d_s = c4*g_s*g_s
        Fy_s = -sign(1.0, v_s)*d_s*ny
        
        u_n = v(i  ,j+1,k  ,1) - u_ref
        v_n = v(i  ,j+1,k  ,2) - v_ref
        w_n = v(i  ,j+1,k  ,3) - w_ref
        g_n = sqrt(u_n*u_n + v_n*v_n + w_n*w_n)
        d_n = c1*g_n*g_n + c2*g_n + c3
        if (g_n < ep) d_n = c4*g_n*g_n
        Fy_n = -sign(1.0, v_n)*d_n*ny
        
        u_b = v(i  ,j  ,k-1,1) - u_ref
        v_b = v(i  ,j  ,k-1,2) - v_ref
        w_b = v(i  ,j  ,k-1,3) - w_ref
        g_b = sqrt(u_b*u_b + v_b*v_b + w_b*w_b)
        d_b = c1*g_b*g_b + c2*g_b + c3
        if (g_b < ep) d_b = c4*g_b*g_b
        Fz_b = -sign(1.0, w_b)*d_b*nz
        
        u_t = v(i  ,j  ,k+1,1) - u_ref
        v_t = v(i  ,j  ,k+1,2) - v_ref
        w_t = v(i  ,j  ,k+1,3) - w_ref
        g_t = sqrt(u_t*u_t + v_t*v_t + w_t*w_t)
        d_t = c1*g_t*g_t + c2*g_t + c3
        if (g_t < ep) d_t = c4*g_t*g_t
        Fz_t = -sign(1.0, w_t)*d_t*nz
        
        re = 0.5*(Fx_p + Fx_e)
        rw = 0.5*(Fx_p + Fx_w)
        rn = 0.5*(Fy_p + Fy_n)
        rs = 0.5*(Fy_p + Fy_s)
        rt = 0.5*(Fz_p + Fz_t)
        rb = 0.5*(Fz_p + Fz_b)
        
        bw = 0.5*( b0 + real(ibits(bd(i-1,j  ,k  ), top_vf, bitw_vf )) * qtz )
        be = 0.5*( b0 + real(ibits(bd(i+1,j  ,k  ), top_vf, bitw_vf )) * qtz )
        bs = 0.5*( b0 + real(ibits(bd(i  ,j-1,k  ), top_vf, bitw_vf )) * qtz )
        bn = 0.5*( b0 + real(ibits(bd(i  ,j+1,k  ), top_vf, bitw_vf )) * qtz )
        bb = 0.5*( b0 + real(ibits(bd(i  ,j  ,k-1), top_vf, bitw_vf )) * qtz )
        bt = 0.5*( b0 + real(ibits(bd(i  ,j  ,k+1), top_vf, bitw_vf )) * qtz )
        src(i,j,k) = -dh * ( be*re - bw*rw + bn*rn - bs*rs + bt*rt - bb*rb ) * b_c

        m = m+1
      end if
    end do
    end do
    end do
    
    flop = flop + real(m*287)

    return
    end subroutine cbc_psrc_hex

!  ************************************************************************************
!> @subroutine cbc_update_vcf_hex (vf, sz, g, st, ed, dt, bd, odr, v00, nv, c, v, flop)
!! @brief 射影プロセスにおいて，圧力損失部によりセルフェイス速度を修正する
!! @param[out] vf セルフェイスの速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dt 時間積分幅
!! @param bd BCindex ID
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param v セルセンターの速度ベクトル n+1
!! @param[out] flop flop count
!<
    subroutine cbc_update_vcf_hex (vf, sz, g, st, ed, dt, bd, odr, v00, nv, c, v, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, m, q
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, dt, u_ref, v_ref, w_ref, nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep
    real                                                        ::  u_e, u_n, u_t, u_p
    real                                                        ::  v_e, v_n, v_t, v_p
    real                                                        ::  w_e, w_n, w_t, w_p
    real                                                        ::  g_e, g_n, g_t, g_p
    real                                                        ::  d_e, d_n, d_t, d_p
    real                                                        ::  Fx_e, Fx_p
    real                                                        ::  Fy_n, Fy_p
    real                                                        ::  Fz_t, Fz_p
    real                                                        ::  be, bn, bt, b0, s_p, s_e, s_n, s_t
    real                                                        ::  re, rn, rt, qtz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf, v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
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
    m = 0
    q = 0
    qtz = 1.0/255.0
    
    ! ループはマイナス1からで候補フェイスをカバー
    do k=st(3)-1,ed(3)
    do j=st(2)-1,ed(2)
    do i=st(1)-1,ed(1)
      idx = bd(i,j,k)
      s_p = real(ibits(idx, State, 1)) ! Fluid=1, Solid=0
      
      if ( s_p == 1 ) then ! 流体セルのみ候補
      
        b0  = real(ibits( idx, top_vf, bitw_vf )) * qtz
      
        u_p = v(i  ,j  ,k  ,1) - u_ref
        v_p = v(i  ,j  ,k  ,2) - v_ref
        w_p = v(i  ,j  ,k  ,3) - w_ref
        g_p = sqrt(u_p*u_p + v_p*v_p + w_p*w_p)
        d_p = c1*g_p*g_p + c2*g_p + c3
        if (g_p < ep) d_p = c4*g_p*g_p
        Fx_p = -sign(1.0, u_p)*d_p*nx
        Fy_p = -sign(1.0, v_p)*d_p*ny
        Fz_p = -sign(1.0, w_p)*d_p*nz
      
        q = q+1
      
        ! X方向
        if ( (ibits(idx, 0, bitw_cmp) == odr) .or. (ibits(bd(i+1,j,k), 0, bitw_cmp) == odr) ) then ! フェイスの両側のどちらかのコンポーネントのエントリがHEXの場合
          u_e = v(i+1,j  ,k  ,1) - u_ref
          v_e = v(i+1,j  ,k  ,2) - v_ref
          w_e = v(i+1,j  ,k  ,3) - w_ref
          g_e = sqrt(u_e*u_e + v_e*v_e + w_e*w_e)
          d_e = c1*g_e*g_e + c2*g_e + c3
          if (g_e < ep) d_e = c4*g_e*g_e
          Fx_e = -sign(1.0, u_e)*d_e*nx
          s_e = real(ibits(bd(i+1,j  ,k  ), State, 1)) ! Fluid=1, Solid=0
        
          re = 0.5*(Fx_p + Fx_e)
          be = 0.5*( b0 + real(ibits(bd(i+1,j  ,k  ), top_vf, bitw_vf )) * qtz )
          vf(i,j,k,1) = vf(i,j,k,1) + dt * be * re * s_p*s_e ! s_p*s_eはセル界面の両側が流体の場合のみ1.0
        
          m = m+1
        endif
      
        ! Y方向
        if ( (ibits(idx, 0, bitw_cmp) == odr) .or. (ibits(bd(i,j+1,k), 0, bitw_cmp) == odr) ) then
          u_n = v(i  ,j+1,k  ,1) - u_ref
          v_n = v(i  ,j+1,k  ,2) - v_ref
          w_n = v(i  ,j+1,k  ,3) - w_ref
          g_n = sqrt(u_n*u_n + v_n*v_n + w_n*w_n)
          d_n = c1*g_n*g_n + c2*g_n + c3
          if (g_n < ep) d_n = c4*g_n*g_n
          Fy_n = -sign(1.0, v_n)*d_n*ny
          s_n = real(ibits(bd(i  ,j+1,k  ), State, 1)) ! Fluid=1, Solid=0
        
          rn = 0.5*(Fy_p + Fy_n)
          bn = 0.5*( b0 + real(ibits(bd(i  ,j+1,k  ), top_vf, bitw_vf )) * qtz )
          vf(i,j,k,2) = vf(i,j,k,2) + dt * bn * rn * s_p*s_n
        
          m = m+1
        endif
      
        ! Z方向
        if ( (ibits(idx, 0, bitw_cmp) == odr) .or. (ibits(bd(i,j,k+1), 0, bitw_cmp) == odr) ) then
          u_t = v(i  ,j  ,k+1,1) - u_ref
          v_t = v(i  ,j  ,k+1,2) - v_ref
          w_t = v(i  ,j  ,k+1,3) - w_ref
          g_t = sqrt(u_t*u_t + v_t*v_t + w_t*w_t)
          d_t = c1*g_t*g_t + c2*g_t + c3
          if (g_t < ep) d_t = c4*g_t*g_t
          Fz_t = -sign(1.0, w_t)*d_t*nz
          s_t = real(ibits(bd(i  ,j  ,k+1), State, 1)) ! Fluid=1, Solid=0
        
          rt = 0.5*(Fz_p + Fz_t)
          bt = 0.5*( b0 + real(ibits(bd(i  ,j  ,k+1), top_vf, bitw_vf )) * qtz )
          vf(i,j,k,3) = vf(i,j,k,3) + dt * bt * rt * s_p*s_t
        
          m = m+1
        endif
      endif

    end do
    end do
    end do
    
    flop = flop + real(q*33 + m*32)

    return
    end subroutine cbc_update_vcf_hex

!  ************************************************************************************
!> @subroutine cbc_update_vcc_hex (v, sz, g, st, ed, dt, bd, odr, v00, nv, c, vm, flop)
!! @brief 射影プロセスにおいて，圧力損失部によりセルセンター速度を修正する
!! @param[out] v セルセンターの速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dt 時間積分幅
!! @param bd BCindex ID
!! @param odr 速度境界条件のエントリ
!! @param v00 参照速度
!! @param nv 法線ベクトル
!! @param c 圧力損失部の係数
!! @param vm ローカルノードにおける無次元通過風速と圧力損失量の総和
!! @param[out] flop flop count
!<
    subroutine cbc_update_vcc_hex (v, sz, g, st, ed, dt, bd, odr, v00, nv, c, vm, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, m
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, dt, r, u_ref, v_ref, w_ref, nx, ny, nz
    real                                                        ::  c1, c2, c3, c4, ep, qtz, a_vel, a_prs, s
    real                                                        ::  u_p, v_p, w_p, g_p, d_p, Fx_p, Fy_p, Fz_p, rt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    real, dimension(6)                                          ::  c
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  nv
    real, dimension(2)                                          ::  vm

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
    m = 0
    qtz = 1.0/255.0
    a_prs = 0.0
    a_vel = 0.0

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      idx = bd(i,j,k)
      if ( ibits(idx, 0, bitw_cmp) == odr ) then ! コンポーネントのエントリがHEXの場合
        r = real(ibits( idx, top_vf, bitw_vf )) * qtz
        rt= r * dt
        u_p = v(i  ,j  ,k  ,1) - u_ref
        v_p = v(i  ,j  ,k  ,2) - v_ref
        w_p = v(i  ,j  ,k  ,3) - w_ref
        g_p = sqrt(u_p*u_p + v_p*v_p + w_p*w_p)
        d_p = c1*g_p*g_p + c2*g_p + c3
        if (g_p < ep) d_p = c4*g_p*g_p
        Fx_p = -sign(1.0, u_p)*d_p*nx
        Fy_p = -sign(1.0, v_p)*d_p*ny
        Fz_p = -sign(1.0, w_p)*d_p*nz
        
        v(i,j,k,1) = v(i,j,k,1) + rt * Fx_p
        v(i,j,k,2) = v(i,j,k,2) + rt * Fy_p
        v(i,j,k,3) = v(i,j,k,3) + rt * Fz_p
        
        s = nx*u_p + ny*v_p + nz*w_p
        a_vel = a_vel + sign(1.0, s)*g_p*r
        a_prs = a_prs - sign(1.0, s)*d_p*r
        
        m = m+1
      end if
    end do
    end do
    end do
    
    vm(1)= a_vel
    vm(2)= a_prs
    flop = flop + real(m*55)

    return
    end subroutine cbc_update_vcc_hex