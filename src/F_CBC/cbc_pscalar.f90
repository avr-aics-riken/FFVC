!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************

!> @file cbc_pscalar.f90
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN
    
!  **************************************************************************************
!> @subroutine cbc_ps_muscl (ws, sz, g, dh, c_scheme, v00, v, t, bv, bh1, bh2, swt, flop)
!! @brief パッシブスカラの移流部分の積分
!! @param[out] ws スカラ値の対流項による増分
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param c_scheme 対流項スキームのモード（1-UWD, 2-center, 3-MUSCL）
!! @param v00 参照速度
!! @param v cell center u^{n+1}
!! @param t 温度
!! @param bv  BCindex V
!! @param bh1 BCindex H1
!! @param bh2 BCindex H2
!! @param swt 固定壁の扱い（0-断熱，1-共役熱移動）
!! @param[out] flop
!<
    subroutine cbc_ps_muscl (ws, sz, g, dh, c_scheme, v00, v, t, bv, bh1, bh2, swt, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, c_scheme, idx, swt, hdx
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
    integer, dimension(3)                                       ::  sz
    real                                                        ::  UPe, UPw, VPn, VPs, WPt, WPb
    real                                                        ::  Up0, Ue1, Uw1
    real                                                        ::  Vp0, Vs1, Vn1
    real                                                        ::  Wp0, Wb1, Wt1
    real                                                        ::  Fp0, Fe1, Fe2, Fw1, Fw2, Fs1, Fs2, Fn1, Fn2, Fb1, Fb2, Ft1, Ft2
    real                                                        ::  ck, u_ref, v_ref, w_ref, dh, dh1, flop, actv
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  a_e, a_w, a_n, a_s, a_t, a_b
    real                                                        ::  d1, d2, d3, d4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4
    real                                                        ::  Fr_r, Fr_l, Fl_r, Fl_l
    real                                                        ::  cr, cl, acr, acl, cnv, ss, b
    real                                                        ::  w_e, w_w, w_n, w_s, w_t, w_b
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  t, ws
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bh1, bh2, bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    ss = 1.0
    ck = 0.0
    b = 0.0
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
            
    if ( c_scheme == 1 ) then      !     1st order upwind
      ss = 0.0
    else if ( c_scheme == 2 ) then !     2nd order central 
      ck = 1.0     
    else if ( c_scheme == 3 ) then !     3rd order MUSCL
      ck = 1.0/3.0
      b  = (3.0-ck)/(1.0-ck)
    endif

    ! /3 + 2 = 26 ! DP 41
    ! loop : 6+1+6+36 + ( 4+36+32+4+17 )*3dir + 2 + 3 = 333
    flop = flop + real(ix)*real(jx)*real(kx)*333.0 + 26.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*333.0 + 41.0 ! DP
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      cnv = 0.0

      ! セル状態 (0-solid / 1-fluid)
      b_p = ibits(bv(i  ,j  ,k  ), State, 1)
      b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      
      ! セル界面のフラグ (0-wall face / 1-fluid)
      w_e = real(b_e1 * b_p)
      w_w = real(b_w1 * b_p)
      w_n = real(b_n1 * b_p)
      w_s = real(b_s1 * b_p)
      w_t = real(b_t1 * b_p)
      w_b = real(b_b1 * b_p) ! real*6 flop
      
      ! 変数のロード
      Fp0 = t(i  ,j  ,k  )
      Fw2 = t(i-2,j  ,k  )
      Fw1 = t(i-1,j  ,k  )
      Fe1 = t(i+1,j  ,k  )
      Fe2 = t(i+2,j  ,k  )
      Fs2 = t(i  ,j-2,k  )
      Fs1 = t(i  ,j-1,k  )
      Fn1 = t(i  ,j+1,k  )
      Fn2 = t(i  ,j+2,k  )
      Fb2 = t(i  ,j  ,k-2)
      Fb1 = t(i  ,j  ,k-1)
      Ft1 = t(i  ,j  ,k+1)
      Ft2 = t(i  ,j  ,k+2)
      
      idx = bh1(i,j,k)
      hdx = bh2(i,j,k)
      actv= real(ibits(hdx, State, 1)) ! 1 flop 対流マスクなのでbh2の状態を参照
      
      ! セル状態 (0-solid / 1-fluid)
      b_w2= ibits(bh1(i-2,j  ,k  ), State, 1)
      b_w1= ibits(bh1(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bh1(i+1,j  ,k  ), State, 1)
      b_e2= ibits(bh1(i+2,j  ,k  ), State, 1)
      b_s2= ibits(bh1(i  ,j-2,k  ), State, 1)
      b_s1= ibits(bh1(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bh1(i  ,j+1,k  ), State, 1)
      b_n2= ibits(bh1(i  ,j+2,k  ), State, 1)
      b_b2= ibits(bh1(i  ,j  ,k-2), State, 1)
      b_b1= ibits(bh1(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bh1(i  ,j  ,k+1), State, 1)
      b_t2= ibits(bh1(i  ,j  ,k+2), State, 1)

      ! 各面のHBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC) 
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(idx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(idx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(idx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(idx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(idx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(idx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0

      ! 断熱マスク
      a_e = real(ibits(hdx, adbtc_E, 1))
      a_w = real(ibits(hdx, adbtc_W, 1))
      a_n = real(ibits(hdx, adbtc_N, 1))
      a_s = real(ibits(hdx, adbtc_S, 1))
      a_t = real(ibits(hdx, adbtc_T, 1))
      a_b = real(ibits(hdx, adbtc_B, 1)) ! real*6 = 6 flop

      Up0 = v(1, i  ,j  ,k  )
      Vp0 = v(2, i  ,j  ,k  )
      Wp0 = v(3, i  ,j  ,k  )
      Uw1 = v(1, i-1,j  ,k  )
      Ue1 = v(1, i+1,j  ,k  )
      Vs1 = v(2, i  ,j-1,k  )
      Vn1 = v(2, i  ,j+1,k  )
      Wb1 = v(3, i  ,j  ,k-1)
      Wt1 = v(3, i  ,j  ,k+1)
      
      ! 界面速度（スタガード位置）> 36 flop
      UPe = 0.5*(Up0+Ue1)*w_e + u_ref*(1.0-w_e)
      UPw = 0.5*(Up0+Uw1)*w_w + u_ref*(1.0-w_w)
      VPn = 0.5*(Vp0+Vn1)*w_n + v_ref*(1.0-w_n)
      VPs = 0.5*(Vp0+Vs1)*w_s + v_ref*(1.0-w_s)
      WPt = 0.5*(Wp0+Wt1)*w_t + w_ref*(1.0-w_t)
      WPb = 0.5*(Wp0+Wb1)*w_b + w_ref*(1.0-w_b)

      ! X方向 ---------------------------------------
      ! 壁面の扱い　共役熱移動の場合は，固体セル内部の値をそのまま使う
      if (swt == 0) then ! (i,j,k)からみて断熱壁の処理
        if ( b_e2 == 0 ) Fe2 = Fe1
        if ( b_e1 == 0 ) Fe1 = Fp0
        if ( b_w1 == 0 ) Fw1 = Fp0
        if ( b_w2 == 0 ) Fw2 = Fw1
      endif
      
      d4 = Fe2 - Fe1
      d3 = Fe1 - Fp0
      d2 = Fp0 - Fw1
      d1 = Fw1 - Fw2 ! 4 flop
      
      include 'muscl.h' ! 36 flop
      
      Fr_r = Fe1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Fr_l = Fp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Fl_r = Fp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Fl_l = Fw1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss ! 32 flop
      
      ! 流束　壁面上で速度ゼロ->対流熱流束がゼロになる >  4 flop
      cr  = UPe - u_ref
      cl  = UPw - u_ref
      acr = abs(cr)
      acl = abs(cl)
      
      ! 境界条件の面の寄与はスキップ > 17 flop
      cnv = 0.5*(cr*(Fr_r+Fr_l) - acr*(Fr_r-Fr_l)) * c_e * a_e &
          - 0.5*(cl*(Fl_r+Fl_l) - acl*(Fl_r-Fl_l)) * c_w * a_w

      ! Y方向 ---------------------------------------
      ! 壁面の扱い
      if (swt == 0) then ! (i,j,k)からみて断熱壁の処理
        if ( b_n2 == 0 ) Fn2 = Fn1
        if ( b_n1 == 0 ) Fn1 = Fp0
        if ( b_s1 == 0 ) Fs1 = Fp0
        if ( b_s2 == 0 ) Fs2 = Fs1
      endif
      
      d4 = Fn2 - Fn1
      d3 = Fn1 - Fp0
      d2 = Fp0 - Fs1
      d1 = Fs1 - Fs2
      
      include 'muscl.h'
      
      Fr_r = Fn1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Fr_l = Fp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Fl_r = Fp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Fl_l = Fs1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      
      cr  = VPn - v_ref
      cl  = VPs - v_ref
      acr = abs(cr)
      acl = abs(cl)
      
      cnv = 0.5*(cr*(Fr_r+Fr_l) - acr*(Fr_r-Fr_l)) * c_n * a_n &
          - 0.5*(cl*(Fl_r+Fl_l) - acl*(Fl_r-Fl_l)) * c_s * a_s + cnv

      ! Z方向 ---------------------------------------
      ! 壁面の扱い
      if (swt == 0) then ! (i,j,k)からみて断熱壁の処理
        if ( b_t2 == 0 ) Ft2 = Ft1
        if ( b_t1 == 0 ) Ft1 = Fp0
        if ( b_b1 == 0 ) Fb1 = Fp0
        if ( b_b2 == 0 ) Fb2 = Fb1
      endif
      
      d4 = Ft2 - Ft1
      d3 = Ft1 - Fp0
      d2 = Fp0 - Fb1
      d1 = Fb1 - Fb2
      
      include 'muscl.h'
      
      Fr_r = Ft1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Fr_l = Fp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Fl_r = Fp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Fl_l = Fb1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      
      cr  = WPt - w_ref
      cl  = WPb - w_ref
      acr = abs(cr)
      acl = abs(cl)
      
      cnv = 0.5*(cr*(Fr_r+Fr_l) - acr*(Fr_r-Fr_l)) * c_t * a_t &
          - 0.5*(cl*(Fl_r+Fl_l) - acl*(Fl_r-Fl_l)) * c_b * a_b + cnv

      ws(i,j,k) = ( -cnv*dh1 ) * actv
    end do
    end do
    end do

    return
    end subroutine cbc_ps_muscl

!  ****************************************************************
!> @subroutine cbc_ps_buoyancy (v, sz, g, dt, gr, rei, t, bd, flop)
!! @brief 浮力項の計算
!! @param[in/out] v 速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param gr Grashof数
!! @param rei Reynolds数の逆数
!! @param t 温度
!! @param bd BCindex ID
!! @param[out] flop
!<
    subroutine cbc_ps_buoyancy (v, sz, g, dt, gr, rei, t, bd, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  dt, gr, rei, dgr, flop
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  t
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dgr = dt*gr*rei*rei
    
    flop = flop + real(ix)*real(jx)*real(kx)*4.0 + 3.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
        v(3,i,j,k) = v(3,i,j,k) + dgr*t(i,j,k) * real(ibits(bd(i,j,k), State, 1))
    end do
    end do
    end do

    return
    end subroutine cbc_ps_buoyancy

!  **************************************************************************
!> @subroutine cbc_ps_diff_ee (t, sz, g, res, dh, dt, pei, qbc, bh, ws, flop)
!! @brief 温度の拡散項の半陰的時間積分（対流項を積分した結果を用いて粘性項を計算）
!! @param t 温度
!! @param sz 配列長
!! @param g ガイドセル長
!! @param res 残差
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param pei ペクレ数の逆数
!! @param qbc 境界条件の熱流束
!! @param bh BCindex H2
!! @param ws 部分段階の温度
!! @param flop 浮動小数演算数
!<
    subroutine cbc_ps_diff_ee (t, sz, g, res, dh, dt, pei, qbc, bh, ws, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  dh, dt, pei, dth1, dth2, res, flop, delta
    real                                                      ::  t_p, t_w, t_e, t_s, t_n, t_b, t_t
    real                                                      ::  g_p, g_w, g_e, g_s, g_n, g_b, g_t
    real                                                      ::       a_w, a_e, a_s, a_n, a_b, a_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  t, ws
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  qbc
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bh

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dth1 = dt/dh
    dth2 = dth1*pei/dh
    res  = 0.0
    
    ! /2 + 1 = 17 flop ! DP 27 flop
    ! loop : 6 + 6 + 1 + 51 = 64 flop
    flop = flop + real(ix)*real(jx)*real(kx)*64.0 + 17.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*64.0 + 27.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bh(i,j,k)
      
      t_p = t(i  , j  , k  )
      t_w = t(i-1, j  , k  )
      t_e = t(i+1, j  , k  )
      t_s = t(i  , j-1, k  )
      t_n = t(i  , j+1, k  )
      t_b = t(i  , j  , k-1)
      t_t = t(i  , j  , k+1)
      
      a_w = real(ibits(idx, adbtc_W, 1))
      a_e = real(ibits(idx, adbtc_E, 1))
      a_s = real(ibits(idx, adbtc_S, 1))
      a_n = real(ibits(idx, adbtc_N, 1))
      a_b = real(ibits(idx, adbtc_B, 1))
      a_t = real(ibits(idx, adbtc_T, 1)) ! real*6 = 6 flop
      
      g_w = real(ibits(idx, gma_W, 1))
      g_e = real(ibits(idx, gma_E, 1))
      g_s = real(ibits(idx, gma_S, 1))
      g_n = real(ibits(idx, gma_N, 1))
      g_b = real(ibits(idx, gma_B, 1))
      g_t = real(ibits(idx, gma_T, 1)) ! real*6 = 6 flop

      g_p = real(ibits(idx, h_diag, 3))  ! diagonal > 1 flop
      
      delta =(dth2*( g_w * a_w * t_w  & ! west  
                   + g_e * a_e * t_e  & ! east  
                   + g_s * a_s * t_s  & ! south 
                   + g_n * a_n * t_n  & ! north 
                   + g_b * a_b * t_b  & ! bottom
                   + g_t * a_t * t_t  & ! top
                   - g_p *       t_p) &
            - dth1*(-(1.0-g_w)*a_w * qbc(1, i-1, j  , k  )  & ! west   gamma
                    +(1.0-g_e)*a_e * qbc(1, i  , j  , k  )  & ! east   gamma
                    -(1.0-g_s)*a_s * qbc(2, i  , j-1, k  )  & ! south  gamma
                    +(1.0-g_n)*a_n * qbc(2, i  , j  , k  )  & ! north  gamma
                    -(1.0-g_b)*a_b * qbc(3, i  , j  , k-1)  & ! bottom gamma
                    +(1.0-g_t)*a_t * qbc(3, i  , j  , k  )  & ! top    gamma
            ) )  * real(ibits(idx, Active,  1))
      t(i,j,k) = ws(i,j,k) + delta
      res = res + delta*delta
    end do
    end do
    end do

    return
    end subroutine cbc_ps_diff_ee

!  **********************************************************
!> @subroutine cbc_hbc_drchlt (t, sz, g, st, ed, bh, odr, tc)
!! @brief 温度指定境界条件を設定するために必要な参照値をセットする
!! @param[in/out] t 温度
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bh BCindex H1
!! @param odr 内部境界処理時の境界条件のエントリ
!! @param tc 指定する値
!<
    subroutine cbc_hbc_drchlt (t, sz, g, st, ed, bh, odr, tc)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  tc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  t
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bh
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      idx = bh(i,j,k)

      if ( 0 /= iand(idx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        if ( ibits(idx, bc_face_E, bitw_5) == odr ) then
          t(i+1,j,k) = tc
        end if
        
        if ( ibits(idx, bc_face_W, bitw_5) == odr ) then
          t(i-1,j,k) = tc
        end if

        if ( ibits(idx, bc_face_N, bitw_5) == odr ) then
          t(i,j+1,k) = tc
        end if
        
        if ( ibits(idx, bc_face_S, bitw_5) == odr ) then
          t(i,j-1,k) = tc
        end if
			
        if ( ibits(idx, bc_face_T, bitw_5) == odr ) then
          t(i,j,k+1) = tc
        end if
        
        if ( ibits(idx, bc_face_B, bitw_5) == odr ) then
          t(i,j,k-1) = tc
        end if

      endif
    end do
    end do
    end do
    
    return
    end subroutine cbc_hbc_drchlt
