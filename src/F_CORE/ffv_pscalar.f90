!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################
    
!> @file   ffv_pscalar.f90
!! @brief  パッシブスカラー計算のルーチン群
!! @author aics
!<

!> ********************************************************************
!! @brief パッシブスカラの移流部分の積分
!! @param [in,out] ws     内部エネルギーの対流項による増分
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     dh     格子幅
!! @param [in]     scheme 対流項スキームのモード（1-UWD, 3-MUSCL）
!! @param [in]     v00    参照速度
!! @param [in]     v      cell center u^{n+1}
!! @param [in]     ie     内部エネルギー
!! @param [in]     bid    Cut ID
!! @param [in]     cdf    BCindex C
!! @param [in[     bh     BCindex B
!! @param [in]     swt    固定壁の扱い（0-断熱，1-共役熱移動）
!! @param [out]    flop   浮動小数点演算数
!<
    subroutine ps_muscl (ws, sz, g, dh, scheme, v00, v, ie, bid, cdf, bh, swt, flop)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, scheme, idx, swt, hdx
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer                                                     ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
    integer, dimension(3)                                       ::  sz
    double precision                                            ::  flop
    real                                                        ::  UPe, UPw, VPn, VPs, WPt, WPb
    real                                                        ::  Up0, Ue1, Uw1
    real                                                        ::  Vp0, Vs1, Vn1
    real                                                        ::  Wp0, Wb1, Wt1
    real                                                        ::  Fp0, Fe1, Fe2, Fw1, Fw2, Fs1, Fs2, Fn1, Fn2, Fb1, Fb2, Ft1, Ft2
    real                                                        ::  ck, u_ref, v_ref, w_ref, dh, dh1, actv
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  a_e, a_w, a_n, a_s, a_t, a_b
    real                                                        ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4
    real                                                        ::  Fr_r, Fr_l, Fl_r, Fl_l
    real                                                        ::  cr, cl, acr, acl, cnv, ss, b, cm1, cm2, ss_4
    real                                                        ::  w_e, w_w, w_n, w_s, w_t, w_b
    real                                                        ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  ie, ws
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  cdf, bh, bid
    
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
            
    if ( scheme == 1 ) then      !     1st order upwind
      ss = 0.0
    else if ( scheme == 3 ) then !     3rd order MUSCL
      ck = 1.0/3.0
      b  = (3.0-ck)/(1.0-ck)
    endif

    ss_4 = 0.25*ss

    cm1 = 1.0 - ck
    cm2 = 1.0 + ck

    ! /3 + 2 = 26 ! DP 41
    ! loop : 6+1+6+36 + ( 4+36+22+4+17 )*3dir + 2 + 3 = 333
    flop = flop + dble(ix)*dble(jx)*dble(kx)*303.0d0 + 26.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*303.0d0 + 41.0d0 ! DP
    

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, ss, ck, b, u_ref, v_ref, w_ref, swt, cm1, cm2, ss_4) &
!$OMP PRIVATE(cnv, idx, hdx, actv) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(Fp0, Fe1, Fe2, Fw1, Fw2, Fs1, Fs2, Fn1, Fn2, Fb1, Fb2, Ft1, Ft2) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(a_e, a_w, a_n, a_s, a_t, a_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(Up0, Ue1, Uw1) &
!$OMP PRIVATE(Vp0, Vs1, Vn1) &
!$OMP PRIVATE(Wp0, Wb1, Wt1) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Fr_r, Fr_l, Fl_r, Fl_l) &
!$OMP PRIVATE(cr, cl, acr, acl)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      cnv = 0.0

      ! セル状態 (0-solid / 1-fluid)
      b_p = ibits(cdf(i  ,j  ,k  ), State, 1)
      b_w2= ibits(cdf(i-2,j  ,k  ), State, 1)
      b_w1= ibits(cdf(i-1,j  ,k  ), State, 1)
      b_e1= ibits(cdf(i+1,j  ,k  ), State, 1)
      b_e2= ibits(cdf(i+2,j  ,k  ), State, 1)
      b_s2= ibits(cdf(i  ,j-2,k  ), State, 1)
      b_s1= ibits(cdf(i  ,j-1,k  ), State, 1)
      b_n1= ibits(cdf(i  ,j+1,k  ), State, 1)
      b_n2= ibits(cdf(i  ,j+2,k  ), State, 1)
      b_b2= ibits(cdf(i  ,j  ,k-2), State, 1)
      b_b1= ibits(cdf(i  ,j  ,k-1), State, 1)
      b_t1= ibits(cdf(i  ,j  ,k+1), State, 1)
      b_t2= ibits(cdf(i  ,j  ,k+2), State, 1)

      ! セル界面のフラグ (0-wall face / 1-fluid)
      w_e = real(b_e1 * b_p)
      w_w = real(b_w1 * b_p)
      w_n = real(b_n1 * b_p)
      w_s = real(b_s1 * b_p)
      w_t = real(b_t1 * b_p)
      w_b = real(b_b1 * b_p) ! real*6 flop
      
      ! 変数のロード
      Fp0 = ie(i  ,j  ,k  )
      Fw2 = ie(i-2,j  ,k  )
      Fw1 = ie(i-1,j  ,k  )
      Fe1 = ie(i+1,j  ,k  )
      Fe2 = ie(i+2,j  ,k  )
      Fs2 = ie(i  ,j-2,k  )
      Fs1 = ie(i  ,j-1,k  )
      Fn1 = ie(i  ,j+1,k  )
      Fn2 = ie(i  ,j+2,k  )
      Fb2 = ie(i  ,j  ,k-2)
      Fb1 = ie(i  ,j  ,k-1)
      Ft1 = ie(i  ,j  ,k+1)
      Ft2 = ie(i  ,j  ,k+2)
      
      idx = cdf(i,j,k)
      hdx = bh(i,j,k)
      actv= real(ibits(hdx, State, 1)) ! 1 flop 対流マスクなのでbhの状態を参照


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

      ! 延びたステンシルの参照先がIBCである場合のスキームの破綻を回避，１次精度におとす
      lmt_w = 1.0
      lmt_e = 1.0
      lmt_s = 1.0
      lmt_n = 1.0
      lmt_b = 1.0
      lmt_t = 1.0

      if ( (ibits(cdf(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bid(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
      if ( (ibits(cdf(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bid(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
      if ( (ibits(cdf(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bid(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
      if ( (ibits(cdf(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bid(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
      if ( (ibits(cdf(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bid(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
      if ( (ibits(cdf(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bid(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0

      ! 外部境界条件の場合
      if ( (i == 1)  .and. (ibits(bid(0   , j   , k   ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
      if ( (i == ix) .and. (ibits(bid(ix+1, j   , k   ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
      if ( (j == 1)  .and. (ibits(bid(i   , 0   , k   ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
      if ( (j == jx) .and. (ibits(bid(i   , jx+1, k   ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
      if ( (k == 1)  .and. (ibits(bid(i   , j   , 0   ), vbc_uwd, 1) == 1) ) lmt_b = 0.0
      if ( (k == kx) .and. (ibits(bid(i   , j   , kx+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0

      Uw1 = v(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)

      Vs1 = v(i  ,j-1,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)

      Wb1 = v(i  ,j  ,k-1, 3)
      Wp0 = v(i  ,j  ,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3)
      
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
      
      dv4 = Fe2 - Fe1
      dv3 = Fe1 - Fp0
      dv2 = Fp0 - Fw1
      dv1 = Fw1 - Fw2 ! 4 flop
      
s4 = sign(1.0, dv4) ! sign is zero flop
s3 = sign(1.0, dv3)
s2 = sign(1.0, dv2)
s1 = sign(1.0, dv1)

g6 = s4 * max(0.0, min( abs(dv4), s4 * b * dv3))
g5 = s3 * max(0.0, min( abs(dv3), s3 * b * dv4))
g4 = s3 * max(0.0, min( abs(dv3), s3 * b * dv2))
g3 = s2 * max(0.0, min( abs(dv2), s2 * b * dv3))
g2 = s2 * max(0.0, min( abs(dv2), s2 * b * dv1))
g1 = s1 * max(0.0, min( abs(dv1), s1 * b * dv2))

      Fr_r = Fe1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Fr_l = Fp0 + (cm1*g3+cm2*g4)*ss_4
      Fl_r = Fp0 - (cm1*g4+cm2*g3)*ss_4
      Fl_l = Fw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w ! 22 flop
      
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
      
      dv4 = Fn2 - Fn1
      dv3 = Fn1 - Fp0
      dv2 = Fp0 - Fs1
      dv1 = Fs1 - Fs2
      
s4 = sign(1.0, dv4) ! sign is zero flop
s3 = sign(1.0, dv3)
s2 = sign(1.0, dv2)
s1 = sign(1.0, dv1)

g6 = s4 * max(0.0, min( abs(dv4), s4 * b * dv3))
g5 = s3 * max(0.0, min( abs(dv3), s3 * b * dv4))
g4 = s3 * max(0.0, min( abs(dv3), s3 * b * dv2))
g3 = s2 * max(0.0, min( abs(dv2), s2 * b * dv3))
g2 = s2 * max(0.0, min( abs(dv2), s2 * b * dv1))
g1 = s1 * max(0.0, min( abs(dv1), s1 * b * dv2))

      Fr_r = Fn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Fr_l = Fp0 + (cm1*g3+cm2*g4)*ss_4
      Fl_r = Fp0 - (cm1*g4+cm2*g3)*ss_4
      Fl_l = Fs1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      
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
      
      dv4 = Ft2 - Ft1
      dv3 = Ft1 - Fp0
      dv2 = Fp0 - Fb1
      dv1 = Fb1 - Fb2
      
s4 = sign(1.0, dv4) ! sign is zero flop
s3 = sign(1.0, dv3)
s2 = sign(1.0, dv2)
s1 = sign(1.0, dv1)

g6 = s4 * max(0.0, min( abs(dv4), s4 * b * dv3))
g5 = s3 * max(0.0, min( abs(dv3), s3 * b * dv4))
g4 = s3 * max(0.0, min( abs(dv3), s3 * b * dv2))
g3 = s2 * max(0.0, min( abs(dv2), s2 * b * dv3))
g2 = s2 * max(0.0, min( abs(dv2), s2 * b * dv1))
g1 = s1 * max(0.0, min( abs(dv1), s1 * b * dv2))

      Fr_r = Ft1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Fr_l = Fp0 + (cm1*g3+cm2*g4)*ss_4
      Fl_r = Fp0 - (cm1*g4+cm2*g3)*ss_4
      Fl_l = Fb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      
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
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine ps_muscl


!> ********************************************************************
!! @brief 浮力項の計算（単相流）
!! @param [in,out] v     速度ベクトル
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     dgr   係数
!! @param [in]     ie    内部エネルギー
!! @param [in]     bd    BCindex B
!! @param [in]     ncompo コンポーネント数
!! @param [in]     mtbl   コンポーネントの物性値
!! @param [in,out] flop  浮動小数点演算数
!! @todo 対象セルは流体だけでよい？　active flag?
!<
    subroutine ps_buoyancy (v, sz, g, dgr, ie, bd, ncompo, mtbl, flop)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, ncompo, l, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  dgr, r, rcp, t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ie
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    real*8, dimension(3, 0:ncompo)                            ::  mtbl

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    r = dgr
    
    flop = flop + dble(ix)*dble(jx)*dble(kx)*12.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, r) &
!$OMP PRIVATE(idx, l, rcp, t)

!$OMP DO SCHEDULE(static)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bd(i,j,k)
      l = ibits(idx, 0, bitw_5)
      rcp = mtbl(1, l) * mtbl(2, l)
      t = ie(i, j, k) / rcp
      v(i,j,k,3) = v(i,j,k,3) + r * t * real(ibits(idx, State, 1))
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine ps_buoyancy

!> ********************************************************************
!! @brief 温度の拡散項の半陰的時間積分（対流項を積分した結果を用いて粘性項を計算）
!! @param [in,out] ie     内部エネルギー
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [out]    res_l2 残差の自乗和
!! @param [in]     dh     格子幅
!! @param [in]     dt     時間積分幅
!! @param [in]     qbc    境界条件の熱流束
!! @param [in]     bh     BCindex B
!! @param [in,out] ws     部分段階の温度
!! @param [in]     ncompo コンポーネント数
!! @param [in]     mtbl   コンポーネントの物性値
!! @param [in]     h_mode mode(0-conjugate heat transfer / 1-othres)
!! @param [in,out] flop   浮動小数演算数
!<
    subroutine ps_diff_ee (ie, sz, g, res_l2, dh, dt, qbc, bh, ws, ncompo, mtbl, h_mode, flop)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx, ncompo, h_mode, hm
    integer                                                   ::  l_p, l_w, l_e, l_s, l_n, l_b, l_t
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, res, res_l2
    real                                                      ::  dh, dt, dth1, dth2, delta, sw
    real                                                      ::  t_p, t_w, t_e, t_s, t_n, t_b, t_t
    real                                                      ::  g_w, g_e, g_s, g_n, g_b, g_t, g_p
    real                                                      ::  a_w, a_e, a_s, a_n, a_b, a_t
    real                                                      ::  lmd_p, lmd_w, lmd_e, lmd_s, lmd_n, lmd_b, lmd_t
    real                                                      ::  rcp_p, rcp_w, rcp_e, rcp_s, rcp_n, rcp_b, rcp_t
    real                                                      ::  tc_w, tc_e, tc_s, tc_n, tc_b, tc_t
    real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ie, ws
    real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  qbc
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bh
    real*8, dimension(3, 0:ncompo)                            ::  mtbl

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dth1 = dt/dh
    dth2 = dth1/dh
    res  = 0.0
    hm = 1
    if ( h_mode == 0) hm = 0

    ! /2 + 1 = 17 flop ! DP 27 flop
    ! loop : 6 + 6 + 1 + 51 = 64 flop
    flop = flop + dble(ix)*dble(jx)*dble(kx)*195.0d0 + 24.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*64.0d0 + 27.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dth1, dth2, hm) &
!$OMP PRIVATE(idx, delta, sw) &
!$OMP PRIVATE(t_p, t_w, t_e, t_s, t_n, t_b, t_t) &
!$OMP PRIVATE(g_w, g_e, g_s, g_n, g_b, g_t, g_p) &
!$OMP PRIVATE(lmd_p, lmd_w, lmd_e, lmd_s, lmd_n, lmd_b, lmd_t) &
!$OMP PRIVATE(rcp_p, rcp_w, rcp_e, rcp_s, rcp_n, rcp_b, rcp_t) &
!$OMP PRIVATE(tc_w, tc_e, tc_s, tc_n, tc_b, tc_t) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t) &
!$OMP PRIVATE(l_p, l_w, l_e, l_s, l_n, l_b, l_t) &
!$OMP PRIVATE(a_w, a_e, a_s, a_n, a_b, a_t)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix

      idx = bh(i,j,k)

      sw = real(ibits(idx, Active,  1))
      if (hm == 0) sw = 1.0 ! conjugate heat transfer
      
      a_w = real(ibits(idx, adbtc_W, 1))
      a_e = real(ibits(idx, adbtc_E, 1))
      a_s = real(ibits(idx, adbtc_S, 1))
      a_n = real(ibits(idx, adbtc_N, 1))
      a_b = real(ibits(idx, adbtc_B, 1))
      a_t = real(ibits(idx, adbtc_T, 1))
      
      g_w = real(ibits(idx, gma_W, 1))
      g_e = real(ibits(idx, gma_E, 1))
      g_s = real(ibits(idx, gma_S, 1))
      g_n = real(ibits(idx, gma_N, 1))
      g_b = real(ibits(idx, gma_B, 1))
      g_t = real(ibits(idx, gma_T, 1))

      c_w = g_w * a_w
      c_e = g_e * a_e
      c_s = g_s * a_s
      c_n = g_n * a_n
      c_b = g_b * a_b
      c_t = g_t * a_t ! 6

      l_p = ibits(idx,               0, bitw_5)
      l_w = ibits(bh(i-1, j  , k  ), 0, bitw_5)
      l_e = ibits(bh(i+1, j  , k  ), 0, bitw_5)
      l_s = ibits(bh(i  , j-1, k  ), 0, bitw_5)
      l_n = ibits(bh(i  , j+1, k  ), 0, bitw_5)
      l_b = ibits(bh(i  , j  , k-1), 0, bitw_5)
      l_t = ibits(bh(i  , j  , k+1), 0, bitw_5)

      rcp_p = mtbl(1, l_p) * mtbl(2, l_p)
      rcp_w = mtbl(1, l_w) * mtbl(2, l_w)
      rcp_e = mtbl(1, l_e) * mtbl(2, l_e)
      rcp_s = mtbl(1, l_s) * mtbl(2, l_s)
      rcp_n = mtbl(1, l_n) * mtbl(2, l_n)
      rcp_b = mtbl(1, l_b) * mtbl(2, l_b)
      rcp_t = mtbl(1, l_t) * mtbl(2, l_t) ! 7

      lmd_p = mtbl(3, l_p)
      lmd_w = mtbl(3, l_w)
      lmd_e = mtbl(3, l_e)
      lmd_s = mtbl(3, l_s)
      lmd_n = mtbl(3, l_n)
      lmd_b = mtbl(3, l_b)
      lmd_t = mtbl(3, l_t)

      t_p = ie(i  , j  , k  ) / rcp_p
      t_w = ie(i-1, j  , k  ) / rcp_w
      t_e = ie(i+1, j  , k  ) / rcp_e
      t_s = ie(i  , j-1, k  ) / rcp_s
      t_n = ie(i  , j+1, k  ) / rcp_n
      t_b = ie(i  , j  , k-1) / rcp_b
      t_t = ie(i  , j  , k+1) / rcp_t ! 7*8 = 56

      tc_w = lmd_p * lmd_w / (lmd_p + lmd_w) * 2.0
      tc_e = lmd_p * lmd_e / (lmd_p + lmd_e) * 2.0
      tc_s = lmd_p * lmd_s / (lmd_p + lmd_s) * 2.0
      tc_n = lmd_p * lmd_n / (lmd_p + lmd_n) * 2.0
      tc_b = lmd_p * lmd_b / (lmd_p + lmd_b) * 2.0
      tc_t = lmd_p * lmd_t / (lmd_p + lmd_t) * 2.0 ! (3+8)*6 = 66

      g_p = c_w * tc_w &
          + c_e * tc_e &
          + c_s * tc_s &
          + c_n * tc_n &
          + c_b * tc_b &
          + c_t * tc_t    ! 11

      delta =(dth2*( c_w * tc_w * t_w  & ! west
                   + c_e * tc_e * t_e  & ! east
                   + c_s * tc_s * t_s  & ! south
                   + c_n * tc_n * t_n  & ! north
                   + c_b * tc_b * t_b  & ! bottom
                   + c_t * tc_t * t_t  & ! top
                   - g_p *        t_p) &
            - dth1*(-(1.0-g_w)*a_w * qbc(1, i, j, k)  & ! west   gamma
                    +(1.0-g_e)*a_e * qbc(2, i, j, k)  & ! east   gamma
                    -(1.0-g_s)*a_s * qbc(3, i, j, k)  & ! south  gamma
                    +(1.0-g_n)*a_n * qbc(4, i, j, k)  & ! north  gamma
                    -(1.0-g_b)*a_b * qbc(5, i, j, k)  & ! bottom gamma
                    +(1.0-g_t)*a_t * qbc(6, i, j, k)  & ! top    gamma
            ) )  * sw                        ! 20 + 26 = 46
      ie(i,j,k) = ws(i,j,k) + delta
      res = res + dble(delta*delta) ! 3
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    res_l2 = res

    return
    end subroutine ps_diff_ee
