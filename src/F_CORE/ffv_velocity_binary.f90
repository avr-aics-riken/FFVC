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

!> @file   ffv_velocity_binary.f90
!! @brief  速度計算のルーチン群（バイナリモデル）
!! @author aics
!<

!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（1-UWD, 3-MUSCL）
!! @param [in]  v00       参照速度
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bv        BCindex C
!! @param [in]  bid       Cut ID
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [out] flop      浮動小数点演算数
!<
    subroutine pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, vf, bv, bid, vcs_coef, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, bix
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
    real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
    real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
    real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
    real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
    real                                                      ::  ck, dh, dh1, dh2, vcs, vcs_coef
    real                                                      ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b, cm1, cm2, ss_4
    real                                                      ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
    real                                                      ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
    real                                                      ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
    real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
    real                                                      ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bid
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1

    ! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
    vcs = vcs_coef
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_ref2 = 2.0*u_ref
    v_ref2 = 2.0*v_ref
    w_ref2 = 2.0*w_ref
    
    ck = 0.0
    b  = 0.0
    ss = 1.0
            
    if ( c_scheme == 1 ) then      !     1st order upwind
      ss = 0.0
    else if ( c_scheme == 3 ) then !     3rd order MUSCL
      ck = 1.0/3.0
      b  = (3.0-ck)/(1.0-ck)
    else
      write(*,*) 'out of scheme selection'
      stop
    endif
    
    ss_4 = 0.25*ss
    
    cm1 = 1.0 - ck
    cm2 = 1.0 + ck

    
! /4 + 12 = 8*4+13  = 45 flops
! /4 + 12 = 13*4+13 = 65 flops ! DP
! ループ内セットアップ   6 + 36 + 2 = 44 flops
! 対流項の各方向　( 7+6+7 + 78*3 + 12 ) * 3dir = 798
! 粘性項 36 + 36 + 12 = 84
! 粘性項置換部　67 ! DP 92
! Total : 45 + 44 + 798 + 84 = 971 ! DP 65 + 44 + 798 + 84 = 991

    flop = flop + dble(ix)*dble(jx)*dble(kx)*798.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*991.0d0 ! DP


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, vcs, b, ck, ss_4, ss, cm1, cm2) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, rei) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bix, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(cr, cl, acr, acl) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t, EX, EY, EZ)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      cnv_u = 0.0
      cnv_v = 0.0
      cnv_w = 0.0

      ! 各軸方向5点の変数ロード
      Ub2 = v(i  ,j  ,k-2, 1)
      Ub1 = v(i  ,j  ,k-1, 1)
      Us2 = v(i  ,j-2,k  , 1)
      Us1 = v(i  ,j-1,k  , 1)
      Uw2 = v(i-2,j  ,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)
      Ue2 = v(i+2,j  ,k  , 1)
      Un1 = v(i  ,j+1,k  , 1)
      Un2 = v(i  ,j+2,k  , 1)
      Ut1 = v(i  ,j  ,k+1, 1)
      Ut2 = v(i  ,j  ,k+2, 1)

      Vb2 = v(i  ,j  ,k-2, 2)
      Vb1 = v(i  ,j  ,k-1, 2)
      Vs2 = v(i  ,j-2,k  , 2)
      Vs1 = v(i  ,j-1,k  , 2)
      Vw2 = v(i-2,j  ,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2)
      Ve2 = v(i+2,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)
      Vn2 = v(i  ,j+2,k  , 2)
      Vt1 = v(i  ,j  ,k+1, 2)
      Vt2 = v(i  ,j  ,k+2, 2)

      Wb2 = v(i  ,j  ,k-2, 3)
      Wb1 = v(i  ,j  ,k-1, 3)
      Ws2 = v(i  ,j-2,k  , 3)
      Ws1 = v(i  ,j-1,k  , 3)
      Ww2 = v(i-2,j  ,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3)
      Wp0 = v(i  ,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3)
      We2 = v(i+2,j  ,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3)
      Wn2 = v(i  ,j+2,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3)
      Wt2 = v(i  ,j  ,k+2, 3)

      bvx = bv(i,j,k)
      bix = bid(i,j,k)

      ! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
      b_p = real(ibits(bvx, State, 1))

      ! 各方向に物体があれば，マスクはゼロ
      ! セル界面のフラグ b_?1={0.0-wall face / 1.0-fluid} <<= bix={0-fluid, ID-wall}
      b_w1 = 1.0
      b_e1 = 1.0
      b_s1 = 1.0
      b_n1 = 1.0
      b_b1 = 1.0
      b_t1 = 1.0
      if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w1 = 0.0
      if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e1 = 0.0
      if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s1 = 0.0
      if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n1 = 0.0
      if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b1 = 0.0
      if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t1 = 0.0


      ! (i,j,k)を基準にした遠い方向なので，隣接セルで判断
      b_w2 = 1.0
      b_e2 = 1.0
      b_s2 = 1.0
      b_n2 = 1.0
      b_b2 = 1.0
      b_t2 = 1.0
      if ( ibits(bid(i-1,j  ,k  ), bc_face_W, bitw_5) /= 0 ) b_w2 = 0.0
      if ( ibits(bid(i+1,j  ,k  ), bc_face_E, bitw_5) /= 0 ) b_e2 = 0.0
      if ( ibits(bid(i  ,j-1,k  ), bc_face_S, bitw_5) /= 0 ) b_s2 = 0.0
      if ( ibits(bid(i  ,j+1,k  ), bc_face_N, bitw_5) /= 0 ) b_n2 = 0.0
      if ( ibits(bid(i  ,j  ,k-1), bc_face_B, bitw_5) /= 0 ) b_b2 = 0.0
      if ( ibits(bid(i  ,j  ,k+1), bc_face_T, bitw_5) /= 0 ) b_t2 = 0.0


      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0


      ! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとす
      lmt_w = 1.0
      lmt_e = 1.0
      lmt_s = 1.0
      lmt_n = 1.0
      lmt_b = 1.0
      lmt_t = 1.0
      if ( ibits(bid(i-1, j  , k  ), vbc_uwd, 1) == 1) lmt_w = 0.0
      if ( ibits(bid(i+1, j  , k  ), vbc_uwd, 1) == 1) lmt_e = 0.0
      if ( ibits(bid(i  , j-1, k  ), vbc_uwd, 1) == 1) lmt_s = 0.0
      if ( ibits(bid(i  , j+1, k  ), vbc_uwd, 1) == 1) lmt_n = 0.0
      if ( ibits(bid(i  , j  , k-1), vbc_uwd, 1) == 1) lmt_b = 0.0
      if ( ibits(bid(i  , j  , k+1), vbc_uwd, 1) == 1) lmt_t = 0.0


      ! 界面速度（スタガード位置） > 24 flops
      UPe = vf(i  , j  , k  ,1)*b_e1 + u_ref*(1.0-b_e1)
      UPw = vf(i-1, j  , k  ,1)*b_w1 + u_ref*(1.0-b_w1)
      VPn = vf(i  , j  , k  ,2)*b_n1 + v_ref*(1.0-b_n1)
      VPs = vf(i  , j-1, k  ,2)*b_s1 + v_ref*(1.0-b_s1)
      WPt = vf(i  , j  , k  ,3)*b_t1 + w_ref*(1.0-b_t1)
      WPb = vf(i  , j  , k-1,3)*b_b1 + w_ref*(1.0-b_b1)


      ! セルセンターからの壁面修正速度 > 3 flops
      uq = u_ref2 - Up0
      vq = v_ref2 - Vp0
      wq = w_ref2 - Wp0
			
      ! X方向 ---------------------------------------
      
      ! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
      if ( (b_e2 == 0.0)  ) then  ! 7 flops
        Ue2 = u_ref2 - v(i+1,j  ,k  , 1)
        Ve2 = v_ref2 - v(i+1,j  ,k  , 2)
        We2 = w_ref2 - v(i+1,j  ,k  , 3)
      endif
      
      if ( b_e1 == 0.0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0.0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
      
      if ( (b_w2 == 0.0)  ) then ! 7 flops
        Uw2 = u_ref2 - v(i-1,j  ,k  , 1)
        Vw2 = v_ref2 - v(i-1,j  ,k  , 2)
        Ww2 = w_ref2 - v(i-1,j  ,k  , 3)
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = UPe - u_ref
      cl  = UPw - u_ref
      acr = abs(cr)
      acl = abs(cl)
      
      dv4 = Ue2-Ue1
      dv3 = Ue1-Up0
      dv2 = Up0-Uw1
      dv1 = Uw1-Uw2

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

      Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_e1
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_w1 ! > 4 + 4 + 36 + 5*4+7*2 = 78 flops

      dv4 = Ve2-Ve1
      dv3 = Ve1-Vp0
      dv2 = Vp0-Vw1
      dv1 = Vw1-Vw2
      
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

      Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_e1
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_w1

      dv4 = We2-We1
      dv3 = We1-Wp0
      dv2 = Wp0-Ww1
      dv1 = Ww1-Ww2
      
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

      Wrr = We1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_e1
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_w1
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
      cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
      cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! > 4*3 = 12 flops
			
      ! Y方向 ---------------------------------------
      
      if ( (b_n2 == 0.0)  ) then
        Un2 = u_ref2 - v(i  ,j+1,k  , 1)
        Vn2 = v_ref2 - v(i  ,j+1,k  , 2)
        Wn2 = w_ref2 - v(i  ,j+1,k  , 3)
      endif
      
      if ( b_n1 == 0.0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      endif
      
      if ( b_s1 == 0.0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      endif
      
      if ( (b_s2 == 0.0)  ) then
        Us2 = u_ref2 - v(i  ,j-1,k  , 1)
        Vs2 = v_ref2 - v(i  ,j-1,k  , 2)
        Ws2 = w_ref2 - v(i  ,j-1,k  , 3)
      endif
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = VPn - v_ref
      cl  = VPs - v_ref
      acr = abs(cr)
      acl = abs(cl)
      
      dv4 = Un2-Un1
      dv3 = Un1-Up0
      dv2 = Up0-Us1
      dv1 = Us1-Us2
      
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

      Urr = Un1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Us1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_n1
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_s1

      dv4 = Vn2-Vn1
      dv3 = Vn1-Vp0
      dv2 = Vp0-Vs1
      dv1 = Vs1-Vs2
      
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

      Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_n1
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_s1

      dv4 = Wn2-Wn1
      dv3 = Wn1-Wp0
      dv2 = Wp0-Ws1
      dv1 = Ws1-Ws2
      
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

      Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_n1
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_s1
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
      cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
      cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
      ! Z方向 ---------------------------------------
      
      ! 壁面の場合の参照速度の修正
      if ( (b_t2 == 0.0)  ) then
        Ut2 = u_ref2 - v(i  ,j  ,k+1, 1)
        Vt2 = v_ref2 - v(i  ,j  ,k+1, 2)
        Wt2 = w_ref2 - v(i  ,j  ,k+1, 3)
      end if
      
      if ( b_t1 == 0.0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0.0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if
      
      if ( (b_b2 == 0.0)  ) then
        Ub2 = u_ref2 - v(i  ,j  ,k-1, 1)
        Vb2 = v_ref2 - v(i  ,j  ,k-1, 2)
        Wb2 = w_ref2 - v(i  ,j  ,k-1, 3)
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = WPt - w_ref
      cl  = WPb - w_ref
      acr = abs(cr)
      acl = abs(cl)

      dv4 = Ut2-Ut1
      dv3 = Ut1-Up0
      dv2 = Up0-Ub1
      dv1 = Ub1-Ub2
      
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

      Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_t1
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_b1

      dv4 = Vt2-Vt1
      dv3 = Vt1-Vp0
      dv2 = Vp0-Vb1
      dv1 = Vb1-Vb2
      
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

      Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_t1
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_b1

      dv4 = Wt2-Wt1
      dv3 = Wt1-Wp0
      dv2 = Wp0-Wb1
      dv1 = Wb1-Wb2
      
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

      Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_t1
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_b1
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
      cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
      cnv_w = cnv_w + fw_r*c_t - fw_l*c_b


     ! 粘性項の計算　セル界面の剪断力を計算し，必要に応じて置換する
      EX = ( ( Ue1 - Up0 ) * c_e &
           - ( Up0 - Uw1 ) * c_w &
           + ( Un1 - Up0 ) * c_n &
           - ( Up0 - Us1 ) * c_s &
           + ( Ut1 - Up0 ) * c_t &
           - ( Up0 - Ub1 ) * c_b ) * dh1 * dh2

      EY = ( ( Ve1 - Vp0 ) * c_e &
           - ( Vp0 - Vw1 ) * c_w &
           + ( Vn1 - Vp0 ) * c_n &
           - ( Vp0 - Vs1 ) * c_s &
           + ( Vt1 - Vp0 ) * c_t &
           - ( Vp0 - Vb1 ) * c_b ) * dh1 * dh2

      EZ = ( ( We1 - Wp0 ) * c_e &
           - ( Wp0 - Ww1 ) * c_w &
           + ( Wn1 - Wp0 ) * c_n &
           - ( Wp0 - Ws1 ) * c_s &
           + ( Wt1 - Wp0 ) * c_t &
           - ( Wp0 - Wb1 ) * c_b ) * dh1 * dh2

      ! 対流項と粘性項の和 > 4*3 = 12 flops
      wv(i,j,k,1) = -cnv_u*dh1 + EX*vcs
      wv(i,j,k,2) = -cnv_v*dh1 + EY*vcs
      wv(i,j,k,3) = -cnv_w*dh1 + EZ*vcs
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine pvec_muscl


!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（1-UWD, 3-MUSCL）
!! @param [in]  v00       参照速度
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bv        BCindex C
!! @param [in]  bid       Cut ID
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [in]  Cs        定数CS
!! @param [in]  imodel    乱流モデル
!! @param [in]  nu        動粘性係数
!! @param [in]  rho       密度
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_muscl_les (wv, sz, g, dh, c_scheme, v00, rei, v, vf, bv, bid, vcs_coef, Cs, imodel, nu, rho, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, bix
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  ck, dh, dh1, dh2, vcs, vcs_coef
real                                                      ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2
real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b, cm1, cm2, ss_4
real                                                      ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
real                                                      ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
real                                                      ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
real                                                      ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(0:3)                                      ::  v00
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bid
integer                                                   ::  imodel
real                                                      ::  Cs, nu, rho, Cw
real                                                      ::  DUDX, DUDY, DUDZ
real                                                      ::  DVDX, DVDY, DVDZ
real                                                      ::  DWDX, DWDY, DWDZ
real                                                      ::  fs, DUDY_w, tauw, utau, yc, yp, up, min_h
real                                                      ::  S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS
real                                                      ::  W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW
real                                                      ::  S_1, S_2, S_3, W_1, W_2, W_3
real                                                      ::  S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d
real                                                      ::  Fcs, E_csm, Q_csm, Sijd2, nut
double precision                                          ::  EPS


ix = sz(1)
jx = sz(2)
kx = sz(3)

dh1 = 1.0/dh
dh2 = rei*dh1
EPS = 1.0d-10
Cw  = 0.325d0

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef

! 参照座標系の速度
u_ref = v00(1)
v_ref = v00(2)
w_ref = v00(3)
u_ref2 = 2.0*u_ref
v_ref2 = 2.0*v_ref
w_ref2 = 2.0*w_ref

ck = 0.0
b  = 0.0
ss = 1.0

if ( c_scheme == 1 ) then      !     1st order upwind
ss = 0.0
else if ( c_scheme == 3 ) then !     3rd order MUSCL
ck = 1.0/3.0
b  = (3.0-ck)/(1.0-ck)
else
write(*,*) 'out of scheme selection'
stop
endif

ss_4 = 0.25*ss

cm1 = 1.0 - ck
cm2 = 1.0 + ck


! /4 + 12 = 8*4+13  = 45 flops
! /4 + 12 = 13*4+13 = 65 flops ! DP
! ループ内セットアップ   6 + 36 + 2 = 44 flops
! 対流項の各方向　( 7+6+7 + 78*3 + 12 ) * 3dir = 798
! 粘性項 36 + 36 + 12 = 84
! 粘性項置換部　67 ! DP 92
! Total : 45 + 44 + 798 + 84 = 971 ! DP 65 + 44 + 798 + 84 = 991

flop = flop + dble(ix)*dble(jx)*dble(kx)*798.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*991.0d0 ! DP


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, vcs, b, ck, ss_4, ss, cm1, cm2) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, rei) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bix, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(cr, cl, acr, acl) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t, EX, EY, EZ) &
!$OMP FIRSTPRIVATE(imodel, EPS, Cs, nu, rho, Cw) &
!$OMP PRIVATE(DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ) &
!$OMP PRIVATE(S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS) &
!$OMP PRIVATE(W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW) &
!$OMP PRIVATE(S_1, S_2, S_3, W_1, W_2, W_3) &
!$OMP PRIVATE(S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d) &
!$OMP PRIVATE(Fcs, E_csm, Q_csm, Sijd2, nut) &
!$OMP PRIVATE(fs, DUDY_w, tauw, utau, yc, yp, up, min_h)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0

! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bvx = bv(i,j,k)
bix = bid(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bvx, State, 1))

! 各方向に物体があれば，マスクはゼロ
! セル界面のフラグ b_?1={0.0-wall face / 1.0-fluid} <<= bix={0-fluid, ID-wall}
b_w1 = 1.0
b_e1 = 1.0
b_s1 = 1.0
b_n1 = 1.0
b_b1 = 1.0
b_t1 = 1.0
if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w1 = 0.0
if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e1 = 0.0
if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s1 = 0.0
if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n1 = 0.0
if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b1 = 0.0
if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t1 = 0.0


! (i,j,k)を基準にした遠い方向なので，隣接セルで判断
b_w2 = 1.0
b_e2 = 1.0
b_s2 = 1.0
b_n2 = 1.0
b_b2 = 1.0
b_t2 = 1.0
if ( ibits(bid(i-1,j  ,k  ), bc_face_W, bitw_5) /= 0 ) b_w2 = 0.0
if ( ibits(bid(i+1,j  ,k  ), bc_face_E, bitw_5) /= 0 ) b_e2 = 0.0
if ( ibits(bid(i  ,j-1,k  ), bc_face_S, bitw_5) /= 0 ) b_s2 = 0.0
if ( ibits(bid(i  ,j+1,k  ), bc_face_N, bitw_5) /= 0 ) b_n2 = 0.0
if ( ibits(bid(i  ,j  ,k-1), bc_face_B, bitw_5) /= 0 ) b_b2 = 0.0
if ( ibits(bid(i  ,j  ,k+1), bc_face_T, bitw_5) /= 0 ) b_t2 = 0.0


! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
c_e = 1.0
c_w = 1.0
c_n = 1.0
c_s = 1.0
c_t = 1.0
c_b = 1.0
if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとす
lmt_w = 1.0
lmt_e = 1.0
lmt_s = 1.0
lmt_n = 1.0
lmt_b = 1.0
lmt_t = 1.0
if ( ibits(bid(i-1, j  , k  ), vbc_uwd, 1) == 1) lmt_w = 0.0
if ( ibits(bid(i+1, j  , k  ), vbc_uwd, 1) == 1) lmt_e = 0.0
if ( ibits(bid(i  , j-1, k  ), vbc_uwd, 1) == 1) lmt_s = 0.0
if ( ibits(bid(i  , j+1, k  ), vbc_uwd, 1) == 1) lmt_n = 0.0
if ( ibits(bid(i  , j  , k-1), vbc_uwd, 1) == 1) lmt_b = 0.0
if ( ibits(bid(i  , j  , k+1), vbc_uwd, 1) == 1) lmt_t = 0.0

!if ( (ibits(bv(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
!if ( (ibits(bv(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
!if ( (ibits(bv(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
!if ( (ibits(bv(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
!if ( (ibits(bv(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
!if ( (ibits(bv(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0


! 界面速度（スタガード位置） > 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1 + u_ref*(1.0-b_e1)
UPw = vf(i-1, j  , k  ,1)*b_w1 + u_ref*(1.0-b_w1)
VPn = vf(i  , j  , k  ,2)*b_n1 + v_ref*(1.0-b_n1)
VPs = vf(i  , j-1, k  ,2)*b_s1 + v_ref*(1.0-b_s1)
WPt = vf(i  , j  , k  ,3)*b_t1 + w_ref*(1.0-b_t1)
WPb = vf(i  , j  , k-1,3)*b_b1 + w_ref*(1.0-b_b1)


! セルセンターからの壁面修正速度 > 3 flops
uq = u_ref2 - Up0
vq = v_ref2 - Vp0
wq = w_ref2 - Wp0

! X方向 ---------------------------------------

! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
if ( (b_e2 == 0.0)  ) then  ! 7 flops
Ue2 = u_ref2 - v(i+1,j  ,k  , 1)
Ve2 = v_ref2 - v(i+1,j  ,k  , 2)
We2 = w_ref2 - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
Ue1 = uq
Ve1 = vq
We1 = wq
endif

if ( b_w1 == 0.0 ) then
Uw1 = uq
Vw1 = vq
Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then ! 7 flops
Uw2 = u_ref2 - v(i-1,j  ,k  , 1)
Vw2 = v_ref2 - v(i-1,j  ,k  , 2)
Ww2 = w_ref2 - v(i-1,j  ,k  , 3)
end if

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = UPe - u_ref
cl  = UPw - u_ref
acr = abs(cr)
acl = abs(cl)

dv4 = Ue2-Ue1
dv3 = Ue1-Up0
dv2 = Up0-Uw1
dv1 = Uw1-Uw2

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

Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_e1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_w1 ! > 4 + 4 + 36 + 5*4+7*2 = 78 flops

dv4 = Ve2-Ve1
dv3 = Ve1-Vp0
dv2 = Vp0-Vw1
dv1 = Vw1-Vw2

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

Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_e1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_w1

dv4 = We2-We1
dv3 = We1-Wp0
dv2 = Wp0-Ww1
dv1 = Ww1-Ww2

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

Wrr = We1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_e1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_w1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! > 4*3 = 12 flops

! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
Un2 = u_ref2 - v(i  ,j+1,k  , 1)
Vn2 = v_ref2 - v(i  ,j+1,k  , 2)
Wn2 = w_ref2 - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
Un1 = uq
Vn1 = vq
Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
Us1 = uq
Vs1 = vq
Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
Us2 = u_ref2 - v(i  ,j-1,k  , 1)
Vs2 = v_ref2 - v(i  ,j-1,k  , 2)
Ws2 = w_ref2 - v(i  ,j-1,k  , 3)
endif

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = VPn - v_ref
cl  = VPs - v_ref
acr = abs(cr)
acl = abs(cl)

dv4 = Un2-Un1
dv3 = Un1-Up0
dv2 = Up0-Us1
dv1 = Us1-Us2

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

Urr = Un1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Us1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_n1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_s1

dv4 = Vn2-Vn1
dv3 = Vn1-Vp0
dv2 = Vp0-Vs1
dv1 = Vs1-Vs2

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

Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_n1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_s1

dv4 = Wn2-Wn1
dv3 = Wn1-Wp0
dv2 = Wp0-Ws1
dv1 = Ws1-Ws2

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

Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_n1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_s1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
cnv_w = cnv_w + fw_r*c_n - fw_l*c_s

! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
Ut2 = u_ref2 - v(i  ,j  ,k+1, 1)
Vt2 = v_ref2 - v(i  ,j  ,k+1, 2)
Wt2 = w_ref2 - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
Ut1 = uq
Vt1 = vq
Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
Ub1 = uq
Vb1 = vq
Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
Ub2 = u_ref2 - v(i  ,j  ,k-1, 1)
Vb2 = v_ref2 - v(i  ,j  ,k-1, 2)
Wb2 = w_ref2 - v(i  ,j  ,k-1, 3)
end if

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = WPt - w_ref
cl  = WPb - w_ref
acr = abs(cr)
acl = abs(cl)

dv4 = Ut2-Ut1
dv3 = Ut1-Up0
dv2 = Up0-Ub1
dv1 = Ub1-Ub2

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

Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_t1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_b1

dv4 = Vt2-Vt1
dv3 = Vt1-Vp0
dv2 = Vp0-Vb1
dv1 = Vb1-Vb2

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

Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_t1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_b1

dv4 = Wt2-Wt1
dv3 = Wt1-Wp0
dv2 = Wp0-Wb1
dv1 = Wb1-Wb2

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

Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_t1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_b1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
cnv_w = cnv_w + fw_r*c_t - fw_l*c_b


! 粘性項の評価
DUDX = ( Ue1 - Uw1 ) * dh1 * 0.5d0
DUDY = ( Un1 - Us1 ) * dh1 * 0.5d0
DUDZ = ( Ut1 - Ub1 ) * dh1 * 0.5d0

DVDX = ( Ve1 - Vw1 ) * dh1 * 0.5d0
DVDY = ( Vn1 - Vs1 ) * dh1 * 0.5d0
DVDZ = ( Vt1 - Vb1 ) * dh1 * 0.5d0

DWDX = ( We1 - Ww1 ) * dh1 * 0.5d0
DWDY = ( We1 - Ws1 ) * dh1 * 0.5d0
DWDZ = ( We1 - Wb1 ) * dh1 * 0.5d0


!----------- Sij (strain rate tensor, 歪み速度テンソル)
S11 = DUDX
S12 = (DVDX + DUDY) * 0.5d0
S13 = (DWDX + DUDZ) * 0.5d0

S21 = S12
S22 = DVDY
S23 = (DWDY + DVDZ) * 0.5d0

S31 = S13
S32 = S23
S33 = DWDZ

SSS = DSQRT( 2.0d0 * (S11*S11 + S22*S22 + S33*S33) &
           + 4.0d0 * (S12*S12 + S23*S23 + S31*S31) )


!----------- Wij (vorticity tensor, 渦度テンソル)
W11 = 0.0d0
W12 = (DVDX - DUDY) * 0.5d0
W13 = (DWDX - DUDZ) * 0.5d0

W21 = (DUDY - DVDX) * 0.5d0
W22 = 0.0d0
W23 = (DWDY - DVDZ) * 0.5d0

W31 = (DUDZ - DWDX) * 0.5d0
W32 = (DVDZ - DWDY) * 0.5d0
W33 = 0.0d0

WWW = DSQRT( 2.0d0 * (W12*W12 + W13*W13 + W21*W21 &
                    + W23*W23 + W31*W31 + W32*W32 ) )


!----------- Smagorinsky model
if (imodel == 1) then
yc     = (j - 0.5d0)*dh
min_h  = min(yc, dh*jx - yc)
DUDY_w = v(i, jx, k, 1) * (dh1*2.0d0)
tauw   = (nu * rho) * abs(DUDY_w)
utau   = sqrt(tauw/rho)
yp     = min_h * utau / nu
up     = Up0 / utau
fs     = 1.0d0 - exp(-yp/26.0d0)
nut = (Cs * fs * dh) * (Cs * fs * dh) * SSS
end if


!----------- CSM model
if (imodel == 2) then
Q_csm   = (WWW * WWW - SSS * SSS) * 0.25d0
E_csm   = (WWW * WWW + SSS * SSS) * 0.25d0
Fcs     = Q_csm / (E_csm + EPS)

nut = (1.0d0/22.0d0) * abs(Fcs) * sqrt(abs(Fcs)) * (1.0d0 - Fcs) * (dh*dh) * SSS

end if


!----------- WALE model
if (imodel == 3) then

S11d = 0.5d0*(DUDX*DUDX + DUDX*DUDX) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S12d = 0.5d0*(DUDY*DUDY + DVDX*DVDX)
S13d = 0.5d0*(DUDZ*DUDZ + DWDX*DWDX)
S21d = S12d
S22d = 0.5d0*(DVDY*DVDY + DVDY*DVDY) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S23d = 0.5d0*(DVDZ*DVDZ + DWDY*DWDY)
S31d = S13d
S32d = S23d
S33d = 0.5d0*(DWDZ*DWDZ + DWDZ*DWDZ) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)

Sijd2 =  S11d*S11d + S21d*S21d + S31d*S31d &
       + S12d*S12d + S22d*S22d + S32d*S32d &
       + S13d*S13d + S23d*S23d + S33d*S33d

nut = (Cw * dh)* (Cw * dh) * (Sijd2)**(3.0d0/2.0d0)  &
    / ( (SSS * SSS)**(5.0d0/2.0d0) + (Sijd2)**(5.0d0/4.0d0) )

end if
!===============


! 粘性項の計算　セル界面の剪断力を計算し，必要に応じて置換する
EX = ( ( Ue1 - Up0 ) * c_e &
     - ( Up0 - Uw1 ) * c_w &
     + ( Un1 - Up0 ) * c_n &
     - ( Up0 - Us1 ) * c_s &
     + ( Ut1 - Up0 ) * c_t &
     - ( Up0 - Ub1 ) * c_b ) * dh1 * dh2

EY = ( ( Ve1 - Vp0 ) * c_e &
     - ( Vp0 - Vw1 ) * c_w &
     + ( Vn1 - Vp0 ) * c_n &
     - ( Vp0 - Vs1 ) * c_s &
     + ( Vt1 - Vp0 ) * c_t &
     - ( Vp0 - Vb1 ) * c_b ) * dh1 * dh2

EZ = ( ( We1 - Wp0 ) * c_e &
     - ( Wp0 - Ww1 ) * c_w &
     + ( Wn1 - Wp0 ) * c_n &
     - ( Wp0 - Ws1 ) * c_s &
     + ( Wt1 - Wp0 ) * c_t &
     - ( Wp0 - Wb1 ) * c_b ) * dh1 * dh2


! 対流項と粘性項の和 > 4*3 = 12 flops
wv(i, j, k, 1) = -cnv_u*dh1 + EX * ( 1.0d0 + nut/nu ) * vcs
wv(i, j, k, 2) = -cnv_v*dh1 + EY * ( 1.0d0 + nut/nu ) * vcs
wv(i, j, k, 3) = -cnv_w*dh1 + EZ * ( 1.0d0 + nut/nu ) * vcs

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pvec_muscl_les


!> ********************************************************************
!! @brief 次ステップのセルセンター，フェイスの速度と発散値を更新
!! @param [out] v    n+1時刻のセルセンター速度ベクトル
!! @param [in]  vf   n+1時刻のセルフェイス速度ベクトル
!! @param [out] div  \sum {u^{n+1}}
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dt   時間積分幅
!! @param [in]  dh   格子幅
!! @param [in]  vc   セルセンター疑似速度ベクトル
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [in]  bv   BCindex C
!! @param [out] flop 浮動小数点演算数
!! @note
!!    - actvのマスクはSPEC_VEL/OUTFLOWの参照セルをマスクしないようにbvを使う
!<
    subroutine update_vec (v, vf, div, sz, g, dt, dh, vc, p, bp, bv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, bvx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  dh, dt, dd, actv
    real                                                      ::  pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt
    real                                                      ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                      ::  Uef, Uwf, Vnf, Vsf, Wtf, Wbf
    real                                                      ::  c1, c2, c3, c4, c5, c6
    real                                                      ::  N_e, N_w, N_n, N_s, N_t, N_b
    real                                                      ::  D_e, D_w, D_n, D_s, D_t, D_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, vf
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, p
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp, bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dd = dt/dh

    flop = flop + dble(ix)*dble(jx)*dble(kx)*106.0 + 8.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(bpx, bvx, actv) &
!$OMP PRIVATE(c1, c2, c3, c4, c5, c6) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(D_e, D_w, D_n, D_s, D_t, D_b) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Uef, Uwf, Vnf, Vsf, Wtf, Wbf) &
!$OMP PRIVATE(pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dd)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      bpx = bp(i,j,k)
      bvx = bv(i,j,k)
      actv = real(ibits(bvx, State,  1))

      c1 = 1.0
      c2 = 1.0
      c3 = 1.0
      c4 = 1.0
      c5 = 1.0
      c6 = 1.0

      ! Neumann条件のとき，0.0 > 6 flop
      ! 物体があればノイマン条件なので，セルフェイスマスクとしても利用
      N_w = real(ibits(bpx, bc_n_W, 1))  ! w
      N_e = real(ibits(bpx, bc_n_E, 1))  ! e
      N_s = real(ibits(bpx, bc_n_S, 1))  ! s
      N_n = real(ibits(bpx, bc_n_N, 1))  ! n
      N_b = real(ibits(bpx, bc_n_B, 1))  ! b
      N_t = real(ibits(bpx, bc_n_T, 1))  ! t

      ! Dirichlet(p=0)条件のとき 0.0 > 6 flop
      D_w = real(ibits(bpx, bc_d_W, 1))  ! w
      D_e = real(ibits(bpx, bc_d_E, 1))  ! e
      D_s = real(ibits(bpx, bc_d_S, 1))  ! s
      D_n = real(ibits(bpx, bc_d_N, 1))  ! n
      D_b = real(ibits(bpx, bc_d_B, 1))  ! b
      D_t = real(ibits(bpx, bc_d_T, 1))  ! t


      ! 疑似ベクトル
      Uw0 = vc(i-1,j  ,k  , 1)
      Up0 = vc(i  ,j  ,k  , 1)
      Ue0 = vc(i+1,j  ,k  , 1)

      Vs0 = vc(i  ,j-1,k  , 2)
      Vp0 = vc(i  ,j  ,k  , 2)
      Vn0 = vc(i  ,j+1,k  , 2)

      Wb0 = vc(i  ,j  ,k-1, 3)
      Wp0 = vc(i  ,j  ,k  , 3)
      Wt0 = vc(i  ,j  ,k+1, 3)

      Uw = 0.5*( Up0 + Uw0 )*N_w ! 18 flop
      Ue = 0.5*( Up0 + Ue0 )*N_e
      Vs = 0.5*( Vp0 + Vs0 )*N_s
      Vn = 0.5*( Vp0 + Vn0 )*N_n
      Wb = 0.5*( Wp0 + Wb0 )*N_b
      Wt = 0.5*( Wp0 + Wt0 )*N_t

      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c1 = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c2 = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c3 = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c4 = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c5 = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c6 = 0.0
      
      ! 圧力勾配 18flop
      pc  = p(i,  j,  k  )
      pxw = (-p(i-1,j  ,k  )*D_w + pc) * N_w
      pxe = ( p(i+1,j  ,k  )*D_e - pc) * N_e
      pys = (-p(i  ,j-1,k  )*D_s + pc) * N_s
      pyn = ( p(i  ,j+1,k  )*D_n - pc) * N_n
      pzb = (-p(i  ,j  ,k-1)*D_b + pc) * N_b
      pzt = ( p(i  ,j  ,k+1)*D_t - pc) * N_t
      px = 0.5*(pxe + pxw)
      py = 0.5*(pyn + pys)
      pz = 0.5*(pzt + pzb)

      ! セルフェイス VBCの寄与と壁面の影響は除外 18flop
      Uwf = (Uw - dd * pxw) * c1 * N_w
      Uef = (Ue - dd * pxe) * c2 * N_e
      Vsf = (Vs - dd * pys) * c3 * N_s
      Vnf = (Vn - dd * pyn) * c4 * N_n
      Wbf = (Wb - dd * pzb) * c5 * N_b
      Wtf = (Wt - dd * pzt) * c6 * N_t

      ! i=1...ix >> vfは0...ixの範囲をカバーするので，通信不要 6flop
      vf(i-1,j  ,k  ,1) = Uwf
      vf(i  ,j  ,k  ,1) = Uef
      vf(i  ,j-1,k  ,2) = Vsf
      vf(i  ,j  ,k  ,2) = Vnf
      vf(i  ,j  ,k-1,3) = Wbf
      vf(i  ,j  ,k  ,3) = Wtf

      div(i,j,k) = (Uef - Uwf + Vnf - Vsf + Wtf - Wbf) * actv ! 6flop
      
      ! セルセンタの速度更新 9flop
      v(i,j,k,1) = ( Up0-px*dd ) * actv
      v(i,j,k,2) = ( Vp0-py*dd ) * actv
      v(i,j,k,3) = ( Wp0-pz*dd ) * actv

    end do
    end do
    end do
!$OMP END DO

!$OMP END PARALLEL

    return
    end subroutine update_vec



!> ********************************************************************
!! @brief 次ステップのセルセンターの速度を更新
!! @param [out] v        n+1時刻の速度ベクトル
!! @param [out] vf       n+1時刻のセルフェイス速度ベクトル
!! @param [out] div      \sum {u^{n+1}}
!! @param [in]  sz       配列長
!! @param [in]  g        ガイドセル長
!! @param [in]  dt       時間積分幅
!! @param [in]  dh       格子幅
!! @param [in]  vc       セルセンター疑似速度ベクトル
!! @param [in]  p        圧力
!! @param [in]  bp       BCindex P
!! @param [in]  bv       BCindex C
!! @param [in]  bid      Cut ID
!! @param [out] flop     浮動小数点演算数
!! @param [in   c_scheme 対流項スキーム
!! @note
!!    - actvのマスクはSPEC_VEL/OUTFLOWの参照セルをマスクしないようにbvを使う
!<
subroutine update_vec4 (v, vf, div, sz, g, dt, dh, vc, p, bp, bv, bid, flop, c_scheme)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, bvx, c_scheme
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
real                                                      ::  dh, dt, dd, actv, ss, c1, c2, sw1, sw2
real                                                      ::  N_e, N_w, N_n, N_s, N_t, N_b
real                                                      ::  D_e, D_w, D_n, D_s, D_t, D_b
real                                                      ::  N_e2, N_w2, N_n2, N_s2, N_t2, N_b2
real                                                      ::  D_e2, D_w2, D_n2, D_s2, D_t2, D_b2
real                                                      ::  Ue2, Uw2, Vn2, Vs2, Wt2, Wb2
real                                                      ::  Ue1, Uw1, Vn1, Vs1, Wt1, Wb1
real                                                      ::  Up0, Vp0, Wp0
real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
real                                                      ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
real                                                      ::  pc, pe1, pe2, pw1, pw2, pn1, pn2, ps1, ps2, pt1, pt2, pb1, pb2
real                                                      ::  dp_e4, dp_w4, dp_n4, dp_s4, dp_t4, dp_b4
real                                                      ::  dp_e2, dp_w2, dp_n2, dp_s2, dp_t2, dp_b2
real                                                      ::  dpx, dpy, dpz, dp_e, dp_w, dp_n, dp_s, dp_t, dp_b
real                                                      ::  Uef, Uwf, Vnf, Vsf, Wtf, Wbf
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, vf
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, p
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp, bv, bid


c1 = 9.0/8.0
c2 = 1.0/24.0

if ( c_scheme == 2 ) then      ! 2nd order accuracy
ss = 0.0
else if ( c_scheme == 4 ) then ! 4th order accuracy
ss = 1.0
else
write(*, *) 'out of scheme selection'
stop
end if

ix = sz(1)
jx = sz(2)
kx = sz(3)
dd = dt/dh

flop = flop + dble(ix)*dble(jx)*dble(kx)*159.0 + 8.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(bpx, bvx, actv) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(D_e, D_w, D_n, D_s, D_t, D_b) &
!$OMP PRIVATE(N_e2, N_w2, N_n2, N_s2, N_t2, N_b2) &
!$OMP PRIVATE(D_e2, D_w2, D_n2, D_s2, D_t2, D_b2) &
!$OMP PRIVATE(Uw2, Uw1, Up0, Ue1, Ue2) &
!$OMP PRIVATE(Vs2, Vs1, Vp0, Vn1, Vn2) &
!$OMP PRIVATE(Wb2, Wb1, Wp0, Wt1, Wt2) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t, sw1, sw2) &
!$OMP PRIVATE(pc, pe1, pe2, pw1, pw2, pn1, pn2, ps1, ps2, pt1, pt2, pb1, pb2) &
!$OMP PRIVATE(dp_e4, dp_w4, dp_n4, dp_s4, dp_t4, dp_b4) &
!$OMP PRIVATE(dp_e2, dp_w2, dp_n2, dp_s2, dp_t2, dp_b2) &
!$OMP PRIVATE(dpx, dpy, dpz, dp_e, dp_w, dp_n, dp_s, dp_t, dp_b) &
!$OMP PRIVATE(Uef, Uwf, Vnf, Vsf, Wtf, Wbf) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dd, ss, c1, c2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
bpx = bp(i,j,k)
bvx = bv(i,j,k)
actv = real(ibits(bvx, State,  1))

c_w = 1.0
c_e = 1.0
c_s = 1.0
c_n = 1.0
c_b = 1.0
c_t = 1.0

! Neumann条件のとき，0.0 > 6 flop
! 物体があればノイマン条件なので，セルフェイスマスクとしても利用
N_w = real(ibits(bpx, bc_n_W, 1))  ! w
N_e = real(ibits(bpx, bc_n_E, 1))  ! e
N_s = real(ibits(bpx, bc_n_S, 1))  ! s
N_n = real(ibits(bpx, bc_n_N, 1))  ! n
N_b = real(ibits(bpx, bc_n_B, 1))  ! b
N_t = real(ibits(bpx, bc_n_T, 1))  ! t

N_w2 = real( ibits(bp(i-1, j,   k  ), bc_n_W, 1) )
N_e2 = real( ibits(bp(i+1, j,   k  ), bc_n_E, 1) )
N_s2 = real( ibits(bp(i,   j-1, k  ), bc_n_S, 1) )
N_n2 = real( ibits(bp(i,   j+1, k  ), bc_n_N, 1) )
N_b2 = real( ibits(bp(i,   j,   k-1), bc_n_B, 1) )
N_t2 = real( ibits(bp(i,   j,   k+1), bc_n_T, 1) )


! Dirichlet(p=0)条件のとき 0.0
D_w = real(ibits(bpx, bc_d_W, 1))  ! w
D_e = real(ibits(bpx, bc_d_E, 1))  ! e
D_s = real(ibits(bpx, bc_d_S, 1))  ! s
D_n = real(ibits(bpx, bc_d_N, 1))  ! n
D_b = real(ibits(bpx, bc_d_B, 1))  ! b
D_t = real(ibits(bpx, bc_d_T, 1))  ! t

D_w2 = real( ibits(bp(i-1, j,   k  ), bc_d_W, 1) )
D_e2 = real( ibits(bp(i+1, j,   k  ), bc_d_E, 1) )
D_s2 = real( ibits(bp(i,   j-1, k  ), bc_d_S, 1) )
D_n2 = real( ibits(bp(i,   j+1, k  ), bc_d_N, 1) )
D_b2 = real( ibits(bp(i,   j,   k-1), bc_d_B, 1) )
D_t2 = real( ibits(bp(i,   j,   k+1), bc_d_T, 1) )


! 疑似ベクトル
Uw2 = vc(i-2,j  ,k  , 1)
Uw1 = vc(i-1,j  ,k  , 1)
Up0 = vc(i  ,j  ,k  , 1)
Ue1 = vc(i+1,j  ,k  , 1)
Ue2 = vc(i+2,j  ,k  , 1)

Vs2 = vc(i  ,j-2,k  , 2)
Vs1 = vc(i  ,j-1,k  , 2)
Vp0 = vc(i  ,j  ,k  , 2)
Vn1 = vc(i  ,j+1,k  , 2)
Vn2 = vc(i  ,j+2,k  , 2)

Wb2 = vc(i  ,j  ,k-2, 3)
Wb1 = vc(i  ,j  ,k-1, 3)
Wp0 = vc(i  ,j  ,k  , 3)
Wt1 = vc(i  ,j  ,k+1, 3)
Wt2 = vc(i  ,j  ,k+2, 3)

Uw = 0.5*( Up0 + Uw1 )*N_w ! 18 flop
Ue = 0.5*( Up0 + Ue1 )*N_e
Vs = 0.5*( Vp0 + Vs1 )*N_s
Vn = 0.5*( Vp0 + Vn1 )*N_n
Wb = 0.5*( Wp0 + Wb1 )*N_b
Wt = 0.5*( Wp0 + Wt1 )*N_t


!      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0

! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避(lmt_?=0.0)、デフォルトは1.0
lmt_w = 1.0
lmt_e = 1.0
lmt_s = 1.0
lmt_n = 1.0
lmt_b = 1.0
lmt_t = 1.0
if ( ibits(bid(i-1, j,   k),   vbc_uwd, 1) == 1) lmt_w = 0.0
if ( ibits(bid(i+1, j,   k),   vbc_uwd, 1) == 1) lmt_e = 0.0
if ( ibits(bid(i,   j-1, k),   vbc_uwd, 1) == 1) lmt_s = 0.0
if ( ibits(bid(i,   j+1, k),   vbc_uwd, 1) == 1) lmt_n = 0.0
if ( ibits(bid(i,   j,   k-1), vbc_uwd, 1) == 1) lmt_b = 0.0
if ( ibits(bid(i,   j,   k+1), vbc_uwd, 1) == 1) lmt_t = 0.0


!---- target cell >> 24 flops
pc = p(i, j, k)

!---- x-direction
pe1 = p(i+1, j, k) * D_e  * N_e
pe2 = p(i+2, j, k) * D_e2 * N_e2
pw1 = p(i-1, j, k) * D_w  * N_w
pw2 = p(i-2, j, k) * D_w2 * N_w2

!---- y-direction
pn1 = p(i, j+1, k) * D_n  * N_n
pn2 = p(i, j+2, k) * D_n2 * N_n2
ps1 = p(i, j-1, k) * D_s  * N_s
ps2 = p(i, j-2, k) * D_s2 * N_s2

!---- z-direction
pt1 = p(i, j, k+1) * D_t  * N_t
pt2 = p(i, j, k+2) * D_t2 * N_t2
pb1 = p(i, j, k-1) * D_b  * N_b
pb2 = p(i, j, k-2) * D_b2 * N_b2

!---- pressure difference for each direction >> 26*3 flops
! N_e = N_w = 0: wall face, 1: fluid face
sw1   = N_e * lmt_e * ss
sw2   = N_w * lmt_w * ss
dp_e2 = pe1 - pc
dp_w2 = pc  - pw1
dp_e4 = c1 * dp_e2 - c2 * (pe2 - pw1)
dp_w4 = c1 * dp_w2 - c2 * (pe1 - pw2)
dp_e  = dp_e4 * sw1 + dp_e2 * (1.0 - sw1)
dp_w  = dp_w4 * sw2 + dp_w2 * (1.0 - sw2)
dpx   = 0.5d0 * ( dp_e * c_e + dp_w * c_w )

! N_n = N_s = 0: wall face, 1: fluid face
sw1   = N_n * lmt_n * ss
sw2   = N_s * lmt_s * ss
dp_n2 = pn1 - pc
dp_s2 = pc  - ps1
dp_n4 = c1 * dp_n2 - c2 * (pn2 - ps1)
dp_s4 = c1 * dp_s2 - c2 * (pn1 - ps2)
dp_n  = dp_n4 * sw1 + dp_n2 * (1.0 - sw1)
dp_s  = dp_s4 * sw2 + dp_s2 * (1.0 - sw2)
dpy   = 0.5d0 * ( dp_n * c_n + dp_s * c_s )

! N_t = N_b = 0: wall face, 1: fluid face
sw1   = N_t * lmt_t * ss
sw2   = N_b * lmt_b * ss
dp_t2 = pt1 - pc
dp_b2 = pc  - pb1
dp_t4 = c1 * dp_t2 - c2 * (pt2 - pb1)
dp_b4 = c1 * dp_b2 - c2 * (pt1 - pb2)
dp_t  = dp_t4 * sw1 + dp_t2 * (1.0 - sw1)
dp_b  = dp_b4 * sw2 + dp_b2 * (1.0 - sw2)
dpz   = 0.5d0 * ( dp_t * c_t + dp_b * c_b )


! セルフェイス VBCの寄与と壁面の影響は除外 24flop
Uwf = (Uw - dd * dp_w) * c_w * N_w
Uef = (Ue - dd * dp_e) * c_e * N_e
Vsf = (Vs - dd * dp_s) * c_s * N_s
Vnf = (Vn - dd * dp_n) * c_n * N_n
Wbf = (Wb - dd * dp_b) * c_b * N_b
Wtf = (Wt - dd * dp_t) * c_t * N_t

! i=1...ix >> vfは0...ixの範囲をカバーするので，通信不要
vf(i-1,j  ,k  ,1) = Uwf
vf(i  ,j  ,k  ,1) = Uef
vf(i  ,j-1,k  ,2) = Vsf
vf(i  ,j  ,k  ,2) = Vnf
vf(i  ,j  ,k-1,3) = Wbf
vf(i  ,j  ,k  ,3) = Wtf

div(i,j,k) = (Uef - Uwf + Vnf - Vsf + Wtf - Wbf) * actv ! 6flop

! セルセンタの速度更新 9flop
v(i, j, k, 1) = ( Up0 - dpx * dd ) * actv
v(i, j, k, 2) = ( Vp0 - dpy * dd ) * actv
v(i, j, k, 3) = ( Wp0 - dpz * dd ) * actv

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine update_vec4


!> ********************************************************************
!! @brief 速度の発散に使う\sum{u^*}を計算する
!! @param [out]    div  速度の和
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     vc   セルセンター疑似ベクトル
!! @param [in]     bv   BCindex C
!! @param [in]     bid  Cut ID
!! @param [in,out] flop 浮動小数点演算数
!<
    subroutine divergence_cc (div, sz, g, vc, bv, bid, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx, bix
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb, actv
    real                                                      ::  Ue0, Uw0, Up0, Vn0, Vs0, Vp0, Wb0, Wt0, Wp0
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bid

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! loop : 1 + 42 + 13
    flop  = flop + dble(ix)*dble(jx)*dble(kx)*37.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(bvx, actv, bix) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue0, Uw0, Up0, Vn0, Vs0, Vp0, Wb0, Wt0, Wp0) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bvx = bv(i,j,k)
      bix = bid(i,j,k)
      actv= real(ibits(bvx, State, 1))
      
      ! 各セルセンター位置の変数ロード
      Uw0 = vc(i-1,j  ,k  , 1)
      Up0 = vc(i  ,j  ,k  , 1)
      Ue0 = vc(i+1,j  ,k  , 1)

      Vs0 = vc(i  ,j-1,k  , 2)
      Vp0 = vc(i  ,j  ,k  , 2)
      Vn0 = vc(i  ,j+1,k  , 2)

      Wb0 = vc(i  ,j  ,k-1, 3)
      Wp0 = vc(i  ,j  ,k  , 3)
      Wt0 = vc(i  ,j  ,k+1, 3)

      ! 0-solid / 1-fluid
b_w = 1.0
b_e = 1.0
b_s = 1.0
b_n = 1.0
b_b = 1.0
b_t = 1.0
if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w = 0.0
if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e = 0.0
if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s = 0.0
if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n = 0.0
if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b = 0.0
if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t = 0.0

      Uw = 0.5*( Up0 + Uw0 )*b_w
      Ue = 0.5*( Up0 + Ue0 )*b_e
      Vs = 0.5*( Vp0 + Vs0 )*b_s
      Vn = 0.5*( Vp0 + Vn0 )*b_n
      Wb = 0.5*( Wp0 + Wb0 )*b_b
      Wt = 0.5*( Wp0 + Wt0 )*b_t

      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
      c_w = 1.0
      c_e = 1.0
      c_s = 1.0
      c_n = 1.0
      c_b = 1.0
      c_t = 1.0

      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      
      ! VBC面の影響をフラグで無効化 >> OBC_SPEC_VEL, OBC_WALL  13flops
      div(i,j,k) = ( Ue*c_e - Uw*c_w + Vn*c_n - Vs*c_s + Wt*c_t - Wb*c_b ) * actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine divergence_cc
    
!> ********************************************************************
!! @brief 疑似ベクトルの時間積分（Euler陽解法）
!! @param [in,out] vc   対流項と粘性項の和 > 疑似ベクトル
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dt   時間積分幅
!! @param [in]     v    速度ベクトル（n-step, collocated）
!! @param [in]     bd   BCindex B
!! @param [in,out] flop 浮動小数点演算数
!! @note ここのマスクはIDのこと，VSPEC, OUTFLOWの増分をキャンセルするため
!<
    subroutine euler_explicit (vc, sz, g, dt, v, bd, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  actv, dt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*8.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(actv) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      actv = dt * real(ibits(bd(i,j,k), State, 1))

      vc(i,j,k,1) = v(i,j,k,1) + vc(i,j,k,1)* actv
      vc(i,j,k,2) = v(i,j,k,2) + vc(i,j,k,2)* actv
      vc(i,j,k,3) = v(i,j,k,3) + vc(i,j,k,3)* actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine euler_explicit


!> ********************************************************************
!! @brief Smagorinsky Modelによる乱流渦粘性係数の計算，減衰関数を併用
!! @param vt 乱流渦粘性係数
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param re レイノルズ数
!! @param cs 須磨五輪スキー定数
!! @param v 速度ベクトル
!! @param bx BC index for V
!! @param vt_range 渦粘性係数の最小値と最大値
!! @param yp_range Y+の最小値と最大値
!! @param v00 参照速度
!! @note
!!    - vtmin, vtmax > vt_range(2)
!!    - ypmin, ypmax > yp_range(2)
!! @todo
!!    - 境界条件は必要か？
!! @note NOCHECK
!<
    subroutine eddy_viscosity (vt, sz, g, dh, re, cs, v, bx, vt_range, yp_range, v00)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, m
    integer, dimension(3)                                       ::  sz
    real, dimension(2)                                          ::  vt_range, yp_range
    real                                                        ::  dh, re, cs, yp, Ut, Ut0
    real                                                        ::  DUDX, DUDY, DUDZ
    real                                                        ::  DVDX, DVDY, DVDZ
    real                                                        ::  DWDX, DWDY, DWDZ
    real                                                        ::  d1, d2, ddd, c1, c2, c3, delta, dd
    real                                                        ::  fs, aaa, Vmag, dis, tw, up1
    real                                                        ::  u_ref, v_ref, w_ref, u1, u2, u3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  vt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vt_range(1) =  1.0e6
    vt_range(2) = -1.0e6
    yp_range(1) =  1.0e6
    yp_range(2) = -1.0e6
    Delta = (dh*dh*dh)**(0.333)
    dd    = 0.5/dh
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      DUDX=(v(i+1,j  ,k  , 1)-v(i-1,j  ,k  , 1)) * dd
      DUDY=(v(i  ,j+1,k  , 1)-v(i  ,j-1,k  , 1)) * dd
      DUDZ=(v(i  ,j  ,k+1, 1)-v(i  ,j  ,k-1, 1)) * dd

      DVDX=(v(i+1,j  ,k  , 2)-v(i-1,j  ,k  , 2)) * dd
      DVDY=(v(i  ,j+1,k  , 2)-v(i  ,j-1,k  , 2)) * dd
      DVDZ=(v(i  ,j  ,k+1, 2)-v(i  ,j  ,k-1, 2)) * dd

      DWDX=(v(i+1,j  ,k  , 3)-v(i-1,j  ,k  , 3)) * dd
      DWDY=(v(i  ,j+1,k  , 3)-v(i  ,j-1,k  , 3)) * dd
      DWDZ=(v(i  ,j  ,k+1, 3)-v(i  ,j  ,k-1, 3)) * dd

      D1  = DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ
      c1  = DUDY+DVDX
      c2  = DVDZ+DWDY
      c3  = DWDX+DUDZ
      D2  = c1*c1 + c2*c2 + c3*c3
      DDD = SQRT(2.0D0*D1 + D2)

!     ------- Y+ & fs--
      fs =1.0
      !AAA= bx(i,j,k) * bx(i,j,k)   &
      !    *bx(i,j,k) * bx(i,j,k)   &
      !    *bx(i,j,k) * bx(i,j,k)
      u1 = v(i,j,k,1) - u_ref
      u2 = v(i,j,k,2) - v_ref
      u3 = v(i,j,k,3) - w_ref
      Vmag = SQRT(u1*u1 + u2*u2 + u3*u3)

      IF ( (AAA < 1.0) .AND. (Vmag >= 0.001) ) THEN
        !DIS = AMIN1(bnd(i,j,k), bnd(i,j,k),    &
        !            bnd(i,j,k), bnd(i,j,k),    &
        !            bnd(i,j,k), bnd(i,j,k) ) * dh
        Tw = Vmag/(DIS*RE)
        Ut0= SQRT(Tw)
        YP = DIS*Ut0*RE
        Up1= Vmag
        
        DO M=1,5
          Ut = Ut0 + (Up1-Ut0*(ALOG(YP)/0.4+5.5)) / (Up1/Ut0+1.0/0.4)
          IF (ABS(Ut-Ut0) <= 1.0E-6) exit
          IF (Ut <= 0.010) Ut=0.0100
          Ut0=Ut
          YP =DIS*Ut0*RE
        end do
        
        fs = 1.0-exp(-YP/25.0)
        yp_range(1) = amin1(yp_range(1), YP)
        yp_range(2) = amax1(yp_range(2), YP)
      ENDIF

      VT(i,j,k) = ((Cs*fs*Delta)**2)*DDD*RE
      vt_range(1) = MIN(vt_range(1), VT(i,j,k))
      vt_range(2) = MAX(vt_range(2), VT(i,j,k))
    end do
    end do
    end do

    do k=1,kx
    do i=1,ix
      VT(i,0   ,k)=VT(i,1 ,k)
      VT(i,jx+1,k)=VT(i,jx,k)
    end do
    end do

    do k=1,kx
    do j=1,jx
      VT(0   ,j,k)=VT(1 ,j,k)
      VT(ix+1,j,k)=VT(ix,j,k)
    end do
    end do
 
    do j=1,jx
    do i=1,ix
      VT(i,j,0   )=VT(i,j,1 )
      VT(i,j,kx+1)=VT(i,j,kx)
    end do
    end do

    return
    end subroutine eddy_viscosity

!> ********************************************************************
!! @brief Adams-Bashforth法による疑似ベクトルの時間積分
!! @param[out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param v 速度ベクトル（n-step, collocated）
!! @param ab 前ステップの対流項（＋粘性項）の計算値
!! @param bd BCindex B
!! @param v00 参照速度
!! @param[out] flop
!! @note NOCHECK
!<
    subroutine ab2 (vc, sz, g, dt, v, ab, bd, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  actv, dt, ab_u, ab_v, ab_w, u_ref, v_ref, w_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v, ab
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + dble(ix*jx*kx) * 27.0d0
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      actv = real(ibits(bd(i,j,k), State, 1))
      
      ab_u = ab(i,j,k,1)
      ab_v = ab(i,j,k,2)
      ab_w = ab(i,j,k,3)
      
      ab(i,j,k,1) = vc(i,j,k,1)
      ab(i,j,k,2) = vc(i,j,k,2)
      ab(i,j,k,3) = vc(i,j,k,3)

      vc(i,j,k,1) = ( v(i,j,k,1) + dt*0.5*( 3.0*vc(i,j,k,1)-ab_u ) )* actv + (1.0-actv)*u_ref
      vc(i,j,k,2) = ( v(i,j,k,2) + dt*0.5*( 3.0*vc(i,j,k,2)-ab_v ) )* actv + (1.0-actv)*v_ref
      vc(i,j,k,3) = ( v(i,j,k,3) + dt*0.5*( 3.0*vc(i,j,k,3)-ab_w ) )* actv + (1.0-actv)*w_ref
    end do
    end do
    end do

    return
    end subroutine ab2
    

!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項 u \frac{\partial u}{\partial x}
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（2-Central_2nd, 4-Central_4th）
!! @param [in]  v00       参照速度
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bv        BCindex C
!! @param [in]  bid       Cut ID
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_central (wv, sz, g, dh, c_scheme, v00, rei, v, vf, bv, bid, vcs_coef, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, bix
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  dh, dh1, dh2, vcs, vcs_coef, ss, sw1, sw2, c1, c2
real                                                      ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2
real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
real                                                      ::  cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  uq, vq, wq
real                                                      ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
real                                                      ::  ufr2, ufl2, vfr2, vfl2, wfr2, wfl2
real                                                      ::  ufr4, ufl4, vfr4, vfl4, wfr4, wfl4
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(0:3)                                      ::  v00
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bid

ix = sz(1)
jx = sz(2)
kx = sz(3)

dh1= 1.0/dh
dh2= rei*dh1

c1 = 9.0/8.0
c2 = 1.0/24.0

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef


if ( c_scheme == 2 ) then      !     2nd
ss = 0.0
else if ( c_scheme == 4 ) then !     4th
ss = 1.0
else
write(*,*) 'out of scheme selection'
stop
endif


! 参照座標系の速度
u_ref = v00(1)
v_ref = v00(2)
w_ref = v00(3)
u_ref2 = 2.0*u_ref
v_ref2 = 2.0*v_ref
w_ref2 = 2.0*w_ref


flop = flop + dble(ix)*dble(jx)*dble(kx)*330.0d0 + 28.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*991.0d0 ! DP


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, vcs, ss, c1, c2) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, rei) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bix, uq, vq, wq, sw1, sw2) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(ufr2, ufl2, vfr2, vfl2, wfr2, wfl2) &
!$OMP PRIVATE(ufr4, ufl4, vfr4, vfl4, wfr4, wfl4) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t, EX, EY, EZ)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0
EX = 0.0
EY = 0.0
EZ = 0.0

! 変数のロード
! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bvx = bv(i,j,k)
bix = bid(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bvx, State, 1))

! 各方向に物体があれば，マスクはゼロ
! セル界面のフラグ b_?1={0.0-wall face / 1.0-fluid} <<= bix={0-fluid, ID-wall}
b_w1 = 1.0
b_e1 = 1.0
b_s1 = 1.0
b_n1 = 1.0
b_b1 = 1.0
b_t1 = 1.0
if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w1 = 0.0
if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e1 = 0.0
if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s1 = 0.0
if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n1 = 0.0
if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b1 = 0.0
if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t1 = 0.0


! (i,j,k)を基準にした遠い方向なので，隣接セルで判断
b_w2 = 1.0
b_e2 = 1.0
b_s2 = 1.0
b_n2 = 1.0
b_b2 = 1.0
b_t2 = 1.0
if ( ibits(bid(i-1,j  ,k  ), bc_face_W, bitw_5) /= 0 ) b_w2 = 0.0
if ( ibits(bid(i+1,j  ,k  ), bc_face_E, bitw_5) /= 0 ) b_e2 = 0.0
if ( ibits(bid(i  ,j-1,k  ), bc_face_S, bitw_5) /= 0 ) b_s2 = 0.0
if ( ibits(bid(i  ,j+1,k  ), bc_face_N, bitw_5) /= 0 ) b_n2 = 0.0
if ( ibits(bid(i  ,j  ,k-1), bc_face_B, bitw_5) /= 0 ) b_b2 = 0.0
if ( ibits(bid(i  ,j  ,k+1), bc_face_T, bitw_5) /= 0 ) b_t2 = 0.0


! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
c_e = 1.0
c_w = 1.0
c_n = 1.0
c_s = 1.0
c_t = 1.0
c_b = 1.0
if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0

! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避(lmt_?=0.0)、デフォルトは1.0
lmt_w = 1.0
lmt_e = 1.0
lmt_s = 1.0
lmt_n = 1.0
lmt_b = 1.0
lmt_t = 1.0
if ( ibits(bid(i-1, j  , k  ), vbc_uwd, 1) == 1) lmt_w = 0.0
if ( ibits(bid(i+1, j  , k  ), vbc_uwd, 1) == 1) lmt_e = 0.0
if ( ibits(bid(i  , j-1, k  ), vbc_uwd, 1) == 1) lmt_s = 0.0
if ( ibits(bid(i  , j+1, k  ), vbc_uwd, 1) == 1) lmt_n = 0.0
if ( ibits(bid(i  , j  , k-1), vbc_uwd, 1) == 1) lmt_b = 0.0
if ( ibits(bid(i  , j  , k+1), vbc_uwd, 1) == 1) lmt_t = 0.0


!if ( (ibits(bv(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
!if ( (ibits(bv(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
!if ( (ibits(bv(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
!if ( (ibits(bv(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
!if ( (ibits(bv(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
!if ( (ibits(bv(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0

! 界面速度（スタガード位置） 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1 + u_ref*(1.0-b_e1)
UPw = vf(i-1, j  , k  ,1)*b_w1 + u_ref*(1.0-b_w1)
VPn = vf(i  , j  , k  ,2)*b_n1 + v_ref*(1.0-b_n1)
VPs = vf(i  , j-1, k  ,2)*b_s1 + v_ref*(1.0-b_s1)
WPt = vf(i  , j  , k  ,3)*b_t1 + w_ref*(1.0-b_t1)
WPb = vf(i  , j  , k-1,3)*b_b1 + w_ref*(1.0-b_b1)

! セルセンターからの壁面修正速度 3 flops
uq = u_ref2 - Up0
vq = v_ref2 - Vp0
wq = w_ref2 - Wp0

! X方向 ---------------------------------------

! 速度指定の場合に参照先として，固体内にテンポラリに与えた値を使う
if ( (b_e2 == 0.0)  ) then
  Ue2 = u_ref2 - v(i+1,j  ,k  , 1)
  Ve2 = v_ref2 - v(i+1,j  ,k  , 2)
  We2 = w_ref2 - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
  Ue1 = uq
  Ve1 = vq
  We1 = wq
endif

if ( b_w1 == 0.0 ) then
  Uw1 = uq
  Vw1 = vq
  Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then
  Uw2 = u_ref2 - v(i-1,j  ,k  , 1)
  Vw2 = v_ref2 - v(i-1,j  ,k  , 2)
  Ww2 = w_ref2 - v(i-1,j  ,k  , 3)
end if

! flux u \frac{\partial u}{\partial x}
! 壁面がある場合　b_e1, b_w1 (0-wall face / 1-fluid)
! vspec/outflowを参照する場合(lmt_e, lmt_w)は二次精度へおとす
! 二次精度の時には、ssで強制
sw1 = b_e1 * lmt_e * ss ! 4 flops
sw2 = b_w1 * lmt_w * ss

ufr2 = Ue1-Up0 ! 20 flops
ufl2 = Up0-Uw1
ufr4 = c1 * ufr2 - c2 * (Ue2-Uw1)
ufl4 = c1 * ufl2 - c2 * (Ue1-Uw2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * dh1
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * dh1

vfr2 = Ve1-Vp0
vfl2 = Vp0-Vw1
vfr4 = c1 * vfr2 - c2 * (Ve2-Vw1)
vfl4 = c1 * vfl2 - c2 * (Ve1-Vw2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * dh1
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * dh1

wfr2 = We1-Wp0
wfl2 = Wp0-Ww1
wfr4 = c1 * wfr2 - c2 * (We2-Ww1)
wfl4 = c1 * wfl2 - c2 * (We1-Ww2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * dh1
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * dh1


! 流束の加算　平均勾配　VBCでない面の寄与のみを評価する 21 flops
cnv_u = cnv_u + 0.5*(Upe * fu_r * c_e + Upw * fu_l * c_w)
cnv_v = cnv_v + 0.5*(Upe * fv_r * c_e + Upw * fv_l * c_w)
cnv_w = cnv_w + 0.5*(Upe * fw_r * c_e + Upw * fw_l * c_w)

! 粘性項の加算 12 flops
EX = EX + (fu_r * c_e - fu_l * c_w)
EY = EY + (fv_r * c_e - fv_l * c_w)
EZ = EZ + (fw_r * c_e - fw_l * c_w)


! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
  Un2 = u_ref2 - v(i  ,j+1,k  , 1)
  Vn2 = v_ref2 - v(i  ,j+1,k  , 2)
  Wn2 = w_ref2 - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
  Un1 = uq
  Vn1 = vq
  Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
  Us1 = uq
  Vs1 = vq
  Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
  Us2 = u_ref2 - v(i  ,j-1,k  , 1)
  Vs2 = v_ref2 - v(i  ,j-1,k  , 2)
  Ws2 = w_ref2 - v(i  ,j-1,k  , 3)
endif

! flux
sw1 = b_n1 * lmt_n * ss
sw2 = b_s1 * lmt_s * ss

ufr2 = Un1-Up0
ufl2 = Up0-Us1
ufr4 = c1 * ufr2 - c2 * (Un2-Us1)
ufl4 = c1 * ufl2 - c2 * (Un1-Us2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * dh1
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * dh1

vfr2 = Vn1-Vp0
vfl2 = Vp0-Vs1
vfr4 = c1 * vfr2 - c2 * (Vn2-Vs1)
vfl4 = c1 * vfl2 - c2 * (Vn1-Vs2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * dh1
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * dh1

wfr2 = Wn1-Wp0
wfl2 = Wp0-Ws1
wfr4 = c1 * wfr2 - c2 * (Wn2-Ws1)
wfl4 = c1 * wfl2 - c2 * (Wn1-Ws2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * dh1
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * dh1


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Vpn * fu_r * c_n + Vps * fu_l * c_s)
cnv_v = cnv_v + 0.5*(Vpn * fv_r * c_n + Vps * fv_l * c_s)
cnv_w = cnv_w + 0.5*(Vpn * fw_r * c_n + Vps * fw_l * c_s)

! 粘性項の加算
EX = EX + (fu_r * c_n - fu_l * c_s)
EY = EY + (fv_r * c_n - fv_l * c_s)
EZ = EZ + (fw_r * c_n - fw_l * c_s)


! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
  Ut2 = u_ref2 - v(i  ,j  ,k+1, 1)
  Vt2 = v_ref2 - v(i  ,j  ,k+1, 2)
  Wt2 = w_ref2 - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
  Ut1 = uq
  Vt1 = vq
  Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
  Ub1 = uq
  Vb1 = vq
  Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
  Ub2 = u_ref2 - v(i  ,j  ,k-1, 1)
  Vb2 = v_ref2 - v(i  ,j  ,k-1, 2)
  Wb2 = w_ref2 - v(i  ,j  ,k-1, 3)
end if

! flux
sw1 = b_t1 * lmt_t * ss
sw2 = b_b1 * lmt_b * ss

ufr2 = Ut1-Up0
ufl2 = Up0-Ub1
ufr4 = c1 * ufr2 - c2 * (Ut2-Ub1)
ufl4 = c1 * ufl2 - c2 * (Ut1-Ub2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * dh1
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * dh1

vfr2 = Vt1-Vp0
vfl2 = Vp0-Vb1
vfr4 = c1 * vfr2 - c2 * (Vt2-Vb1)
vfl4 = c1 * vfl2 - c2 * (Vt1-Vb2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * dh1
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * dh1

wfr2 = Wt1-Wp0
wfl2 = Wp0-Wb1
wfr4 = c1 * wfr2 - c2 * (Wt2-Wb1)
wfl4 = c1 * wfl2 - c2 * (Wt1-Wb2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * dh1
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * dh1


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Wpt * fu_r * c_t + Wpb * fu_l * c_b)
cnv_v = cnv_v + 0.5*(Wpt * fv_r * c_t + Wpb * fv_l * c_b)
cnv_w = cnv_w + 0.5*(Wpt * fw_r * c_t + Wpb * fw_l * c_b)

! 粘性項の加算
EX = EX + (fu_r * c_t - fu_l * c_b)
EY = EY + (fv_r * c_t - fv_l * c_b)
EZ = EZ + (fw_r * c_t - fw_l * c_b)


! 対流項と粘性項の和 > 12 flops
wv(i,j,k,1) = -cnv_u + EX * vcs * dh2
wv(i,j,k,2) = -cnv_v + EY * vcs * dh2
wv(i,j,k,3) = -cnv_w + EZ * vcs * dh2
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pvec_central


!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項 u \frac{\partial u}{\partial x}
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（2-Central_2nd, 4-Central_4th）
!! @param [in]  v00       参照速度
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bv        BCindex C
!! @param [in]  bid       Cut ID
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [in]  Cs        定数CS
!! @param [in]  imodel    乱流モデル
!! @param [in]  nu        動粘性係数
!! @param [in]  rho       密度
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_central_les (wv, sz, g, dh, c_scheme, v00, rei, v, vf, bv, bid, vcs_coef, Cs, imodel, nu, rho, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, bix
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  dh, dh1, dh2, vcs, vcs_coef, ss, sw1, sw2, c1, c2
real                                                      ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2
real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
real                                                      ::  cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  uq, vq, wq
real                                                      ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
real                                                      ::  ufr2, ufl2, vfr2, vfl2, wfr2, wfl2
real                                                      ::  ufr4, ufl4, vfr4, vfl4, wfr4, wfl4
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(0:3)                                      ::  v00
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bid
integer                                                   :: imodel
real                                                      ::  Cs, nu, rho, Cw
real                                                      ::  DUDX, DUDY, DUDZ
real                                                      ::  DVDX, DVDY, DVDZ
real                                                      ::  DWDX, DWDY, DWDZ
real                                                      ::  fs, DUDY_w, tauw, utau, yc, yp, up, min_h
real                                                      ::  S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS
real                                                      ::  W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW
real                                                      ::  S_1, S_2, S_3, W_1, W_2, W_3
real                                                      ::  S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d
real                                                      ::  Fcs, E_csm, Q_csm, Sijd2, nut
double precision                                          :: EPS

ix = sz(1)
jx = sz(2)
kx = sz(3)

dh1 = 1.0/dh
dh2 = rei*dh1
EPS = 1.0d-10
Cw  = 0.325d0

c1 = 9.0/8.0
c2 = 1.0/24.0

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef


if ( c_scheme == 2 ) then      !     2nd
ss = 0.0
else if ( c_scheme == 4 ) then !     4th
ss = 1.0
else
write(*,*) 'out of scheme selection'
stop
endif


! 参照座標系の速度
u_ref = v00(1)
v_ref = v00(2)
w_ref = v00(3)
u_ref2 = 2.0*u_ref
v_ref2 = 2.0*v_ref
w_ref2 = 2.0*w_ref


flop = flop + dble(ix)*dble(jx)*dble(kx)*330.0d0 + 28.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*991.0d0 ! DP


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, vcs, ss, c1, c2) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, rei) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bix, uq, vq, wq, sw1, sw2) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(ufr2, ufl2, vfr2, vfl2, wfr2, wfl2) &
!$OMP PRIVATE(ufr4, ufl4, vfr4, vfl4, wfr4, wfl4) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t, EX, EY, EZ) &
!$OMP FIRSTPRIVATE(imodel, EPS, Cs, nu, rho, Cw) &
!$OMP PRIVATE(DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ) &
!$OMP PRIVATE(S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS) &
!$OMP PRIVATE(W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW) &
!$OMP PRIVATE(S_1, S_2, S_3, W_1, W_2, W_3) &
!$OMP PRIVATE(S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d) &
!$OMP PRIVATE(Fcs, E_csm, Q_csm, Sijd2, nut) &
!$OMP PRIVATE(fs, DUDY_w, tauw, utau, yc, yp, up, min_h)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0
EX = 0.0
EY = 0.0
EZ = 0.0

! 変数のロード
! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bvx = bv(i,j,k)
bix = bid(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bvx, State, 1))

! 各方向に物体があれば，マスクはゼロ
! セル界面のフラグ b_?1={0.0-wall face / 1.0-fluid} <<= bix={0-fluid, ID-wall}
b_w1 = 1.0
b_e1 = 1.0
b_s1 = 1.0
b_n1 = 1.0
b_b1 = 1.0
b_t1 = 1.0
if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w1 = 0.0
if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e1 = 0.0
if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s1 = 0.0
if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n1 = 0.0
if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b1 = 0.0
if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t1 = 0.0


! (i,j,k)を基準にした遠い方向なので，隣接セルで判断
b_w2 = 1.0
b_e2 = 1.0
b_s2 = 1.0
b_n2 = 1.0
b_b2 = 1.0
b_t2 = 1.0
if ( ibits(bid(i-1,j  ,k  ), bc_face_W, bitw_5) /= 0 ) b_w2 = 0.0
if ( ibits(bid(i+1,j  ,k  ), bc_face_E, bitw_5) /= 0 ) b_e2 = 0.0
if ( ibits(bid(i  ,j-1,k  ), bc_face_S, bitw_5) /= 0 ) b_s2 = 0.0
if ( ibits(bid(i  ,j+1,k  ), bc_face_N, bitw_5) /= 0 ) b_n2 = 0.0
if ( ibits(bid(i  ,j  ,k-1), bc_face_B, bitw_5) /= 0 ) b_b2 = 0.0
if ( ibits(bid(i  ,j  ,k+1), bc_face_T, bitw_5) /= 0 ) b_t2 = 0.0


! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
c_e = 1.0
c_w = 1.0
c_n = 1.0
c_s = 1.0
c_t = 1.0
c_b = 1.0
if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0

! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避(lmt_?=0.0)、デフォルトは1.0
lmt_w = 1.0
lmt_e = 1.0
lmt_s = 1.0
lmt_n = 1.0
lmt_b = 1.0
lmt_t = 1.0
if ( ibits(bid(i-1, j  , k  ), vbc_uwd, 1) == 1) lmt_w = 0.0
if ( ibits(bid(i+1, j  , k  ), vbc_uwd, 1) == 1) lmt_e = 0.0
if ( ibits(bid(i  , j-1, k  ), vbc_uwd, 1) == 1) lmt_s = 0.0
if ( ibits(bid(i  , j+1, k  ), vbc_uwd, 1) == 1) lmt_n = 0.0
if ( ibits(bid(i  , j  , k-1), vbc_uwd, 1) == 1) lmt_b = 0.0
if ( ibits(bid(i  , j  , k+1), vbc_uwd, 1) == 1) lmt_t = 0.0

!if ( (ibits(bv(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
!if ( (ibits(bv(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
!if ( (ibits(bv(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
!if ( (ibits(bv(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
!if ( (ibits(bv(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
!if ( (ibits(bv(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0

! 界面速度（スタガード位置） 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1 + u_ref*(1.0-b_e1)
UPw = vf(i-1, j  , k  ,1)*b_w1 + u_ref*(1.0-b_w1)
VPn = vf(i  , j  , k  ,2)*b_n1 + v_ref*(1.0-b_n1)
VPs = vf(i  , j-1, k  ,2)*b_s1 + v_ref*(1.0-b_s1)
WPt = vf(i  , j  , k  ,3)*b_t1 + w_ref*(1.0-b_t1)
WPb = vf(i  , j  , k-1,3)*b_b1 + w_ref*(1.0-b_b1)

! セルセンターからの壁面修正速度 3 flops
uq = u_ref2 - Up0
vq = v_ref2 - Vp0
wq = w_ref2 - Wp0

! X方向 ---------------------------------------

! 速度指定の場合に参照先として，固体内にテンポラリに与えた値を使う
if ( (b_e2 == 0.0)  ) then
Ue2 = u_ref2 - v(i+1,j  ,k  , 1)
Ve2 = v_ref2 - v(i+1,j  ,k  , 2)
We2 = w_ref2 - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
Ue1 = uq
Ve1 = vq
We1 = wq
endif

if ( b_w1 == 0.0 ) then
Uw1 = uq
Vw1 = vq
Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then
Uw2 = u_ref2 - v(i-1,j  ,k  , 1)
Vw2 = v_ref2 - v(i-1,j  ,k  , 2)
Ww2 = w_ref2 - v(i-1,j  ,k  , 3)
end if

! flux u \frac{\partial u}{\partial x}
! 壁面がある場合　b_e1, b_w1 (0-wall face / 1-fluid)
! vspec/outflowを参照する場合(lmt_e, lmt_w)は二次精度へおとす
! 二次精度の時には、ssで強制
sw1 = b_e1 * lmt_e * ss ! 4 flops
sw2 = b_w1 * lmt_w * ss

ufr2 = Ue1-Up0 ! 20 flops
ufl2 = Up0-Uw1
ufr4 = c1 * ufr2 - c2 * (Ue2-Uw1)
ufl4 = c1 * ufl2 - c2 * (Ue1-Uw2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * dh1
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * dh1
DUDX = 0.5d0 * ( fu_r * c_e + fu_l * c_w )


vfr2 = Ve1-Vp0
vfl2 = Vp0-Vw1
vfr4 = c1 * vfr2 - c2 * (Ve2-Vw1)
vfl4 = c1 * vfl2 - c2 * (Ve1-Vw2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * dh1
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * dh1
DVDX = 0.5d0 * ( fv_r * c_e + fv_l * c_w )


wfr2 = We1-Wp0
wfl2 = Wp0-Ww1
wfr4 = c1 * wfr2 - c2 * (We2-Ww1)
wfl4 = c1 * wfl2 - c2 * (We1-Ww2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * dh1
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * dh1
DWDX = 0.5d0 * ( fw_r * c_e + fw_l * c_w )


! 流束の加算　平均勾配　VBCでない面の寄与のみを評価する 21 flops
cnv_u = cnv_u + 0.5*(Upe * fu_r * c_e + Upw * fu_l * c_w)
cnv_v = cnv_v + 0.5*(Upe * fv_r * c_e + Upw * fv_l * c_w)
cnv_w = cnv_w + 0.5*(Upe * fw_r * c_e + Upw * fw_l * c_w)

! 粘性項の加算 12 flops
EX = EX + (fu_r * c_e - fu_l * c_w)
EY = EY + (fv_r * c_e - fv_l * c_w)
EZ = EZ + (fw_r * c_e - fw_l * c_w)


! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
Un2 = u_ref2 - v(i  ,j+1,k  , 1)
Vn2 = v_ref2 - v(i  ,j+1,k  , 2)
Wn2 = w_ref2 - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
Un1 = uq
Vn1 = vq
Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
Us1 = uq
Vs1 = vq
Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
Us2 = u_ref2 - v(i  ,j-1,k  , 1)
Vs2 = v_ref2 - v(i  ,j-1,k  , 2)
Ws2 = w_ref2 - v(i  ,j-1,k  , 3)
endif

! flux
sw1 = b_n1 * lmt_n * ss
sw2 = b_s1 * lmt_s * ss

ufr2 = Un1-Up0
ufl2 = Up0-Us1
ufr4 = c1 * ufr2 - c2 * (Un2-Us1)
ufl4 = c1 * ufl2 - c2 * (Un1-Us2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * dh1
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * dh1
DUDY = 0.5d0 * ( fu_r * c_n + fu_l * c_s )


vfr2 = Vn1-Vp0
vfl2 = Vp0-Vs1
vfr4 = c1 * vfr2 - c2 * (Vn2-Vs1)
vfl4 = c1 * vfl2 - c2 * (Vn1-Vs2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * dh1
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * dh1
DVDY = 0.5d0 * ( fv_r * c_n + fv_l * c_s )


wfr2 = Wn1-Wp0
wfl2 = Wp0-Ws1
wfr4 = c1 * wfr2 - c2 * (Wn2-Ws1)
wfl4 = c1 * wfl2 - c2 * (Wn1-Ws2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * dh1
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * dh1
DWDY = 0.5d0 * ( fw_r * c_n + fw_l * c_s )


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Vpn * fu_r * c_n + Vps * fu_l * c_s)
cnv_v = cnv_v + 0.5*(Vpn * fv_r * c_n + Vps * fv_l * c_s)
cnv_w = cnv_w + 0.5*(Vpn * fw_r * c_n + Vps * fw_l * c_s)

! 粘性項の加算
EX = EX + (fu_r * c_n - fu_l * c_s)
EY = EY + (fv_r * c_n - fv_l * c_s)
EZ = EZ + (fw_r * c_n - fw_l * c_s)


! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
Ut2 = u_ref2 - v(i  ,j  ,k+1, 1)
Vt2 = v_ref2 - v(i  ,j  ,k+1, 2)
Wt2 = w_ref2 - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
Ut1 = uq
Vt1 = vq
Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
Ub1 = uq
Vb1 = vq
Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
Ub2 = u_ref2 - v(i  ,j  ,k-1, 1)
Vb2 = v_ref2 - v(i  ,j  ,k-1, 2)
Wb2 = w_ref2 - v(i  ,j  ,k-1, 3)
end if

! flux
sw1 = b_t1 * lmt_t * ss
sw2 = b_b1 * lmt_b * ss

ufr2 = Ut1-Up0
ufl2 = Up0-Ub1
ufr4 = c1 * ufr2 - c2 * (Ut2-Ub1)
ufl4 = c1 * ufl2 - c2 * (Ut1-Ub2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * dh1
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * dh1
DUDZ = 0.5d0 * ( fu_r * c_t + fu_l * c_b )


vfr2 = Vt1-Vp0
vfl2 = Vp0-Vb1
vfr4 = c1 * vfr2 - c2 * (Vt2-Vb1)
vfl4 = c1 * vfl2 - c2 * (Vt1-Vb2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * dh1
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * dh1
DVDZ = 0.5d0 * ( fv_r * c_t + fv_l * c_b )


wfr2 = Wt1-Wp0
wfl2 = Wp0-Wb1
wfr4 = c1 * wfr2 - c2 * (Wt2-Wb1)
wfl4 = c1 * wfl2 - c2 * (Wt1-Wb2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * dh1
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * dh1
DWDZ = 0.5d0 * ( fw_r * c_t + fw_l * c_b )


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Wpt * fu_r * c_t + Wpb * fu_l * c_b)
cnv_v = cnv_v + 0.5*(Wpt * fv_r * c_t + Wpb * fv_l * c_b)
cnv_w = cnv_w + 0.5*(Wpt * fw_r * c_t + Wpb * fw_l * c_b)



!----------- Sij (strain rate tensor, 歪み速度テンソル)
S11 = DUDX
S12 = (DVDX + DUDY) * 0.5d0
S13 = (DWDX + DUDZ) * 0.5d0

S21 = S12
S22 = DVDY
S23 = (DWDY + DVDZ) * 0.5d0

S31 = S13
S32 = S23
S33 = DWDZ

SSS = DSQRT( 2.0d0 * (S11*S11 + S22*S22 + S33*S33) &
           + 4.0d0 * (S12*S12 + S23*S23 + S31*S31) )


!----------- Wij (vorticity tensor, 渦度テンソル)
W11 = 0.0d0
W12 = (DVDX - DUDY) * 0.5d0
W13 = (DWDX - DUDZ) * 0.5d0

W21 = (DUDY - DVDX) * 0.5d0
W22 = 0.0d0
W23 = (DWDY - DVDZ) * 0.5d0

W31 = (DUDZ - DWDX) * 0.5d0
W32 = (DVDZ - DWDY) * 0.5d0
W33 = 0.0d0

WWW = DSQRT( 2.0d0 * (W12*W12 + W13*W13 + W21*W21 &
                    + W23*W23 + W31*W31 + W32*W32 ) )


!----------- Smagorinsky model
if (imodel == 1) then
yc     = (j - 0.5d0)*dh
min_h  = min(yc, dh*jx - yc)
DUDY_w = v(i, jx, k, 1) * (dh1*2.0d0)
tauw   = (nu * rho) * abs(DUDY_w)
utau   = sqrt(tauw/rho)
yp     = min_h * utau / nu
up     = Up0 / utau
fs     = 1.0d0 - exp(-yp/26.0d0)

nut = (Cs * fs * dh) * (Cs * fs * dh) * SSS
end if


!----------- CSM model
if (imodel == 2) then
Q_csm = (WWW * WWW - SSS * SSS) * 0.25d0
E_csm = (WWW * WWW + SSS * SSS) * 0.25d0
Fcs   = Q_csm / (E_csm + EPS)

nut = (1.0d0/22.0d0) * abs(Fcs) * sqrt(abs(Fcs)) * (1.0d0 - Fcs) * (dh*dh) * SSS

end if


!----------- WALE model
if (imodel == 3) then

S11d = 0.5d0*(DUDX*DUDX + DUDX*DUDX) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S12d = 0.5d0*(DUDY*DUDY + DVDX*DVDX)
S13d = 0.5d0*(DUDZ*DUDZ + DWDX*DWDX)
S21d = S12d
S22d = 0.5d0*(DVDY*DVDY + DVDY*DVDY) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S23d = 0.5d0*(DVDZ*DVDZ + DWDY*DWDY)
S31d = S13d
S32d = S23d
S33d = 0.5d0*(DWDZ*DWDZ + DWDZ*DWDZ) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)

Sijd2 =  S11d*S11d + S21d*S21d + S31d*S31d &
       + S12d*S12d + S22d*S22d + S32d*S32d &
       + S13d*S13d + S23d*S23d + S33d*S33d

nut = (Cw * dh)* (Cw * dh) * (Sijd2)**(3.0d0/2.0d0)  &
    / ( (SSS * SSS)**(5.0d0/2.0d0) + (Sijd2)**(5.0d0/4.0d0) )

end if


! 粘性項の加算
EX = EX + (fu_r * c_t - fu_l * c_b)
EY = EY + (fv_r * c_t - fv_l * c_b)
EZ = EZ + (fw_r * c_t - fw_l * c_b)


! 対流項と粘性項の和 > 12 flops
wv(i, j, k, 1) = -cnv_u + EX * ( 1.0d0 + nut/nu ) * vcs * dh2
wv(i, j, k, 2) = -cnv_v + EY * ( 1.0d0 + nut/nu ) * vcs * dh2
wv(i, j, k, 3) = -cnv_w + EZ * ( 1.0d0 + nut/nu ) * vcs * dh2

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pvec_central_les

!> ********************************************************************
!! @brief 4次精度の打ち切り誤差を持つポアソン反復の1ステップ目のソース項
!! @param [out]    rhs   ソース項
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     dt    時間積分幅
!! @param [in]     dh    格子幅
!! @param [in]     u_sum \sum{u}
!! @param [in,out] flop  浮動小数点演算数
!<
subroutine src_trnc (rhs, sz, g, dt, dh, u_sum, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dh, dt, dth
real                                                      ::  bw2, bw1, be2, be1, bs2, bs1, bn2, bn1, bb2, bb1, bt2, bt1, b0
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  rhs, u_sum

ix = sz(1)
jx = sz(2)
kx = sz(3)

dth = 1.0 / (dt*dh)

flop  = flop + dble(ix)*dble(jx)*dble(kx)*24.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(bw2, bw1, be2, be1, bs2, bs1, bn2, bn1, bb2, bb1, bt2, bt1, b0) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dth)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

b0  = u_sum(i  , j  , k  )

bw2 = u_sum(i-2, j  , k  )
bw1 = u_sum(i-1, j  , k  )
be1 = u_sum(i+1, j  , k  )
be2 = u_sum(i+2, j  , k  )

bs2 = u_sum(i  , j-2, k  )
bs1 = u_sum(i  , j-1, k  )
bn1 = u_sum(i  , j+1, k  )
bn2 = u_sum(i  , j+2, k  )

bb2 = u_sum(i  , j  , k-2)
bb1 = u_sum(i  , j  , k-1)
bt1 = u_sum(i  , j  , k+1)
bt2 = u_sum(i  , j  , k+2)

rhs(i,j,k) = ( bw2 - 4.0*bw1 + 6.0*b0 - 4.0*be1 + be2 &
+   bs2 - 4.0*bs1 + 6.0*b0 - 4.0*bn1 + bn2 &
+   bb2 - 4.0*bb1 + 6.0*b0 - 4.0*bt1 + bt2 ) * dth

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine src_trnc


!> ********************************************************************
!! @brief 4次精度の打ち切り誤差を持つポアソン反復の1ステップ目のソース項
!! @param [out]    rhs  ソース項
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dt   時間積分幅
!! @param [in]     dh   格子幅
!! @param [in]     u_sum \sum{u}
!! @param [in]     b    h^2/dt \nabla{u}
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine src_1st (rhs, sz, g, dt, dh, u_sum, b, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dh, dt, dth
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  rhs, u_sum, b

ix = sz(1)
jx = sz(2)
kx = sz(3)

dth = 0.25 * dh / dt

flop  = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0 + 9.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dth)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

rhs(i,j,k) = b(i,j,k) - &
( u_sum(i+1,j  ,k  ) + u_sum(i-1,j  ,k  ) &
+ u_sum(i  ,j+1,k  ) + u_sum(i  ,j-1,k  ) &
+ u_sum(i  ,j  ,k+1) + u_sum(i  ,j  ,k-1) &
- 6.0 * u_sum(i,j,k) ) * dth

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine src_1st


!> ********************************************************************
!! @brief 4次精度の打ち切り誤差を持つポアソン反復の2ステップ目のソース項
!! @param [out]    rhs  ソース項
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dh   格子幅
!! @param [in]     b    h^2/dt \nabla{u}
!! @param [in]     psi  1ステップ目の反復結果
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine src_2nd (rhs, sz, g, dh, b, psi, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dh, h2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  rhs, b, psi

ix = sz(1)
jx = sz(2)
kx = sz(3)

h2 = 0.25*dh*dh

flop  = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, h2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

rhs(i,j,k) = b(i,j,k) - h2 * psi(i,j,k)

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine src_2nd

