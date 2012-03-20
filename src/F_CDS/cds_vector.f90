!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************
!
!> @file cds_vector.f90
!> @brief vector routines for CDS
!> @author keno, FSI Team, VCAD, RIKEN

!  ********************************************************************************************
!> @subroutine cds_pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, bv, bp, v_mode, cut, flop)
!! @brief 疑似ベクトルの計算
!! @param wv 仮のベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param c_scheme 対流項スキームのモード（1-UWD, 2-center 3-MUSCL）
!! @param v00 参照速度
!! @param rei レイノルズ数の逆数
!! @param v 速度ベクトル（n-step, collocated）
!! @param bv BCindex V
!! @param bp BCindex P
!! @param v_mode 粘性項のモード (0=対流項のみ, 1=対流項と粘性項，2=粘性項は壁法則)
!! @param cut カット情報(float)
!! @param[out] flop
!<
    subroutine cds_pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, bv, bp, v_mode, cut, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, c_scheme, bpx, bvx, v_mode
    integer, dimension(3)                                       ::  sz
    real                                                        ::  dh, dh1, dh2, ck, cnv_u, cnv_v, cnv_w, b, flop
    real                                                        ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, rei, ss, vcs
    real                                                        ::  UPe, UPw, VPn, VPs, WPt, WPb
    real                                                        ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
    real                                                        ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
    real                                                        ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
    real                                                        ::  Ue1_t, Ue2_t, Uw1_t, Uw2_t, Us1_t, Us2_t, Un1_t, Un2_t, Ub1_t, Ub2_t, Ut1_t, Ut2_t
    real                                                        ::  Ve1_t, Ve2_t, Vw1_t, Vw2_t, Vs1_t, Vs2_t, Vn1_t, Vn2_t, Vb1_t, Vb2_t, Vt1_t, Vt2_t
    real                                                        ::  We1_t, We2_t, Ww1_t, Ww2_t, Ws1_t, Ws2_t, Wn1_t, Wn2_t, Wb1_t, Wb2_t, Wt1_t, Wt2_t
    real                                                        ::  cp_e, cp_w, cp_n, cp_s, cp_t, cp_b, ce_e, cw_w, cn_n, cs_s, ct_t, cb_b
    real                                                        ::  hw, he, hs, hn, hb, ht, hww, hee, hss, hnn, hbb, htt
    real                                                        ::  d1, d2, d3, d4, s1, s2, s3, s4, g1, g2, g3, g4, g5, g6
    real                                                        ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  tmp1, tmp2, tmp3, tmp4, r_tmp1, r_tmp2, r_tmp3, r_tmp4, cm1, cm2, ss_4
    real                                                        ::  u_L, u_R, v_L, v_R, w_L, w_R, cr, cl, acr, acl
    real                                                        ::  uu_e, uu_w, uu_s, uu_n, uu_b, uu_t
    real                                                        ::  vv_e, vv_w, vv_s, vv_n, vv_b, vv_t
    real                                                        ::  ww_e, ww_w, ww_s, ww_n, ww_b, ww_t
		real                                                        ::  EX, EY, EZ, delta_x, delta_y, delta_z
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv
    real*4, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    real(4), dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  cvf
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv, bp
    real, dimension(0:3)                                        ::  v00
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
		dh2= rei*dh1*dh1
    
    ! 粘性項のマスク
    if ( v_mode == 0 ) then ! 粘性項は計算しない
      vcs = 0.0
    else
      vcs = 1.0
    endif
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_ref2 = 1.5*u_ref
    v_ref2 = 1.5*v_ref
    w_ref2 = 1.5*w_ref
    
		ck = 0.0
    b  = 0.0
    ss = 1.0
    
    ! スキーム精度の切り替え
    if ( c_scheme == 1 ) then      !     1st order upwind
      ss = 0.0
    else if ( c_scheme == 2 ) then !     2nd order central 
      ck = 1.0     
    else if ( c_scheme == 3 ) then !     3rd order MUSCL
      ck = 1.0/3.0
      b  = (3.0-ck)/(1.0-ck)
    endif
    
    ss_4 = 0.25*ss
    
    cm1 = 1.0 - ck
    cm2 = 1.0 + ck
    
    flop = flop + real(ix)*real(jx)*real(kx)*1203.0 + 34.0
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, vcs, b, ck, ss_4, ss, cm1, cm2) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bpx) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(cp_e, cp_w, cp_n, cp_s, cp_t, cp_b, ce_e, cw_w, cn_n, cs_s, ct_t, cb_b) &
!$OMP PRIVATE(hw, he, hs, hn, hb, ht, hww, hee, hss, hnn, hbb, htt) &
!$OMP PRIVATE(tmp1, tmp2, tmp3, tmp4, r_tmp1, r_tmp2, r_tmp3, r_tmp4) &
!$OMP PRIVATE(u_L, u_R, v_L, v_R, w_L, w_R, cr, cl, acr, acl) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(Ue1_t, Ue2_t, Uw1_t, Uw2_t, Us1_t, Us2_t, Un1_t, Un2_t, Ub1_t, Ub2_t, Ut1_t, Ut2_t) &
!$OMP PRIVATE(Ve1_t, Ve2_t, Vw1_t, Vw2_t, Vs1_t, Vs2_t, Vn1_t, Vn2_t, Vb1_t, Vb2_t, Vt1_t, Vt2_t) &
!$OMP PRIVATE(We1_t, We2_t, Ww1_t, Ww2_t, Ws1_t, Ws2_t, Wn1_t, Wn2_t, Wb1_t, Wb2_t, Wt1_t, Wt2_t) &
!$OMP PRIVATE(d1, d2, d3, d4, s1, s2, s3, s4, g1, g2, g3, g4, g5, g6) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(uu_e, uu_w, uu_s, uu_n, uu_b, uu_t) &
!$OMP PRIVATE(vv_e, vv_w, vv_s, vv_n, vv_b, vv_t) &
!$OMP PRIVATE(ww_e, ww_w, ww_s, ww_n, ww_b, ww_t)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      cnv_u = 0.0
      cnv_v = 0.0
      cnv_w = 0.0
      
      ! 変数のロード
      include '../F_CBC/load_var_stencil5.h'
            
      bvx = bv(i,j,k)
      bpx = bp(i,j,k)
      
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
      
      ! カット情報
			cw_w = cut(1,i-1,j  ,k  ) ! d_{i-1}^-
      cp_w = cut(1,i  ,j  ,k  ) ! d_{i}^-
      cp_e = cut(2,i  ,j  ,k  ) ! d_{i}^+
      ce_e = cut(2,i+1,j  ,k  ) ! d_{i+1}^+
			
			cs_s = cut(3,i  ,j-1,k  ) ! d_{j-1}^-
      cp_s = cut(3,i  ,j  ,k  ) ! d_{j}^-
      cp_n = cut(4,i  ,j  ,k  ) ! d_{j}^+
      cn_n = cut(4,i  ,j+1,k  ) ! d_{j+1}^+

			cb_b = cut(5,i  ,j  ,k-1) ! d_{k-1}^-
      cp_b = cut(5,i  ,j  ,k  ) ! d_{k}^-
      cp_t = cut(6,i  ,j  ,k  ) ! d_{k}^+
      ct_t = cut(6,i  ,j  ,k+1) ! d_{k+1}^+
      
      ! 外挿時の分母 >  12*9 flop
      hww= 1.0/(0.5 + cw_w) ! for (i-2)
      hw = 1.0/(0.5 + cp_w) ! for (i-1)
      he = 1.0/(0.5 + cp_e) ! for (i+1)
      hee= 1.0/(0.5 + ce_e) ! for (i+2)
      
      hss= 1.0/(0.5 + cs_s) ! for (j-2)
      hs = 1.0/(0.5 + cp_s) ! for (j-1)
      hn = 1.0/(0.5 + cp_n) ! for (j+1)
      hnn= 1.0/(0.5 + cn_n) ! for (j+2)
      
      hbb= 1.0/(0.5 + cb_b) ! for (k-2)
      hb = 1.0/(0.5 + cp_b) ! for (k-1)
      ht = 1.0/(0.5 + cp_t) ! for (k+1)
      htt= 1.0/(0.5 + ct_t) ! for (k+2)
			
      ! X方向  > 4 + 12 + 84 + 232 = 332 flop ---------------------------------------
      
      ! flag >  4 flop
      tmp1 = 1.0
      tmp2 = 1.0
      tmp3 = 1.0
      tmp4 = 1.0
      if ( cw_w < 1.0 ) tmp1=0.0
      if ( cp_w < 1.0 ) tmp2=0.0
      if ( cp_e < 1.0 ) tmp3=0.0
      if ( ce_e < 1.0 ) tmp4=0.0
      r_tmp1 = 1.0-tmp1
      r_tmp2 = 1.0-tmp2
      r_tmp3 = 1.0-tmp3
      r_tmp4 = 1.0-tmp4
      
      ! 流体部分のセルフェイス位置の値を内挿 >  12 flop
      u_L = 0.5*(Up0+Uw1)
      v_L = 0.5*(Vp0+Vw1)
      w_L = 0.5*(Wp0+Ww1)
      u_R = 0.5*(Up0+Ue1)
      v_R = 0.5*(Vp0+Ve1)
      w_R = 0.5*(Wp0+We1)
      
      UPw = u_L
      UPe = u_R
      
      ! カットがある場合の参照セルの値の計算  >  21*4=84 flop
      ! d_{i-1}^-
      Uw2_t = (u_ref2 + (cw_w-1.0)*u_L )*hww
      Vw2_t = (v_ref2 + (cw_w-1.0)*v_L )*hww
      Ww2_t = (w_ref2 + (cw_w-1.0)*w_L )*hww
      Uw2 = tmp1*Uw2 + r_tmp1*Uw2_t
      Vw2 = tmp1*Vw2 + r_tmp1*Vw2_t
      Ww2 = tmp1*Ww2 + r_tmp1*Ww2_t
      
      ! d_{i}^-
      Uw1_t = (u_ref2 + (cp_w-1.0)*u_R )*hw
      Vw1_t = (v_ref2 + (cp_w-1.0)*v_R )*hw
      Ww1_t = (w_ref2 + (cp_w-1.0)*w_R )*hw
      Uw1 = tmp2*Uw1 + r_tmp2*Uw1_t
      Vw1 = tmp2*Vw1 + r_tmp2*Vw1_t
      Ww1 = tmp2*Ww1 + r_tmp2*Ww1_t
      
      ! d_{i}^+
      Ue1_t = (u_ref2 + (cp_e-1.0)*u_L )*he
      Ve1_t = (v_ref2 + (cp_e-1.0)*v_L )*he
      We1_t = (w_ref2 + (cp_e-1.0)*w_L )*he
      Ue1 = tmp3*Ue1 + r_tmp3*Ue1_t
      Ve1 = tmp3*Ve1 + r_tmp3*Ve1_t
      We1 = tmp3*We1 + r_tmp3*We1_t
      
      ! d_{i+1}^+
      Ue2_t = (u_ref2 + (ce_e-1.0)*u_R )*hee
      Ve2_t = (v_ref2 + (ce_e-1.0)*v_R )*hee
      We2_t = (w_ref2 + (ce_e-1.0)*w_R )*hee
      Ue2 = tmp4*Ue2 + r_tmp4*Ue2_t
      Ve2 = tmp4*Ve2 + r_tmp4*Ve2_t
      We2 = tmp4*We2 + r_tmp4*We2_t
      
      ! 界面の左右の状態を計算  > 4+ (4+36+20+16)*3 = 232 flop
      cr  = UPe - u_ref
      cl  = UPw - u_ref
      acr = abs(cr)
      acl = abs(cl)
      
      d4 = Ue2-Ue1
      d3 = Ue1-Up0
      d2 = Up0-Uw1
      d1 = Uw1-Uw2
      
      include '../F_CBC/muscl.h'
      
      Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4
      cnv_u = 0.5*(cr*(Urr+Url) - acr*(Urr-Url))*c_e &
            - 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull))*c_w + cnv_u
      
      d4 = Ve2-Ve1
      d3 = Ve1-Vp0
      d2 = Vp0-Vw1
      d1 = Vw1-Vw2
      
      include '../F_CBC/muscl.h'
      
      Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4
      cnv_v = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl))*c_e &
            - 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll))*c_w + cnv_v
            
      d4 = We2-We1
      d3 = We1-Wp0
      d2 = Wp0-Ww1
      d1 = Ww1-Ww2
      
      include '../F_CBC/muscl.h'
      
      Wrr = We1 - (cm1*g6+cm2*g5)*ss_4
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4
      cnv_w = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl))*c_e &
            - 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll))*c_w + cnv_w
			
      ! Y方向 ---------------------------------------
      
      tmp1 = 1.0
      tmp2 = 1.0
      tmp3 = 1.0
      tmp4 = 1.0
      if ( cs_s < 1.0 ) tmp1=0.0
      if ( cp_s < 1.0 ) tmp2=0.0
      if ( cp_n < 1.0 ) tmp3=0.0
      if ( cn_n < 1.0 ) tmp4=0.0
      r_tmp1 = 1.0-tmp1
      r_tmp2 = 1.0-tmp2
      r_tmp3 = 1.0-tmp3
      r_tmp4 = 1.0-tmp4
      
      ! セルフェイス位置の値を内挿
      u_L = 0.5*(Up0+Us1)
      v_L = 0.5*(Vp0+Vs1)
      w_L = 0.5*(Wp0+Ws1)
      u_R = 0.5*(Up0+Un1)
      v_R = 0.5*(Vp0+Vn1)
      w_R = 0.5*(Wp0+Wn1)
      
      VPs = v_L
      VPn = v_R
      
      ! カットがある場合の参照セルの値の計算
      ! d_{j-1}^-
      Us2_t = (u_ref2 + (cs_s-1.0)*u_L )*hss
      Vs2_t = (v_ref2 + (cs_s-1.0)*v_L )*hss
      Ws2_t = (w_ref2 + (cs_s-1.0)*w_L )*hss
      Us2 = tmp1*Us2 + r_tmp1*Us2_t
      Vs2 = tmp1*Vs2 + r_tmp1*Vs2_t
      Ws2 = tmp1*Ws2 + r_tmp1*Ws2_t
      
      ! d_{j}^-
      Us1_t = (u_ref2 + (cp_s-1.0)*u_R )*hs
      Vs1_t = (v_ref2 + (cp_s-1.0)*v_R )*hs
      Ws1_t = (w_ref2 + (cp_s-1.0)*w_R )*hs
      Us1 = tmp2*Us1 + r_tmp2*Us1_t
      Vs1 = tmp2*Vs1 + r_tmp2*Vs1_t
      Ws1 = tmp2*Ws1 + r_tmp2*Ws1_t
      
      ! d_{j}^+
      Un1_t = (u_ref2 + (cp_n-1.0)*u_L )*hn
      Vn1_t = (v_ref2 + (cp_n-1.0)*v_L )*hn
      Wn1_t = (w_ref2 + (cp_n-1.0)*w_L )*hn
      Un1 = tmp3*Un1 + r_tmp3*Un1_t
      Vn1 = tmp3*Vn1 + r_tmp3*Vn1_t
      Wn1 = tmp3*Wn1 + r_tmp3*Wn1_t
      
      ! d_{j+1}^+
      Un2_t = (u_ref2 + (cn_n-1.0)*u_R )*hnn
      Vn2_t = (v_ref2 + (cn_n-1.0)*v_R )*hnn
      Wn2_t = (w_ref2 + (cn_n-1.0)*w_R )*hnn
      Un2 = tmp4*Un2 + r_tmp4*Un2_t
      Vn2 = tmp4*Vn2 + r_tmp4*Vn2_t
      Wn2 = tmp4*Wn2 + r_tmp4*Wn2_t
      
      ! 界面の左右の状態を計算
      cr  = VPn - v_ref
      cl  = VPs - v_ref
      acr = abs(cr)
      acl = abs(cl)
      
      d4 = Un2-Un1
      d3 = Un1-Up0
      d2 = Up0-Us1
      d1 = Us1-Us2
      
      include '../F_CBC/muscl.h'
      
      Urr = Un1 - (cm1*g6+cm2*g5)*ss_4
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Us1 + (cm1*g1+cm2*g2)*ss_4
      cnv_u = 0.5*(cr*(Urr+Url) - acr*(Urr-Url))*c_n &
            - 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull))*c_s + cnv_u
      
      d4 = Vn2-Vn1
      d3 = Vn1-Vp0
      d2 = Vp0-Vs1
      d1 = Vs1-Vs2
      
      include '../F_CBC/muscl.h'
      
      Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4
      cnv_v = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl))*c_n &
            - 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll))*c_s + cnv_v

      d4 = Wn2-Wn1
      d3 = Wn1-Wp0
      d2 = Wp0-Ws1
      d1 = Ws1-Ws2
      
      include '../F_CBC/muscl.h'
      
      Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4
      cnv_w = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl))*c_n &
            - 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll))*c_s + cnv_w
			
      ! Z方向 ---------------------------------------
      
      tmp1 = 1.0
      tmp2 = 1.0
      tmp3 = 1.0
      tmp4 = 1.0
      if ( cb_b < 1.0 ) tmp1=0.0
      if ( cp_b < 1.0 ) tmp2=0.0
      if ( cp_t < 1.0 ) tmp3=0.0
      if ( ct_t < 1.0 ) tmp4=0.0
      r_tmp1 = 1.0-tmp1
      r_tmp2 = 1.0-tmp2
      r_tmp3 = 1.0-tmp3
      r_tmp4 = 1.0-tmp4
      
      ! セルフェイス位置の値を内挿
      u_L = 0.5*(Up0+Ub1)
      v_L = 0.5*(Vp0+Vb1)
      w_L = 0.5*(Wp0+Wb1)
      u_R = 0.5*(Up0+Ut1)
      v_R = 0.5*(Vp0+Vt1)
      w_R = 0.5*(Wp0+Wt1)
      
      WPb = w_L
      WPt = w_R
      
      ! カットがある場合の参照セルの値の計算
      ! d_{k-1}^-
      Ub2_t = (u_ref2 + (cb_b-1.0)*u_L )*hbb
      Vb2_t = (v_ref2 + (cb_b-1.0)*v_L )*hbb
      Wb2_t = (w_ref2 + (cb_b-1.0)*w_L )*hbb
      Ub2 = tmp1*Ub2 + r_tmp1*Ub2_t
      Vb2 = tmp1*Vb2 + r_tmp1*Vb2_t
      Wb2 = tmp1*Wb2 + r_tmp1*Wb2_t
      
      ! d_{k}^-
      Ub1_t = (u_ref2 + (cp_b-1.0)*u_R )*hb
      Vb1_t = (v_ref2 + (cp_b-1.0)*v_R )*hb
      Wb1_t = (w_ref2 + (cp_b-1.0)*w_R )*hb
      Ub1 = tmp2*Ub1 + r_tmp2*Ub1_t
      Vb1 = tmp2*Vb1 + r_tmp2*Vb1_t
      Wb1 = tmp2*Wb1 + r_tmp2*Wb1_t
      
      ! d_{k}^+
      Ut1_t = (u_ref2 + (cp_t-1.0)*u_L )*ht
      Vt1_t = (v_ref2 + (cp_t-1.0)*v_L )*ht
      Wt1_t = (w_ref2 + (cp_t-1.0)*w_L )*ht
      Ut1 = tmp3*Ut1 + r_tmp3*Ut1_t
      Vt1 = tmp3*Vt1 + r_tmp3*Vt1_t
      Wt1 = tmp3*Wt1 + r_tmp3*Wt1_t
      
      ! d_{k+1}^+
      Ut2_t = (u_ref2 + (ct_t-1.0)*u_R )*htt
      Vt2_t = (v_ref2 + (ct_t-1.0)*v_R )*htt
      Wt2_t = (w_ref2 + (ct_t-1.0)*w_R )*htt
      Ut2 = tmp4*Ut2 + r_tmp4*Ut2_t
      Vt2 = tmp4*Vt2 + r_tmp4*Vt2_t
      Wt2 = tmp4*Wt2 + r_tmp4*Wt2_t
      
      ! 界面の左右の状態を計算
      cr  = WPt - w_ref
      cl  = WPb - w_ref
      acr = abs(cr)
      acl = abs(cl)

      d4 = Ut2-Ut1
      d3 = Ut1-Up0
      d2 = Up0-Ub1
      d1 = Ub1-Ub2
      
      include '../F_CBC/muscl.h'
      
      Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4
      cnv_u = 0.5*(cr*(Urr+Url) - acr*(Urr-Url))*c_t &
            - 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull))*c_b + cnv_u

      d4 = Vt2-Vt1
      d3 = Vt1-Vp0
      d2 = Vp0-Vb1
      d1 = Vb1-Vb2
      
      include '../F_CBC/muscl.h'
      
      Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4
      cnv_v = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl))*c_t &
            - 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll))*c_b + cnv_v

      d4 = Wt2-Wt1
      d3 = Wt1-Wp0
      d2 = Wp0-Wb1
      d1 = Wb1-Wb2
      
      include '../F_CBC/muscl.h'
      
      Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4
      cnv_w = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl))*c_t &
            - 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll))*c_b + cnv_w

      
      ! 粘性項の計算  >  18 + 27 + 15*3 + 9 = 99 flop -----------------------------------------------
      uu_e = ( Ue1 - Up0 )
      uu_w = ( Up0 - Uw1 )
      uu_n = ( Un1 - Up0 )
      uu_s = ( Up0 - Us1 )
      uu_t = ( Ut1 - Up0 )
      uu_b = ( Up0 - Ub1 )
      
      vv_e = ( Ve1 - Vp0 )
      vv_w = ( Vp0 - Vw1 )
      vv_n = ( Vn1 - Vp0 )
      vv_s = ( Vp0 - Vs1 )
      vv_t = ( Vt1 - Vp0 )
      vv_b = ( Vp0 - Vb1 )
      
      ww_e = ( We1 - Wp0 )
      ww_w = ( Wp0 - Ww1 )
      ww_n = ( Wn1 - Wp0 )
      ww_s = ( Wp0 - Ws1 )
      ww_t = ( Wt1 - Wp0 )
      ww_b = ( Wp0 - Wb1 )
      
      delta_x = 2.0/(cp_e + cp_w)
      delta_y = 2.0/(cp_n + cp_s)
      delta_z = 2.0/(cp_t + cp_b)
      
      Ex = ( (uu_e * c_e - uu_w * c_w) * delta_x &
           + (uu_n * c_n - uu_s * c_s) * delta_y &
           + (uu_t * c_t - uu_b * c_b) * delta_z &
           ) * dh2
           
      Ey = ( (vv_e * c_e - vv_w * c_w) * delta_x &
           + (vv_n * c_n - vv_s * c_s) * delta_y &
           + (vv_t * c_t - vv_b * c_b) * delta_z &
           ) * dh2
           
      Ez = ( (ww_e * c_e - ww_w * c_w) * delta_x &
           + (ww_n * c_n - ww_s * c_s) * delta_y &
           + (ww_t * c_t - ww_b * c_b) * delta_z &
           ) * dh2
			
      ! 対流項と粘性項の和
      wv(1,i,j,k) = -cnv_u*dh1 + EX*vcs
      wv(2,i,j,k) = -cnv_v*dh1 + EY*vcs
      wv(3,i,j,k) = -cnv_w*dh1 + EZ*vcs
      
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

		return
    end subroutine cds_pvec_muscl

!  ***************************************************************************************
!> @subroutine cds_update_vec (v, div, sz, g, dt, dh, vc, p, bp, bv, cut, v00, coef, flop)
!! @brief 次のステップのセルセンターの速度を更新し，発散値を計算する
!! @param[out] v n+1時刻のセルセンターの速度ベクトル
!! @param div 発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param dh 格子幅
!! @param vc セルセンタの疑似速度ベクトル
!! @param p 圧力
!! @param bp BCindex P
!! @param bv BCindex V
!! @param cut カット情報(float)
!! @param v00 参照速度
!! @param coef 係数
!! @param[out] flop flop count
!<
    subroutine cds_update_vec (v, div, sz, g, dt, dh, vc, p, bp, bv, cut, v00, coef, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bpx, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  dh, dt, dd, flop, coef, actv, r_actv
    real                                                        ::  pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue_f, Uw_f, Vn_f, Vs_f, Wt_f, Wb_f
    real                                                        ::  Ue_c, Uw_c, Vn_c, Vs_c, Wt_c, Wb_c
    real                                                        ::  d1, d2, d3, d4, d5, d6
    real                                                        ::  q1, q2, q3, q4, q5, q6
    real                                                        ::  c1, c2, c3, c4, c5, c6
    real                                                        ::  h1, h2, h3, h4, h5, h6
    real                                                        ::  N_e, N_w, N_n, N_s, N_t, N_b
    real                                                        ::  u_ref2, v_ref2, w_ref2
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  p
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vc, v
    real, dimension(0:3)                                        ::  v00
    real*4, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp, bv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dd = dt/dh
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_ref2 = 1.5*u_ref
    v_ref2 = 1.5*v_ref
    w_ref2 = 1.5*w_ref
    
    flop = flop + real(ix)*real(jx)*real(kx)*133.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bpx = bp(i,j,k)
      bvx = bv(i,j,k)
      actv = real(ibits(bvx, State, 1))
      r_actv = 1.0 - actv
      
      c1 = 1.0
      c2 = 1.0
      c3 = 1.0
      c4 = 1.0
      c5 = 1.0
      c6 = 1.0
      
      ! Neumann条件のとき，0.0
      N_w = real(ibits(bpx, bc_n_W, 1))  ! w
      N_e = real(ibits(bpx, bc_n_E, 1))  ! e
      N_s = real(ibits(bpx, bc_n_S, 1))  ! s
      N_n = real(ibits(bpx, bc_n_N, 1))  ! n
      N_b = real(ibits(bpx, bc_n_B, 1))  ! b
      N_t = real(ibits(bpx, bc_n_T, 1))  ! t
      
      Up0 = vc(i  ,j  ,k  ,1)
      Vp0 = vc(i  ,j  ,k  ,2)
      Wp0 = vc(i  ,j  ,k  ,3)
      Uw0 = vc(i-1,j  ,k  ,1)
      Ue0 = vc(i+1,j  ,k  ,1)
      Vs0 = vc(i  ,j-1,k  ,2)
      Vn0 = vc(i  ,j+1,k  ,2)
      Wb0 = vc(i  ,j  ,k-1,3)
      Wt0 = vc(i  ,j  ,k+1,3)
      
      ! 距離情報 [0, 1]
      d1 = cut(1,i,j,k) ! d_{i}^-
      d2 = cut(2,i,j,k) ! d_{i}^+
      d3 = cut(3,i,j,k) ! d_{j}^-
      d4 = cut(4,i,j,k) ! d_{j}^+
      d5 = cut(5,i,j,k) ! d_{k}^-
      d6 = cut(6,i,j,k) ! d_{k}^+
      
      ! 交点がある場合-> q=0.0
      q1 = 1.0
      q2 = 1.0
      q3 = 1.0
      q4 = 1.0
      q5 = 1.0
      q6 = 1.0
      if (d1 < 1.0) q1 = 0.0 ! q1 = aint(d1) 
      if (d2 < 1.0) q2 = 0.0 ! q2 = aint(d2)
      if (d3 < 1.0) q3 = 0.0 ! q3 = aint(d3)
      if (d4 < 1.0) q4 = 0.0 ! q4 = aint(d4)
      if (d5 < 1.0) q5 = 0.0 ! q5 = aint(d5)
      if (d6 < 1.0) q6 = 0.0 ! q6 = aint(d6)
      
      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c1 = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c2 = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c3 = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c4 = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c5 = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c6 = 0.0
      
      h1 = 1.0/(d1 + 0.5)
      h2 = 1.0/(d2 + 0.5)
      h3 = 1.0/(d3 + 0.5)
      h4 = 1.0/(d4 + 0.5)
      h5 = 1.0/(d5 + 0.5)
      h6 = 1.0/(d6 + 0.5)
      
      ! 圧力勾配
      pc  = p(i,  j,  k  )
      pxw = (pc - p(i-1,j  ,k  )) * N_w
      pxe = (p(i+1,j  ,k  ) - pc) * N_e
      pys = (pc - p(i  ,j-1,k  )) * N_s
      pyn = (p(i  ,j+1,k  ) - pc) * N_n
      pzb = (pc - p(i  ,j  ,k-1)) * N_b
      pzt = (p(i  ,j  ,k+1) - pc) * N_t
      px = 0.5*(pxe + pxw)
      py = 0.5*(pyn + pys)
      pz = 0.5*(pzt + pzb)
      
      ! 両側流体の場合のセルフェイス速度
      Uw_f = (Uw0 + Up0) * 0.50
      Ue_f = (Ue0 + Up0) * 0.50
      Vs_f = (Vs0 + Vp0) * 0.50
      Vn_f = (Vn0 + Vp0) * 0.50
      Wt_f = (Wt0 + Wp0) * 0.50
      Wb_f = (Wb0 + Wp0) * 0.50
      
      ! 交点がある場合の補正 stable
      Uw_c = (u_ref2 + (d1-1.0)*Ue_f )*h1
      Ue_c = (u_ref2 + (d2-1.0)*Uw_f )*h2
      Vs_c = (v_ref2 + (d3-1.0)*Vn_f )*h3
      Vn_c = (v_ref2 + (d4-1.0)*Vs_f )*h4
      Wb_c = (w_ref2 + (d5-1.0)*Wt_f )*h5
      Wt_c = (w_ref2 + (d6-1.0)*Wb_f )*h6
      
      ! カットがある場合のcompactな壁面修正
      !Uw_c = (0.5*u_ref + (d1-0.5)*Up0 )/d1 
      !Ue_c = (0.5*u_ref + (d2-0.5)*Up0 )/d2
      !Vs_c = (0.5*v_ref + (d3-0.5)*Vp0 )/d3
      !Vn_c = (0.5*v_ref + (d4-0.5)*Vp0 )/d4
      !Wb_c = (0.5*w_ref + (d5-0.5)*Wp0 )/d5
      !Wt_c = (0.5*w_ref + (d6-0.5)*Wp0 )/d6
      
      ! セルフェイス速度 q=0.0(交点), 1.0(Fluid)
      Uw = q1*Uw_f + (1.0-q1)*Uw_c - dd * pxw
      Ue = q2*Ue_f + (1.0-q2)*Ue_c - dd * pxe
      Vs = q3*Vs_f + (1.0-q3)*Vs_c - dd * pys
      Vn = q4*Vn_f + (1.0-q4)*Vn_c - dd * pyn
      Wb = q5*Wb_f + (1.0-q5)*Wb_c - dd * pzb
      Wt = q6*Wt_f + (1.0-q6)*Wt_c - dd * pzt
      
      ! 発散値 VBCの寄与は除外
      div(i,j,k) = (Ue * c2 - Uw * c1 &
                  + Vn * c4 - Vs * c3 &
                  + Wt * c6 - Wb * c5) * coef * actv
      
      ! セルセンタの速度更新
      v(i,j,k,1) = ( Up0-px*dd )*actv + r_actv*u_ref
      v(i,j,k,2) = ( Vp0-py*dd )*actv + r_actv*v_ref
      v(i,j,k,3) = ( Wp0-pz*dd )*actv + r_actv*w_ref
    end do
    end do
    end do
    
    return 
    end subroutine cds_update_vec

!  *************************************************************
!> @subroutine cds_div (div, sz, g, coef, v, bv, cut, v00, flop)
!! @brief 速度の発散を計算する
!! @param div 速度の発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param coef 係数
!! @param v 疑似ベクトル
!! @param bv BCindex V
!! @param cut カット情報(float)
!! @param v00 参照速度
!! @param[out] flop flop count
!<
    subroutine cds_div (div, sz, g, coef, v, bv, cut, v00, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                        ::  Ue1, Uw1, Vn1, Vs1, Wt1, Wb1
    real                                                        ::  c1, c2, c3, c4, c5, c6
    real                                                        ::  d1, d2, d3, d4, d5, d6
    real                                                        ::  h1, h2, h3, h4, h5, h6
    real                                                        ::  q1, q2, q3, q4, q5, q6
    real                                                        ::  u_ref, v_ref, w_ref, coef, flop, actv
    real                                                        ::  u_ref2, v_ref2, w_ref2
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
		real, dimension(0:3)                                        ::  v00
    real*4, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_ref2 = 1.5*u_ref
    v_ref2 = 1.5*v_ref
    w_ref2 = 1.5*w_ref
    
    flop = flop + real(ix*jx*kx)*73.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bvx = bv(i,j,k)
      actv= real(ibits(bvx, State, 1))
      
      Up0 = v(i  ,j  ,k  ,1)
      Vp0 = v(i  ,j  ,k  ,2)
      Wp0 = v(i  ,j  ,k  ,3)
      Uw0 = v(i-1,j  ,k  ,1)
      Ue0 = v(i+1,j  ,k  ,1)
      Vs0 = v(i  ,j-1,k  ,2)
      Vn0 = v(i  ,j+1,k  ,2)
      Wb0 = v(i  ,j  ,k-1,3)
      Wt0 = v(i  ,j  ,k+1,3)
      
      d1 = cut(1,i,j,k) ! d_{i}^-
      d2 = cut(2,i,j,k) ! d_{i}^+
      d3 = cut(3,i,j,k) ! d_{j}^-
      d4 = cut(4,i,j,k) ! d_{j}^+
      d5 = cut(5,i,j,k) ! d_{k}^-
      d6 = cut(6,i,j,k) ! d_{k}^+
      
      ! 交点がある場合，切り捨て -> q?=0.0
      q1 = aint(d1)
      q2 = aint(d2)
      q3 = aint(d3)
      q4 = aint(d4)
      q5 = aint(d5)
      q6 = aint(d6)
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c = 1.0(Normal) / 0.0(BC)
      c1 = 1.0
      c2 = 1.0
      c3 = 1.0
      c4 = 1.0
      c5 = 1.0
      c6 = 1.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c1 = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c2 = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c3 = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c4 = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c5 = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c6 = 0.0
      
      h1 = 1.0/(d1 + 0.5)
      h2 = 1.0/(d2 + 0.5)
      h3 = 1.0/(d3 + 0.5)
      h4 = 1.0/(d4 + 0.5)
      h5 = 1.0/(d5 + 0.5)
      h6 = 1.0/(d6 + 0.5)
      
      !Uw1_t = (u_ref + (d1-1.0)*Up0 )/d1
      !Ue1_t = (u_ref + (d2-1.0)*Up0 )/d2
      !Vs1_t = (v_ref + (d3-1.0)*Vp0 )/d3
      !Vn1_t = (v_ref + (d4-1.0)*Vp0 )/d4
      !Wb1_t = (w_ref + (d5-1.0)*Wp0 )/d5
      !Wt1_t = (w_ref + (d6-1.0)*Wp0 )/d6
      
      !Uw0 = q1*Uw0 + (1.0-q1)*Uw1_t
      !Ue0 = q2*Ue0 + (1.0-q2)*Ue1_t
      !Vs0 = q3*Vs0 + (1.0-q3)*Vs1_t
      !Vn0 = q4*Vn0 + (1.0-q4)*Vn1_t
      !Wb0 = q5*Wb0 + (1.0-q5)*Wb1_t
      !Wt0 = q6*Wt0 + (1.0-q6)*Wt1_t
      
      ! 流体セルの場合
      Uw1 = (Uw0 + Up0) * 0.50
      Ue1 = (Ue0 + Up0) * 0.50
      Vs1 = (Vs0 + Vp0) * 0.50
      Vn1 = (Vn0 + Vp0) * 0.50
      Wb1 = (Wb0 + Wp0) * 0.50
      Wt1 = (Wt0 + Wp0) * 0.50
      
      ! 交点がある場合の補正 stable
      Uw = (u_ref2 + (d1-1.0)*Ue1 )*h1
      Ue = (u_ref2 + (d2-1.0)*Uw1 )*h2
      Vs = (v_ref2 + (d3-1.0)*Vn1 )*h3
      Vn = (v_ref2 + (d4-1.0)*Vs1 )*h4
      Wb = (w_ref2 + (d5-1.0)*Wt1 )*h5
      Wt = (w_ref2 + (d6-1.0)*Wb1 )*h6
      
      ! 交点がある場合の補正 compact
      !Uw = (0.5*u_ref + (d1-0.5)*Up0 )/d1
      !Ue = (0.5*u_ref + (d2-0.5)*Up0 )/d2
      !Vs = (0.5*v_ref + (d3-0.5)*Vp0 )/d3
      !Vn = (0.5*v_ref + (d4-0.5)*Vp0 )/d4
      !Wb = (0.5*w_ref + (d5-0.5)*Wp0 )/d5
      !Wt = (0.5*w_ref + (d6-0.5)*Wp0 )/d6
      
      ! 壁面の補正済みの値
      Uw = q1*Uw1 + (1.0-q1)*Uw
      Ue = q2*Ue1 + (1.0-q2)*Ue
      Vs = q3*Vs1 + (1.0-q3)*Vs
      Vn = q4*Vn1 + (1.0-q4)*Vn
      Wb = q5*Wb1 + (1.0-q5)*Wb
      Wt = q6*Wt1 + (1.0-q6)*Wt
      
      ! VBCの場合には寄与をキャンセル
      div(i,j,k) = ( Ue * c2 - Uw * c1 + Vn * c4 - Vs * c3 + Wt * c6 - Wb * c5 ) * coef * actv
    end do
    end do
    end do

    return
    end subroutine cds_div

!  *****************************************************************************************
!> @subroutine cds_eddy_viscosity (vt, sz, g, dh, re, cs, v, cut, vt_range, yp_range, v00)
!! @brief Smagorinsky Modelによる乱流渦粘性係数の計算，減衰関数を併用
!! @param vt 乱流渦粘性係数
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param re レイノルズ数
!! @param cs 須磨五輪スキー定数
!! @param v 速度ベクトル
!! @param cut カット情報
!! @param vt_range 渦粘性係数の最小値と最大値
!! @param yp_range Y+の最小値と最大値
!! @param v00 参照速度
!! @date 2009/02/17
!! @note
!!    - vtmin, vtmax > vt_range(2)
!!    - ypmin, ypmax > yp_range(2)
!! @todo
!!    - 境界条件は必要か？
!<
    subroutine cds_eddy_viscosity (vt, sz, g, dh, re, cs, v, cut, vt_range, yp_range, v00)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, m
    integer, dimension(3)                                       ::  sz
    real, dimension(2)                                          ::  vt_range, yp_range
    real                                                        ::  dh, re, cs, yp, Ut, Ut0
    real                                                        ::  dudx, dudy, dudz
    real                                                        ::  dvdx, dvdy, dvdz
    real                                                        ::  dwdx, dwdy, dwdz
    real                                                        ::  d1, d2, ddd, c1, c2, c3, delta, dd
    real                                                        ::  fs, aaa, Vmag, dis, tw, up1
    real                                                        ::  u_ref, v_ref, w_ref, u1, u2, u3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  vt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(0:3)                                        ::  v00
    real, dimension(6, 1-g:sz(1)+1, 1-g:sz(2)+1, 1-g:sz(3)+1)   ::  cut

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
      DUDX=(v(i+1,j,k,1)-v(i-1,j,k,1)) * dd
      DUDY=(v(i,j+1,k,1)-v(i,j-1,k,1)) * dd
      DUDZ=(v(i,j,k+1,1)-v(i,j,k-1,1)) * dd

      DVDX=(v(i+1,j,k,2)-v(i-1,j,k,2)) * dd
      DVDY=(v(i,j+1,k,2)-v(i,j-1,k,2)) * dd
      DVDZ=(v(i,j,k+1,2)-v(i,j,k-1,2)) * dd

      DWDX=(v(i+1,j,k,3)-v(i-1,j,k,3)) * dd
      DWDY=(v(i,j+1,k,3)-v(i,j-1,k,3)) * dd
      DWDZ=(v(i,j,k+1,3)-v(i,j,k-1,3)) * dd

      D1  = DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ
      c1  = DUDY+DVDX
      c2  = DVDZ+DWDY
      c3  = DWDX+DUDZ
      D2  = c1*c1 + c2*c2 + c3*c3
      DDD = SQRT(2.0D0*D1 + D2)

!     ------- Y+ & fs--
      fs =1.0
      AAA= cut(1,i,j,k) * cut(2,i,j,k)   &
          *cut(3,i,j,k) * cut(4,i,j,k)   &
          *cut(5,i,j,k) * cut(6,i,j,k)
      u1 = v(i,j,k,1) - u_ref
      u2 = v(i,j,k,2) - v_ref
      u3 = v(i,j,k,3) - w_ref
      Vmag = SQRT(u1*u1 + u2*u2 + u3*u3)

      IF ( (AAA < 1.0) .AND. (Vmag >= 0.001) ) THEN
        DIS = AMIN1(cut(1,i,j,k), cut(2,i,j,k),    &
                    cut(3,i,j,k), cut(4,i,j,k),    &
                    cut(5,i,j,k), cut(6,i,j,k) ) * dh
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
    end subroutine cds_eddy_viscosity
