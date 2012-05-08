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
    real                                                        ::  Uw_r, Ue_r, Us_r, Un_r, Ub_r, Ut_r
    real                                                        ::  Vw_r, Ve_r, Vs_r, Vn_r, Vb_r, Vt_r
    real                                                        ::  Ww_r, We_r, Ws_r, Wn_r, Wb_r, Wt_r
    real                                                        ::  Uw_t, Ue_t, Us_t, Un_t, Ub_t, Ut_t
    real                                                        ::  Vw_t, Ve_t, Vs_t, Vn_t, Vb_t, Vt_t
    real                                                        ::  Ww_t, We_t, Ws_t, Wn_t, Wb_t, Wt_t
    real                                                        ::  Uw_c, Ue_c, Us_c, Un_c, Ub_c, Ut_c
    real                                                        ::  Vw_c, Ve_c, Vs_c, Vn_c, Vb_c, Vt_c
    real                                                        ::  Ww_c, We_c, Ws_c, Wn_c, Wb_c, Wt_c
    real                                                        ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
    real                                                        ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
    real                                                        ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
    real                                                        ::  Ue1_t, Ue2_t, Uw1_t, Uw2_t, Us1_t, Us2_t, Un1_t, Un2_t, Ub1_t, Ub2_t, Ut1_t, Ut2_t
    real                                                        ::  Ve1_t, Ve2_t, Vw1_t, Vw2_t, Vs1_t, Vs2_t, Vn1_t, Vn2_t, Vb1_t, Vb2_t, Vt1_t, Vt2_t
    real                                                        ::  We1_t, We2_t, Ww1_t, Ww2_t, Ws1_t, Ws2_t, Wn1_t, Wn2_t, Wb1_t, Wb2_t, Wt1_t, Wt2_t
    real                                                        ::  cw, ce, cs, cn, cb, ct
    real                                                        ::  ww, we, ws, wn, wb, wt
    real                                                        ::  dw, de, ds, dn, db, dt, dww, dee, dss, dnn, dbb, dtt
    real                                                        ::  hw, he, hs, hn, hb, ht, hww, hee, hss, hnn, hbb, htt
    real                                                        ::  qw, qe, qs, qn, qb, qt, qww, qee, qss, qnn, qbb, qtt
    real                                                        ::  rw, re, rs, rn, rb, rt, rww, ree, rss, rnn, rbb, rtt
    real                                                        ::  r_qw, r_qe, r_qs, r_qn, r_qb, r_qt, r_qww, r_qee, r_qss, r_qnn, r_qbb, r_qtt
    real                                                        ::  dv1, dv2, dv3, dv4, s1, s2, s3, s4, g1, g2, g3, g4, g5, g6
    real                                                        ::  Uer, Uel, Uwr, Uwl, Ver, Vel, Vwr, Vwl, Wer, Wel, Wwr, Wwl
    real                                                        ::  Unr, Unl, Usr, Usl, Vnr, Vnl, Vsr, Vsl, Wnr, Wnl, Wsr, Wsl
    real                                                        ::  Utr, Utl, Ubr, Ubl, Vtr, Vtl, Vbr, Vbl, Wtr, Wtl, Wbr, Wbl
    real                                                        ::  ufw, ufe, vfs, vfn, wfb, wft
    real                                                        ::  afw, afe, afs, afn, afb, aft
    real                                                        ::  cm1, cm2, ss_4
    real                                                        ::  du_e, du_w, du_s, du_n, du_b, du_t
    real                                                        ::  dv_e, dv_w, dv_s, dv_n, dv_b, dv_t
    real                                                        ::  dw_e, dw_w, dw_s, dw_n, dw_b, dw_t
		real                                                        ::  EX, EY, EZ
    real                                                        ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, wv
    real*4, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
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
    
    flop = flop + real(ix)*real(jx)*real(kx)*1221.0 + 34.0
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, vcs, b, ck, ss_4, ss, cm1, cm2) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bpx) &
!$OMP PRIVATE(cw, ce, cs, cn, cb, ct) &
!$OMP PRIVATE(ww, we, ws, wn, wb, wt) &
!$OMP PRIVATE(dw, de, ds, dn, db, dt, dww, dee, dss, dnn, dbb, dtt) &
!$OMP PRIVATE(hw, he, hs, hn, hb, ht, hww, hee, hss, hnn, hbb, htt) &
!$OMP PRIVATE(qw, qe, qs, qn, qb, qt, qww, qee, qss, qnn, qbb, qtt) &
!$OMP PRIVATE(rw, re, rs, rn, rb, rt, rww, ree, rss, rnn, rbb, rtt) &
!$OMP PRIVATE(r_qw, r_qe, r_qs, r_qn, r_qb, r_qt, r_qww, r_qee, r_qss, r_qnn, r_qbb, r_qtt) &
!$OMP PRIVATE(ufw, ufe, vfs, vfn, wfb, wft) &
!$OMP PRIVATE(afw, afe, afs, afn, afb, aft) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(Uw_r, Ue_r, Us_r, Un_r, Ub_r, Ut_r) &
!$OMP PRIVATE(Vw_r, Ve_r, Vs_r, Vn_r, Vb_r, Vt_r) &
!$OMP PRIVATE(Ww_r, We_r, Ws_r, Wn_r, Wb_r, Wt_r) &
!$OMP PRIVATE(Uw_t, Ue_t, Us_t, Un_t, Ub_t, Ut_t) &
!$OMP PRIVATE(Vw_t, Ve_t, Vs_t, Vn_t, Vb_t, Vt_t) &
!$OMP PRIVATE(Ww_t, We_t, Ws_t, Wn_t, Wb_t, Wt_t) &
!$OMP PRIVATE(Uw_c, Ue_c, Us_c, Un_c, Ub_c, Ut_c) &
!$OMP PRIVATE(Vw_c, Ve_c, Vs_c, Vn_c, Vb_c, Vt_c) &
!$OMP PRIVATE(Ww_c, We_c, Ws_c, Wn_c, Wb_c, Wt_c) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(Ue1_t, Ue2_t, Uw1_t, Uw2_t, Us1_t, Us2_t, Un1_t, Un2_t, Ub1_t, Ub2_t, Ut1_t, Ut2_t) &
!$OMP PRIVATE(Ve1_t, Ve2_t, Vw1_t, Vw2_t, Vs1_t, Vs2_t, Vn1_t, Vn2_t, Vb1_t, Vb2_t, Vt1_t, Vt2_t) &
!$OMP PRIVATE(We1_t, We2_t, Ww1_t, Ww2_t, Ws1_t, Ws2_t, Wn1_t, Wn2_t, Wb1_t, Wb2_t, Wt1_t, Wt2_t) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, s1, s2, s3, s4, g1, g2, g3, g4, g5, g6) &
!$OMP PRIVATE(Uer, Uel, Uwr, Uwl, Ver, Vel, Vwr, Vwl, Wer, Wel, Wwr, Wwl) &
!$OMP PRIVATE(Unr, Unl, Usr, Usl, Vnr, Vnl, Vsr, Vsl, Wnr, Wnl, Wsr, Wsl) &
!$OMP PRIVATE(Utr, Utl, Ubr, Ubl, Vtr, Vtl, Vbr, Vbl, Wtr, Wtl, Wbr, Wbl) &
!$OMP PRIVATE(du_e, du_w, du_s, du_n, du_b, du_t) &
!$OMP PRIVATE(dv_e, dv_w, dv_s, dv_n, dv_b, dv_t) &
!$OMP PRIVATE(dw_e, dw_w, dw_s, dw_n, dw_b, dw_t) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t, EX, EY, EZ)

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
      cnv_u = 0.0
      cnv_v = 0.0
      cnv_w = 0.0
      
      ! 変数のロード
      include '../F_CBC/load_var_stencil5.h'
            
      bvx = bv(i,j,k)
      bpx = bp(i,j,k)
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> ce = 1.0(Normal) / 0.0(BC) 
      cw = 1.0
      ce = 1.0
      cs = 1.0
      cn = 1.0
      cb = 1.0
      ct = 1.0

      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) cw = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) ce = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) cs = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) cn = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) cb = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) ct = 0.0

      ! 延びたステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとす
      lmt_w = 1.0
      lmt_e = 1.0
      lmt_s = 1.0
      lmt_n = 1.0
      lmt_b = 1.0
      lmt_t = 1.0

      if ( (ibits(bv(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
      if ( (ibits(bv(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
      if ( (ibits(bv(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
      if ( (ibits(bv(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
      if ( (ibits(bv(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
      if ( (ibits(bv(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0
      
      ! 外部境界条件の場合
      if ( (i == 1)  .and. (ibits(bp(0   , j   , k   ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
      if ( (i == ix) .and. (ibits(bp(ix+1, j   , k   ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
      if ( (j == 1)  .and. (ibits(bp(i   , 0   , k   ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
      if ( (j == jx) .and. (ibits(bp(i   , jx+1, k   ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
      if ( (k == 1)  .and. (ibits(bp(i   , j   , 0   ), vbc_uwd, 1) == 1) ) lmt_b = 0.0
      if ( (k == kx) .and. (ibits(bp(i   , j   , kx+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0
      
      ! カット情報
			dww= cut(1,i-1,j  ,k  ) ! d_{i-1}^-
      dw = cut(1,i  ,j  ,k  ) ! d_{i}^-
      de = cut(2,i  ,j  ,k  ) ! d_{i}^+
      dee= cut(2,i+1,j  ,k  ) ! d_{i+1}^+
			
			dss= cut(3,i  ,j-1,k  ) ! d_{j-1}^-
      ds = cut(3,i  ,j  ,k  ) ! d_{j}^-
      dn = cut(4,i  ,j  ,k  ) ! d_{j}^+
      dnn= cut(4,i  ,j+1,k  ) ! d_{j+1}^+

			dbb= cut(5,i  ,j  ,k-1) ! d_{k-1}^-
      db = cut(5,i  ,j  ,k  ) ! d_{k}^-
      dt = cut(6,i  ,j  ,k  ) ! d_{k}^+
      dtt= cut(6,i  ,j  ,k+1) ! d_{k+1}^+
      
      ! 外挿時の分母 >  12*9 flop
      hww= 1.0 / (0.5 + dww) ! for (i-2)
      hw = 1.0 / (0.5 + dw ) ! for (i-1)
      he = 1.0 / (0.5 + de ) ! for (i+1)
      hee= 1.0 / (0.5 + dee) ! for (i+2)
      
      hss= 1.0 / (0.5 + dss) ! for (j-2)
      hs = 1.0 / (0.5 + ds ) ! for (j-1)
      hn = 1.0 / (0.5 + dn ) ! for (j+1)
      hnn= 1.0 / (0.5 + dnn) ! for (j+2)
      
      hbb= 1.0 / (0.5 + dbb) ! for (k-2)
      hb = 1.0 / (0.5 + db ) ! for (k-1)
      ht = 1.0 / (0.5 + dt ) ! for (k+1)
      htt= 1.0 / (0.5 + dtt) ! for (k+2)
      
      ! フラグと係数 > 18 flop
      qw = 1.0
      qe = 1.0
      qs = 1.0
      qn = 1.0
      qb = 1.0
      qt = 1.0
      qww= 1.0
      qee= 1.0
      qss= 1.0
      qnn= 1.0
      qbb= 1.0
      qtt= 1.0

      if (dw  < 1.0) qw = 0.0 
      if (de  < 1.0) qe = 0.0
      if (ds  < 1.0) qs = 0.0
      if (dn  < 1.0) qn = 0.0
      if (db  < 1.0) qb = 0.0
      if (dt  < 1.0) qt = 0.0
      if (dww < 1.0) qww= 0.0
      if (dee < 1.0) qee= 0.0
      if (dss < 1.0) qss= 0.0
      if (dnn < 1.0) qnn= 0.0
      if (dbb < 1.0) qbb= 0.0
      if (dtt < 1.0) qtt= 0.0

      r_qw = 1.0 - qw
      r_qe = 1.0 - qe
      r_qs = 1.0 - qs
      r_qn = 1.0 - qn
      r_qb = 1.0 - qb
      r_qt = 1.0 - qt
      r_qww= 1.0 - qww
      r_qee= 1.0 - qee
      r_qss= 1.0 - qss
      r_qnn= 1.0 - qnn
      r_qbb= 1.0 - qbb
      r_qtt= 1.0 - qtt

      ww = dw - 0.5
      we = de - 0.5
      ws = ds - 0.5
      wn = dn - 0.5
      wb = db - 0.5
      wt = dt - 0.5

      ! X方向  > 12 + 20 + 76 + 238 = 346 flop ---------------------------------------

      ! カットがない場合のセルフェイス速度　MUSCL内挿で参照する >  12 flop
      Uw_t = 0.5 * (Up0 + Uw1)
      Vw_t = 0.5 * (Vp0 + Vw1)
      Ww_t = 0.5 * (Wp0 + Ww1)
      Ue_t = 0.5 * (Up0 + Ue1)
      Ve_t = 0.5 * (Vp0 + Ve1)
      We_t = 0.5 * (Wp0 + We1)
      
      ! セルフェイスの移流速度 >  20 flop
      Uw_r = Uw_t * qw + r_qw * u_ref
      Ue_r = Ue_t * qe + r_qe * u_ref
      Uw_c = hw * ( ww * Ue_r + u_ref )
      Ue_c = he * ( we * Uw_r + u_ref )
      UPw  = ( Uw_t * qw + r_qw * Uw_c ) * qw
      UPe  = ( Ue_t * qe + r_qe * Ue_c ) * qe

      ! カットがある場合の参照セルの値の計算  >  19*4=76 flop
      ! u_{i-1}
      rw = dw - 1.0
      Uw1_t = (u_ref2 + rw * Ue_t) * hw
      Vw1_t = (v_ref2 + rw * Ve_t) * hw
      Ww1_t = (w_ref2 + rw * We_t) * hw
      Uw1 = qw * Uw1 + r_qw * Uw1_t
      Vw1 = qw * Vw1 + r_qw * Vw1_t
      Ww1 = qw * Ww1 + r_qw * Ww1_t
      
      ! u_{i+1}
      re = de - 1.0
      Ue1_t = (u_ref2 + re * Uw_t) * he
      Ve1_t = (v_ref2 + re * Vw_t) * he
      We1_t = (w_ref2 + re * Ww_t) * he
      Ue1 = qe * Ue1 + r_qe * Ue1_t
      Ve1 = qe * Ve1 + r_qe * Ve1_t
      We1 = qe * We1 + r_qe * We1_t

      ! u_{i-2}
      rww = dww - 1.0
      Uw2_t = (u_ref2 + rww * Uw_t) * hww
      Vw2_t = (v_ref2 + rww * Vw_t) * hww
      Ww2_t = (w_ref2 + rww * Ww_t) * hww
      Uw2 = qww * Uw2 + r_qww * Uw2_t
      Vw2 = qww * Vw2 + r_qww * Vw2_t
      Ww2 = qww * Ww2 + r_qww * Ww2_t

      ! u_{i+2}
      ree = dee - 1.0
      Ue2_t = (u_ref2 + ree * Ue_t) * hee
      Ve2_t = (v_ref2 + ree * Ve_t) * hee
      We2_t = (w_ref2 + ree * We_t) * hee
      Ue2 = qee * Ue2 + r_qee * Ue2_t
      Ve2 = qee * Ve2 + r_qee * Ve2_t
      We2 = qee * We2 + r_qee * We2_t
      
      ! 界面の左右の状態を計算  > 4 + (4+36+22+16)*3 = 238 flop
      ufw  = UPw - u_ref
      ufe  = UPe - u_ref
      afw = abs(ufw)
      afe = abs(ufe)
      
      dv4 = Ue2 - Ue1
      dv3 = Ue1 - Up0
      dv2 = Up0 - Uw1
      dv1 = Uw1 - Uw2
      
      include '../F_CBC/muscl.h'
      
      Uer = Ue1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_e
      Uel = Up0 + (cm1 * g3 + cm2 * g4) * ss_4
      Uwr = Up0 - (cm1 * g4 + cm2 * g3) * ss_4
      Uwl = Uw1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_w
      cnv_u = 0.5 * (ufe * (Uer + Uel) - afe * (Uer - Uel)) * ce &
            - 0.5 * (ufw * (Uwr + Uwl) - afw * (Uwr - Uwl)) * cw + cnv_u
      
      dv4 = Ve2 - Ve1
      dv3 = Ve1 - Vp0
      dv2 = Vp0 - Vw1
      dv1 = Vw1 - Vw2
      
      include '../F_CBC/muscl.h'
      
      Ver = Ve1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_e
      Vel = Vp0 + (cm1 * g3 + cm2 * g4) * ss_4
      Vwr = Vp0 - (cm1 * g4 + cm2 * g3) * ss_4
      Vwl = Vw1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_w
      cnv_v = 0.5 * (ufe * (Ver + Vel) - afe * (Ver - Vel)) * ce &
            - 0.5 * (ufw * (Vwr + Vwl) - afw * (Vwr - Vwl)) * cw + cnv_v
            
      dv4 = We2 - We1
      dv3 = We1 - Wp0
      dv2 = Wp0 - Ww1
      dv1 = Ww1 - Ww2
      
      include '../F_CBC/muscl.h'
      
      Wer = We1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_e
      Wel = Wp0 + (cm1 * g3 + cm2 * g4) * ss_4
      Wwr = Wp0 - (cm1 * g4 + cm2 * g3) * ss_4
      Wwl = Ww1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_w
      cnv_w = 0.5 * (ufe * (Wer + Wel) - afe * (Wer - Wel)) * ce &
            - 0.5 * (ufw * (Wwr + Wwl) - afw * (Wwr - Wwl)) * cw + cnv_w
			
      ! Y方向 ---------------------------------------
      
      ! セルフェイス位置
      Us_t = 0.5 * (Up0 + Us1)
      Vs_t = 0.5 * (Vp0 + Vs1)
      Ws_t = 0.5 * (Wp0 + Ws1)
      Un_t = 0.5 * (Up0 + Un1)
      Vn_t = 0.5 * (Vp0 + Vn1)
      Wn_t = 0.5 * (Wp0 + Wn1)
      
      ! セルフェイスの移流速度
      Vs_r = Vs_t * qs + r_qs * v_ref
      Vn_r = Vn_t * qn + r_qn * v_ref
      Vs_c = hs * ( ws * Vn_r + v_ref )
      Vn_c = hn * ( wn * Vs_r + v_ref )
      VPs  = ( Vs_t * qs + r_qs * Vs_c ) * qs
      VPn  = ( Vn_t * qn + r_qn * Vn_c ) * qn
      
      ! カットがある場合の参照セルの値の計算
      ! u_{j-1}
      rs = ds - 1.0
      Us1_t = (u_ref2 + rs * Un_t ) * hs
      Vs1_t = (v_ref2 + rs * Vn_t ) * hs
      Ws1_t = (w_ref2 + rs * Wn_t ) * hs
      Us1 = qs * Us1 + r_qs * Us1_t
      Vs1 = qs * Vs1 + r_qs * Vs1_t
      Ws1 = qs * Ws1 + r_qs * Ws1_t
      
      ! u_{j+1}
      rn = dn - 1.0
      Un1_t = (u_ref2 + rn * Us_t ) * hn
      Vn1_t = (v_ref2 + rn * Vs_t ) * hn
      Wn1_t = (w_ref2 + rn * Ws_t ) * hn
      Un1 = qn * Un1 + r_qn * Un1_t
      Vn1 = qn * Vn1 + r_qn * Vn1_t
      Wn1 = qn * Wn1 + r_qn * Wn1_t
    
      ! u_{j-2}
      rss = dss - 1.0
      Us2_t = (u_ref2 + rss * Us_t ) * hss
      Vs2_t = (v_ref2 + rss * Vs_t ) * hss
      Ws2_t = (w_ref2 + rss * Ws_t ) * hss
      Us2 = qss * Us2 + r_qss * Us2_t
      Vs2 = qss * Vs2 + r_qss * Vs2_t
      Ws2 = qss * Ws2 + r_qss * Ws2_t

      ! u_{j+2}
      rnn = dnn - 1.0
      Un2_t = (u_ref2 + rnn * Un_t ) * hnn
      Vn2_t = (v_ref2 + rnn * Vn_t ) * hnn
      Wn2_t = (w_ref2 + rnn * Wn_t ) * hnn
      Un2 = qnn * Un2 + r_qnn * Un2_t
      Vn2 = qnn * Vn2 + r_qnn * Vn2_t
      Wn2 = qnn * Wn2 + r_qnn * Wn2_t
      
      ! 界面の左右の状態を計算
      vfs  = VPs - v_ref
      vfn  = VPn - v_ref
      afs = abs(vfs)
      afn = abs(vfn)
      
      dv4 = Un2 - Un1
      dv3 = Un1 - Up0
      dv2 = Up0 - Us1
      dv1 = Us1 - Us2
      
      include '../F_CBC/muscl.h'
      
      Unr = Un1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_n
      Unl = Up0 + (cm1 * g3 + cm2 * g4) * ss_4
      Usr = Up0 - (cm1 * g4 + cm2 * g3) * ss_4
      Usl = Us1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_s
      cnv_u = 0.5 * (vfn * (Unr + Unl) - afn * (Unr - Unl)) * cn &
            - 0.5 * (vfs * (Usr + Usl) - afs * (Usr - Usl)) * cs + cnv_u
      
      dv4 = Vn2 - Vn1
      dv3 = Vn1 - Vp0
      dv2 = Vp0 - Vs1
      dv1 = Vs1 - Vs2
      
      include '../F_CBC/muscl.h'
      
      Vnr = Vn1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_n
      Vnl = Vp0 + (cm1 * g3 + cm2 * g4) * ss_4
      Vsr = Vp0 - (cm1 * g4 + cm2 * g3) * ss_4
      Vsl = Vs1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_s
      cnv_v = 0.5 * (vfn * (Vnr + Vnl) - afn * (Vnr - Vnl)) * cn &
            - 0.5 * (vfs * (Vsr + Vsl) - afs * (Vsr - Vsl)) * cs + cnv_v

      dv4 = Wn2 - Wn1
      dv3 = Wn1 - Wp0
      dv2 = Wp0 - Ws1
      dv1 = Ws1 - Ws2
      
      include '../F_CBC/muscl.h'
      
      Wnr = Wn1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_n
      Wnl = Wp0 + (cm1 * g3 + cm2 * g4) * ss_4
      Wsr = Wp0 - (cm1 * g4 + cm2 * g3) * ss_4
      Wsl = Ws1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_s
      cnv_w = 0.5 * (vfn * (Wnr + Wnl) - afn * (Wnr - Wnl)) * cn &
            - 0.5 * (vfs * (Wsr + Wsl) - afs * (Wsr - Wsl)) * cs + cnv_w
			
      ! Z方向 ---------------------------------------
      
      ! セルフェイス位置の値を内挿
      Ub_t = 0.5 * (Up0 + Ub1)
      Vb_t = 0.5 * (Vp0 + Vb1)
      Wb_t = 0.5 * (Wp0 + Wb1)
      Ut_t = 0.5 * (Up0 + Ut1)
      Vt_t = 0.5 * (Vp0 + Vt1)
      Wt_t = 0.5 * (Wp0 + Wt1)
      
      ! セルフェイスの移流速度
      Wb_r = Wb_t * qb + r_qb * w_ref
      Wt_r = Wt_t * qt + r_qt * w_ref
      Wb_c = hb * ( wb * Wt_r + w_ref )
      Wt_c = ht * ( wt * Wb_r + w_ref )
      WPb  = ( Wb_t * qb + r_qb * Wb_c ) * qb
      WPt  = ( Wt_t * qt + r_qt * Wt_c ) * qt
      
      ! カットがある場合の参照セルの値の計算
      ! u_{k-1}
      rb = db - 1.0
      Ub1_t = (u_ref2 + rb * Ut_t ) * hb
      Vb1_t = (v_ref2 + rb * Vt_t ) * hb
      Wb1_t = (w_ref2 + rb * Wt_t ) * hb
      Ub1 = qb * Ub1 + r_qb * Ub1_t
      Vb1 = qb * Vb1 + r_qb * Vb1_t
      Wb1 = qb * Wb1 + r_qb * Wb1_t
      
      ! u_{k+1}
      rt = dt - 1.0
      Ut1_t = (u_ref2 + rt * Ub_t ) * ht
      Vt1_t = (v_ref2 + rt * Vb_t ) * ht
      Wt1_t = (w_ref2 + rt * Wb_t ) * ht
      Ut1 = qt * Ut1 + r_qt * Ut1_t
      Vt1 = qt * Vt1 + r_qt * Vt1_t
      Wt1 = qt * Wt1 + r_qt * Wt1_t
      
      ! u_{k-2}
      rbb = dbb - 1.0
      Ub2_t = (u_ref2 + rbb * Ub_t ) * hbb
      Vb2_t = (v_ref2 + rbb * Vb_t ) * hbb
      Wb2_t = (w_ref2 + rbb * Wb_t ) * hbb
      Ub2 = qbb * Ub2 + r_qbb * Ub2_t
      Vb2 = qbb * Vb2 + r_qbb * Vb2_t
      Wb2 = qbb * Wb2 + r_qbb * Wb2_t

      ! u_{k+2}
      rtt = dtt - 1.0
      Ut2_t = (u_ref2 + rtt * Ut_t ) * htt
      Vt2_t = (v_ref2 + rtt * Vt_t ) * htt
      Wt2_t = (w_ref2 + rtt * Wt_t ) * htt
      Ut2 = qtt * Ut2 + r_qtt * Ut2_t
      Vt2 = qtt * Vt2 + r_qtt * Vt2_t
      Wt2 = qtt * Wt2 + r_qtt * Wt2_t
      
      ! 界面の左右の状態を計算
      wfb  = WPb - w_ref
      wft  = WPt - w_ref
      afb = abs(wfb)
      aft = abs(wft)

      dv4 = Ut2 - Ut1
      dv3 = Ut1 - Up0
      dv2 = Up0 - Ub1
      dv1 = Ub1 - Ub2
      
      include '../F_CBC/muscl.h'
      
      Utr = Ut1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_t
      Utl = Up0 + (cm1 * g3 + cm2 * g4) * ss_4
      Ubr = Up0 - (cm1 * g4 + cm2 * g3) * ss_4
      Ubl = Ub1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_b
      cnv_u = 0.5 * (wft * (Utr + Utl) - aft * (Utr - Utl)) * ct &
            - 0.5 * (wfb * (Ubr + Ubl) - afb * (Ubr - Ubl)) * cb + cnv_u

      dv4 = Vt2 - Vt1
      dv3 = Vt1 - Vp0
      dv2 = Vp0 - Vb1
      dv1 = Vb1 - Vb2
      
      include '../F_CBC/muscl.h'
      
      Vtr = Vt1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_t
      Vtl = Vp0 + (cm1 * g3 + cm2 * g4) * ss_4
      Vbr = Vp0 - (cm1 * g4 + cm2 * g3) * ss_4
      Vbl = Vb1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_b
      cnv_v = 0.5 * (wft * (Vtr + Vtl) - aft * (Vtr - Vtl)) * ct &
            - 0.5 * (wfb * (Vbr + Vbl) - afb * (Vbr - Vbl)) * cb + cnv_v

      dv4 = Wt2 - Wt1
      dv3 = Wt1 - Wp0
      dv2 = Wp0 - Wb1
      dv1 = Wb1 - Wb2
      
      include '../F_CBC/muscl.h'
      
      Wtr = Wt1 - (cm1 * g6 + cm2 * g5) * ss_4 * lmt_t
      Wtl = Wp0 + (cm1 * g3 + cm2 * g4) * ss_4
      Wbr = Wp0 - (cm1 * g4 + cm2 * g3) * ss_4
      Wbl = Wb1 + (cm1 * g1 + cm2 * g2) * ss_4 * lmt_b
      cnv_w = 0.5 * (wft * (Wtr + Wtl) - aft * (Wtr - Wtl)) * ct &
            - 0.5 * (wfb * (Wbr + Wbl) - afb * (Wbr - Wbl)) * cb + cnv_w

      
      ! 粘性項の計算  >  18 + 10*3 + 9 = 57 flop -----------------------------------------------
      du_e = ( Ue1 - Up0 )
      du_w = ( Up0 - Uw1 )
      du_n = ( Un1 - Up0 )
      du_s = ( Up0 - Us1 )
      du_t = ( Ut1 - Up0 )
      du_b = ( Up0 - Ub1 )
      
      dv_e = ( Ve1 - Vp0 )
      dv_w = ( Vp0 - Vw1 )
      dv_n = ( Vn1 - Vp0 )
      dv_s = ( Vp0 - Vs1 )
      dv_t = ( Vt1 - Vp0 )
      dv_b = ( Vp0 - Vb1 )
      
      dw_e = ( We1 - Wp0 )
      dw_w = ( Wp0 - Ww1 )
      dw_n = ( Wn1 - Wp0 )
      dw_s = ( Wp0 - Ws1 )
      dw_t = ( Wt1 - Wp0 )
      dw_b = ( Wp0 - Wb1 )
      
      Ex = ( (du_e * ce - du_w * cw)  &
           + (du_n * cn - du_s * cs)  &
           + (du_t * ct - du_b * cb)  &
           ) * dh2
           
      Ey = ( (dv_e * ce - dv_w * cw)  &
           + (dv_n * cn - dv_s * cs)  &
           + (dv_t * ct - dv_b * cb)  &
           ) * dh2
           
      Ez = ( (dw_e * ce - dw_w * cw)  &
           + (dw_n * cn - dw_s * cs)  &
           + (dw_t * ct - dw_b * cb)  &
           ) * dh2
			
      ! 対流項と粘性項の和
      wv(1,i,j,k) = -cnv_u * dh1 + EX * vcs
      wv(2,i,j,k) = -cnv_v * dh1 + EY * vcs
      wv(3,i,j,k) = -cnv_w * dh1 + EZ * vcs
      
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

		return
    end subroutine cds_pvec_muscl

!  **************************************************************************************
!> @subroutine cds_update_vec (v, div, sz, g, delta_t, dh, vc, p, bp, bv, cut, v00, flop)
!! @brief 次のステップのセルセンターの速度を更新し，発散値を計算する
!! @param[out] v n+1時刻のセルセンターの速度ベクトル
!! @param div 発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param delta_t 時間積分幅
!! @param dh 格子幅
!! @param vc セルセンタの疑似速度ベクトル
!! @param p 圧力
!! @param bp BCindex P
!! @param bv BCindex V
!! @param cut カット情報(float)
!! @param v00 参照速度
!! @param[out] flop flop count
!<
    subroutine cds_update_vec (v, div, sz, g, delta_t, dh, vc, p, bp, bv, cut, v00, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bpx, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  dh, delta_t, dd, flop, coef, actv, r_actv
    real                                                        ::  pc, gpx, gpy, gpz
    real                                                        ::  gpw_r, gpe_r, gps_r, gpn_r, gpb_r, gpt_r
    real                                                        ::  gpw_c, gpe_c, gps_c, gpn_c, gpb_c, gpt_c
    real                                                        ::  gpw, gpe, gps, gpn, gpb, gpt
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                        ::  Ue_f, Uw_f, Vn_f, Vs_f, Wt_f, Wb_f
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  Ue_r, Uw_r, Vn_r, Vs_r, Wt_r, Wb_r
    real                                                        ::  Ue_c, Uw_c, Vn_c, Vs_c, Wt_c, Wb_c
    real                                                        ::  dw, de, ds, dn, db, dt
    real                                                        ::  ww, we, ws, wn, wb, wt
    real                                                        ::  qw, qe, qs, qn, qb, qt
    real                                                        ::  r_qw, r_qe, r_qs, r_qn, r_qb, r_qt
    real                                                        ::  cw, ce, cs, cn, cb, ct
    real                                                        ::  hw, he, hs, hn, hb, ht
    real                                                        ::  N_e, N_w, N_n, N_s, N_t, N_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  p
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, vc
    real*4, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp, bv
    real, dimension(0:3)                                        ::  v00
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dd = delta_t/dh
    coef = dh/delta_t
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    flop = flop + real(ix)*real(jx)*real(kx)*228.0 + 16.0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dd, u_ref, v_ref, w_ref, coef) &
!$OMP PRIVATE(bpx, bvx, actv, r_actv) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(dw, de, ds, dn, db, dt) &
!$OMP PRIVATE(ww, we, ws, wn, wb, wt) &
!$OMP PRIVATE(qw, qe, qs, qn, qb, qt) &
!$OMP PRIVATE(r_qw, r_qe, r_qs, r_qn, r_qb, r_qt) &
!$OMP PRIVATE(cw, ce, cs, cn, cb, ct) &
!$OMP PRIVATE(hw, he, hs, hn, hb, ht) &
!$OMP PRIVATE(gpw_r, gpe_r, gps_r, gpn_r, gpb_r, gpt_r) &
!$OMP PRIVATE(gpw_c, gpe_c, gps_c, gpn_c, gpb_c, gpt_c) &
!$OMP PRIVATE(gpw, gpe, gps, gpn, gpb, gpt) &
!$OMP PRIVATE(pc, gpx, gpy, gpz) &
!$OMP PRIVATE(Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t) &
!$OMP PRIVATE(Ue_r, Uw_r, Vn_r, Vs_r, Wt_r, Wb_r) &
!$OMP PRIVATE(Ue_c, Uw_c, Vn_c, Vs_c, Wt_c, Wb_c) &
!$OMP PRIVATE(Ue_f, Uw_f, Vn_f, Vs_f, Wt_f, Wb_f)

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
      bpx = bp(i,j,k)
      bvx = bv(i,j,k)
      actv = real(ibits(bvx, State, 1))
      r_actv = 1.0 - actv
      
      ! Neumann条件のとき，0.0
      N_w = real(ibits(bpx, bc_n_W, 1))  ! w
      N_e = real(ibits(bpx, bc_n_E, 1))  ! e
      N_s = real(ibits(bpx, bc_n_S, 1))  ! s
      N_n = real(ibits(bpx, bc_n_N, 1))  ! n
      N_b = real(ibits(bpx, bc_n_B, 1))  ! b
      N_t = real(ibits(bpx, bc_n_T, 1))  ! t
      
      Up0 = vc(1, i  ,j  ,k  )
      Vp0 = vc(2, i  ,j  ,k  )
      Wp0 = vc(3, i  ,j  ,k  )
      Uw0 = vc(1, i-1,j  ,k  )
      Ue0 = vc(1, i+1,j  ,k  )
      Vs0 = vc(2, i  ,j-1,k  )
      Vn0 = vc(2, i  ,j+1,k  )
      Wb0 = vc(3, i  ,j  ,k-1)
      Wt0 = vc(3, i  ,j  ,k+1)
      
      ! 距離情報 [0, 1]
      dw = cut(1,i,j,k) ! d_{i}^-
      de = cut(2,i,j,k) ! d_{i}^+
      ds = cut(3,i,j,k) ! d_{j}^-
      dn = cut(4,i,j,k) ! d_{j}^+
      db = cut(5,i,j,k) ! d_{k}^-
      dt = cut(6,i,j,k) ! d_{k}^+

      ww = dw - 0.5
      we = de - 0.5
      ws = ds - 0.5
      wn = dn - 0.5
      wb = db - 0.5
      wt = dt - 0.5
      
      ! if cut -> q=0.0
      qw = 1.0
      qe = 1.0
      qs = 1.0
      qn = 1.0
      qb = 1.0
      qt = 1.0

      if (dw < 1.0) qw = 0.0
      if (de < 1.0) qe = 0.0
      if (ds < 1.0) qs = 0.0
      if (dn < 1.0) qn = 0.0
      if (db < 1.0) qb = 0.0
      if (dt < 1.0) qt = 0.0

      r_qw = 1.0 - qw
      r_qe = 1.0 - qe
      r_qs = 1.0 - qs
      r_qn = 1.0 - qn
      r_qb = 1.0 - qb
      r_qt = 1.0 - qt
      
      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
      cw = 1.0
      ce = 1.0
      cs = 1.0
      cn = 1.0
      cb = 1.0
      ct = 1.0

      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) cw = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) ce = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) cs = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) cn = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) cb = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) ct = 0.0
      
      hw = 1.0 / (dw + 0.5)
      he = 1.0 / (de + 0.5)
      hs = 1.0 / (ds + 0.5)
      hn = 1.0 / (dn + 0.5)
      hb = 1.0 / (db + 0.5)
      ht = 1.0 / (dt + 0.5)
      
      ! Reference pressure gradient
      pc    = p(i,  j,  k  )
      gpw_r = (pc - p(i-1,j  ,k  )) * N_w
      gpe_r = (p(i+1,j  ,k  ) - pc) * N_e
      gps_r = (pc - p(i  ,j-1,k  )) * N_s
      gpn_r = (p(i  ,j+1,k  ) - pc) * N_n
      gpb_r = (pc - p(i  ,j  ,k-1)) * N_b
      gpt_r = (p(i  ,j  ,k+1) - pc) * N_t
      
      ! Correction of pressure gradient by wall effect
      gpw_c = hw * ww * gpe_r
      gpe_c = he * we * gpw_r
      gps_c = hs * ws * gpn_r
      gpn_c = hn * wn * gps_r
      gpb_c = hb * wb * gpt_r
      gpt_c = ht * wt * gpb_r
      
      ! Cell face gradient
      gpw = gpw_r * qw + r_qw * gpw_c
      gpe = gpe_r * qe + r_qe * gpe_c
      gps = gps_r * qs + r_qs * gps_c
      gpn = gpn_r * qn + r_qn * gpn_c
      gpb = gpb_r * qb + r_qb * gpb_c
      gpt = gpt_r * qt + r_qt * gpt_c
      
      ! Cell center gradient
      gpx = 0.5 * (gpe + gpw)
      gpy = 0.5 * (gpn + gps)
      gpz = 0.5 * (gpt + gpb)
      
      ! セルフェイス速度　カットなし
      Uw_t = 0.5 * (Uw0 + Up0)
      Ue_t = 0.5 * (Ue0 + Up0)
      Vs_t = 0.5 * (Vs0 + Vp0)
      Vn_t = 0.5 * (Vn0 + Vp0)
      Wb_t = 0.5 * (Wb0 + Wp0)
      Wt_t = 0.5 * (Wt0 + Wp0)

      ! reference side at cell face
      Uw_r = Uw_t * qw + r_qw * u_ref
      Ue_r = Ue_t * qe + r_qe * u_ref
      Vs_r = Vs_t * qs + r_qs * v_ref
      Vn_r = Vn_t * qn + r_qn * v_ref
      Wb_r = Wb_t * qb + r_qb * w_ref
      Wt_r = Wt_t * qt + r_qt * w_ref
      
      ! Correction of cell face velocity to be interpolated
      Uw_c = hw * ( ww * Ue_r + u_ref )
      Ue_c = he * ( we * Uw_r + u_ref )
      Vs_c = hs * ( ws * Vn_r + v_ref )
      Vn_c = hn * ( wn * Vs_r + v_ref )
      Wb_c = hb * ( wb * Wt_r + w_ref )
      Wt_c = ht * ( wt * Wb_r + w_ref )
      
      ! Cell face velocity
      Uw_f = Uw_t * qw + r_qw * Uw_c
      Ue_f = Ue_t * qe + r_qe * Ue_c
      Vs_f = Vs_t * qs + r_qs * Vs_c
      Vn_f = Vn_t * qn + r_qn * Vn_c
      Wb_f = Wb_t * qb + r_qb * Wb_c
      Wt_f = Wt_t * qt + r_qt * Wt_c
      
      ! 発散値 VBCの寄与は除外
      div(i,j,k) = ((Ue_f - dd * gpe) * ce &
                   -(Uw_f - dd * gpw) * cw &
                   +(Vn_f - dd * gpn) * cn &
                   -(Vs_f - dd * gps) * cs &
                   +(Wt_f - dd * gpt) * ct &
                   -(Wb_f - dd * gpb) * cb ) * coef * actv
      
      ! セルセンタの速度更新
      v(1,i,j,k) = ( Up0 - gpx * dd ) * actv + r_actv * u_ref
      v(2,i,j,k) = ( Vp0 - gpy * dd ) * actv + r_actv * v_ref
      v(3,i,j,k) = ( Wp0 - gpz * dd ) * actv + r_actv * w_ref

    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
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
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                        ::  Ue_f, Uw_f, Vn_f, Vs_f, Wt_f, Wb_f
    real                                                        ::  Ue_r, Uw_r, Vn_r, Vs_r, Wt_r, Wb_r
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  Ue_c, Uw_c, Vn_c, Vs_c, Wt_c, Wb_c
    real                                                        ::  cw, ce, cs, cn, cb, ct
    real                                                        ::  dw, de, ds, dn, db, dt
    real                                                        ::  hw, he, hs, hn, hb, ht
    real                                                        ::  ww, we, ws, wn, wb, wt
    real                                                        ::  qw, qe, qs, qn, qb, qt
    real                                                        ::  r_qw, r_qe, r_qs, r_qn, r_qb, r_qt
    real                                                        ::  coef, flop, actv, r_actv, u_ref, v_ref, w_ref
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real*4, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(0:3)                                        ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    flop = flop + real(ix)*real(jx)*real(kx)*152.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef, u_ref, v_ref, w_ref) &
!$OMP PRIVATE(bvx, actv, r_actv) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(dw, de, ds, dn, db, dt) &
!$OMP PRIVATE(ww, we, ws, wn, wb, wt) &
!$OMP PRIVATE(qw, qe, qs, qn, qb, qt) &
!$OMP PRIVATE(r_qw, r_qe, r_qs, r_qn, r_qb, r_qt) &
!$OMP PRIVATE(cw, ce, cs, cn, cb, ct) &
!$OMP PRIVATE(hw, he, hs, hn, hb, ht) &
!$OMP PRIVATE(Ue_r, Uw_r, Vn_r, Vs_r, Wt_r, Wb_r) &
!$OMP PRIVATE(Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t) &
!$OMP PRIVATE(Ue_c, Uw_c, Vn_c, Vs_c, Wt_c, Wb_c) &
!$OMP PRIVATE(Ue_f, Uw_f, Vn_f, Vs_f, Wt_f, Wb_f)

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
      bvx = bv(i,j,k)
      actv= real(ibits(bvx, State, 1))
      
      Up0 = v(1, i  ,j  ,k  )
      Vp0 = v(2, i  ,j  ,k  )
      Wp0 = v(3, i  ,j  ,k  )
      Uw0 = v(1, i-1,j  ,k  )
      Ue0 = v(1, i+1,j  ,k  )
      Vs0 = v(2, i  ,j-1,k  )
      Vn0 = v(2, i  ,j+1,k  )
      Wb0 = v(3, i  ,j  ,k-1)
      Wt0 = v(3, i  ,j  ,k+1)
      
      dw = cut(1,i,j,k) ! d_{i}^-
      de = cut(2,i,j,k) ! d_{i}^+
      ds = cut(3,i,j,k) ! d_{j}^-
      dn = cut(4,i,j,k) ! d_{j}^+
      db = cut(5,i,j,k) ! d_{k}^-
      dt = cut(6,i,j,k) ! d_{k}^+

      ww = dw - 0.5
      we = de - 0.5
      ws = ds - 0.5
      wn = dn - 0.5
      wb = db - 0.5
      wt = dt - 0.5
      
      ! if cut -> q=0.0
      qw = 1.0
      qe = 1.0
      qs = 1.0
      qn = 1.0
      qb = 1.0
      qt = 1.0

      if (dw < 1.0) qw = 0.0 
      if (de < 1.0) qe = 0.0
      if (ds < 1.0) qs = 0.0
      if (dn < 1.0) qn = 0.0
      if (db < 1.0) qb = 0.0
      if (dt < 1.0) qt = 0.0

      r_qw = 1.0 - qw
      r_qe = 1.0 - qe
      r_qs = 1.0 - qs
      r_qn = 1.0 - qn
      r_qb = 1.0 - qb
      r_qt = 1.0 - qt
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c = 1.0(Normal) / 0.0(BC)
      cw = 1.0
      ce = 1.0
      cs = 1.0
      cn = 1.0
      cb = 1.0
      ct = 1.0

      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) cw = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) ce = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) cs = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) cn = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) cb = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) ct = 0.0
      
      hw = 1.0 / (dw + 0.5)
      he = 1.0 / (de + 0.5)
      hs = 1.0 / (ds + 0.5)
      hn = 1.0 / (dn + 0.5)
      hb = 1.0 / (db + 0.5)
      ht = 1.0 / (dt + 0.5)

      ! セルフェイス速度　カットなし
      Uw_t = 0.5 * (Uw0 + Up0)
      Ue_t = 0.5 * (Ue0 + Up0)
      Vs_t = 0.5 * (Vs0 + Vp0)
      Vn_t = 0.5 * (Vn0 + Vp0)
      Wb_t = 0.5 * (Wb0 + Wp0)
      Wt_t = 0.5 * (Wt0 + Wp0)
      
      ! reference side at cell face
      Uw_r = Uw_t * qw + r_qw * u_ref
      Ue_r = Ue_t * qe + r_qe * u_ref
      Vs_r = Vs_t * qs + r_qs * v_ref
      Vn_r = Vn_t * qn + r_qn * v_ref
      Wb_r = Wb_t * qb + r_qb * w_ref
      Wt_r = Wt_t * qt + r_qt * w_ref

      ! Correction of cell face velocity to be interpolated
      Uw_c = hw * ( ww * Ue_r + u_ref )
      Ue_c = he * ( we * Uw_r + u_ref )
      Vs_c = hs * ( ws * Vn_r + v_ref )
      Vn_c = hn * ( wn * Vs_r + v_ref )
      Wb_c = hb * ( wb * Wt_r + w_ref )
      Wt_c = ht * ( wt * Wb_r + w_ref )

      ! Cell face
      Uw_f = Uw_t * qw + r_qw * Uw_c
      Ue_f = Ue_t * qe + r_qe * Ue_c
      Vs_f = Vs_t * qs + r_qs * Vs_c
      Vn_f = Vn_t * qn + r_qn * Vn_c
      Wb_f = Wb_t * qb + r_qb * Wb_c
      Wt_f = Wt_t * qt + r_qt * Wt_c
      
      ! VBCの場合には寄与をキャンセル
      div(i,j,k) = ( (Ue_f - u_ref) * ce  &
                   - (Uw_f - u_ref) * cw  &
                   + (Vn_f - v_ref) * cn  &
                   - (Vs_f - v_ref) * cs  &
                   + (Wt_f - w_ref) * ct  &
                   - (Wb_f - w_ref) * cb ) * coef * actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

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
